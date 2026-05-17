import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from itertools import product, permutations, combinations, combinations_with_replacement
from tqdm import tqdm
import networkx as nx
from joblib import Parallel, delayed # Added for parallelization

def _process_algo(algo, outDir, EdgeDict, numEdges, refFB, refFF, refMI):
    '''
    Worker function: processes a single algorithm and returns its motif ratios.
 
    Isolated so that joblib can dispatch it independently for each algorithm.
    All arguments are either primitives or thread-safe read-only structures,
    which is why we use prefer="threads" in the Parallel call — no pickling
    needed, and NetworkX releases the GIL during its C-level graph traversals.
 
    :param algo:     One entry from inputSettings.algorithms (list/tuple).
    :param outDir:   Output directory string.
    :param EdgeDict: Dict of all possible directed edges (read-only).
    :param numEdges: Number of edges in the reference network.
    :param refFB:    Reference feedback-loop count (already floored to 1).
    :param refFF:    Reference feedforward-loop count (already floored to 1).
    :param refMI:    Reference mutual-interaction count (already floored to 1).
 
    :returns: (algo_name, fbl_ratio, ffl_ratio, mi_ratio)
    '''
    # PPCOR and PIDC produce undirected scores; skip them (same as original).
    if algo[0] in ('PPCOR', 'PIDC'):
        return None
 
    rankedPath = Path(outDir + '/' + algo[0] + '/rankedEdges.csv')
 
    if not rankedPath.exists():
        print(rankedPath, ' does not exist. Skipping...')
        return (algo[0], 0, 0, 0)
 
    predDF = pd.read_csv(rankedPath, sep='\t', header=0, index_col=None)
    predDF = predDF.loc[predDF['Gene1'] != predDF['Gene2']]
    predDF.drop_duplicates(keep='first', inplace=True)
    predDF.reset_index(drop=True, inplace=True)
 
    if predDF.shape[0] == 0:
        return (algo[0], 0, 0, 0)
 
    # Keep top-k edges by edge weight (same threshold logic as original).
    predDF.EdgeWeight = predDF.EdgeWeight.round(6).abs()
    maxk = min(predDF.shape[0], numEdges)
    edgeWeightTopk = predDF.iloc[maxk - 1].EdgeWeight
    nonZeroMin = np.nanmin(predDF.EdgeWeight.replace(0, np.nan).values)
    bestVal = max(nonZeroMin, edgeWeightTopk)
    newDF = predDF.loc[predDF['EdgeWeight'] >= bestVal]
 
    # Build the predicted directed graph.
    predGraph = nx.DiGraph()
    for key, (u, v) in ((k, k.split('|')) for k in EdgeDict):
        if len(newDF.loc[(newDF['Gene1'] == u) & (newDF['Gene2'] == v)]) > 0:
            predGraph.add_edge(u, v)
 
    fb, ff, mi = getNetProp(predGraph)
    return (algo[0], fb / refFB, ff / refFF, mi / refMI)
 
 
def Motifs(datasetDict, inputSettings, n_jobs=-1):
    '''
    Computes ratios of the counts of various network motifs for each algorithm
    for a given dataset. The ratios are computed by dividing the counts of
    various network motifs in the predicted top-k network to their respective
    values in the reference network.
 
    Parallelization strategy
    ------------------------
    Each algorithm's ranked-edge file is processed independently inside
    `_process_algo`. We dispatch those calls in parallel with joblib using
    threads (prefer="threads") so that:
      - NetworkX objects don't need to be pickled (no multiprocessing overhead).
      - NetworkX's C extensions release the GIL, giving real concurrency for
        the graph-traversal bottleneck (cycle detection, path enumeration).
 
    :param datasetDict:   A dictionary containing the dataset name and path to
                          the reference network.
    :type datasetDict: dict
 
    :param inputSettings: An object of class :class:`BLEval.InputSettings`.
    :type inputSettings: :class:`BLEval.InputSettings`
 
    :param n_jobs: Number of parallel jobs. -1 uses all available cores.
                   Passed directly to joblib.Parallel.
    :type n_jobs: int
 
    :returns:
        - FBL: A dataframe containing ratios of three-node feedback loop motifs
        - FFL: A dataframe containing ratios of three-node feedforward loop motifs
        - MI:  A dataframe containing ratios of two-node mutual interaction motifs
    '''
 
    # --- Reference network ---
    trueEdgesDF = pd.read_csv(
        str(inputSettings.datadir) + '/' + datasetDict['name'] + '/' + datasetDict['trueEdges'],
        sep=',', header=0, index_col=None
    )
 
    possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:, ['Gene1', 'Gene2']]), r=2))
    EdgeDict = {'|'.join(p): 0 for p in possibleEdges}
 
    refGraph = nx.DiGraph()
    for key in EdgeDict:
        u, v = key.split('|')
        if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == u) & (trueEdgesDF['Gene2'] == v)]) > 0:
            refGraph.add_edge(u, v)
 
    numEdges = len(refGraph.edges())
    refFB, refFF, refMI = getNetProp(refGraph)
 
    # Floor at 1 to avoid division-by-zero when a motif is absent in the
    # reference network.
    refFB = refFB or 1
    refFF = refFF or 1
    refMI = refMI or 1
 
    outDir = "outputs/" + str(inputSettings.datadir).split("inputs/")[1] + '/' + datasetDict['name']
 
    # --- Parallel dispatch over algorithms ---
    # tqdm wraps the input list so the progress bar still works.
    results = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(_process_algo)(algo, outDir, EdgeDict, numEdges, refFB, refFF, refMI)
        for algo in tqdm(inputSettings.algorithms, total=len(inputSettings.algorithms), unit=" Algorithms")
    )
 
    # --- Collect results ---
    dataDict = {'FBL': {}, 'FFL': {}, 'Mutual': {}}
    for result in results:
        if result is None:          # PPCOR / PIDC were skipped
            continue
        name, fbl, ffl, mi = result
        dataDict['FBL'][name]    = fbl
        dataDict['FFL'][name]    = ffl
        dataDict['Mutual'][name] = mi
 
    dataDF = pd.DataFrame(dataDict)
    return dataDF['FBL'], dataDF['FFL'], dataDF['Mutual']

    
def getNetProp(inGraph):
    '''
    A helper function to compute
    counts of various network motifs.
    

    :param inGraph: An graph object of class :class:`networkx.DiGraph`.
    :type inGraph: :obj:networkx.DiGraph

    :returns:
        - A value corresponding to the number of three-node feedback loops
        - A value corresponding to the number of three-node feedforward loops
        - A value corresponding to the number of two-node mutual interaction


    '''

    # number of weakly connected components in 
    # reference network
    # numCC = len(list(nx.weakly_connected_components(inGraph)))
    
    # number of feedback loop 
    # in reference network
    allCyc = nx.simple_cycles(inGraph)
    cycSet = set()
    for cyc in allCyc:
        if len(cyc) == 3:
            cycSet.add(frozenset(cyc))
    
    numFB = len(cycSet)
    
    # number of feedfwd loops
    # in reference network
    allPaths = []
    allPathsSet = set()   
    for u,v in inGraph.edges():
        allPaths = nx.all_simple_paths(inGraph, u, v, cutoff=2)
        for p in allPaths:
            if len(p) > 2:
                allPathsSet.add(frozenset(p))
                
    numFF= len(allPathsSet)
    
    
    # number of mutual interactions
    numMI = 0.0
    for u,v in inGraph.edges():
        if (v,u) in inGraph.edges():
            numMI += 0.5

    return numFB, numFF, numMI
