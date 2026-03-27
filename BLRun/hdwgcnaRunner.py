import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for hdWGCNA.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    if not RunnerObj.inputDir.joinpath("hdWGCNA").exists():
        print("Input folder for hdWGCNA does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("hdWGCNA").mkdir(exist_ok = False)
    
    if not RunnerObj.inputDir.joinpath("hdWGCNA/ExpressionData.csv").exists():
        # input data
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        
        # Write .csv file
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath("hdWGCNA/ExpressionData.csv"),
                             sep = ',', header  = True, index = True)

def run(RunnerObj):
    '''
    Function to run hdWGCNA algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/hdWGCNA/ExpressionData.csv"
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/hdWGCNA/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath =  "data/" + str(outDir) + 'outFile.txt'

    run_tf = str(RunnerObj.params['run_tf'])
    nThreads = str(RunnerObj.params['nThreads'])
    nVariableFeatures = str(RunnerObj.params['nVariableFeatures'])
    clustering_resolution = str(RunnerObj.params['clustering_resolution'])
    metacells_k = str(RunnerObj.params['metacells_k'])
    min_cells_per_metacell = str(RunnerObj.params['min_cells_per_metacell'])
    max_shared_cells = str(RunnerObj.params['max_shared_cells'])
    networkType = str(RunnerObj.params['networkType'])
    soft_power = str(RunnerObj.params['soft_power'])
    delete_grey_module = str(RunnerObj.params['delete_grey_module'])

    #cmdToRun = ' '.join(['docker run --rm -v', 
    #                     str(Path.cwd())+':/data/',  # <-- add the missing space and slash
    #                     'hdWGCNA:base',          # <-- image name as a separate argument
    #                     '/bin/sh -c \"time -v -o',
    #                     "data/" + str(outDir) + 'time.txt', 
    #                     'Rscript runhdWGCNA.R',
    #                     inputPath, 
    #                     outPath, '\"'])

    cmdToRun = ' '.join(['/software/singularity/3.11.4//bin/singularity exec --cleanenv -B',
                        str(Path.cwd()) + ':/data',
                        str(Path.cwd() / 'Algorithms/TF_hdWGCNA/tf_hdwgcna_1.0.sif'),
                        #'/bin/sh -c "echo OK > /data/'+ str(outDir) +'test2.txt"'
                        '/bin/sh -c \"time -v -o',
                        '/data/' + str(outDir) + 'time.txt',
                        'Rscript /data/Algorithms/TF_hdWGCNA/runTF_hdWGCNA.R',
                        inputPath, outPath, 
                        run_tf, nThreads, nVariableFeatures, clustering_resolution, metacells_k, min_cells_per_metacell, max_shared_cells, networkType, soft_power, delete_grey_module, '\"']) # for singularity

    print(cmdToRun)
    os.system(cmdToRun)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from hdWGCNA

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/hdWGCNA/"
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
    else:
        print("Output file found, proceeding in formating file...")
    
    print("Output directory for hdWGCNA:", outDir)

    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = '\t', header = 0)

    # Formats the outFile if necessary
    required_cols = ['Gene1', 'Gene2', 'EdgeWeight']

    if not all(col in OutDF.columns for col in required_cols):
        outFile = OutDF[["tf", "gene", "Cor"]].copy()
        outFile.columns = ['Gene1', 'Gene2', 'EdgeWeight']
    else:
        outFile = OutDF[['Gene1', 'Gene2', 'EdgeWeight']].copy()
    
    outFile.to_csv(f'{outDir}/rankedEdges.csv', sep='\t', index=False)
