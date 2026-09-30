import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for GENIE3.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    if not RunnerObj.inputDir.joinpath("GENIE3").exists():
        print("Input folder for GENIE3 does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("GENIE3").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("GENIE3/ExpressionData.csv").exists():
        # input data
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)

        # Write .csv file
        ExpressionData.T.to_csv(RunnerObj.inputDir.joinpath("GENIE3/ExpressionData.csv"),
                             sep = '\t', header  = True, index = True)
    
def run(RunnerObj):
    '''
    Function to run GENIE3 algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    inputPath = "/data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/GENIE3/ExpressionData.csv"
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/GENIE3/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath = "/data/" +  str(outDir) + 'outFile.txt'
    
    # cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/data/ --expose=41269', 
    #                      'grnbeeline/arboreto:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time.txt', 'python runArboreto.py --algo=GENIE3',
    #                      '--inFile='+inputPath, '--outFile='+outPath, '\"'])

    n_cpus = os.environ.get("LSB_DJOB_NUMPROC", "1")
    os.environ["SINGULARITYENV_LSB_DJOB_NUMPROC"] = str(n_cpus)

    cmdToRun = (
            "/software/singularity/3.11.4/bin/singularity exec --cleanenv "
            f"-B {str(Path.cwd())}:/data "
            f"{str(Path.cwd())}/Algorithms/ARBORETO/arboreto_1.0.sif "
            f'/bin/sh -c \"'
            f"time -v -o /data/{outDir}time.txt "
            f"python /data/Algorithms/ARBORETO/runArboreto.py "
            f"--algo=GENIE3 "
            f"--inFile={inputPath} "
            f"--outFile={outPath}"
            "\""
        )

    print(cmdToRun)
    os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from GENIE3.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/GENIE3/"

        
    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = '\t', header = 0)
    
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
    
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for idx, row in OutDF.iterrows():
        outFile.write('\t'.join([row['TF'],row['target'],str(row['importance'])])+'\n')
    outFile.close()
    
