import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for INFERELATOR.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    if not RunnerObj.inputDir.joinpath("TF_hdWGCNA").exists():
        print("Input folder for TF_hdWGCNA does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("TF_hdWGCNA").mkdir(exist_ok = False)
    
    if not RunnerObj.inputDir.joinpath("TF_hdWGCNA/ExpressionData.csv").exists():
        # input data
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)

        # Write .csv file
        ExpressionData.T.to_csv(RunnerObj.inputDir.joinpath("TF_hdWGCNA/ExpressionData.csv"),
                             sep = '\t', header  = True, index = True)
        

def run(RunnerObj):
    '''
    Function to run INFERELATOR algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/TF_hdWGCNA/ExpressionData.csv"
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/TF_hdWGCNA/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath =  "data/" + str(outDir) + 'outFile.txt'
    cmdToRun = ' '.join([#'docker run --rm -v',
                         #str(Path.cwd()) + ':/data',
                         #'tf_hdwgcna:base',
                         #'/bin/sh -c \"time -v -o',
                         'docker run --rm -v', str(Path.cwd())+':/data/ tf_hdwgcna:base /bin/sh -c \"time -v -o',
                         "data/" + str(outDir) + 'time.txt', 
                         'Rscript runTF_hdWGCNA.R',
                         inputPath, 
                         outPath, '\"'])
    print(cmdToRun)
    os.system(cmdToRun)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from INFERELATOR.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/TF_hdWGCNA/"
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
    else:
        print("Output file found, proceeding in formating file...")
    
    print("Output directory for TF_hdWGCNA:", outDir)

    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = '\t', header = 0)

    # Formats the outFile if necessary
    print("Formats the outFile if necessary...")
