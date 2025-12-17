import os
import pandas as pd
from pathlib import Path
import numpy as np

env = os.environ.copy() # added for debugging

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for PPCOR.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("PPCOR").exists():
        print("Input folder for PPCOR does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("PPCOR").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("PPCOR/ExpressionData.csv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        
        newExpressionData = ExpressionData.copy()
        
        # Write .csv file
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath("PPCOR/ExpressionData.csv"),
                             sep = ',', header  = True, index = True)
    
def run(RunnerObj):
    '''
    Function to run PPCOR algorithm
    '''
    inputPath = "/data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/PPCOR/ExpressionData.csv"
    

    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/PPCOR/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath = "/data/" +  str(outDir) + 'outFile.txt'
    print("The outPath is:...")
    print(outPath)

    #print(RunnerObj.container) #julieth

    #cmdToRun = ' '.join(['docker run --rm -v', 
    #                    str(Path.cwd())+':/data/ grnbeeline/ppcor:base /bin/sh -c \"time -v -o',
    #                    "data/" + str(outDir) + 'time.txt', 
    #                    'Rscript runPPCOR.R',
    #                    inputPath, outPath, '\"']) #for docker
    #cmdToRun = ' '.join(['pwd && cd Algorithms/PPCOR && ls -la']) #check current folder and content

    cmdToRun = ' '.join(['/software/singularity/3.11.4//bin/singularity exec --cleanenv -B',
                        str(Path.cwd()) + ':/data',
                        str(Path.cwd() / 'Algorithms/PPCOR/ppcor_1.0.sif'),
                        #'/bin/sh -c "echo OK > /data/'+ str(outDir) +'test2.txt"'
                        '/bin/sh -c \"time -v -o',
                        '/data/' + str(outDir) + 'time.txt',
                        'Rscript /data/Algorithms/PPCOR/runPPCOR.R',
                        inputPath, outPath, '\"']) # for singularity

    print("The command to run is:...")
    print(cmdToRun.strip())
    os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from PPCOR.
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/PPCOR/"
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
        
    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = '\t', header = 0)
    # edges with significant p-value
    part1 = OutDF.loc[OutDF['pValue'] <= float(RunnerObj.params['pVal'])]
    part1 = part1.assign(absCorVal = part1['corVal'].abs())
    # edges without significant p-value
    part2 = OutDF.loc[OutDF['pValue'] > float(RunnerObj.params['pVal'])]
    
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for idx, row in part1.sort_values('absCorVal', ascending = False).iterrows():
        outFile.write('\t'.join([row['Gene1'],row['Gene2'],str(row['corVal'])])+'\n')
    
    for idx, row in part2.iterrows():
        outFile.write('\t'.join([row['Gene1'],row['Gene2'],str(0)])+'\n')
    outFile.close()
    
