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
    if not RunnerObj.inputDir.joinpath("INFERELATOR").exists():
        print("Input folder for INFERELATOR does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("INFERELATOR").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("INFERELATOR/ExpressionData.tsv").exists():
        # input data
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)

        # Write .csv file
        ExpressionData.T.to_csv(RunnerObj.inputDir.joinpath("INFERELATOR/ExpressionData.tsv"),
                             sep = '\t', header  = True, index = True)
        
    if not RunnerObj.inputDir.joinpath("INFERELATOR/refNetwork.tsv").exists():
        refNetworkData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges),
                                     header = 0, index_col = None)
        
        # Convert simbols to numbers
        refNetworkData['Type_binary'] = [1 if i == '+' else -1 if i == '-' else 0 for i in refNetworkData.Type]
        # file to long format
        short_df = refNetworkData.pivot_table(index='Gene2', columns='Gene1', values='Type_binary').reset_index() #TF,gene,type of conection
        short_df = short_df.fillna(0)
        short_df = short_df.rename_axis(None, axis=1)
        # pass from float to int
        short_df.iloc[:, 1:] = short_df.iloc[:, 1:].fillna(0).astype(int)
        short_df = short_df.rename(columns={short_df.columns[0]: "Gene1/Gene2"})
        
        # Write reference network data in INFERELATOR folder
        short_df.to_csv(RunnerObj.inputDir.joinpath("INFERELATOR/refNetwork.tsv"),
                        sep="\t", index=False, header=True)
    
    if not RunnerObj.inputDir.joinpath("INFERELATOR/refNetwork_TF.tsv").exists():
        refNetworkData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges),
                                     header = 0, index_col = None)
        # Select transcription factors
        refNetworkData_TF = refNetworkData.Gene1.unique()
        refNetworkData_TF = pd.DataFrame(refNetworkData_TF, columns=['Gene1'])

        # Write TF network data in INFERELATOR folder
        refNetworkData_TF.to_csv(RunnerObj.inputDir.joinpath("INFERELATOR/refNetwork_TF.tsv"),
                             sep = '\t', index=False, header=False)
    
    
def run(RunnerObj):
    '''
    Function to run INFERELATOR algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    
    #Input directory
    inputDir= "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1]

    # Get path for ExpressionData.csv generated in INFERELATOR folder for certain type of network in inputs
    #tambien llamado en otros cuadernos inputPath
    expressionDataPath = f"/{inputDir}/INFERELATOR/ExpressionData.tsv"

    # Get path for refNetwork.tsv generated in INFERELATOR folder for certain type of network in inputs
    refNetworkPath = f"/{inputDir}/INFERELATOR/refNetwork.tsv"
    
    # Get path for refNetwork_TF.tsv generated in INFERELATOR folder for certain type of network in inputs
    refNetworkTFPath = f"/{inputDir}/INFERELATOR/refNetwork_TF.tsv"

    # make output dirs if they do not exist:
    outDir = "data/outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/INFERELATOR/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath = outDir + 'figure_5a_shuffled/random_seed_42/network.tsv.gz'
    
    # Command to run INFERELATOR
    #cmdToRun = ' '.join([
    #    'docker run --rm -v',
    #    str(Path.cwd()) + ':/data',
    #    'inferelator:base',
    #    '/bin/sh -c \"time -v -o',
    #    'data/' + str(outDir) + 'time.txt',
    #    'conda run -n inferenv python runINFERELATOR.py',
    #    '--input_dir=' + inputDir,
    #    '--expression_file=' + expressionDataPath,
    #    '--ref_net_file=' + refNetworkPath,
    #    '--ref_net_tf_names=' + refNetworkTFPath,
    #    '--priori_net_file=' + refNetworkPath,
    #    '--output_dir=' + outDir,
    #    '--out_file=' + outPath, '\"'])
    #    #'> data/' + str(outDir) + 'console_output.txt 2>&1\"'])

    cmdToRun = (
        "/software/singularity/3.11.4/bin/singularity exec --cleanenv "
        f"-B {str(Path.cwd())}:/data "
        f"{str(Path.cwd())}/Algorithms/INFERELATOR/inferelator_1.0.sif "
        f'/bin/sh -c \"'
        f"time -v -o /{outDir}time.txt "
        f"conda run -n inferenv python "
        f"/data/Algorithms/INFERELATOR/runINFERELATOR.py "
        f"--input_dir={inputDir} "
        f"--expression_file={expressionDataPath} "
        f"--ref_net_file={refNetworkPath} "
        f"--ref_net_tf_names={refNetworkTFPath} "
        f"--priori_net_file={refNetworkPath} "
        f"--output_dir={outDir}"
        "\""
    )
                        
    print(cmdToRun)
    os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from INFERELATOR.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    
    # Quit if output directory does not exist
        # Build base output path
    input_subdir = str(RunnerObj.inputDir).split("inputs/")[1]
    outDir = Path("outputs") / input_subdir / "INFERELATOR"
    print("Output directory for INFERELATOR:", outDir)

    # Create full nested folder structure
    finalDir = outDir / "figure_5a_shuffled" / "random_seed_42"
    if not finalDir.exists():
        print("Creating full folder path for INFERELATOR output...")
        finalDir.mkdir(parents=True, exist_ok=False)  # creates all intermediate dirs

    print("Final output directory for INFERELATOR:", finalDir)

    out_file = finalDir / 'network.tsv.gz'

    if not out_file.exists():
        print(f"{out_file} does not exist, skipping...")
        return
    
    else:
        print("Output file found, proceeding in formating file...")
        # Read output
        OutDF = pd.read_csv(out_file, sep = '\t', index_col=None) #to read network file


        #sort if necessary
        #OutDF.sort_values(by="EdgeWeight", ascending=False, inplace=True)

        # Formats the outFile into a ranked edgelist tab separated file
        outFile = OutDF[["regulator", "target", "model_coefficient"]]
        outFile.columns = ['Gene1','Gene2','EdgeWeight']