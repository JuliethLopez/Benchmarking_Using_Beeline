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
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        
        # read reference network data
        refNetworkData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges),
                                    header = 0, index_col = None)
    
        # Convert simbols to numbers
        #refNetworkData['Type_binary'] = [1 if i == '+' else -1 if i == '-' else 0 for i in refNetworkData.Type]
        
        # Get all unique TFs (The "Sources" / Columns)
        all_tfs = refNetworkData['Gene1'].unique()

        # Get all unique Genes from Expression data (The "Targets" / Rows)
        all_genes = ExpressionData.index

        # Create matrix with the correct dimensions and fill it
        prior_df = pd.DataFrame(0, index=all_genes, columns=all_tfs)
        for _, row in refNetworkData.iterrows():
            # row['Gene1'] is the TF (Column)
            # row['Gene2'] is the Target (Row)
            # Convert '+' to 1, '-' to -1
            val = 1 if row['Type'] == '+' else -1
            prior_df.at[row['Gene2'], row['Gene1']] = val

        prior_df.index.name = "Gene1/Gene2"
        prior_df = prior_df.reset_index()
        prior_df.iloc[:, 1:] = prior_df.iloc[:, 1:].fillna(0).astype(int)
        # prior_df = prior_df.rename(columns={prior_df.columns[0]: "Gene1/Gene2"})
        
        # # file to long format
        # short_df = refNetworkData.pivot_table(index='Gene2', columns='Gene1', values='Type_binary').reset_index() #TF,gene,type of conection
        # short_df = short_df.fillna(0)
        # short_df = short_df.rename_axis(None, axis=1)
        # # pass from float to int
        # short_df.iloc[:, 1:] = short_df.iloc[:, 1:].fillna(0).astype(int)
        # short_df = short_df.rename(columns={short_df.columns[0]: "Gene1/Gene2"})
        
        # Write reference network data in INFERELATOR folder
        prior_df.to_csv(RunnerObj.inputDir.joinpath("INFERELATOR/refNetwork.tsv"),
                        sep="\t", index=False, header=True)


    if not RunnerObj.inputDir.joinpath("INFERELATOR/priori.tsv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        refNetworkData_TF = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.params["ref_net_tf_list"]),
                                        header = None, index_col = None)
        
        if RunnerObj.params['use_priors'] == False:
            # Create an empty matrix (all zeros)
            # Rows = All Genes in experiment, Columns = All TFs
            tf_list = refNetworkData_TF.iloc[:, 0].tolist()

            # Use prior_df name to keep consistency with the saving line
            prior_df = pd.DataFrame(0, index=ExpressionData.index, columns=tf_list)
            prior_df.index.name = "Gene1/Gene2"
            prior_df = prior_df.reset_index()
        else:
            # read reference network data
            refNetworkData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges),
                                     header = 0, index_col = None)
        
            # Convert simbols to numbers
            #refNetworkData['Type_binary'] = [1 if i == '+' else -1 if i == '-' else 0 for i in refNetworkData.Type]
            
            # Get all unique TFs (The "Sources" / Columns)
            all_tfs = refNetworkData['Gene1'].unique()

            # Get all unique Genes from Expression data (The "Targets" / Rows)
            all_genes = ExpressionData.index

            # Create matrix with the correct dimensions and fill it
            prior_df = pd.DataFrame(0, index=all_genes, columns=all_tfs)
            for _, row in refNetworkData.iterrows():
                # row['Gene1'] is the TF (Column)
                # row['Gene2'] is the Target (Row)
                # Convert '+' to 1, '-' to -1
                val = 1 if row['Type'] == '+' else -1
                prior_df.at[row['Gene2'], row['Gene1']] = val

            prior_df.index.name = "Gene1/Gene2"
            prior_df = prior_df.reset_index()
            prior_df.iloc[:, 1:] = prior_df.iloc[:, 1:].fillna(0).astype(int)
            # prior_df = prior_df.rename(columns={prior_df.columns[0]: "Gene1/Gene2"})
            
            # # file to long format
            # short_df = refNetworkData.pivot_table(index='Gene2', columns='Gene1', values='Type_binary').reset_index() #TF,gene,type of conection
            # short_df = short_df.fillna(0)
            # short_df = short_df.rename_axis(None, axis=1)
            # # pass from float to int
            # short_df.iloc[:, 1:] = short_df.iloc[:, 1:].fillna(0).astype(int)
            # short_df = short_df.rename(columns={short_df.columns[0]: "Gene1/Gene2"})
            
        # Write reference network data in INFERELATOR folder
        prior_df.to_csv(RunnerObj.inputDir.joinpath("INFERELATOR/priori.tsv"),
                        sep="\t", index=False, header=True)
    
    if not RunnerObj.inputDir.joinpath("INFERELATOR/list_TF.tsv").exists():
        refNetworkData_TF = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.params["ref_net_tf_list"]),
                                        header = None, index_col = None)
        # Select transcription factors
        #refNetworkData_TF = refNetworkData.Gene1.unique()
        refNetworkData_TF = pd.DataFrame(refNetworkData_TF.iloc[:, 0])

        # Write TF network data in INFERELATOR folder
        refNetworkData_TF.to_csv(RunnerObj.inputDir.joinpath("INFERELATOR/list_TF.tsv"),
                            sep = '\t', index=False, header=False)

    
def run(RunnerObj):
    '''
    Function to run INFERELATOR algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    
    #Input directory
    #container_inputDir = f"/data/inputs/{input_subdir}"
    inputDir= "/data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1]

    # Get path for ExpressionData.csv generated in INFERELATOR folder for certain type of network in inputs
    #tambien llamado en otros cuadernos inputPath
    expressionDataPath = f"{inputDir}/INFERELATOR/ExpressionData.tsv"

    # Get path for refNetwork.tsv generated in INFERELATOR folder for certain type of network in inputs
    refNetworkPath = f"{inputDir}/INFERELATOR/refNetwork.tsv"

    # Get path for priori.tsv generated in INFERELATOR folder for certain type of network in inputs
    prioriPath = f"{inputDir}/INFERELATOR/priori.tsv"
    
    # Get path for list_TF.tsv generated in INFERELATOR folder for certain type of network in inputs
    refNetworkTFPath = f"{inputDir}/INFERELATOR/list_TF.tsv"

    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/INFERELATOR/"
    os.makedirs(outDir, exist_ok = True)
    
    #output_dir = "/" + outDir

    num_bootstraps = int(RunnerObj.params['num_bootstraps'])
    cross_val = str(RunnerObj.params['cross_val'])
    count_minimum = float(RunnerObj.params['count_minimum'])
    random_seed_min = int(RunnerObj.params['random_seed_min'])
    random_seed_selection = int(RunnerObj.params['random_seed_selection'])
    random_seed_max = int(RunnerObj.params['random_seed_max'])

    folder_results = f"random_seed_{0}" #{random_seed_selection}"

    #outPath =  outDir + "not_imputation_not_prior/" + folder_results + "/network.tsv.gz"

    # Create the nested directory structure for output
    nested_outDir = Path(outDir) / "not_imputation_not_prior" / folder_results
    nested_outDir.mkdir(parents=True, exist_ok=True)

    outPath = str(nested_outDir / "network.tsv.gz")


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
    #    '--priori_net_file=' + prioriPath,
    #    '--output_dir=' + outDir,
    #    '--out_file=' + outPath, '\"'])
    #    #'> data/' + str(outDir) + 'console_output.txt 2>&1\"'])

    n_cpus = os.environ.get("LSB_DJOB_NUMPROC", "1")
    os.environ["SINGULARITYENV_LSB_DJOB_NUMPROC"] = str(n_cpus)

    cmdToRun = (
        "/software/singularity/3.11.4/bin/singularity exec --cleanenv "
        f"-B {str(Path.cwd())}:/data "
        f"{str(Path.cwd())}/Algorithms/INFERELATOR/inferelator_1.0.sif "
        f'/bin/sh -c \"'
        f"time -v -o /data/{outDir}time.txt "
        f"conda run -n inferenv python "
        f"/data/Algorithms/INFERELATOR/runINFERELATOR.py "
        f"--input_dir={inputDir} "
        f"--expression_file={expressionDataPath} "
        f"--ref_net_file={refNetworkPath} "
        f"--ref_net_tf_names={refNetworkTFPath} "
        f"--priori_net_file={prioriPath} "
        f"--output_dir=/data/{outDir} "
        f"--out_file=/data/{outPath} "
        f"--num_bootstraps={num_bootstraps} "
        f"--cross_val={cross_val} "
        f"--count_minimum={count_minimum} "
        f"--random_seed_min={random_seed_min} "
        f"--random_seed_max={random_seed_max} "
        f"--random_seed_selection={random_seed_selection}"
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
    outDir = str(Path("outputs") / input_subdir / "INFERELATOR")
    print("Output directory for INFERELATOR:", outDir)

    # Create full nested folder structure
    #random_seed_selection = int(RunnerObj.params['random_seed_selection'])
    #folder_results = f"random_seed_{random_seed_selection}"

    #read the final network file from INFERELATOR
    final_network_file = pd.read_csv(f'{outDir}/not_imputation_not_prior/random_seed_0/network.tsv.gz', sep='\t', index_col=None)

    #finalDir = Path(outDir) / "not_imputation_not_prior" / folder_results
    finalDir = Path(outDir) / "not_imputation_not_prior/random_seed_0"
    if not finalDir.exists():
        print("Creating full folder path for INFERELATOR output...")
        finalDir.mkdir(parents=True, exist_ok=False)  # creates all intermediate dirs

    print("Final output directory for INFERELATOR:", finalDir)

    out_file = Path(finalDir) / "network.tsv.gz"

    print("Outfile path:", out_file)

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
        outFile.to_csv(f'{outDir}/rankedEdges.csv', sep='\t', index=False)
