# Link to paper: https://doi.org/10.1093/bioinformatics/btac117
# Tutorial: https://inferelator.readthedocs.io/en/latest/

# Load modules
import os
import argparse
import pandas as pd #to load read GSD dataset
import numpy as np
import sys
from pathlib import Path
#print(sys.path)
sys.path.append('INFERELATOR') #to add a path to search for the requested module

from inferelator import inferelator_workflow, inferelator_verbose_level, MPControl, CrossValidationManager

# Set verbosity level to "Talky"
inferelator_verbose_level(1)

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Run scSGL algorithm.') 
    parser.add_argument('--input_dir', #DATA_DIR
        help='Path to input data directory containing expression matrix and metadata files')
    parser.add_argument('--output_dir', #OUTPUT_DIR
        help='Path to output directory where results will be saved')
    parser.add_argument('--expression_file', #EXPRESSION_FILE_NAME
        help='Path to ExpressionData file')
    parser.add_argument('--ref_net_file', #GOLD STANDARD, COULD BE USED THE NETWORK_PRIOR
        help='Path to refNetwork file')
    parser.add_argument('--gene_metadata_file_name',
        help='File with gene metadata: example SystematicName, GeneName, etc.')
    parser.add_argument('--metadata_columns', default=[],
        help='Metadata columns to extract from expression matrix')
    parser.add_argument('--ref_net_tf_names', #TF_NAMES
        help='Path to the names of TF from refNetwork')
    parser.add_argument('--priori_net_file', #PRIORI_NET IF KNOWN, IF NOT CREATE A DUMMI FILE WITH THE SAME FORMAT AS REFNETWORK BUT WITH ALL ZEROS, THEN THE ALGORITHM WILL NOT USE PRIOR INFORMATION
        help='Path to priors file, if available. If not, set to None or refNetwork file')
    parser.add_argument('--out_file',
        help='Path to output file')
    parser.add_argument('--num_bootstraps', type=int, default=2,
        help='Number of bootstraps for cross-validation')
    parser.add_argument('--cross_val',
        help='To do crossvalidation, True or False')
    parser.add_argument('--count_minimum', type=float, default=0.05,
        help='Minimum count threshold for filtering genes')
    parser.add_argument('--random_seed_min', type=int, default=42,
        help='Minimum random seed for grid search')
    parser.add_argument('--random_seed_max', type=int, default=43,
        help='Maximum random seed for grid search')
    parser.add_argument('--random_seed_selection', type=int, default=42,
        help='Random seed selection for results folder naming')

    return parser

def parse_arguments():
    parser = get_parser()
    opts = parser.parse_args()

    return opts

def set_up_workflow(args):

    opts = parse_arguments()
    
    # Start Multiprocessing Engine
    n_cores_local = int(os.environ.get("LSB_DJOB_NUMPROC", 1))
    print(f"Detected CPUs from LSF: {n_cores_local}")
    
    local_engine = True

    # Multiprocessing uses the pathos implementation of multiprocessing (with dill instead of cPickle)
    # This is suited for a single computer, but will likely be too slow for the example here
    if __name__ == '__main__' and local_engine:
        MPControl.set_multiprocess_engine("multiprocessing")  # local thread-based better for singularity container
        MPControl.client.processes = n_cores_local
        MPControl.connect()

    wkf = inferelator_workflow(regression="bbsr", workflow="single-cell")

    wkf.set_file_paths(input_dir=opts.input_dir,
                       output_dir=opts.output_dir,
                       expression_matrix_file=opts.expression_file,
                       gene_metadata_file=opts.gene_metadata_file_name,
                       gold_standard_file=None, #opts.ref_net_file,
                       priors_file=opts.priori_net_file,
                       tf_names_file=opts.ref_net_tf_names) # in this case is possible to use the reference list due to is used in the beggining to create sintetic network
    # Set metadata columns
    #wkf.set_file_properties(extract_metadata_from_expression_matrix=True,
    #                        expression_matrix_metadata=METADATA_COLUMNS,
    #                        expression_matrix_columns_are_genes=True,
    #                        gene_list_index="SystematicName")
    wkf.set_crossvalidation_parameters(split_gold_standard_for_crossvalidation=opts.cross_val, #here we use all gold standard for evaluation
                                       cv_split_ratio=0.5)
    wkf.set_run_parameters(num_bootstraps=opts.num_bootstraps)
    wkf.set_count_minimum(opts.count_minimum)
    wkf.add_preprocess_step("log2")
    return wkf


def cross_validation(wkf):
    opts = parse_arguments()
    
    cv_wrap = CrossValidationManager(wkf)
    cv_wrap.add_gridsearch_parameter('random_seed', list(range(opts.random_seed_min, opts.random_seed_max))) #this affects bbotstraping resampling and priorshuffling if done
    return cv_wrap

def main(args):
    import os
    print(f"Running in CWD: {os.getcwd()}")
    print("using args:")
    print(args)
    #print(os.system("ls -lah"))
    #print(os.system("ls /"))
    #print(os.system("find / -name \"ExpressionData.tsv\" 2>/dev/null"))
    opts = parse_arguments()
    print(opts)

    # create the output dir if the output directory does not exist
    if not os.path.exists(opts.output_dir):
        os.makedirs(opts.output_dir)

    # Create workflow and run without imputation
    # without prior information, whithout shuffle prior information
    # with cross validation evaluation using gold standard
    worker = set_up_workflow(args)
    #worker.set_shuffle_parameters(shuffle_prior_axis=0) # shuffle prior and alter biological meaning
    worker.set_network_data_flags(use_no_gold_standard=True) #not to use gold standard for evaluation, use all gold standard for evaluation
    worker.set_tfa(tfa_driver=False) #not use prior information, then not calculate transcription factor activity
    worker.append_to_path('output_dir', 'not_imputation_not_prior')

    # run
    if opts.cross_val == True:
        # run method
        cross_validation(worker).run()
    else:
        # run method
        worker.run()
    
    # select and save results
    if opts.cross_val == True:
        # select random seed
        if opts.random_seed_selection==0:
            print("Random seed is cero... Selecting the best one...")
            # Read crossvalidation performance results
            cv_results = pd.read_csv(f"{opts.output_dir}not_imputation_not_prior/crossvalidation_performance.tsv", sep='\t')
            cv_results = cv_results.sort_values(by=["MCC","F1"], ascending=[False,False])
            
            # select the best random seed
            random_seed_selection = int(cv_results.iloc[0]["random_seed"])
            print(f"Best random seed: {random_seed_selection}")
            folder_results = f"random_seed_{random_seed_selection}"
        else:
            print(f"Using random seed: {opts.random_seed_selection}")
            folder_results = f"random_seed_{opts.random_seed_selection}"

        # save the folder with selected random seed
        final_folder = f"{opts.output_dir}/not_imputation_not_prior/{folder_results}"
        Path(final_folder).mkdir(parents=True, exist_ok=True) #check if folder/path exist, if not creste them

    else:
        # define final folder
        final_folder = f"{opts.output_dir}not_imputation_not_prior"
    
    network_output_name = f'{final_folder}/network.tsv.gz'
    network_output = pd.read_csv(network_output_name, sep='\t', index_col=None)

    # Select only required columns
    output = network_output[["regulator", "target", "model_coefficient"]]

    # Write final output
    output.to_csv(opts.out_file, sep='\t', index=False)
    #output.to_csv(f'{final_folder}/network.tsv.gz', sep='\t', index=False)

if __name__ == "__main__":
    main(sys.argv)