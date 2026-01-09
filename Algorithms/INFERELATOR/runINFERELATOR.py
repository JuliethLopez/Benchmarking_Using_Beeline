# Run pcor using spearman's correlation as mentioned in the PNI paper
#pcorResults=  pcor(x= t(as.matrix(inputExpr)), method = "spearman")

# Write output to a file
#DF = data.frame(Gene1 = geneNames[c(row(pcorResults$estimate))], Gene2 = geneNames[c(col(pcorResults$estimate))]
#                , corVal = c(pcorResults$estimate), pValue =  c(pcorResults$p.value))
#outDF <- DF[order(DF$corVal, decreasing=TRUE), ]
#write.table(outDF, outFile, sep = "\t", quote = FALSE, row.names = FALSE)

# Donde se define outDF como argumento de ppcor??

##########################################################################################
# Link to paper: https://doi.org/10.1093/bioinformatics/btac117
# Tutorial: https://inferelator.readthedocs.io/en/latest/

# Load modules
import os
import argparse
import pandas as pd #to load read GSD dataset
import numpy as np
import sys
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
    parser.add_argument('--ref_net_file', #NETWORK_PRIOR
        help='Path to refNetwork file')
    parser.add_argument('--gene_metadata_file_name',
        help='File with gene metadata: example SystematicName, GeneName, etc.')
    parser.add_argument('--metadata_columns', default=[],
        help='Metadata columns to extract from expression matrix')
    parser.add_argument('--ref_net_tf_names', #TF_NAMES
        help='Path to the names of TF from refNetwork')
    parser.add_argument('--priori_net_file', #PRIORI_NET
        help='Path to priors file, if available. If not, set to None or refNetwork file')
    parser.add_argument('--out_file',
        help='Path to output file')
    parser.add_argument('--num_bootstraps', type=int, default=2,
        help='Number of bootstraps for cross-validation')
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
    n_cores_local = 8
    local_engine = True

    # Multiprocessing uses the pathos implementation of multiprocessing (with dill instead of cPickle)
    # This is suited for a single computer, but will likely be too slow for the example here
    if __name__ == '__main__' and local_engine:
        MPControl.set_multiprocess_engine("multiprocessing")
        MPControl.client.processes = n_cores_local
        MPControl.connect()

    wkf = inferelator_workflow(regression="bbsr", workflow="single-cell")

    wkf.set_file_paths(input_dir=opts.input_dir,
                       output_dir=opts.output_dir,
                       expression_matrix_file=opts.expression_file,
                       gene_metadata_file=opts.gene_metadata_file_name,
                       gold_standard_file=opts.ref_net_file,
                       priors_file=opts.priori_net_file,
                       tf_names_file=opts.ref_net_tf_names) # in this case is possible to use the reference list due to is used in the beggining to create sintetic network
    # Set metadata columns
    #wkf.set_file_properties(extract_metadata_from_expression_matrix=True,
    #                        expression_matrix_metadata=METADATA_COLUMNS,
    #                        expression_matrix_columns_are_genes=True,
    #                        gene_list_index="SystematicName")
    wkf.set_crossvalidation_parameters(split_gold_standard_for_crossvalidation=False, #here we use all gold standard for evaluation
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
    worker.set_tfa(tfa_driver=False) #not to use prior information
    worker.append_to_path('output_dir', 'not_imputation_not_prior')

    cross_validation(worker).run()

    folder_results = f"random_seed_{random_seed_selection}"

    network_output_name = f'{opts.output_dir}/not_imputation_not_prior/{folder_results}/network.tsv.gz'
    network_output = pd.read_csv(network_output_name, sep='\t', index_col=None)

    # Select only required columns
    output = network_output[["regulator", "target", "model_coefficient"]]

    # Write final output
    output.to_csv(opts.out_file, sep='\t', index=False)


if __name__ == "__main__":
    main(sys.argv)