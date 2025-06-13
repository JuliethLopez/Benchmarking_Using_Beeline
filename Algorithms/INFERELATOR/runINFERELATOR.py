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
                       tf_names_file=opts.ref_net_tf_names)
    #wkf.set_file_properties(extract_metadata_from_expression_matrix=True,
    #                        expression_matrix_metadata=METADATA_COLUMNS,
    #                        expression_matrix_columns_are_genes=True,
    #                        gene_list_index="SystematicName")
    wkf.set_crossvalidation_parameters(split_gold_standard_for_crossvalidation=True,
                                       cv_split_ratio=0.5)
    wkf.set_run_parameters(num_bootstraps=2)
    wkf.set_count_minimum(0.05)
    wkf.add_preprocess_step("log2")
    return wkf


def set_up_fig5a(wkf):
    cv_wrap = CrossValidationManager(wkf)
    cv_wrap.add_gridsearch_parameter('random_seed', list(range(42, 43)))
    return cv_wrap

def main(args):
    #import os
    #print(f"Running in CWD: {os.getcwd()}")
    opts = parse_arguments()

    # create the output dir if the output directory does not exist
    if not os.path.exists(opts.output_dir):
        os.makedirs(opts.output_dir)

    # Figure 5A: Shuffled Priors
    worker = set_up_workflow(args)
    worker.set_shuffle_parameters(shuffle_prior_axis=0)
    worker.append_to_path('output_dir', 'figure_5a_shuffled')

    set_up_fig5a(worker).run()

    network_output_name = f'{opts.outdir}/figure_5a_shuffled/random_seed_42/network.tsv.gz'
    network_output = pd.read_csv(network_output_name, sep='\t', index_col=None)
    
    # select only Gene1/regulator Gene2/target and EdgeWeight/model_coefficient
    output = network_output[["regulator", "target", "model_coefficient"]]
    #Each edge is also associated with a weight, which is either positive or negative depending on the sign of the edge.
    output.to_csv(opts.out_file, sep='\t', index=False)
    
    G.to_csv(opts.out_file, index = False, sep = '\t')  #to write the output file


if __name__ == "__main__":
    main(sys.argv)