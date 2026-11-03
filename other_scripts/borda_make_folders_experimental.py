# Borda order
import os
import pandas as pd
import numpy as np
import shutil

# Define rute
all_datasets = ['ESC', 'hHep']

for dataset in all_datasets:
    input_folder = f"/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/experimental"
    os.getcwd()
    os.chdir(input_folder)

    # Read borda file
    borda_df = pd.read_csv(f"{dataset}/{dataset}_Chip_seq/Borda.csv")

    # copy borda original file and folder results to other folders for evaluate reference networks
    reference_folders = ['hTFtargets','ns_Chip_seq','RegNetwork','STRING_2026','STRING_reg_2026']
    for reference_folder in reference_folders:
        #os.makedirs(f"{dataset}/{dataset}_{reference_folder}",exist_ok=True)
        shutil.copy(f"{dataset}/{dataset}_Chip_seq/Borda.csv", f"{dataset}/{dataset}_{reference_folder}/Borda.csv")
        
    # Separate borda in different files
    borda_methods = ['BORDA','mBORDA','sBORDA','smBORDA']
    # Create folders and save the corresponding files
    for method in borda_methods:
        print(method)
        method_folder = f"{dataset}/{dataset}_Chip_seq/{method}"
        print(method_folder)
        os.makedirs(method_folder,exist_ok=True)
        method_df = borda_df[['Gene1','Gene2',method]]
        method_df.columns = ['Gene1','Gene2','EdgeWeight']
        # save
        method_df.to_csv(f'{method_folder}/rankedEdges.csv', sep='\t', index=False)
        
        for reference_folder in reference_folders:
            os.makedirs(f"{dataset}/{dataset}_{reference_folder}/{method}",exist_ok=True)
            shutil.copy(f"{method_folder}/rankedEdges.csv", f"{dataset}/{dataset}_{reference_folder}/{method}/rankedEdges.csv")
