# Borda order
import os
import pandas as pd
import numpy as np

# Define rute
all_datasets = ['GSD','HSC','mCAD','VSC'] # 

for dataset in all_datasets:
    input_folder = f"/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/Curated/{dataset}"
    os.getcwd()
    os.chdir(input_folder)

    for ncell in ['200','1000','2000']:
        for rep in range(1,101):
            sub_folder = f"{dataset}-{ncell}-{rep}"
            print(sub_folder)
            # Read borda file
            borda_df = pd.read_csv(f"{sub_folder}/Borda.csv")
            # Separate borda in different files
            borda_methods = ['BORDA','mBORDA','sBORDA','smBORDA']
            # Create folders and save the corresponding files
            for method in borda_methods:
                method_folder = f"{sub_folder}/{method}"
                print(method_folder)
                os.makedirs(method_folder,exist_ok=True)
                method_df = borda_df[['Gene1','Gene2',method]]
                method_df.columns = ['Gene1','Gene2','EdgeWeight']
                # this filter is applied because the network is full connected
                # take the 20% higher connections. As ranks are used, better connections have lower values
                #threshold = np.quantile(method_df['EdgeWeight'], 0.2, axis=1)
                #method_df['EdgeWeight'] = [0 if x > threshold else x for x in method_df['EdgeWeight']]
                # save
                method_df.to_csv(f'{method_folder}/rankedEdges.csv', sep='\t', index=False)