import os
import shutil
import pandas as pd

methods = ['GSD']
n_cells = ['200', '1000', '2000']
reps = range(1, 101)

base_dir = "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/Curated"

for method in methods:
    for n_cell in n_cells:
        for rep in reps:
            try:
                base_path = os.path.join(
                    base_dir,
                    method,
                    f"{method}-{n_cell}-{rep}",
                    "INFERELATOR",
                    "not_imputation_not_prior"
                )

                indir = os.path.join(base_path, "crossvalidation_performance.tsv")

                cv_results = pd.read_csv(indir, sep='\t')
                cv_results = cv_results.sort_values(by=["MCC", "F1"], ascending=[False, False])

                random_seed_selection = int(cv_results.iloc[0]["random_seed"])
                print(f"{method}-{n_cell}-{rep} -> Best random seed: {random_seed_selection}")

                best_seed_folder = os.path.join(
                    base_path,
                    f"random_seed_{random_seed_selection}"
                )

                # copy files
                if os.path.exists(best_seed_folder):
                    for file_name in os.listdir(best_seed_folder):
                        src = os.path.join(best_seed_folder, file_name)
                        dst = os.path.join(base_path, file_name)

                        if os.path.isfile(src):
                            shutil.copy2(src, dst)
                else:
                    print(f"Warning: {best_seed_folder} does not exist")

                # delete all random_seed folders
                for item in os.listdir(base_path):
                    item_path = os.path.join(base_path, item)

                    if item.startswith("random_seed") and os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                        print(f"Deleted: {item_path}")

            except Exception as e:
                print(f"Error in {method}-{n_cell}-{rep}: {e}")


## Filter WGCNA
setwd("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline")
datasets <- c("GSD", "HSC", "mCAD", "VSC")
cell_number <- c(2000, 1000, 200)

for (name in 1:length(datasets)) {
  for (ncell in 1:length(cell_number)) {
    for (rep in 1:100) {
      current_cells <- cell_number[ncell]
      
      path_input <- paste0(
        "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/Curated/", 
        datasets[name], "/", datasets[name], "-", current_cells, "-", rep, "/hdWGCNA/outFile.txt"
      )

      path_save <- paste0(
        "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/Curated/", 
        datasets[name], "/", datasets[name], "-", current_cells, "-", rep, "/hdWGCNA/rankedEdges.csv"
      )
      
      if (file.exists(path_input)) {
        #print(paste("Processing:", path))
        
        data <- read.delim(path_input, header = TRUE, sep = "\t")
        
        data$EdgeWeight[data$EdgeWeight < 0.2] <- 0
        
        write.table(data, path_save, sep = "\t", row.names = FALSE, quote = FALSE)
        
      } else {
        warning(paste("File missing, skipping:", path_input))
      }
    }
  }
}
