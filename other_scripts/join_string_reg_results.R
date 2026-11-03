# combine string and regulatory results

dataset <- "ESC" #hHep

# read string results
whole_results <- read.csv(paste0("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/experimental/", dataset, "/", dataset, "-EPr-directed.csv"), header = TRUE)
order_whole_results <- whole_results$X
# read regulatory results
regulatory_results <- read.csv(paste0("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/experimental_reg/", dataset, "/", dataset, "-EPr-directed.csv"), header = TRUE)
# combine results and keep order of whole_results
combined_results <- merge(whole_results, regulatory_results, by = "X", all = TRUE)
combined_results <- combined_results[order(match(combined_results$X, order_whole_results)), ]
# write results
write.csv(combined_results, paste0("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/experimental/", dataset, "/", dataset, "-EPr-directed-combined.csv"), row.names = FALSE, quote = FALSE)

# combine for EPr-undirected
# read string results
whole_results <- read.csv(paste0("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/experimental/", dataset, "/", dataset, "-EPr-undirected.csv"), header = TRUE)
# read regulatory results
regulatory_results <- read.csv(paste0("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/experimental_reg/", dataset, "/", dataset, "-EPr-undirected.csv"), header = TRUE)
# combine results and keep order of whole_results
combined_results <- merge(whole_results, regulatory_results, by = "X", all = TRUE)
combined_results <- combined_results[order(match(combined_results$X, order_whole_results)), ]
# write results
write.csv(combined_results, paste0("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/experimental/", dataset, "/", dataset, "-EPr-undirected-combined.csv"), row.names = FALSE, quote = FALSE)
