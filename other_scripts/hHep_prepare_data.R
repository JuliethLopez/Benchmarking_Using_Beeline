# libraries
library(tidyr)
library(dplyr)

# read data
data = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/BEELINE-data/inputs/scRNA-Seq/hHep/ExpressionData.csv", row.names = 1)
head(data)
dim(data) # 7465 5520

# read genes variability
genes_variability <- read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/GeneOrdering.csv",
          header=TRUE, row.names=1)
summary(genes_variability)
head(genes_variability)

# select HVG, filter by p value and order by variance
genes_variability_filt = genes_variability[genes_variability["VGAMpValue"] < 0.01, ]
dim(genes_variability_filt) # 7210 2
genes_variability_filt = genes_variability_filt %>% arrange(desc(Variance))
gene_selected = rownames(genes_variability_filt)[0:2000] # 100 more variable genes

# list of TF, ideally around 2% (1639/78000)
tf_all = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/human-tfs.csv")
head(tf_all)
length(tf_all)

# TF in data
tf_in_data = rownames(data)[rownames(data) %in% tf_all$TF]
print(paste0("The total of TF in whole data is: ",length(tf_in_data)," (",round(length(tf_in_data)/nrow(data)*100, 2),")"))

tf_in_selected = gene_selected[gene_selected %in% tf_all$TF]
length(tf_in_selected) # 159
print(paste0("The total of TF in data selected is: ",length(tf_in_selected)," (",round(length(tf_in_selected)/length(gene_selected)*100, 2),")"))
write.table(tf_in_selected, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/list_TF.csv", 
          col.names=FALSE, row.names=FALSE, quote=FALSE)

# filter data
data_filtered = data[row.names(data) %in% gene_selected,]
dim(data_filtered) # 2000 425
head(data_filtered)
write.csv(data_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/ExpressionData.csv", 
          col.names=TRUE, row.names=TRUE, quote=FALSE)

# read reference
# reference 1
chip_reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/Networks/human/HepG2-ChIP-seq-network.csv", header=TRUE)
dim(chip_reference) # 342863      2
head(chip_reference)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
chip_reference = chip_reference %>% mutate(Type = "+")
chip_reference_filtered = chip_reference[(chip_reference$Gene1 %in% gene_selected) & (chip_reference$Gene2 %in% gene_selected),]
dim(chip_reference_filtered) # 20990 3

# save results
write.csv(chip_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/HepG2-ChIP-seq-network.csv", row.names = FALSE, quote = FALSE)

# reference 2
ns_chip_reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/Networks/human/Non-specific-ChIP-seq-network.csv", header=TRUE)
dim(ns_chip_reference) # 386293      2
head(ns_chip_reference)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
ns_chip_reference = ns_chip_reference %>% mutate(Type = "+")
ns_chip_reference_filtered = ns_chip_reference[(ns_chip_reference$Gene1 %in% gene_selected) & (ns_chip_reference$Gene2 %in% gene_selected),]
dim(ns_chip_reference_filtered) # 1648 3

# save results
write.csv(ns_chip_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/Non-specific-ChIP-seq-network.csv", row.names = FALSE, quote = FALSE)

# reference 3
s_chip_reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/Networks/human/STRING-network.csv", header=TRUE)
dim(s_chip_reference) # 198285      2
head(s_chip_reference)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
s_chip_reference = s_chip_reference %>% mutate(Type = "+")
s_chip_reference_filtered = s_chip_reference[(s_chip_reference$Gene1 %in% gene_selected) & (s_chip_reference$Gene2 %in% gene_selected),]
dim(s_chip_reference_filtered) # 5488 3

# save results
write.csv(s_chip_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/STRING-network.csv", row.names = FALSE, quote = FALSE)
