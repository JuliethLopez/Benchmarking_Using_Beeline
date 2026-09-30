# libraries
library(tidyr)
library(dplyr)

# read data
data = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/BEELINE-data/inputs/scRNA-Seq/hHep/ExpressionData.csv", row.names = 1)
head(data)
dim(data) # 11515   425

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
#"The total of TF in whole data is: 1115 (9.68)"

tf_in_selected = gene_selected[gene_selected %in% tf_all$TF]
length(tf_in_selected) # 159
print(paste0("The total of TF in data selected is: ",length(tf_in_selected)," (",round(length(tf_in_selected)/length(gene_selected)*100, 2),")"))
#"The total of TF in data selected is: 159 (7.95)"
write.table(tf_in_selected, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/list_TF.csv", 
          col.names=FALSE, row.names=FALSE, quote=FALSE)

# filter data
data_filtered = data[row.names(data) %in% gene_selected,]
dim(data_filtered) # 2000 425
head(data_filtered)
write.csv(data_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/ExpressionData.csv", 
          col.names=TRUE, row.names=TRUE, quote=FALSE)

### read old reference

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
s_reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/Networks/human/STRING-network.csv", header=TRUE)
dim(s_reference) # 198285      2
head(s_reference)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
s_reference = s_reference %>% mutate(Type = "+")
s_reference_filtered = s_reference[(s_reference$Gene1 %in% gene_selected) & (s_reference$Gene2 %in% gene_selected),]
dim(s_reference_filtered) # 5488 3

# save results
write.csv(s_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/STRING-network.csv", row.names = FALSE, quote = FALSE)

# reference 4
#string_genes <- read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/STRING_reference_2026.csv", header = TRUE)
string_reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/STRING_reference_2026.csv", header=TRUE)
dim(string_reference) # 13715404        3
c(1:nrow(string_reference))[is.na(string_reference[,1])]
string_reference[11417595,]


string_reference[,1]==string_reference[,2]
head(string_reference)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
string_reference = string_reference %>% mutate(Type = "+")
string_reference_filtered = string_reference[(string_reference$Gene1 %in% gene_selected) & (string_reference$Gene2 %in% gene_selected),]
dim(string_reference_filtered) # 5488 3

# save results
write.csv(string_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/STRING-network.csv", row.names = FALSE, quote = FALSE)

# reference 5
regnetwork_reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/Networks/human/STRING-network.csv", header=TRUE)
dim(regnetwork_reference) # 198285      2
head(regnetwork_reference)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
regnetwork_reference = regnetwork_reference %>% mutate(Type = "+")
regnetwork_reference_filtered = regnetwork_reference[(regnetwork_reference$Gene1 %in% gene_selected) & (regnetwork_reference$Gene2 %in% gene_selected),]
dim(regnetwork_reference_filtered) # 5488 3

# save results
write.csv(regnetwork_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/STRING-network.csv", row.names = FALSE, quote = FALSE)

# reference 6
htftarget_reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/Networks/human/STRING-network.csv", header=TRUE)
dim(htftarget_reference) # 198285      2
head(htftarget_reference)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
htftarget_reference = htftarget_reference %>% mutate(Type = "+")
htftarget_reference_filtered = htftarget_reference[(htftarget_reference$Gene1 %in% gene_selected) & (htftarget_reference$Gene2 %in% gene_selected),]
dim(htftarget_reference_filtered) # 5488 3

# save results
write.csv(htftarget_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/STRING-network.csv", row.names = FALSE, quote = FALSE)


### read new reference

# reference 1
reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/Networks/human/HepG2-ChIP-seq-network.csv", header=TRUE)
dim(reference) # 342862      2
head(reference)

#Remove repeated, auto loops
reference <- reference %>% filter(Gene1 != Gene2)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference$Type = "+"

reference_filtered = reference[(reference$Gene1 %in% gene_selected) & (reference$Gene2 %in% gene_selected),]
dim(reference_filtered) # 7004 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/hHep_Chip_seq/HepG2-ChIP-seq-network.csv", row.names = FALSE, quote = FALSE)

# reference 2
reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/hTFtargets_TF-Target-liver.csv", header=TRUE)
dim(reference) # 157768      3
head(reference)

#Remove repeated, auto loops
reference <- reference %>% filter(Gene1 != Gene2)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference$Type = "+"

reference_filtered = reference[(reference$Gene1 %in% gene_selected) & (reference$Gene2 %in% gene_selected),]
dim(reference_filtered) # 2432 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/hHep_hTFtargets/hTFtargets_TF-Target-liver.csv", row.names = FALSE, quote = FALSE)

# reference 3
reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/RegNetwork_human_TF_Target.csv", header=TRUE)
dim(reference) # 1112200      3
head(reference)

#Remove repeated, auto loops
reference <- reference %>% filter(Gene1 != Gene2)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference$Type = "+"

reference_filtered = reference[(reference$Gene1 %in% gene_selected) & (reference$Gene2 %in% gene_selected),]
dim(reference_filtered) # 10328 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/hHep_RegNetwork/RegNetwork_human_TF_Target.csv", row.names = FALSE, quote = FALSE)


# reference 4
reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/STRING_reference_2026.csv", header=TRUE)
dim(reference) # 13715404      3
head(reference)

#Remove repeated, auto loops
reference <- reference %>% filter(Gene1 != Gene2)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference$Type = "+"

reference_filtered = reference[(reference$Gene1 %in% gene_selected) & (reference$Gene2 %in% gene_selected),]
dim(reference_filtered) # 281688 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/hHep_STRING_2026/STRING_reference_2026.csv", row.names = FALSE, quote = FALSE)

# reference 5
reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/Networks/human/Non-specific-ChIP-seq-network.csv", header=TRUE)
dim(reference) # 386293      3
head(reference)

#Remove repeated, auto loops
reference <- reference %>% filter(Gene1 != Gene2)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference$Type = "+"

reference_filtered = reference[(reference$Gene1 %in% gene_selected) & (reference$Gene2 %in% gene_selected),]
dim(reference_filtered) # 1612 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/hHep/hHep_ns_Chip_seq/Non-specific-ChIP-seq-network.csv", row.names = FALSE, quote = FALSE)

