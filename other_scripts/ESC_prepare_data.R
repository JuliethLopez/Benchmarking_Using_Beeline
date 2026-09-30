# libraries
library(tidyr)
library(dplyr)

# read data
data = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/expression_data/normalized/han_GSE107552.csv.gz", row.names = 1)
head(data)
dim(data) # 7465 5520

# calculate genes variability
# genes_variability = apply(data, 1, var)
# names(genes_variability) = rownames(data)
# summary(genes_variability)
# genes_variability = genes_variability[order(genes_variability, decreasing=TRUE)]
# head(genes_variability)
# write.table(genes_variability, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/list_genes_variability.csv", 
#           col.names=FALSE, row.names=TRUE, quote=FALSE)

genes_variability <- read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/list_genes_variability.csv",
          header=FALSE, row.names=1)
head(genes_variability)
tail(genes_variability)

# select HVG
gene_selected = rownames(genes_variability)[0:2000] # 100 more variable genes

# list of TF, ideally around 2% (1639/78000)
tf_all = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/human-tfs.csv")
head(tf_all)
dim(tf_all) #1563 2

# TF in data
tf_in_data = rownames(data)[rownames(data) %in% tf_all$TF]
print(paste0("The total of TF in whole data is: ",length(tf_in_data)," (",round(length(tf_in_data)/nrow(data)*100, 2),")"))

tf_in_selected = gene_selected[gene_selected %in% tf_all$TF]
length(tf_in_selected) # 247
print(paste0("The total of TF in data selected is: ",length(tf_in_selected)," (",round(length(tf_in_selected)/length(gene_selected)*100, 2),")"))
#write.table(tf_in_selected, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/list_TF.csv", colnames=FALSE, row.names=FALSE, quote=FALSE)

# filter data
data_filtered = data[row.names(data) %in% gene_selected,]
dim(data_filtered) # 2000 5520
head(data_filtered)
#write.csv(data_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/han_GSE107552.csv", col.names=TRUE, row.names=TRUE, quote=FALSE)

# read data filtered and cehck NA or ceros
data_filtered = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/han_GSE107552.csv", row.names=1)
head(data_filtered)
dim(data_filtered) # 2000 5520
any(is.na(data_filtered)) # There are NOT NA values in the data

### read old reference
# ref 1
reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_chipunion.txt", sep="\t", header=FALSE)
colnames(reference) = c("Gene1", "Gene2")
dim(reference) # 780470      2
head(reference)

# separate rows with multiple genes in Gene2
reference_long <- reference %>%
    separate_rows(Gene2, sep = ",")
dim(reference_long) # 780477      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference_long$Type = "+"

reference_filtered = reference_long[(reference_long$Gene1 %in% gene_selected) & (reference_long$Gene2 %in% gene_selected),]
dim(reference_filtered) # 20990 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_chipunion.csv", row.names = FALSE, quote = FALSE)


# ref 2
reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_chipunion_KDUnion_intersect.txt", sep="\t", header=FALSE)
colnames(reference) = c("Gene1", "Gene2")
dim(reference) # 780470      2
head(reference)

# separate rows with multiple genes in Gene2
reference_long <- reference %>%
separate_rows(Gene2, sep = ",")
dim(reference_long) # 780477      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference_long$Type = "+"

reference_filtered = reference_long[(reference_long$Gene1 %in% gene_selected) & (reference_long$Gene2 %in% gene_selected),]
dim(reference_filtered) # 846 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_chipunion_KDUnion_intersect.csv", row.names = FALSE, quote = FALSE)

# reference 3
reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_encodechip.txt", sep="\t", header=FALSE)
colnames(reference) = c("Gene1", "Gene2")
dim(reference) # 780470      2
head(reference)

# separate rows with multiple genes in Gene2
reference_long <- reference %>%
separate_rows(Gene2, sep = ",")
dim(reference_long) # 780477      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference_long$Type = "+"

reference_filtered = reference_long[(reference_long$Gene1 %in% gene_selected) & (reference_long$Gene2 %in% gene_selected),]
dim(reference_filtered) # 846 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_encodechip.csv", row.names = FALSE, quote = FALSE)

# reference 4
reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_escapechip.txt", sep="\t", header=FALSE)
colnames(reference) = c("Gene1", "Gene2")
dim(reference) # 780470      2
head(reference)

# separate rows with multiple genes in Gene2
reference_long <- reference %>%
separate_rows(Gene2, sep = ",")
dim(reference_long) # 107262      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference_long$Type = "+"

reference_filtered = reference_long[(reference_long$Gene1 %in% gene_selected) & (reference_long$Gene2 %in% gene_selected),]
dim(reference_filtered) # 4833 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_escapechip.csv", row.names = FALSE, quote = FALSE)

# reference 5
reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_KDUnion.txt", sep="\t", header=FALSE)
colnames(reference) = c("Gene1", "Gene2")
dim(reference) # 780470      2
head(reference)

# separate rows with multiple genes in Gene2
reference_long <- reference %>%
separate_rows(Gene2, sep = ",")
dim(reference_long) # 201805      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference_long$Type = "+"

reference_filtered = reference_long[(reference_long$Gene1 %in% gene_selected) & (reference_long$Gene2 %in% gene_selected),]
dim(reference_filtered) # 4521 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_KDUnion.csv", row.names = FALSE, quote = FALSE)

# reference 6
reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_lit.txt", sep="\t", header=FALSE)
colnames(reference) = c("Gene1", "Gene2")
dim(reference) # 780470      2
head(reference)

# separate rows with multiple genes in Gene2
reference_long <- reference %>%
separate_rows(Gene2, sep = ",")
dim(reference_long) # 330      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference_long$Type = "+"

reference_filtered = reference_long[(reference_long$Gene1 %in% gene_selected) & (reference_long$Gene2 %in% gene_selected),]
dim(reference_filtered) # 330 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_lit.csv", row.names = FALSE, quote = FALSE)

# reference 7
reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_logof.txt", sep="\t", header=FALSE)
colnames(reference) = c("Gene1", "Gene2")
dim(reference) # 780470      2
head(reference)

# separate rows with multiple genes in Gene2
reference_long <- reference %>%
separate_rows(Gene2, sep = ",")
dim(reference_long) # 330      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference_long$Type = "+"

reference_filtered = reference_long[(reference_long$Gene1 %in% gene_selected) & (reference_long$Gene2 %in% gene_selected),]
dim(reference_filtered) # 330 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_logof.csv", row.names = FALSE, quote = FALSE)

# reference 8
reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_nishiyama.txt", sep="\t", header=FALSE)
colnames(reference) = c("Gene1", "Gene2")
dim(reference) # 780470      2
head(reference)

# separate rows with multiple genes in Gene2
reference_long <- reference %>%
separate_rows(Gene2, sep = ",")
dim(reference_long) # 330      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference_long$Type = "+"

reference_filtered = reference_long[(reference_long$Gene1 %in% gene_selected) & (reference_long$Gene2 %in% gene_selected),]
dim(reference_filtered) # 330 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_nishiyama.csv", row.names = FALSE, quote = FALSE)

### read new reference
# reference 1
reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/Networks/human/hESC-ChIP-seq-network.csv", header=TRUE)
dim(reference) # 436563      2
head(reference)

#Remove repeated, auto loops
reference <- reference %>% filter(Gene1 != Gene2)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference$Type = "+"

reference_filtered = reference[(reference$Gene1 %in% gene_selected) & (reference$Gene2 %in% gene_selected),]
dim(reference_filtered) # 17923 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC-ChIP-seq-network.csv", row.names = FALSE, quote = FALSE)

# reference 2
reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/hTFtargets_TF-Target-embryo.csv", header=TRUE)
dim(reference) # 144125      3
head(reference)

#Remove repeated, auto loops
reference <- reference %>% filter(Gene1 != Gene2)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference$Type = "+"

reference_filtered = reference[(reference$Gene1 %in% gene_selected) & (reference$Gene2 %in% gene_selected),]
dim(reference_filtered) # 4913 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/ESC_hTFtargets/hTFtargets_TF-Target-embryo.csv", row.names = FALSE, quote = FALSE)

# reference 3
reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/RegNetwork_human_TF_Target.csv", header=TRUE)
dim(reference) # 1112200      3
head(reference)

#Remove repeated, auto loops
reference <- reference %>% filter(Gene1 != Gene2)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference$Type = "+"

reference_filtered = reference[(reference$Gene1 %in% gene_selected) & (reference$Gene2 %in% gene_selected),]
dim(reference_filtered) # 21679 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/ESC_RegNetwork/RegNetwork_human_TF_Target.csv", row.names = FALSE, quote = FALSE)


# reference 4
reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/STRING_reference_2026.csv", header=TRUE)
dim(reference) # 13715404      3
head(reference)

#Remove repeated, auto loops
reference <- reference %>% filter(Gene1 != Gene2)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference$Type = "+"

reference_filtered = reference[(reference$Gene1 %in% gene_selected) & (reference$Gene2 %in% gene_selected),]
dim(reference_filtered) # 379166 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/ESC_STRING_2026/STRING_reference_2026.csv", row.names = FALSE, quote = FALSE)

# reference 5
reference = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/Networks/human/Non-specific-ChIP-seq-network.csv", header=TRUE)
dim(reference) # 386293      3
head(reference)

#Remove repeated, auto loops
reference <- reference %>% filter(Gene1 != Gene2)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
reference$Type = "+"

reference_filtered = reference[(reference$Gene1 %in% gene_selected) & (reference$Gene2 %in% gene_selected),]
dim(reference_filtered) # 4343 3

# save results
write.csv(reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/ESC_ns_Chip_seq/Non-specific-ChIP-seq-network.csv", row.names = FALSE, quote = FALSE)
