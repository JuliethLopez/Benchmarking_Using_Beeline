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
length(tf_all)

# TF in data
tf_in_data = rownames(data)[rownames(data) %in% tf_all$TF]
print(paste0("The total of TF in whole data is: ",length(tf_in_data)," (",round(length(tf_in_data)/nrow(data)*100, 2),")"))

tf_in_selected = gene_selected[gene_selected %in% tf_all$TF]
length(tf_in_selected) # 247
print(paste0("The total of TF in data selected is: ",length(tf_in_selected)," (",round(length(tf_in_selected)/length(gene_selected)*100, 2),")"))
write.table(tf_in_selected, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/list_TF.csv", 
          col.names=FALSE, row.names=FALSE, quote=FALSE)

# filter data
data_filtered = data[row.names(data) %in% gene_selected,]
dim(data_filtered) # 2000 5520
head(data_filtered)
write.csv(data_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/han_GSE107552.csv", 
          col.names=TRUE, row.names=TRUE, quote=FALSE)

# read data filtered and cehck NA or ceros
data_filtered = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/han_GSE107552.csv", row.names=1)
head(data_filtered)
dim(data_filtered) # 2000 5520
any(is.na(data_filtered)) # There are NOT NA values in the data

# read reference
# ref 1
chip_reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_chipunion.txt", sep="\t", header=FALSE)
colnames(chip_reference) = c("Gene1", "Gene2")
dim(chip_reference) # 780470      2
head(chip_reference)

# separate rows with multiple genes in Gene2
chip_reference_clean <- chip_reference %>%
separate_rows(Gene2, sep = ",")
dim(chip_reference_clean) # 780477      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
chip_reference_clean$Type = "+"

chip_reference_filtered = chip_reference_clean[(chip_reference_clean$Gene1 %in% gene_selected) & (chip_reference_clean$Gene2 %in% gene_selected),]
dim(chip_reference_filtered) # 20990 3

# save results
write.csv(chip_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_chipunion.csv", row.names = FALSE, quote = FALSE)

# ref 2
chip_reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_chipunion_KDUnion_intersect.txt", sep="\t", header=FALSE)
colnames(chip_reference) = c("Gene1", "Gene2")
dim(chip_reference) # 780470      2
head(chip_reference)

# separate rows with multiple genes in Gene2
chip_reference_clean <- chip_reference %>%
separate_rows(Gene2, sep = ",")
dim(chip_reference_clean) # 780477      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
chip_reference_clean$Type = "+"

chip_reference_filtered = chip_reference_clean[(chip_reference_clean$Gene1 %in% gene_selected) & (chip_reference_clean$Gene2 %in% gene_selected),]
dim(chip_reference_filtered) # 846 3

# save results
write.csv(chip_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_chipunion_KDUnion_intersect.csv", row.names = FALSE, quote = FALSE)

# reference 3
chip_reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_encodechip.txt", sep="\t", header=FALSE)
colnames(chip_reference) = c("Gene1", "Gene2")
dim(chip_reference) # 780470      2
head(chip_reference)

# separate rows with multiple genes in Gene2
chip_reference_clean <- chip_reference %>%
separate_rows(Gene2, sep = ",")
dim(chip_reference_clean) # 780477      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
chip_reference_clean$Type = "+"

chip_reference_filtered = chip_reference_clean[(chip_reference_clean$Gene1 %in% gene_selected) & (chip_reference_clean$Gene2 %in% gene_selected),]
dim(chip_reference_filtered) # 846 3

# save results
write.csv(chip_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_encodechip.csv", row.names = FALSE, quote = FALSE)

# reference 4
chip_reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_escapechip.txt", sep="\t", header=FALSE)
colnames(chip_reference) = c("Gene1", "Gene2")
dim(chip_reference) # 780470      2
head(chip_reference)

# separate rows with multiple genes in Gene2
chip_reference_clean <- chip_reference %>%
separate_rows(Gene2, sep = ",")
dim(chip_reference_clean) # 107262      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
chip_reference_clean$Type = "+"

chip_reference_filtered = chip_reference_clean[(chip_reference_clean$Gene1 %in% gene_selected) & (chip_reference_clean$Gene2 %in% gene_selected),]
dim(chip_reference_filtered) # 4833 3

# save results
write.csv(chip_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_escapechip.csv", row.names = FALSE, quote = FALSE)

# reference 5
chip_reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_KDUnion.txt", sep="\t", header=FALSE)
colnames(chip_reference) = c("Gene1", "Gene2")
dim(chip_reference) # 780470      2
head(chip_reference)

# separate rows with multiple genes in Gene2
chip_reference_clean <- chip_reference %>%
separate_rows(Gene2, sep = ",")
dim(chip_reference_clean) # 201805      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
chip_reference_clean$Type = "+"

chip_reference_filtered = chip_reference_clean[(chip_reference_clean$Gene1 %in% gene_selected) & (chip_reference_clean$Gene2 %in% gene_selected),]
dim(chip_reference_filtered) # 4521 3

# save results
write.csv(chip_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_KDUnion.csv", row.names = FALSE, quote = FALSE)

# reference 6
chip_reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_lit.txt", sep="\t", header=FALSE)
colnames(chip_reference) = c("Gene1", "Gene2")
dim(chip_reference) # 780470      2
head(chip_reference)

# separate rows with multiple genes in Gene2
chip_reference_clean <- chip_reference %>%
separate_rows(Gene2, sep = ",")
dim(chip_reference_clean) # 330      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
chip_reference_clean$Type = "+"

chip_reference_filtered = chip_reference_clean[(chip_reference_clean$Gene1 %in% gene_selected) & (chip_reference_clean$Gene2 %in% gene_selected),]
dim(chip_reference_filtered) # 330 3

# save results
write.csv(chip_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_lit.csv", row.names = FALSE, quote = FALSE)

# reference 7
chip_reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_logof.txt", sep="\t", header=FALSE)
colnames(chip_reference) = c("Gene1", "Gene2")
dim(chip_reference) # 780470      2
head(chip_reference)

# separate rows with multiple genes in Gene2
chip_reference_clean <- chip_reference %>%
separate_rows(Gene2, sep = ",")
dim(chip_reference_clean) # 330      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
chip_reference_clean$Type = "+"

chip_reference_filtered = chip_reference_clean[(chip_reference_clean$Gene1 %in% gene_selected) & (chip_reference_clean$Gene2 %in% gene_selected),]
dim(chip_reference_filtered) # 330 3

# save results
write.csv(chip_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_logof.csv", row.names = FALSE, quote = FALSE)

# reference 8
chip_reference = read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090/gold_standards/hESC/hESC_nishiyama.txt", sep="\t", header=FALSE)
colnames(chip_reference) = c("Gene1", "Gene2")
dim(chip_reference) # 780470      2
head(chip_reference)

# separate rows with multiple genes in Gene2
chip_reference_clean <- chip_reference %>%
separate_rows(Gene2, sep = ",")
dim(chip_reference_clean) # 330      2

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
chip_reference_clean$Type = "+"

chip_reference_filtered = chip_reference_clean[(chip_reference_clean$Gene1 %in% gene_selected) & (chip_reference_clean$Gene2 %in% gene_selected),]
dim(chip_reference_filtered) # 330 3

# save results
write.csv(chip_reference_filtered, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_nishiyama.csv", row.names = FALSE, quote = FALSE)
