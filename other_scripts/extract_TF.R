# read reference TF lists
AnimalTFDB <- read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/Homo_sapiens_TF.txt", sep="\t")
head(AnimalTFDB)
BeelineTF <- read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/human-tfs.csv")
head(BeelineTF)

# read data
data <- read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/han_GSE107552.csv", row.names=1)
data[1:3,1:3]

# revise which genes are TF
genes = rownames(data)
length(genes)

list_AnimalTF = genes[genes %in% AnimalTFDB$GeneSymbol]
head(list_AnimalTF)
list_TF = genes[genes %in% BeelineTF$TF]
head(list_TF)
length(list_TF)

#save results
write.csv(list_TF, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/list_TF.csv")
