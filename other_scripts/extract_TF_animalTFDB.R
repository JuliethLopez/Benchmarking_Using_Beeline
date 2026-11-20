# It were explored different reference TF lists
# AnimalTFDB and BeelineTF

# read reference TF lists
AnimalTFDB <- read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/Homo_sapiens_TF.txt", sep="\t")
head(AnimalTFDB)
dim(AnimalTFDB) # 1659 6
BeelineTF <- read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/human-tfs.csv")
head(BeelineTF)
dim(BeelineTF) # 1563 2

# read data
data <- read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/han_GSE107552.csv", row.names=1)
dim(data) # 200

# revise which genes are TF
genes = rownames(data)
length(genes)

list_AnimalTF = genes[genes %in% AnimalTFDB$GeneSymbol]
dim(list_AnimalTF) # 0  of 200 genes
list_TF = genes[genes %in% BeelineTF$TF]
length(list_TF) # 21 of 200 genes
head(list_TF)

#save results
write.csv(list_TF, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/list_TF.csv")
