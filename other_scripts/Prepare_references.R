# This script prepares the references
# hTFtargets. Here tissue specific regulations are extracted
# RegNetwork. Here positive interactions are added if those are unknown
# STRING and STRING regulatory interactions. Here gene names are added to the protein ids

#library
library(biomaRt)
library(readr)
library(dplyr)

library(uniprotREST)
library(dplyr)
library(tidyr)

setwd("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth")

# hTFtargets
#########################################
data_tf <- read_table("hTFtargets_TF-Target-information.txt")

# filter by tissue
data_tf_liver <- data_tf[grep("liver",data_tf$tissue), c("TF","target")]
colnames(data_tf_liver) <- c("Gene1","Gene2")
data_tf_liver$Type <- "+"
write.csv(data_tf_liver, "hTFtargets_TF-Target-liver.csv", row.names = FALSE, quote = FALSE)

data_tf_embryo <- data_tf[grep("embryo",data_tf$tissue), c("TF","target")]
colnames(data_tf_embryo) <- c("Gene1","Gene2")
data_tf_embryo$Type <- "+"
write.csv(data_tf_embryo, "hTFtargets_TF-Target-embryo.csv", row.names = FALSE, quote = FALSE)

# RegNetwork
#########################################
regnet <- read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/RegNetwork_human_TF_Target.txt", header = FALSE)
head(regnet)
regnet_filt <- regnet[,c(1,3)]
colnames(regnet_filt) <- c("Gene1","Gene2")
regnet_filt$Type <- "+"
write.csv(regnet_filt, "RegNetwork_human_TF_Target.csv", row.names = FALSE, quote = FALSE)

# STRING
#########################################
#url <- "https://stringdb-downloads.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz"

#file_name <- "STRING_9606.protein.links.full.v12.0.txt.gz"

#download.file(url = url, destfile = file_name)

string <- read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/STRING_9606.protein.links.full.v12.0.txt.gz", header = TRUE)
head(string)
dim(string) #13715404       16

string_filt <- string[,c(1,2)]
string_filt$Type <- "+"

rm(string)

# remove the "9606." prefix from the protein ids
string_filt$protein1 <- sub("^9606\\.", "", string_filt$protein1)
string_filt$protein2 <- sub("^9606\\.", "", string_filt$protein2)
protein_ids <- unique(unique(string_filt$protein1), unique(string_filt$protein2))
length(protein_ids) #19622

# get gene names from ensembl protein ids
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

annot <- getBM(
  attributes = c(
    "ensembl_peptide_id", #unique code for protein
    "ensembl_gene_id", #unique code for gene
    "external_gene_name", #gene name
    #"description",
    "gene_biotype" #protein-coding, immunoglobuling variable chain
  ),
  filters = "ensembl_peptide_id",
  values = protein_ids,
  mart = mart
)

dim(annot) #19168     4

#write.csv(annot, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/STRING_ensembl_prot_gene_2026.csv", row.names = FALSE, quote = FALSE)
annot <- read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/STRING_ensembl_prot_gene_2026.csv", header = TRUE)

# filter annot table
table(annot$gene_biotype)

annot_pc <- annot %>% filter(gene_biotype == "protein_coding") %>% 
  select(ensembl_peptide_id, external_gene_name)
dim(annot_pc) #19167     2

not_annot <- protein_ids[!(protein_ids %in% annot_pc$ensembl_peptide_id)]
length(unique(annot_pc$external_gene_name)) #18916, there is more than one gene that codify for the same protein
length(unique(annot_pc$ensembl_peptide_id)) #19167. Lacking 453

# get gene names from uniprot
res <- uniprot_map(
    ids = protein_ids,
    from = "Ensembl_Protein",
    to = "UniProtKB"
)

res_reviewed = subset(res, Reviewed == "reviewed")
dim(res_reviewed) #18348 8
head(res_reviewed)

# select proteins with only one gene related
res_reviewed_unique = res_reviewed[grepl("\\s",res_reviewed$`Gene.Names`)==FALSE,]
dim(res_reviewed_unique) #5418 8
head(res_reviewed_unique)

table(not_annot %in% res_reviewed_unique$From) #we can add 3 proteins IDs more

res_reviewed_with_annot <- res_reviewed_unique %>% 
  filter(From %in% not_annot) %>%
  mutate(ensembl_peptide_id = From,
    external_gene_name = Gene.Names) %>%
  distinct(ensembl_peptide_id, .keep_all = TRUE) %>%
  select(ensembl_peptide_id, external_gene_name)

dim(res_reviewed_with_annot)

annot_rbind = annot_pc %>% rbind(res_reviewed_with_annot)

head(annot_rbind)
dim(annot_rbind) #19170 2

# replace protein ids with gene names
head(string_filt)
string_genes <- string_filt %>% 
  left_join(annot_rbind, by = c("protein1" = "ensembl_peptide_id")) %>%
  select("Gene1" = "external_gene_name", "protein2", "Type") %>%
  left_join(annot_rbind, by = c("protein2" = "ensembl_peptide_id")) %>%
  select("Gene1", "Gene2" = "external_gene_name", "Type")

head(string_genes)
dim(string_genes)

write.csv(string_genes, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/STRING_reference_2026.csv", row.names = FALSE, quote = FALSE)


# STRING regulatory interactions
#########################################
string <- read.table("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/9606.protein.regulatory.links.full.v12.5.txt.gz", header = TRUE)
head(string)
dim(string) #165336     13

summary(string$combined_score) #min 42 transfered form any technique, I will keep all

# select direction if evidence is 0.5 higher than the other
string$Fold <- string$positive_transferred/string$negative_transferred
string$Type <- ifelse(string$Fold>=1.5, "+", ifelse(string$Fold<=0.5, "-", "0"))

# check if there are proteins that are both positive and negative
string[string$positive_transferred>0 & string$negative_transferred>0, c("Type","Fold","positive_transferred", "negative_transferred", "protein1", "protein2")] %>% head()
# The NA is 0/0, I will replace them with 0
string$Type[is.na(string$Type)] <- "0"

pos=sum(string$Type=="+", na.rm = TRUE) # 18568
neg=sum(string$Type=="-", na.rm = TRUE) # 6959
cer=sum(string$Type=="0", na.rm = TRUE) # 139809

(pos+neg)/nrow(string) #only 15% of the interactions are positive or negative taking double evidence for each one

# select only protein1 and protein2 and type
string_filt <- string[,c("protein1", "protein2", "Type")]
string_filt <- string_filt[string_filt$Type!="0",]
head(string_filt)
dim(string_filt) # 25527     3
rm(string)

# read regulatory events
regulatory_events <- read_tsv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/9606.protein.regulatory.events.v12.5.txt.gz")
head(regulatory_events[,c(1,2,5)])
sum(regulatory_events$regulation_of_gene_expression!=0)
colnames(regulatory_events)
dim(regulatory_events) # 134611 13

# how many regulatory events of gene expression are present in the string_filt table?
string_filt <- string_filt %>% left_join(regulatory_events[,c(1,2,5)], by = c("protein1" = "protein1", "protein2" = "protein2")) %>%
  filter(regulation_of_gene_expression!=0) %>% 
  filter(!is.na(regulation_of_gene_expression))
dim(string_filt) # 10396 4
head(string_filt)


# remove the "9606." prefix from the protein ids
string_filt$protein1 <- sub("^9606\\.", "", string_filt$protein1)
string_filt$protein2 <- sub("^9606\\.", "", string_filt$protein2)
protein_ids <- unique(unique(string_filt$protein1), unique(string_filt$protein2))
length(protein_ids) #1590

# get gene names from ensembl protein ids
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

annot <- getBM(
  attributes = c(
    "ensembl_peptide_id", #unique code for protein
    "ensembl_gene_id", #unique code for gene
    "external_gene_name", #gene name
    #"description",
    "gene_biotype" #protein-coding, immunoglobuling variable chain
  ),
  filters = "ensembl_peptide_id",
  values = protein_ids,
  mart = mart
)

dim(annot) #19168     4

#write.csv(annot, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/STRING_ensembl_prot_gene_2026.csv", row.names = FALSE, quote = FALSE)
annot <- read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/STRING_ensembl_prot_gene_2026.csv", header = TRUE)

# filter annot table
table(annot$gene_biotype)

annot_pc <- annot %>% filter(gene_biotype == "protein_coding") %>% 
  select(ensembl_peptide_id, external_gene_name)
dim(annot_pc) #19167     2

# check if there are proteins that are not annotated
not_annot <- protein_ids[!(protein_ids %in% annot_pc$ensembl_peptide_id)]
length(unique(annot_pc$external_gene_name)) #18915, there is more than one gene that codify for the same protein
length(unique(annot_pc$ensembl_peptide_id)) #19167. Lacking 19167

# get gene names from uniprot
res <- uniprot_map(
    ids = protein_ids,
    from = "Ensembl_Protein",
    to = "UniProtKB"
)

res_reviewed = subset(res, Reviewed == "reviewed")
dim(res_reviewed) #1523 8
head(res_reviewed)

# select proteins with only one gene related
res_reviewed_unique = res_reviewed[grepl("\\s",res_reviewed$`Gene.Names`)==FALSE,]
dim(res_reviewed_unique) #382 8
head(res_reviewed_unique)

table(not_annot %in% res_reviewed_unique$From) #we can add 1 proteins IDs more

res_reviewed_with_annot <- res_reviewed_unique %>% 
  filter(From %in% not_annot) %>%
  mutate(ensembl_peptide_id = From,
    external_gene_name = Gene.Names) %>%
  distinct(ensembl_peptide_id, .keep_all = TRUE) %>%
  select(ensembl_peptide_id, external_gene_name)

dim(res_reviewed_with_annot) #1 2

annot_rbind = annot_pc %>% rbind(res_reviewed_with_annot)

head(annot_rbind)
dim(annot_rbind) #19168 2

# replace protein ids with gene names
head(string_filt)
string_genes <- string_filt %>% 
  left_join(annot_rbind, by = c("protein1" = "ensembl_peptide_id")) %>%
  select("Gene1" = "external_gene_name", "protein2", "Type") %>%
  left_join(annot_rbind, by = c("protein2" = "ensembl_peptide_id")) %>%
  select("Gene1", "Gene2" = "external_gene_name", "Type")

head(string_genes)
dim(string_genes) #10396 3

write.csv(string_genes, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/STRING_regulatory_reference_2026.csv", row.names = FALSE, quote = FALSE)
