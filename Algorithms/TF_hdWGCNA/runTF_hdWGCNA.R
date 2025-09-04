############################
# load hdWGCNA usual libraries
# single-cell analysis package
library(Seurat)

# plotting and data science packages:
library(tidyverse)
library(cowplot)
library(patchwork)
library(magrittr)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# for gene regulatory network
library(JASPAR2020)
library(JASPAR2024)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(xgboost)
library(RSQLite)

print("hdWGCNA libraries loaded successfully")


# set global options
# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)


# Read input arguments
args <- commandArgs(trailingOnly = T)
inFile <- args[1]
outFile <-  args[2]

print(inFile)
print(outFile)
outdir <- file.path(dirname(outFile), "plots")
if(!dir.exists(outdir)) dir.create(outdir)
print(paste0("outdir:", outdir))

# Add TOM folder
outdir_TOM <- file.path(dirname(outFile), "TOM")
if(!dir.exists(outdir_TOM)) dir.create(outdir_TOM)

# input expression data
inputExpr <- read.table(inFile, sep=",", header = 1, row.names = 1)
print(inputExpr[1:5,1:5])

inputExpr <- as.matrix(inputExpr)

#geneNames <- rownames(inputExpr)
#rownames(inputExpr) <- c(geneNames)

# read as seurat
seurat_object <- CreateSeuratObject(counts = inputExpr, project = "GSD", min.cells = 0, min.features = 0)
seurat_object

# preprocesing
seurat_object@assays$RNA$data = seurat_object@assays$RNA$counts
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object)
ggsave(file.path(outdir, "VariableFeatures.png"), plot1)

#all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = VariableFeatures(seurat_object))

seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object))
plot2 <- DimPlot(seurat_object, reduction = "pca") + NoLegend()
ggsave(file.path(outdir, "PCA.png"), plot2)

seurat_object <- RunUMAP(seurat_object, dims = 1:10)
plot3 <- DimPlot(seurat_object, reduction = "umap")
ggsave(file.path(outdir, "UMAP.png"), plot3)

# clustering
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

plot4 <- DimPlot(seurat_object, reduction = "umap", group.by = "RNA_snn_res.0.5")
ggsave(file.path(outdir, "UMAP_clustering.png"), plot4)

# genes expresados en al menos el 5% de las células
seurat_obj <- SetupForWGCNA(
  seurat_object,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)


# metacells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("RNA_snn_res.0.5"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 5, # nearest-neighbors parameter
  max_shared = 2, # maximum number of shared cells between two metacells
  ident.group = 'RNA_snn_res.0.5' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)
cell_types <- unique(seurat_obj@misc$tutorial$wgcna_metacell_obj@meta.data$RNA_snn_res.0.5)

seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))

plot_pca <- DimPlotMetacells(seurat_obj, reduction="pca")
ggsave(filename = file.path(outdir,"pca_metacells.png") , plot=plot_pca)

#not necessary
#seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='orig.ident')

#run dimensional reduction on metacells
seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='pca', dims=1:15)

umap_mc <- DimPlotMetacells(seurat_obj, group.by='RNA_snn_res.0.5') + umap_theme() + ggtitle("Cell Type Metacells")
ggsave(file.path(outdir,"UMAP_metacells.png"),umap_mc)

# with metacells
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = cell_types, # the name of the group of interest in the group.by column
  group.by='RNA_snn_res.0.5', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)


# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
png(file.path(outdir, "SoftPower.png"))
wrap_plots(plot_list, ncol=2)
dev.off()

power_table <- GetPowerTable(seurat_obj)
write_csv2(power_table, file = file.path(outdir,"power_table.csv"))

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 10,
  tom_name = 'GSD2_non_metacell', # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE
)

# get TOM (tarda mucho. Es necesario despues, pero por ahora se puede seguir sin extraerla del objeto seurat)
# mirar bien donde se guarda
#TOM <- GetTOM(seurat_obj)

png(file.path(outdir, "Dendogram_GSD2_non_metacell.png"))
PlotDendrogram(seurat_obj, main='INH hdWGCNA Dendrogram')
dev.off()

# Identify TFs in genes promoter regions (algorithm MotifScan)
# 1. list of position weight matrices for each TF

JASPAR2024 <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))

pfm_core <- TFBSTools::getMatrixSet(
  x = sq24,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# run the motif scan using: TF and binding sites from curated data base, our exiting data (genes + TF)
# Extract promoter sequencing from EnsDb
# Scan those sequences to detect motives from pfm
# delever matric with moptif per gene and targets per TF
seurat_obj <- MotifScan(
  seurat_obj,
  species_genome = 'hg38',
  pfm = pfm_core,
  EnsDb = EnsDb.Hsapiens.v86
)


# COnstruct the nextwork TF_hdWGCNA
# https://smorabit.github.io/hdWGCNA/articles/tf_network.html

# get the motif df
motif_df <- GetMotifs(seurat_obj)
write_csv2(motif_df, file = file.path(outdir,"motifs_df_non_metacell.csv"))

# keep all TFs, and then remove all genes from the grey module
tf_genes <- unique(motif_df$gene_name)
modules <- GetModules(seurat_obj)
write_csv2(modules, file = file.path(outdir,"modules.csv"))

# grey modules contain genes that did not clustered properly to other modules, so can be ignored
nongrey_genes <- subset(modules, module != 'grey') %>% .$gene_name
genes_use <- c(tf_genes, nongrey_genes)

# update the gene list and re-run SetDatExpr
seurat_obj_filter <- SetWGCNAGenes(seurat_obj, genes_use)
seurat_obj_filter <- SetDatExpr(seurat_obj_filter) # group.by = 'cell_type', group_name='EX'

# define model params:
model_params <- list(
  objective = 'reg:squarederror',
  max_depth = 1,
  eta = 0.1,
  nthread=16,
  alpha=0.5
)

# construct the TF network
# Gain: contribution of each feature to the model. Highest percentage means important feature to predict the label used for the training (only available for tree models)
# Cover: metric of the number of observation related to this feature 
# Weight/Frequency: percentage representing the relative number of times a feature have been taken into trees.
# Cor: correlation between feature (TF) and gene expression
seurat_obj_filter <- ConstructTFNetwork(seurat_obj_filter, model_params) 

# extract results
results <- GetTFNetwork(seurat_obj_filter)
# Write output to a file
write.table(results, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
