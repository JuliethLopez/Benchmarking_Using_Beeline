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

# Read input arguments
args <- commandArgs(trailingOnly = T)
inFile <- paste0("/",args[1])
outFile <-  paste0("/",args[2])
run_tf <- as.logical(args[3])
nThreads <- as.numeric(args[4])
nVariableFeatures <- as.numeric(args[5])
clustering_resolution <- as.numeric(args[6])
metacells_k <- as.numeric(args[7])
min_cells_per_metacell <- as.numeric(args[8])
max_shared_cells <- as.numeric(args[9])
networkType <- args[10]
soft_power <- as.numeric(args[11])
delete_grey_module <- as.logical(args[12])

print(paste0("inFile:", inFile))
print(paste0("outFile:", outFile))
print(paste0("run_tf:", run_tf))
print(paste0("nThreads:", nThreads))
print(paste0("nVariableFeatures:", nVariableFeatures))
print(paste0("clustering_resolution:", clustering_resolution))
print(paste0("metacells_k:", metacells_k))
print(paste0("min_cells_per_metacell:", min_cells_per_metacell))
print(paste0("max_shared_cells:", max_shared_cells))
print(paste0("networkType:", networkType))
print(paste0("soft_power:", soft_power))
print(paste0("Delete_grey_modules:", delete_grey_module))

# optionally enable multithreading
enableWGCNAThreads(nThreads = nThreads)

outdir <- file.path(dirname(outFile), "plots")
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
print(paste0("outdir:", outdir))

# Add TOM folder
outdir_TOM <- file.path(dirname(outFile), "TOM")
if(!dir.exists(outdir_TOM)) dir.create(outdir_TOM, recursive = TRUE)
print(paste0("TOM folder:", outdir_TOM))

# Remove old Tom files to ensure not loading of old results
files <- c(
  Sys.glob(file.path(outdir_TOM,"blockwiseTOM*.RData")),
  Sys.glob(file.path(outdir_TOM,"individualTOM*.RData"))
)

if (length(files) > 0) {
  message("Removing old WGCNA TOM files")
  file.remove(files)
}

# input expression data
#print("Prove locations!")
#getwd()
#list.files("/data/inputs/example/GSD/TF_hdWGCNA")
#list.files("data/inputs/example/GSD/TF_hdWGCNA")

# ---------- Basic WGCNA-----------

inputExpr <- read.csv(inFile,row.names = 1,check.names = FALSE) # first column = gene names
inputExpr <- as.matrix(inputExpr)
print("Gene expression data successfully readed")
print(class(inputExpr[1:5,1:5]))
print(dim(inputExpr))# genes x cells

n_cells <- ncol(inputExpr)
print(paste0("The number of cells is: ",n_cells))

geneNames <- rownames(inputExpr)
rownames(inputExpr) <- c(geneNames)
print("gene names set")
print(head(geneNames)[1:5])

# read as seurat
seurat_object <- CreateSeuratObject(counts = inputExpr, project = "GSD", min.cells = 0, min.features = 0)
seurat_object

# preprocesing
seurat_object@assays$RNA$data = seurat_object@assays$RNA$counts
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = nVariableFeatures)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object)
ggsave(file.path(outdir, "VariableFeatures.png"), plot1)

#all.genes <- rownames(seurat_object)
variable_features <- VariableFeatures(seurat_object)
if (length(variable_features) < nVariableFeatures) {
  print(paste0("The number of variable features found (", length(variable_features), ") is less than nVariableFeatures (", nVariableFeatures, "). Using all available variable features."))
  variable_features <- rownames(seurat_object@assays$RNA$data)
}
#sacle data
seurat_object <- ScaleData(seurat_object, features = variable_features)

# Run PCA
seurat_object <- RunPCA(seurat_object, features = variable_features)
plot2 <- DimPlot(seurat_object, reduction = "pca") + NoLegend()
ggsave(file.path(outdir, "PCA.png"), plot2)

# UMAP
n_pcs <- 10
n_pcs_available <- ncol(seurat_object@reductions$pca@cell.embeddings)
print(paste0("Number of PCs in the single cell seurat object: ", n_pcs_available))
if(n_pcs_available < n_pcs){
  print(paste0("The number of PCs available (", n_pcs_available, ") is less than ", n_pcs,". Then it will be used all available PCs."))
  n_pcs <- n_pcs_available
}
seurat_object <- RunUMAP(seurat_object, dims = 1:n_pcs)
plot3 <- DimPlot(seurat_object, reduction = "umap")
ggsave(file.path(outdir, "UMAP.png"), plot3)

# clustering
seurat_object <- FindNeighbors(seurat_object, dims = 1:n_pcs)
seurat_object <- FindClusters(seurat_object, resolution = clustering_resolution)

metacells_variable_group_by <- paste0("RNA_snn_res.", clustering_resolution)

print("Saving UMAP plot with clusters")
plot4 <- DimPlot(seurat_object, reduction = "umap", group.by = metacells_variable_group_by)
ggsave(file.path(outdir, "UMAP_clustering.png"), plot4)


# se podría filtrar genes expresados en al menos el 5% de las células, pero para las pruebas con redes simuladas no funciona
print(paste0("All variable features: ", length(variable_features)))
seurat_obj <- SetupForWGCNA(
  seurat_object,
  gene_select = "fraction", # the gene selection approach, custom
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  #features = variable_features,
  wgcna_name = "hdWGCNA" # the name of the hdWGCNA experiment
)


# metacells
if(n_cells > 5000){
  print("Creating metacells per cluster")
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c(metacells_variable_group_by), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'pca', # select the dimensionality reduction to perform KNN on
    k = metacells_k, # nearest-neighbors parameter. 25 by default in the fuction
    max_shared = max_shared_cells, # maximum number of shared cells between two metacells. 15 by default in the function
    min_cells = min_cells_per_metacell, # minimum number of cells in group to construct metacell. 100 by default in the function
    ident.group = metacells_variable_group_by # set the Idents of the metacell seurat object
  )

  # normalize metacell expression matrix:
  seurat_obj <- NormalizeMetacells(seurat_obj)
  cell_types <- unique(seurat_obj@misc$hdWGCNA$wgcna_metacell_obj@meta.data[[metacells_variable_group_by]])

  seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
  seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))

  # look for metacells object
  print("Object with metacells calculated:")
  print(seurat_obj@misc$hdWGCNA$wgcna_metacell_obj)
  #saveRDS(seurat_obj, file.path(outdir,"metacells.RDS"))

  print("Saving PCA plot with metacells")
  plot_pca <- DimPlotMetacells(seurat_obj, reduction="pca")
  ggsave(filename = file.path(outdir,"pca_metacells.png") , plot=plot_pca)

  #not necessary
  #seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='orig.ident')

  print("Metacells object summary:")
  print(seurat_obj)
  print(seurat_obj@reductions$pca)

  #run dimensional reduction on metacells
  n_pcs_metacells <- 15
  n_pcs_metacells_available <- ncol(seurat_obj@reductions$pca@cell.embeddings)
  print(paste0("Number of PCs in the metacells seurat object: ", n_pcs_metacells_available))
  if(n_pcs_metacells_available < n_pcs_metacells){
    print(paste0("The number of PCs available (", n_pcs_metacells_available, ") is less than ", n_pcs_metacells,". Then it will be used all available PCs."))
    n_pcs_metacells <- n_pcs_metacells_available
  }
  # Run UMAP on metacells
  # print("Running UMAP on metacells")
  # X <- Embeddings(seurat_obj, reduction = "pca")
  # dim(X)
  # n_dim <- ncol(X)
  # if(n_dim < 30){
  #   print(paste0("The number of PCs (", n_dim, ") is less than 30. Setting n.neighbors to min(5,", n_dim - 1,")"))
  # }
  seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='pca', dims=1:n_pcs_metacells) #   n.neighbors = min(5, n_dim - 1))

  umap_mc <- DimPlotMetacells(seurat_obj, group.by=metacells_variable_group_by) + umap_theme() + ggtitle("Cell Type Metacells")
  ggsave(file.path(outdir,"UMAP_metacells.png"),umap_mc)

  # COnfig matrix of co-expression with metacells
  print(paste("Metadata group by:", metacells_variable_group_by))
  print(head(seurat_obj@meta.data[metacells_variable_group_by]))
  print(paste("Group names:", cell_types))
  print(cell_types)
  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = cell_types, # the name of the group of interest in the group.by column
    group.by = metacells_variable_group_by, # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = 'RNA', # using RNA assay
    layer = 'data' # using normalized data
  )

  # # take a look at the metadata
  # wgcna_name <- seurat_obj@misc$active_wgcna
  # seurat_obj@misc[[wgcna_name]]$datExpr
} else {
  paste0("The number of cells (", n_cells, ") is not higher than 6000. Metacells may not be reliable with low cell numbers, so metacells will not be calculated. Using full Seurat object.")
  
  seurat_obj$dummy_group <- "all"

  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = "all",
    group.by = "dummy_group",
    assay = "RNA",
    layer = "data",
    use_metacells = FALSE
  )
}

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = networkType # it is possible to use "signed", "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
png(file.path(outdir, "SoftPower.png"))
wrap_plots(plot_list, ncol=2)
dev.off()

power_table <- GetPowerTable(seurat_obj)
write_csv2(power_table, file = file.path(outdir,"power_table.csv"))

# select first high soft power value
scale_free_topology = power_table$SFT.R.sq
names(scale_free_topology) = power_table$Power

select_soft_power <- function(lista) {
  # if there is only one value
  if (length(lista) == 1) {
    return(as.numeric(names(lista)[1]))
  }
  else {
    # for several values take the minimum with free scale score higher than 0.8, if not, select the first maximum
    if (any(lista > 0.8)){
      valid = lista[lista > 0.8]
      idx = which.min(valid)
      return(as.numeric(names(valid)[idx]))
    } else {
      # first maximum
      for (i in 2:length(lista)) {
        if (lista[i] < lista[i - 1]) {
          soft_power = names(lista)[i - 1]
          return(as.numeric(soft_power))
        }
      }
      # if the value never decreases, use the maximum of all
      soft_power = names(lista)[lista==max(lista)]
      return(as.numeric(soft_power))
    }
  }  
}

if (soft_power == 0) {
  print("Soft_power is cero... Calculating soft power...")
  
  soft_power = as.numeric(select_soft_power(scale_free_topology))
  print("Soft power set to first maximum value:")
  print(soft_power)
} else {
  soft_power_value = scale_free_topology[as.character(soft_power)]
  print(paste0("Soft power set to: ", soft_power, " with scale free topology value of: ", soft_power_value))
  soft_power_computed = as.numeric(select_soft_power(scale_free_topology))
  soft_power_computed_value = scale_free_topology[as.character(soft_power_computed)]
  print(paste0("Soft power computed value: ", soft_power_computed, " with scale free topology value of: ", soft_power_computed_value))
  if (soft_power_value!=soft_power_computed_value){
    print("Warning: Computed soft power value differs from specified value. Computed soft value will be used for network construction")
    soft_power = soft_power_computed
  }
}

# Construct co-expression network. This function calls blockwiseConsensusModules, which uses the parameter minModuleSize for minimum module size for module detection. This value by default is 20, however given the small number of genes we have in total, it is required a smaller value.
print("Min module size set to:")
if (length(variable_features) < nVariableFeatures){
  min_module_size = floor(length(variable_features)/2)
  
} else {min_module_size = 20}
print(min_module_size)

# Try different soft_power
soft_power_try <- soft_power
success <- FALSE

while (soft_power_try >= 1 && !success) {

  message(paste("Trying soft power:", soft_power_try))

  tryCatch({

    seurat_obj <- ConstructNetwork(
      seurat_obj,
      soft_power = soft_power_try,
      tom_name = 'metacells', # name of the topoligical overlap matrix written to disk
      overwrite_tom = TRUE,
      tom_outdir = outdir_TOM,
      minModuleSize = min_module_size
    )

    success <- TRUE
    message(paste("Success with soft power:", soft_power_try))

  }, error = function(e) {

    message(paste("Failed with soft power:", soft_power_try))
    message("Error:", e$message)

    soft_power_try <- soft_power_try - 1
  })
}


# PLot dendrogram
png(file.path(outdir, "Dendogram_tom_metacells.png"))
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
dev.off()

# save modules
modules <- GetModules(seurat_obj)
write_csv2(modules, file = file.path(outdir,"modules.csv"))

# grey modules contain genes that did not clustered properly to other modules, so can be ignored
nongrey_genes <- subset(modules, module != 'grey') %>% .$gene_name

# ---------- Transcription factor regulatory network analysis WGCNA-----------

if(run_tf==TRUE){
  print("Running transcription factor extension...")
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
  # delever metric with moptif per gene and targets per TF
  seurat_obj <- MotifScan(
    seurat_obj,
    species_genome = 'hg38',
    pfm = pfm_core,
    EnsDb = EnsDb.Hsapiens.v86
  )


  # Construct the nextwork TF_hdWGCNA
  # https://smorabit.github.io/hdWGCNA/articles/tf_network.html

  # get the motif df
  motif_df <- GetMotifs(seurat_obj)
  write_csv2(motif_df, file = file.path(outdir,"motifs_df.csv"))

  # keep all TFs, and then remove all genes from the grey module
  tf_genes <- unique(motif_df$gene_name)
  
  if(delete_grey_module == TRUE){
    genes_use <- c(tf_genes, nongrey_genes)
  } else {
    genes_use <- tf_genes
  }


  # update the gene list and re-run SetDatExpr
  seurat_obj_filter <- SetWGCNAGenes(seurat_obj, genes_use)
  seurat_obj_filter <- SetDatExpr(seurat_obj_filter) # Also can be added group.by = 'cell_type', group_name='EX'

  # define model params, this are used for XGBoost
  model_params <- list(
    objective = 'reg:squarederror',
    max_depth = 1,
    eta = 0.1,
    nthread=nThreads,
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
  #results <- results[,c("tf", "gene", "Cor")]
  #colnames(results) = c('Gene1','Gene2','EdgeWeight')
} else {
  # gene names
  gene_name <- modules$gene_name
  #nongrey_genes <- modules[modules$module!="grey","gene_name"]

  # load TOM matrix
  tom_path <- file.path(outdir_TOM, "metacells_TOM.rda")
  #tom_path2 <- file.path(outdir_TOM, "metacells_block.1.rda")

  # repeat {
  #   if (file.exists(tom_path1)) {
  #     tom_path <- tom_path1
  #     break
  #   } else {
  #     if (file.exists(tom_path2)) {
  #       tom_path <- tom_path2
  #       break
  #     }
  #   }
  # }

  print(paste0("Waiting for TOM file: ", tom_path))

  repeat {
    if (file.exists(tom_path)) {
      s1 <- file.info(tom_path)$size
      Sys.sleep(3)
      s2 <- file.info(tom_path)$size
      if (!is.na(s1) && s1 == s2) break
    } else {
      Sys.sleep(3)
      }
    }

  print(paste0("Loading TOM matrix file from: ", tom_path))
  load(tom_path)
  #TOM <- GetTOM(seurat_obj) # do not works well because the folder in container do not match with real location
  print("TOM matrix loaded successfully")
  
  # TOM as matrix
  TOM <- 1 - as.matrix(consTomDS) # as this is a distance object "TOM dissimilarity", is it necessary to tranform the data to obtain node connections
  rownames(TOM) <- gene_name
  colnames(TOM) <- gene_name

  # filter grey modules
  TOM <- TOM[nongrey_genes, nongrey_genes]
  edges <- which(TOM>0, arr.ind=TRUE) #filter according to edge strength. As this could overstimate the result of this method in the benchmarking, only used the minimum value, 0
  #edges <- edges[edges[,1] < edges[,2], ] #delete duplicated in not directed networks, not applicated to be fair with other similar methods, as PPCOR
  
  # create dataframe
  results <- data.frame(
    Gene1 = rownames(TOM)[edges[,1]],
    Gene2 = colnames(TOM)[edges[,2]],
    EdgeWeight = TOM[edges]
  )
  results_filt <- results[order(results$EdgeWeight),]
  #print("Keep more connected edges")
  #threshold <- 0.6**soft_power_try
  #results_filt$EdgeWeight[results_filt$EdgeWeight < threshold] = 0
  #threshold <- quantile(results$EdgeWeight, probs = 0.8)
  #results_filt$EdgeWeight[results_filt$EdgeWeight < threshold] = 0
}

# --------------Write output to a file-----------------

print("This are the results GetTFNetwork")
print(head(results_filt))

write.table(results_filt, outFile, sep = "\t", quote = FALSE, row.names = FALSE)