# library
library(data.table)

# set directory
setwd("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/Curated/GSD")

# iterate for each folder/sample
folders = list.dirs(recursive=FALSE, full.names=FALSE)

for (folder in folders) {
    folder_sample = paste0("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/Curated/GSD/",folder)
    new_folder = paste0(folder_sample,"/hdWGCNA")
    dir.create(new_folder)

    # read modules
    modules <- read.csv2(paste0(folder_sample,"/TF_hdWGCNA/plots/modules.csv"))
    genes <- modules$gene_name
    non_grey_module <- modules[modules$module!="grey","gene_name"]

    # name of file consTomDS
    load(paste0(folder_sample,"/TF_hdWGCNA/TOM/tom_metacells_TOM.rda"))
    TOM <- as.matrix(consTomDS)
    edges <- which(TOM>0.1, arr.ind=TRUE)
    rownames(TOM) <- genes
    colnames(TOM) <- genes

    # filter grey modules
    TOM <- TOM[non_grey_module, non_grey_module]

    edges <- data.frame(
        Gene1 = rownames(TOM)[row(TOM)[upper.tri(TOM)]],
        Gene2 = colnames(TOM)[col(TOM)[upper.tri(TOM)]],
        EdgeWeight = TOM[upper.tri(TOM)]
    )

    dim(edges)

    # melt dataset
    #edges <- as.data.frame(as.table(TOM_module))
    #colnames(edges) <- c("Gene1","Gene2","EdgeWeight")

    #edges <- edges[edges$Gene1 != edges$Gene2,]
    edges <- edges[edges$EdgeWeight > 0.1, ]

