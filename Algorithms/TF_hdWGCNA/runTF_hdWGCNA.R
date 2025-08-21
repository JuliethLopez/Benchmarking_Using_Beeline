library(JASPAR2020)
library(JASPAR2024)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(xgboost)
library(RSQLite)

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

print("hdWGCNA libraries loaded successfully")

# Read input arguments
args <- commandArgs(trailingOnly = T)
inFile <- args[1]
outFile <-  args[2]

# input expression data
inputExpr <- read.table(inFile, sep=",", header = 1, row.names = 1)

print(head(inputExpr))
