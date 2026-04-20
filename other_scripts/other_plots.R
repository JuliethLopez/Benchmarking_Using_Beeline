library(ggplot2)
library(tidyr)
library(dplyr)

# ----------------------------time plot--------------------------------
# read csv of times
times <- list()
datasets <- c("GSD", "HSC", "mCAD", "VSC")
cell_number <- c(2000, 1000, 200)

for (name in 1:length(datasets)){
  path <- paste0("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/Curated/", datasets[name], "/", datasets[name], "-Times.csv")
  data <- read.csv(path)
  
  method_names <- data[,1] 
  
  for (n in cell_number){
    data.colnames <- colnames(data)[grepl(paste0(".",n,"."), colnames(data))]
    times_cell_number <- data[,data.colnames]
    
    res <- rowMeans(times_cell_number, na.rm = TRUE)
    names(res) <- method_names
    times[[datasets[name]]][[as.character(n)]] <- res
  }
}

#transforma to a dataframa
plot_df <- data.frame()

for (ds in names(times)) {
  for (cn in names(times[[ds]])) {
    temp <- data.frame(
      Dataset = ds,
      CellNumber = cn,
      Method = names(times[[ds]][[cn]]),
      Runtime = as.numeric(times[[ds]][[cn]])
    )
    plot_df <- rbind(plot_df, temp)
  }
}

# Ensure CellNumber is a factor for discrete coloring
plot_df$CellNumber <- factor(plot_df$CellNumber, levels = c("200", "1000", "2000"))

# graphic
times_plot <- ggplot(plot_df, aes(x = CellNumber, y = Runtime, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  # Create a grid (facet) for each method
  facet_wrap(~ Dataset, scales = "free_y") + 
  theme_bw() +
  labs(
    title = "Method Runtimes Across Datasets",
    subtitle = "Faceted by Method, colored by Cell Number",
    x = "Dataset",
    y = "Average Execution Time (seconds)",
    color = "Cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )

ggsave(times_plot, filename = "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs_graphics/Curated/Method_Runtimes.png", width = 12, height = 8)


# ----------------------------motifs plot--------------------------------
# read csv of times
motifs_list <- list()
datasets <- c("GSD", "HSC", "mCAD", "VSC")
motifs <- c("FBL", "FFL", "MI")
cell_number <- c(2000, 1000, 200)

for (name in 1:length(datasets)){
    for (mtf in motifs){
        path <- paste0("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/Curated/", datasets[name], "/", datasets[name], "-NetworkMotifs-", mtf,".csv")
        data <- read.csv(path)
    
        method_names <- data[,1] 
    
        for (n in cell_number){
            data.colnames <- colnames(data)[grepl(paste0(".",n,"."), colnames(data))]
            motifs_cell_number <- data[,data.colnames]
            
            res <- rowMeans(motifs_cell_number, na.rm = TRUE)
            names(res) <- method_names
            motifs_list[[datasets[name]]][[paste(mtf,as.character(n),sep="_")]] <- res
        }
    }
}

#transforma to a dataframa
plot_df <- data.frame()

for (ds in names(motifs_list)) {
  for (cn in names(motifs_list[[ds]])) {
    temp <- data.frame(
      Dataset = ds,
      Motif_cellnumber = cn,
      Method = names(motifs_list[[ds]][[cn]]),
      Motifs_score = as.numeric(motifs_list[[ds]][[cn]])
    )
    plot_df <- rbind(plot_df, temp)
  }
}

#create columns for motifs and number of cells
plot_df$Motif_type <- sapply(as.character(plot_df$Motif_cellnumber), function(x) strsplit(x, "_")[[1]][1])
plot_df$cell_number <- sapply(as.character(plot_df$Motif_cellnumber), function(x) strsplit(x, "_")[[1]][2])
plot_df$motif_method <- paste(plot_df$Motif_type, plot_df$Method, sep = "_")

# Ensure cell_number is a factor for discrete coloring
plot_df$cell_number <- factor(plot_df$cell_number, levels = c("200", "1000", "2000"))

# categorical variables as factors for correct ordering
#method_order <- c("hdWGCNA", "INFERELATOR", "GRNBOOST2", "BORDA", "mBORDA", "sBORDA", "smBORDA")
method_order <- c("INFERELATOR", "GRNBOOST2")
plot_df$Method <- factor(plot_df$Method, levels = method_order)

# Heatmap of motifs scores
plot_df_non_borda <- subset(plot_df, Method %in% method_order)
plot_df_non_borda$method_data <- paste(plot_df_non_borda$Method,plot_df_non_borda$Dataset, sep = "_")

# plot
p <- ggplot(plot_df_non_borda, aes(x = cell_number, y = method_data, fill = Motifs_score)) +
  geom_tile(color = "white") +  # color adds a border between tiles
  # creates the grouping: cell_number on top, split by Motif_type below
  facet_grid(. ~ Motif_type, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 1, name = "Raw Score") +
  theme_minimal() +
  labs(
    title = "Heatmap of Motif Scores",
    subtitle = "Grouped by Cell Number and Motif Type",
    x = "Motif Type",
    y = "Method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.2, "lines") # Adjusts spacing between cell_number blocks
  )

ggsave(p, filename = paste0("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs_graphics/Curated/Heatmap_Motifs.jpeg"), width = 12, height = 8)


# ----------------------------memory plot--------------------------------
# read csv of times
memory <- list()
datasets <- c("GSD", "HSC", "mCAD", "VSC")
cell_number <- c(2000, 1000, 200)

for (name in 1:length(datasets)){
  path <- paste0("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/Curated/", datasets[name], "/", datasets[name], "-Memory.csv")
  data <- read.csv(path)
  
  method_names <- data[,1] 
  
  for (n in cell_number){
    data.colnames <- colnames(data)[grepl(paste0(".",n,"."), colnames(data))]
    memory_cell_number <- data[,data.colnames]
    
    res <- rowMeans(memory_cell_number, na.rm = TRUE)
    names(res) <- method_names
    memory[[datasets[name]]][[as.character(n)]] <- res
  }
}

#transforma to a dataframa
plot_df <- data.frame()

for (ds in names(memory)) {
  for (cn in names(memory[[ds]])) {
    temp <- data.frame(
      Dataset = ds,
      CellNumber = cn,
      Method = names(memory[[ds]][[cn]]),
      Runtime = as.numeric(memory[[ds]][[cn]])
    )
    plot_df <- rbind(plot_df, temp)
  }
}

# Ensure factor for discrete coloring
plot_df$CellNumber <- factor(plot_df$CellNumber, levels = c("200", "1000", "2000"))

# categorical variables as factors for correct ordering
method_order <- c("hdWGCNA", "INFERELATOR", "GRNBOOST2", "PIDC", "PPCOR")
plot_df$Method <- factor(plot_df$Method, levels = method_order)

# Heatmap of motifs scores
plot_df_non_borda <- subset(plot_df, Method %in% method_order)
plot_df_non_borda$Runtime_GB <- plot_df_non_borda$Runtime/1000000 # convert to GB

# plot
memory_plot <- ggplot(plot_df_non_borda, aes(x = CellNumber, y = Runtime_GB, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  # Create a grid (facet) for each method
  facet_wrap(~ Dataset, scales = "free_y") + 
  theme_bw() +
  labs(
    title = "Method Memory Usage Across Datasets",
    subtitle = "Faceted by Method, colored by Cell Number",
    x = "Dataset",
    y = "Average Memory Usage (GB)",
    color = "Cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )

ggsave(memory_plot, filename = "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs_graphics/Curated/Method_Memory.png", width = 12, height = 8)
