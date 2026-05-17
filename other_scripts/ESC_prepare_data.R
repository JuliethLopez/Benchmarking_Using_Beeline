# libraries
library(tidyr)
library(dplyr)

# read data
data = read.csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_chipunion.txt", sep="\t")
head(data)
dim(data)

# add colnames
colnames(data) = c("Gene1", "Gene2")

# separate rows with multiple genes in Gene2
data_clean <- data %>%
  separate_rows(Gene2, sep = ",")
dim(data_clean)

# add type of interaction, all will be positive because it is not known if the interaction is for activation or repression
data_clean$Type = "+"

# save results
write.csv(data_clean, "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/inputs/experimental/ESC/hESC_chipunion.csv", row.names = FALSE, quote = FALSE)
