##############################################################################
# Script information                                                      
# Title: immune_correlation
# Author: Minyu Zhou
# Date: 2024-12-27
# Description: None
##############################################################################

### === Spearman correlation analysis ===
library(dplyr)
library(writexl)
library(corrplot)

# Read metadata and microbiome data
meta_data <- read.table("meta.txt", header = TRUE, sep = "\t", row.names = 1)
microbiome_data <- read.table("microbiome.txt", header = TRUE, sep = "\t", row.names = 1)

# Optional: transpose microbiome data if samples are in columns (default format)
if (!all(rownames(meta_data) %in% rownames(microbiome_data))) {
  message("Transposing microbiome_data to match meta_data structure (samples as rows)...")
  microbiome_data <- t(microbiome_data)
}

# Align samples between metadata and microbiome data
common_samples <- intersect(colnames(meta_data), rownames(microbiome_data_t))
meta_data_aligned <- meta_data[, common_samples, drop = FALSE]
microbiome_data_aligned <- microbiome_data_t[common_samples, , drop = FALSE]

# Initialize result table
results <- data.frame(MetaVariable = character(),
                      Microbe = character(),
                      Correlation = numeric(),
                      PValue = numeric(),
                      stringsAsFactors = FALSE)

# Perform Spearman correlation between each meta variable and each microbe
for (meta_var in rownames(meta_data_aligned)) {
  for (microbe in colnames(microbiome_data_aligned)) {
    x <- as.numeric(meta_data_aligned[meta_var, ])
    y <- as.numeric(microbiome_data_aligned[, microbe])
    
    if (!all(is.na(x)) && !all(is.na(y))) {
      cor_test <- cor.test(x, y, method = "spearman")
      results <- rbind(results, data.frame(
        MetaVariable = meta_var,
        Microbe = microbe,
        Correlation = cor_test$estimate,
        PValue = cor_test$p.value
      ))
    }
  }
}

# Save results to Excel
write_xlsx(results, path = "spearman_correlation_results.xlsx")


### === Visualization ===

# Read correlation results (or reuse the 'results' object directly)
data <- results  # If already in memory, no need to reload

# Create correlation matrix
meta_vars <- unique(data$MetaVariable)
microbes <- unique(data$Microbe)
correlation_matrix <- matrix(NA, nrow = length(meta_vars), ncol = length(microbes),
                             dimnames = list(meta_vars, microbes))

# Fill the matrix
for (i in 1:nrow(data)) {
  row_name <- data$MetaVariable[i]
  col_name <- data$Microbe[i]
  if (!is.na(col_name)) {
    correlation_matrix[row_name, col_name] <- data$Correlation[i]
  }
}

# Replace NA values with 0
correlation_matrix[is.na(correlation_matrix)] <- 0

# Remove microbes (columns) with all 0 values
correlation_matrix <- correlation_matrix[, colSums(correlation_matrix != 0) > 0]

# Custom color palette
col_pal <- colorRampPalette(c("#F7B115", "white", "#D6101E"))(200)

# Draw correlation heatmap
corrplot(
  correlation_matrix,
  is.corr = FALSE,
  method = "circle",
  col = col_pal,
  tl.col = "black",
  tl.srt = 45,
  na.label = " ",
  cl.ratio = 0.2,
  cl.align.text = "c",
  cl.cex = 0.8,
  cl.pos = "r",
  cl.length = 5
)
