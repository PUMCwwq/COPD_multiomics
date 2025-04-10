##############################################################################
# Script information                                                      
# Title: network_analysis
# Author: Minyu Zhou
# Date: 2024-11-22
# Description: None
##############################################################################

# Load genus-level data (species × samples format)
genus_data <- read.table("genus.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

# Convert to numeric matrix
genus_data <- as.matrix(genus_data)
mode(genus_data) <- "numeric"

# === Compute Spearman correlation ===
library(Hmisc)
spearman_res <- rcorr(t(genus_data), type = "spearman")
cor_spearman <- spearman_res$r
pvals_spearman <- spearman_res$P

# Filter by threshold: |correlation| > 0.3 and p < 0.05
cor_spearman[abs(cor_spearman) <= 0.3 | pvals_spearman >= 0.05] <- 0
diag(cor_spearman) <- 0  # Remove self-correlation

# Save final adjacency matrix
write.table(cor_spearman, "spearman_network.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

# === Create Gephi-compatible edge list ===
library(igraph)

adj_matrix <- read.table("spearman_network.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
adj_matrix <- as.matrix(adj_matrix)
mode(adj_matrix) <- "numeric"

# Create igraph network
net <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
E(net)$weight <- abs(E(net)$weight)

# Extract and export edge list
edges <- as.data.frame(as_edgelist(net))
colnames(edges) <- c("Source", "Target")
edges$Weight <- E(net)$weight
write.csv(edges, "gephi_edges.csv", row.names = FALSE, quote = FALSE)

# === Compute and export node attributes ===
nodes <- data.frame(ID = V(net)$name)
nodes$Degree <- degree(net)
nodes$Betweenness <- betweenness(net, weights = E(net)$weight)
nodes$Closeness <- closeness(net, weights = E(net)$weight)
write.csv(nodes, "gephi_nodes.csv", row.names = FALSE, quote = FALSE)
