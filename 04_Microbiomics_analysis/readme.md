# Microbiomics analysis
We first analyzed the relative abundance of microbial communities at the phylum and genus levels across three subtypes. We then assessed the contribution of environmental and host-related factors to the variation in microbiome composition using adonis R² analysis. (adonis_analysis.R)

Next, we constructed a microbial co-occurrence network based on Spearman correlation coefficients. We applied the Louvain algorithm in Gephi (v0.10.1) to detect modular structures within the network, network complexity and key taxa were examined.(network_analysis.R)
Functional prediction of these key taxa was performed using PICRUSt2. The specific command used was:

picrust2_pipeline.py -s filtered_rep-seqs.fna -i abs_filtered.biom -o picrust2_results


Finally, we explored the correlation between the abundance of microbial genera and immune-related cell indices derived from routine blood tests. The results were visualized using dot plots. (immune_correlation.R)
