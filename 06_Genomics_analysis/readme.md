# Genomics analysis
PRS analysis was accomplished using pruning and thresholding (“P+T”) method. The variant effect alleles and their effect size to COPD were accessed from the results of Global Biobank Meta-analysis Initiative (GBMI) multi-ancestry meta-analysis. For the regulatory network identification, GWAS-based subtype enrichment was carried out first. Then, Genes with V2G score > 0.05 were selected to constitute regulatory networks of each subtype. Gene-gene interactions were considered if there was at least one evidence supporting interaction or co-expression between genes.

The code (Metabolomics_analysis.R) includes 2 parts:
1. Polygenic Risk Score (PRS) associated analysis
2. Regulatory network identification
