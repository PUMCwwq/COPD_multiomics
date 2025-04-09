# Radiomics analysis
Differential radiomics features among subtypes were defined as fold change (FC) > 1.2 or < 0.8 with a false discovery rate (FDR) < 0.05, based on the two-tailed Wilcoxon rank sum test. Principal components analysis (PCA) was achieved using the Factorextra R package (v1.0.7). (Radiomics_analysis.R)

The code (Radiomics_analysis.R) includes 3 parts:
1. Differential radiomics feature analysis, and heatmap display
2. PCA analysis
3. Heatmap of emphysema, calcification, and volume associated features at level of whole lung
