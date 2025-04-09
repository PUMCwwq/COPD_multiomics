# Metabolomics analysis
Differential metabolites analysis was conducted using the two-tailed Wilcoxon rank sum test. A metabolite was considered significantly altered if it exhibited FC > 1.2 or < 0.8 with FDR < 0.05. KEGG pathway enrichment analysis was done using clusterProfiler R package (v4.14.4). Orthogonal partial least-squares-discriminant analysis (OPLS-DA) was carried out to determine metabolites capable of distinguishing subtypes using ropls R package (v1.6.2). (Metabolomics_analysis.R)

The code (Metabolomics_analysis.R) includes 5 parts:
1. Differential metabolomics feature analysis, and heatmap display
2. PCA analysis
3. KEGG pathway enrichment analysis
4. Spearman’s correlation analysis
5. OPLS-DA analysis
