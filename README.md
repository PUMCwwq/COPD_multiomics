# Deep multi-omics profiling reveals three distinct molecular subtypes of chronic obstructive pulmonary disease in a unique biomass-exposed Chinese population
This repository contains the analysis pipeline and code for the study "Deep multi-omics profiling reveals three distinct molecular subtypes of chronic obstructive pulmonary disease in a unique biomass-exposed Chinese population"
# Data
Raw sequencing data have been deposited in the [Genome Sequence Archive] (https://ngdc.cncb.ac.cn/) under accession code PRJCA037919.
# Pipeline
# Data preprocessing
We recruited COPD patients from five medical centers in Guizhou, China. CT, nasal swab, plasma and blood were collected to acquire radiomics, microbiomics, metabolomics, and genomics data for further analysis. See details in 01_Data_preprocessing.
# Multi-omics integration and clustering
Multi-omics integration and clustering were performed. See details in 02_Integration_clustering.
# Radiomics analysis
CT imagines were processed using FACT Medical Imaging System (v1.2.0).  extracted include emphysema proportion, lung density, and fissure integrity. Radiomics characteristics across subtypes were profiled. See details in 03_Radiomics_analysis.
# Microbiomics (16S rRNA) analysis
Bacteria DNA was extracted and sequenced for 16S rRNA gene profiling. Data were processed with QIIME2 (v2022.2) and Silva database for taxonomic annotation. Analysis pipeline included decontamination via field blanks and DADA2-based ASV inference. See details in 04_Microbiomics_analysis.
# Metabolomics analysis
Untargeted metabolomics was detected via LC-MS/MS. Data processing utilized XCMS (R-based) and BiotreeDB v3.0 for annotation. Key analyses included differential metabolites identification, OPLS-DA, and KEGG pathway enrichment. See details in 05_Metabolomics_analysis.
# Genomics analysis
The genetic susceptibility across subtypes were calculated based on polygenic risk scores (PRS). See details in 06_Genomics_analysis.
