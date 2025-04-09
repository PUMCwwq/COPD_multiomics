# Data acquirement and preprocessing
Patients with a prior diagnose of COPD were recruited from five medical centers of Guizhou province, southwest China. Nasal swabs, blood samples, chest inspiratory phase CT scans, and baseline information were collected.
#Baseline information
Baseline information including demographic data was collected from health assessment questionnaire. Laboratory tests were performed in the hospitals’ laboratory department.
#Radiomics
The chest inspiratory-phase CT scan was with a 1 mm slice thickness, and DICOM-format images were exported. Patient identity information was removed using de-identification software (Disk), and the images were processed using an imaging segmentation software (FACT Medical Imaging System, v1.2.0) to acquire features including emphysema proportion, lung fissure integrity, and lung density.
#Microbiomics
Nasal swabs were collected by the same trained operator and stored in a non-inactivated preservation solution at -80°C for microbiomics (16S rRNA). Sterile swabs moved around the collection room served as field blanks to recognize environmental and reagent contaminants, these blanks were processed concurrently with nasal swab samples under identical conditions, and the resulting ASV data were filtered accordingly to remove contamination.
#Metabolomics
The plasma was used for non-targeted metabolomics. LC-MS/MS analyses for polar metabolites were conducted using an UHPLC system (Vanquish, Thermo Fisher Scientific) with a Waters ACQUITY UPLC BEH Amide (2.1 mm × 50 mm, 1.7 μm) coupled to Orbitrap Exploris 120 mass spectrometer (Orbitrap MS, Thermo). Metabolite identification was accomplished using the R package and BiotreeDB (v3.0). Relative abundance of metabolites was used for further analysis by aligning with QC samples.
#Genomics acquirement
Genotyping was carried out by WeGene Clinical Laboratory (Shenzhen) utilizing the Illumina Infinium Chinese Genotyping Array BeadChip (Illumina WeGene V3 Arrays, covering approximately 700,000 variants).
#Preprocessing
Radiomics, metabolomics, and microbiomics data, represented as continuous variables, were normalized collectively. The top 20% highly variable features were selected for further integration and clustering. Genomics data, treated as categorical variables, were processed separately. (Data_preprocessing.R)

The code (Data_preprocessing.R) includes 3 parts:
1. KNN method for imputing missing values of radiomics.
2. Continuous variables were normalized collectively.
3. Categorical variables - SNPs associated with COPD or lung function were filtered.
