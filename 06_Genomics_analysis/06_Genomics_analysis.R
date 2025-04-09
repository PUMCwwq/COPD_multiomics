##############################################################################
# Script information                                                      
# Title: Genomics_analysis
# Date: 2025-3-1
##############################################################################

##### 1.Polygenic Risk Score (PRS) associated analysis
library(tidyverse)
library(readxl)
library(ggsci)
library(viridis)
library(RColorBrewer)
# Input: PRS score and subtupe information
df_data <- read.delim("PRS.txt",header = T,row.names = 1)
df_data$CLUSTER_MATCH<-as.character(df_data$CLUSTER_MATCH)
# display the box plot
p1<-ggplot(df_data,aes(x=CLUSTER_MATCH, y=df_data$PRS.0.001, color=CLUSTER_MATCH, fill=CLUSTER_MATCH))+
  geom_boxplot()
p2 <- p1 +
  geom_jitter(width = 0.2, size = 0.9)+
  scale_color_viridis(end = 0.8, discrete = T)
p3 <- p2 +
  scale_fill_viridis(end = 0.8, discrete = T, alpha = 0.5)
p3
p3 + 
  scale_y_continuous(limits = c(-40,-30))+    # Agatston (0,1850), calcification (0,610), aera (14,1200)
  guides(color="none")+
  theme_clean()


##### 2. Regulatory network identification
### 2.1. GWAS-based subtype enrichment
# Chi-square test
library(knitr)
# Input: frequencies of COPD-associated risk alleles from SNPs associated with COPD or lung function previously identified
# odds ratios were calculated based on the actual and expected frequencies of COPD-associated risk alleles
df = read.table('risk_allele.txt', header = TRUE)
kable(df)
out <- data.frame()
for (i in 1:nrow(df)){
  t <- chisq.test(matrix(as.vector(t(df[i, 4:9])), nrow=2, ncol=3),correct = TRUE) 
  #t <- fisher.test(matrix(as.vector(t(df[i, 1:4])), nrow=2, ncol=2))
  d <- df[i, ]
  d$p.value <- t$p.value
  d$OR <- t$estimate[[1]]
  d$OR.lower95 <- t$conf.int[1]
  d$OR.upper95 <- t$conf.int[2]
  out <- rbind(out, d)
}
kable(out)
write.csv(out,file="out_risk_allele.csv")
# For each subtype, SNPs with OR > 1.2 and two-tailed Chi-square test P < 0.05 were defined as subtype-enriched SNPs (filtered in excel).

### 2.2. Functional genomics information
# the variant-to-gene (V2G) approach was employed via the online platform Open Target Genetics (https://platform.opentargets.org/?ref=blog.opentargets.org).
# The “variant” referred to the subtype-enriched SNPs identified in step 2.1.GWAS-based subtype enrichment.
# Genes with V2G score > 0.05 were selected to constitute regulatory networks of each subtype (performed in excel).

### 2.3. Gene-gene interaction
# Gene-gene interactions were considered if there was at least one evidence supporting interaction or co-expression between genes.
# Gene-gene interactions were acquired from STRING (https://cn.string-db.org/).
# The network plot was drawed according to genes in the regulatory networks manually.

### 2.4. Enrichment analysis of genes in regulatory networks
# Metascape was utilized (https://metascape.org/gp/index.html#/main/step1). Pathway database included Reactome, WikiPathway, and Gene Ontology.
# Genes contained in the regulatory network of each subtype were used as the input for the Metascape
# Enriched pathways and P values were acquired from Metascape. 

