##############################################################################
# Script information                                                      
# Title: Metabolomics_analysis
# Date: 2025-3-1
##############################################################################

##### 1.Differential metabolomics feature analysis, and heatmap display
### Wilcoxon rank sum test
# Input: metabolomics features
# This example is the comparision between CS and SS&RS subtypes.
df <- read.table("metabolomics.txt", 
                 header=T, 
                 row.names=1, 
                 sep="\t") 

wilcox_data <- data.frame(meta=colnames(df)[ 6:ncol(df)]) 
for(i in 6:ncol(df)){ 
  print(i) 
  wilcox_data[i- 5, 2] <- wilcox.test(df[,i] ~ CS,data = df,exact = FALSE)[[ "p.value"]] } 
wilcox_data$fdr <- p.adjust(wilcox_data$V2,method = "fdr") 
#Calculate fdr # Calculate logFC
data_A <- df[c(which(df$C3 == "1")),] 
data_B <- df[c(which(df$C3 == "2")),] 
wilcox_data$A_AVE <- colMeans(data_A[ 6:ncol(data_A)])
wilcox_data$B_AVE <- colMeans(data_B[ 6:ncol(data_B)])
wilcox_data$minus <- colMeans(data_A[ 6:ncol(data_A)]) - colMeans(data_B[ 6:ncol(data_B)])
wilcox_data$fc <- colMeans(data_A[ 6:ncol(data_A)]) / colMeans(data_B[ 6:ncol(data_B)])
write.csv(wilcox_data,file="wilcox_out_CS_metabolomics.csv")
# The subsequent filtering of FC > 1.2 or FC < 0.8 is done in Excel.
# The code for comparisions between SS and RS&CS, and the comparision between RS and SS&CS were similar to the comparison above.

### Differential metabolomics feature heatmap display
install.packages("pheatmap")
library(pheatmap)
# Input: Differential metabolomics features
data<-read.delim("heatmap_for_metabo.txt",header = T,row.names = 1)
#Note that the number of columns is the same as the number of samples
Groups<-read.table("cluster.txt",header = T,row.names = 1)
cluster<-as.data.frame(Groups)
data1<- log10(data+1) 
data1<-scale(data1)
p<-pheatmap(data1, cluster_rows = T, cluster_cols = F,   
            show_rownames = F,     
            show_colnames = F,      
            scale= "row",            
            border=F,             
            legend = F,
            treeheight_row = 0,
            color = colorRampPalette(c("#1D6A3B","#1D6A3B","#1D6A3B","#1D6A3B","#1D6A3B","#1D6A3B","#1D6A3B","#1D6A3B","#599455","white","#6c2a1a","#6c2a1a","#6c2a1a","#6c2a1a","#6c2a1a","#6c2a1a","#6c2a1a","#6c2a1a","#6c2a1a"))(100)) 
p 


##### 2.PCA analysis
library(ropls)
# PCA was performed based on metabolomics; SS and RS were combined into one group when calculating the 95% confidence interval
data<-read.table("metabo_for_PCA.txt",header = TRUE,sep = "\t")
data$CLUSTER_MATCH<-as.character(data$CLUSTER_MATCH)
head(data)
#prcomp
com1 <- prcomp(data, center = TRUE,scale. = TRUE)
summary(com1)
#PCA
df1<-com1$x
head(df1)
df1<-data.frame(df1,data$CLUSTER_MATCH)
head(df1)
#The variance contribution rate of principal components was extracted and the axis title was generated.
summ<-summary(com1)
xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
p2<-ggplot(data = df1,aes(x=PC1,y=PC2,color=data.CLUSTER_MATCH))+
  stat_ellipse(aes(fill=data.CLUSTER_MATCH),
               type = "norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)
p2+scale_fill_manual(values = c("purple","orange","blue","pink"))+
  scale_colour_manual(values = c("purple","orange","blue","pink"))
write.csv(p2$data,file = "PCA_metabo_SS&RS_CS.cs


##### 3.KEGG pathway enrichment analysis
# load the hsa_kegg annotation file acquired from https://www.kegg.jp/kegg 
pathway_compounds <- readRDS("hsa_kegg_annotation.rds")
# Background metabolites  (Input: all untargeted metabolites detected in this study)
allkeggid<- read.delim("all_metabolite.txt",header = T)
total <- right_join(pathway_compounds, allkeggid, by = "ID")
total$id <- total$ID
total<-total[,2:3]
colnames(total)<-c("pathway","ID")

# Enrichment results  Input: differential metabolites identified for each subtype previously
# Elevated and attenuated metabolites were processed separately

# CS_Elevated_SS&RS_attenuated 
diffkeggid <-read.delim("kegg_CS_up.txt",header = T)
x_CS <- clusterProfiler::enricher(gene = diffkeggid$ID,TERM2GENE = total,minGSSize = 1,pvalueCutoff = 1,qvalueCutoff = 1)
write.csv(x_CS@result,file = "kegg_out_CS_up.csv")

# SS&RS_Elevated_CS_attenuated 
diffkeggid <-read.delim("kegg_SS&RS_up.txt",header = T)
x_SS&RS <- clusterProfiler::enricher(gene = diffkeggid$ID,TERM2GENE = total,minGSSize = 1,pvalueCutoff = 1,qvalueCutoff = 1)
write.csv(x_SS&RS@result,file = "kegg_out_SS&RS_up.csv")


##### 4.Spearman’s correlation analysis
### correlation between lipids to get lipid patterns, and correlation between lipids and baseline information (including clinical indicators)
library(corrplot)
library(psych)

# Impute missing values of baseline information using the KNN method (Input: baseline information and differential lipids identified)
# load the baseline information
data <- read.delim("baseline_lipids.txt",quote = "", 
                   row.names = 1, 
                   stringsAsFactors = FALSE)
# KNN
set.seed(123)
filled_dataset <- knnImputation(data, k = 10)

#Standardization does not affect the calculation of correlation coefficients, but it allows the data to follow a structure with a mean of 0 and a standard deviation of 1, ensuring homoscedasticity.
sam_metabo <- scale(filled_dataset)

# Symmetric matrix of correlation coefficients between all variables
corr_matrix <- corr.test(sam_metabo, method = 'spearman')
corr_matrix$r    # matrix of correlation rho
corr_matrix$p    # matrix of P
# download the rho and P
write.table(corr_matrix$r, 'corr_SAM_sigmetabo_R.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(corr_matrix$p, 'corr_SAM_sigmetabo_P.txt', sep = '\t', col.names = NA, quote = FALSE)
# Based on the previous step, filter for p < 0.05
PP<-corr_matrix$p
PP[PP >= 0.05] <- -1
PP[PP < 0.05 & PP>=0 ] <- 1
PP[PP == -1] <- 0
PP_filter <- RR * PP
write.table(PP_filter, 'corr_SAM_sigmetabo_P_0.05.txt', sep = '\t', col.names = NA, quote = FALSE)
# Based on the previous step, further filter for p < 0.05 and |R| >= 0.2.
RR[abs(RR) < 0.2 ] <- 0
PP_filter_RR_filter <- RR * PP
write.table(PP_filter_RR_filter, 'corr_SAM_3OMICS_p_0.05_r_0.02_new.txt', sep = '\t', col.names = NA, quote = FALSE)

### display the lipid pattern using heatmap
# Input: the correlation rho matrix of differential lipids from correlation analysis above
data<-read.delim("correlation_lipids.txt",header = T,row.names = 1)
Groups<-read.table("group_cor.txt",header = T,row.names = 1)
cluster<-as.data.frame(Groups)
# 9 color blocks in segments ranging from rho = -0.6 to NA to 1.
p<-pheatmap(data, cluster_rows = F, cluster_cols = F,   
            show_rownames = T,      
            show_colnames = F,     
            annotation_col = annotation_c,  
            scale= "none",  
            border=F,   
            treeheight_row = 0,
            legend = T,
            color = colorRampPalette(c("#61AACF","#98CADD","#EAEFF6","white","#F9EFEF","#E9C6C6","#DA9599","#C16E71","#C16E71"))(9))


##### 5.OPLS-DA analysis
set.seed(123)
# Input: differential lipids identified across subtypes
meta_pls<-read.table("lipid_OPLSDA.txt",header = TRUE,row.names = 1,sep = "\t")
# Input: the subtype information which differential lipids belong to
meta_hint<-read.table("OPLSDA_GROUP.txt",header = TRUE,row.names = 1,sep = "\t")
CLUSTER <- meta_hint$CLUSTER
CLUSTER <-as.factor(CLUSTER)
# OPLS-DA
meta_pls.plsda <- opls(meta_pls, CLUSTER, predI = 1,scale = "standard")
# VIP
vipVn <- meta_pls.plsda@vipVn  # getVipVn()
vipVn_select <- vipVn[vipVn > 1] 
write.csv(vipVn_select, VIP_score.csv")


##### 6. Heatmap of subtype elevated lipids









