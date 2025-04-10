##############################################################################
# Script information                                                      
# Title: Radiomics_analysis
# Date: 2024-8-1
##############################################################################

##### 1. Differential radiomics feature analysis, and heatmap display
### Wilcoxon rank sum test
# Input: radiomics features. 
# This example is the comparision between CS and SS&RS subtypes.
df <- read.table("CT.txt", 
                 header=T, 
                 row.names=1, 
                 sep="\t") 
wilcox_data <- data.frame(meta=colnames(df)[ 6:ncol(df)]) 
for(i in 6:ncol(df)){ 
  print(i) 
  wilcox_data[i- 5, 2] <- wilcox.test(df[,i] ~ CS,data = df,exact = FALSE)[[ "p.value"]] } 
wilcox_data$fdr <- p.adjust(wilcox_data$V2,method = "fdr") 
# calculate fdr and logFC  (rename CS as "1", rename SS&RS as "2")
data_A <- df[c(which(df$CS == "1")),] 
data_B <- df[c(which(df$CS == "2")),] 
wilcox_data$A_AVE <- colMeans(data_A[ 6:ncol(data_A)])
wilcox_data$B_AVE <- colMeans(data_B[ 6:ncol(data_B)])
wilcox_data$minus <- colMeans(data_A[ 6:ncol(data_A)]) - colMeans(data_B[ 6:ncol(data_B)])
wilcox_data$fc <- colMeans(data_A[ 6:ncol(data_A)]) / colMeans(data_B[ 6:ncol(data_B)])
write.csv(wilcox_data,file="wilcox_out_CS.csv")
# The subsequent filtering of FC > 1.2 or FC < 0.8 is done in Excel.
# The code for comparisions between SS and RS&CS, and the comparision between RS and SS&CS were similar to the comparison above.

### Differential radiomics feature heatmap display
install.packages("pheatmap")
library(pheatmap)
# Input: differential radiomics features of three subtypes
data<-read.delim("differential_CT.txt",header = T,row.names = 1)
#Note that the number of columns is the same as the number of samples
Groups<-read.table("cluster.txt",header = T,row.names = 1)
cluster<-Groups[,1]
cluster<-as.data.frame(cluster)
data1<-scale(data)
p<-pheatmap(data1, cluster_rows = T, cluster_cols = F,    
            show_rownames = F,       
            show_colnames = F,      
            scale= "row",            
            border=F,              
            legend = F,
            treeheight_row = 0,
            color = colorRampPalette(c("#F9CC76","#F9CC76","#F9CC76","#F9CC76","#F9CC76","#F9CC76","#F9CC76","#F9CC76","white","#7388C1","#7388C1","#7388C1","#7388C1","#7388C1","#7388C1","#7388C1","#7388C1"))(100)) 
p 


##### 2. PCA analysis
library(ropls)
# Input: Radiomics features
data<-read.table("CT_for_PCA.txt",header = TRUE,sep = "\t")
data$CLUSTER_MATCH<-as.character(data$CLUSTER_MATCH)
head(data)
#prcomp
com1 <- prcomp(data, center = TRUE,scale. = TRUE)
summary(com1)
#PCA
#PC score；
df1<-com1$x
head(df1)
df1<-data.frame(df1,data$CLUSTER_MATCH)
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
write.csv(p2$data,file = "PCA_CT.csv")


###### 3.Heatmap of emphysema, calcification, and volume associated features at level of whole lung
install.packages("pheatmap")
library(pheatmap)
# Input: Differential emphysema, calcification, and volume associated features across subtypes
data<-read.delim("Lung_index_heatmap_v2.txt",header = T,row.names = 1)
Groups<-read.table("cluster.txt",header = T,row.names = 1)
cluster<-as.data.frame(Groups)
data1<-log10(data)
data1<-scale(data1)
data1<-t(data1)
p<-pheatmap(data1, cluster_rows = F, cluster_cols = F,     
            show_rownames = T,       
            show_colnames = F,      
            annotation_col = annotation_c,  
            scale= "row",             
            border=F,               
            #legend = F,
            treeheight_row = 5,
            color = colorRampPalette(c("#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","white","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779"))(100))
p 
