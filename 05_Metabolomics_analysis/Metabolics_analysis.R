rm(list=ls())
######################################### INPUT_NUM  rank sum test
setwd("/NASdata_2/wuwq/COPD/OMIC_ALL_SUM")
df <- read.table("NUM_test.txt", 
                 header=T, 
                 row.names=1, 
                 sep="\t") 

wilcox_data <- data.frame(meta=colnames(df)[ 6:ncol(df)]) 
for(i in 6:ncol(df)){ 
  print(i) 
  wilcox_data[i- 5, 2] <- wilcox.test(df[,i] ~ C3,data = df,exact = FALSE)[[ "p.value"]] } 
wilcox_data$fdr <- p.adjust(wilcox_data$V2,method = "fdr") 
#Calculate fdr # Calculate logFC
data_A <- df[c(which(df$C3 == "1")),] 
data_B <- df[c(which(df$C3 == "2")),] 
wilcox_data$A_AVE <- colMeans(data_A[ 6:ncol(data_A)])
wilcox_data$B_AVE <- colMeans(data_B[ 6:ncol(data_B)])
wilcox_data$minus <- colMeans(data_A[ 6:ncol(data_A)]) - colMeans(data_B[ 6:ncol(data_B)])
wilcox_data$fc <- colMeans(data_A[ 6:ncol(data_A)]) / colMeans(data_B[ 6:ncol(data_B)])
write.csv(wilcox_data,file="wilcox_out_C3_NUM.csv")

##########  heatmap
install.packages("pheatmap")
library(pheatmap)
setwd("/NASdata_2/wuwq/COPD/OMIC_ALL_SUM")
data<-read.delim("heatmap_metabo_newrank_2.txt",header = T,row.names = 1)
#Note that the number of columns is the same as the number of samples
Groups<-read.table("groups_new_rank_2.txt",header = T,row.names = 1)
cluster<-Groups[,1]
cluster<-as.data.frame(cluster)
annotation_c<-cluster
annotation_c
rownames(annotation_c)<-rownames(Groups)

data1<- log10(data+1) 
data1<-scale(data1)

p<-pheatmap(data1, cluster_rows = T, cluster_cols = F,   
            show_rownames = F,     
            show_colnames = F,      
            #annotation_col = annotation_c,  
            scale= "row",            
            border=F,             
            legend = F,
            treeheight_row = 0,
            color = colorRampPalette(c("#1D6A3B","#1D6A3B","#1D6A3B","#1D6A3B","#1D6A3B","#1D6A3B","#1D6A3B","#1D6A3B","#599455","white","#6c2a1a","#6c2a1a","#6c2a1a","#6c2a1a","#6c2a1a","#6c2a1a","#6c2a1a","#6c2a1a","#6c2a1a"))(100)) 
p 



### PCA      
rm(list=ls())
dev.off()
dev.new()
library(ropls)

### PCA was performed by combining SS/RS and CS
data<-read.table("/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/Metabo/metabo_input_408_C12C3_PCA.txt",header = TRUE,sep = "\t")
data$CLUSTER_MATCH<-as.character(data$CLUSTER_MATCH)
head(data)
# scale. = TRUE
#prcomp函数 
com1 <- prcomp(data[,3:410], center = TRUE,scale. = TRUE)
summary(com1)
#PCA
#PC score；
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
write.csv(p2$data,file = "/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/Metabo/PCA_metabo_408_C12&C3.cs

rm(list=ls())
############################### Mann–Whitney U-test
setwd("/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/Metabo")

df <- read.table("metabo_all_ANOVA.txt", 
                 header=T, 
                 row.names=1, 
                 sep="\t")

wilcox_data <- data.frame(meta=colnames(df)[ 5:ncol(df)]) 
for(i in 5:ncol(df)){ 
  print(i) 
  wilcox_data[i- 4, 2] <- wilcox.test(df[,i] ~ C3,data = df,exact = FALSE)[[ "p.value"]] } 
wilcox_data$fdr <- p.adjust(wilcox_data$V2,method = "fdr") 
#计算fdr #计算logFC 
data_A <- df[c(which(df$C3 == "1")),] 
data_B <- df[c(which(df$C3 == "2")),] 
wilcox_data$A_AVE <- colMeans(data_A[ 5:ncol(data_A)])
wilcox_data$B_AVE <- colMeans(data_B[ 5:ncol(data_B)])
wilcox_data$minus <- colMeans(data_A[ 5:ncol(data_A)]) - colMeans(data_B[ 5:ncol(data_B)])
wilcox_data$fc <- colMeans(data_A[ 5:ncol(data_A)]) / colMeans(data_B[ 5:ncol(data_B)])
write.csv(wilcox_data,file="wilcox_out_C3_metabo_all.csv")

