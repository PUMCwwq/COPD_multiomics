##################################################### CT
rm(list=ls())
install.packages("pheatmap")
library(pheatmap)
#################################################### Subtype overall heatmap Heatmap
setwd("/NASdata_2/wuwq/COPD/OMIC_ALL_SUM")
data<-read.delim("heatmap_CT_newrank_2.txt",header = T,row.names = 1)
#Note that the number of columns is the same as the number of samples
Groups<-read.table("groups_new_rank_2.txt",header = T,row.names = 1)
cluster<-Groups[,1]
cluster<-as.data.frame(cluster)
annotation_c<-cluster
annotation_c
rownames(annotation_c)<-rownames(Groups)

data1<-scale(data)

p<-pheatmap(data1, cluster_rows = T, cluster_cols = F,    
            show_rownames = F,       
            show_colnames = F,      
            #annotation_col = annotation_c,  
            scale= "row",            
            border=F,              
            legend = F,
            treeheight_row = 0,
            color = colorRampPalette(c("#F9CC76","#F9CC76","#F9CC76","#F9CC76","#F9CC76","#F9CC76","#F9CC76","#F9CC76","white","#7388C1","#7388C1","#7388C1","#7388C1","#7388C1","#7388C1","#7388C1","#7388C1"))(100)) 
p 


########################################################
######################################################## 
rm(list=ls())
dev.off()
dev.new()
library(ropls)

data<-read.table("/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/CT/CT_155_PCA.txt",header = TRUE,sep = "\t")
data$CLUSTER_MATCH<-as.character(data$CLUSTER_MATCH)
head(data)
# scale. = TRUE   indicates normalization of data before analysis
#prcomp
com1 <- prcomp(data[,6:160], center = TRUE,scale. = TRUE)
summary(com1)
#PCA
#PC score；
df1<-com1$x
head(df1)
#Column 5 of the iris dataset was merged
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
write.csv(p2$data,file = "/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/CT/PCA_CT_155_C1_C2_C3.csv")



#######
######Heat map of the LAA series -Whole lung goes to see, except for the differential LAA series and the calcification series, plus lung volume
rm(list=ls())
install.packages("pheatmap")
library(pheatmap)
########################### 
setwd("/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/CT")
data<-read.delim("Lung_index_heatmap_v2.txt",header = T,row.names = 1)
#Note that the number of columns is the same as the number of samples
Groups<-read.table("groups_rank_for_CT.txt",header = T,row.names = 1)
cluster<-Groups[,1]
cluster<-as.data.frame(cluster)
annotation_c<-cluster
annotation_c
rownames(annotation_c)<-rownames(Groups)


######原来scale指的是默认对列标准化，此处转置是为了对行scale以区分出列的差别
data1<-log10(data)
data1<-scale(data1)
data1<-t(data1)
#write.csv(data1,file="scaled_whole_CT_index.csv")
p<-pheatmap(data1, cluster_rows = F, cluster_cols = F,     #行（不）聚类，列(不)聚类
            show_rownames = T,       #不显示行名
            show_colnames = F,      #显示列明 angle_row="15"，行名旋转15度，列明相似
            annotation_col = annotation_c,  #对列进行注释即对列进行分组
            scale= "row",             #将数据按行进行标准化
            border=F,               #不显示小格子边框
            #legend = F,
            treeheight_row = 5,
            #color = colorRampPalette(c("#D2A6C7","#D2A6C7","#D2A6C7","#D2A6C7","#D2A6C7","#D2A6C7","#D2A6C7","#D2A6C7","white","#83C75D","#83C75D","#83C75D","#83C75D","#83C75D","#83C75D","#83C75D","#83C75D"))(100)) 
            #color = colorRampPalette(c("#349979","#349979","#349979","#349979","#a9c287","#a9c287","#a9c287","#a9c287","#cce982","#cce982","#cce982","#cce982","white","#F9CC76","#F9CC76","#F9CC76","#F9CC76","#dda717","#dda717","#dda717","#dda717","#ea6800","#ea6800","#ea6800","#ea6800"))(50)) 
            #color = colorRampPalette(c("#349979","#349979","#349979","#349979","#a9c287","#a9c287","#a9c287","#a9c287","#cce982","#cce982","#cce982","#cce982","white","#ffae49","#ffae49","#ffae49","#ffae49","#dda717","#dda717","#dda717","#dda717","#ea6800","#ea6800","#ea6800","#ea6800"))(50)) 
            #color = colorRampPalette(c("#349979","#349979","#349979","#349979","#cce982","#cce982","#cce982","#cce982","white","#fec739","#fec739","#fec739","#fec739","#ea6800","#ea6800","#ea6800","#ea6800"))(50)) 
            #color = colorRampPalette(c("#3d71a3","#3d71a3","#3d71a3","#3d71a3","#64c8c8","#64c8c8","#64c8c8","#64c8c8","white","#fec739","#fec739","#fec739","#fec739","#ea6800","#ea6800","#ea6800","#ea6800"))(50)) 
            #color = colorRampPalette(c("#349979","#349979","#349979","#349979","#cce982","#cce982","#cce982","#cce982","white","#3d71a3","#3d71a3","#3d71a3","#3d71a3","#3e266b","#3e266b","#3e266b","#3e266b"))(50)) 
            #color = colorRampPalette(c("#3d71a3","#3d71a3","#3d71a3","#3d71a3","#3d71a3","#3d71a3","#3d71a3","#3d71a3","white","#cce982","#cce982","#cce982","#cce982","#cce982","#cce982","#cce982","#cce982"))(50)) 
            #color = colorRampPalette(c("#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","white","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779"))(50))
            #color = colorRampPalette(c("#cce0ea","#cce0ea","#cce0ea","#cce0ea","#cce0ea","#cce0ea","#cce0ea","#cce0ea","white","#f5d7c0","#f5d7c0","#f5d7c0","#f5d7c0","#f5d7c0","#f5d7c0","#f5d7c0","#f5d7c0"))(50))
            #color = colorRampPalette(c("#feb779","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779","white","#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","#7570b3"))(50))
            #color = colorRampPalette(c("#98cadd","#98cadd","#98cadd","#98cadd","#98cadd","#98cadd","#98cadd","#98cadd","white","#da9599","#da9599","#da9599","#da9599","#da9599","#da9599","#da9599","#da9599"))(100))
            #color = colorRampPalette(c("#61aacf","#61aacf","#61aacf","#61aacf","#61aacf","#61aacf","#61aacf","#61aacf","white","#c16e71","#c16e71","#c16e71","#c16e71","#c16e71","#c16e71","#c16e71","#c16e71"))(100))  
            #v2#color = colorRampPalette(c("#61aacf","#61aacf","#61aacf","#61aacf","#61aacf","#61aacf","#61aacf","#61aacf","white","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779"))(100))  
            color = colorRampPalette(c("#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","#7570b3","white","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779","#feb779"))(100))
p 