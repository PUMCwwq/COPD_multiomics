######################################################### Integration
#load data
iris1<-read.table("/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/SNF_num.txt" ,header = T,sep = "\t")
iris2<-read.table("/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/SNF_integer.txt" ,header = T,sep = "\t")
rownames(iris1) <- iris1[,1]
iris1<-as.matrix(iris1)
iris1=iris1[,2:675]
rownames(iris2) <- iris2[,1]
iris2<-as.matrix(iris2)
iris2=iris2[,2:181]
#Continuous variables were converted to numerical values
iris1_vec=as.numeric(as.vector(iris1))
iris1_1 <- matrix(iris1_vec, nrow = nrow(iris1))
rownames(iris1_1)<-rownames(iris1)
colnames(iris1_1)<-colnames(iris1)
iris1_1<-scale(iris1_1)
#Convert the category variable to an integer
iris2_vec=as.integer(as.vector(iris2))
iris2_1 <- matrix(iris2_vec, nrow = nrow(iris2))
rownames(iris2_1)<-rownames(iris2)
colnames(iris2_1)<-colnames(iris2)
#combination_FINALLY
iris3L <- list()
iris3L[1]=list(iris1_1)
iris3L[2]=list(iris2_1)

# Set the other parameters
K = 20 # number of neighbours
alpha = 0.5 # hyperparameter in affinityMatrix
T = 20 # number of iterations of SNF

#Calculate distance
##### gower Distance _3 numerical omics
gower_dist_3o <- daisy(iris1_1, metric = "gower")
gower_dist_3o<-as.matrix(gower_dist_3o)
write.csv(gower_dist_3o,file="/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/NUM_matrix.csv")
##### gower Distance_1 taxonomics
gower_dist_1o <- daisy(iris2_1, metric = "gower")
gower_dist_1o<-as.matrix(gower_dist_1o)
write.csv(gower_dist_1o,file="/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/SNP_matrix.csv")
##########dist_iris3L
dist_iris3 <- list()
dist_iris3[1]=list(gower_dist_3o)
dist_iris3[2]=list(gower_dist_1o)
# Construct the similarity graphs
affinity_iris3L = lapply(dist_iris3, function(x) affinityMatrix(x, K, alpha))
write.csv(affinity_iris3L,file="/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/affinity_matrix.csv")
W = SNF(affinity_iris3L, K, T)
write.csv(W,file="/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/W_matrix.csv")

###function loading
#SpectralClustering
SpectralClustering <- function(affinity, K, type=3) {
  d = rowSums(affinity)
  d[d == 0] = .Machine$double.eps
  D = diag(d)
  L = D - affinity
  if (type == 1) {
    NL = L
  } else if (type == 2) {
    Di = diag(1 / d)
    NL = Di %*% L
  } else if(type == 3) {
    Di = diag(1 / sqrt(d))
    NL = Di %*% L %*% Di
  }
  eig = eigen(NL)
  res = sort(abs(eig$values),index.return = TRUE)
  U = eig$vectors[,res$ix[1:K]]
  normalize <- function(x) x / sqrt(sum(x^2))
  if (type == 3) {
    U = t(apply(U,1,normalize))
  }
  eigDiscrete = discretisation(U)
  eigDiscrete = eigDiscrete$discrete
  labels = apply(eigDiscrete,1,which.max)
  
  clusterSize = apply(as.matrix(1:K), 1, function(x) sum(x==labels))
  if (min(clusterSize) <= 5) {
    labels = kmeans(eigDiscrete, K, 1000, 20)$cluster
  }
  
  return(labels)
}

#discretisation
discretisation <- function(eigenVectors) {
  
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  
  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = discretisationEigenVectorData(eigenVectors %*% R)
    
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
    
  }
  
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}

#discretisationEigenVectorData
discretisationEigenVectorData <- function(eigenVector) {
  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1
  
  return(Y)
  
}

###Spectral clustering results
clustering = SpectralClustering(W,3)
clustering <- as.matrix(clustering)
write.csv(clustering, file="/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/CLUSTER.csv")
write.csv(clustering, file="/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/CLUSTER1.csv")
sil <- silhouette(clustering, W)
fviz_silhouette(sil)


################################################# W matrix display
rm(list=ls())
dev.off()
dev.new()
BiocManager::install("pheatmap")
library(pheatmap)
rm(list=ls())
############################################### All omics affinity matrices are presented
#########################################################Annotated heat map
data<-read.delim("/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/Affinity_heatmap.txt",header = T,row.names = 1)
#Note that the number of columns is the same as the number of samples
Groups<-read.table("/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/groups_affinity.txt")
rownames(Groups)<-Groups$V1
Groups<-Groups[2:319,]
cluster<-Groups[,2]
cluster<-as.data.frame(cluster)
annotation_c<-cluster
annotation_c
rownames(annotation_c)<-rownames(Groups)


data1<- log10(data)  
data1<-scale(data1)
data1<- log10(data1+10)  
data1<-scale(data1)
data1<- log10(data1+10)  
data1<-scale(data1)
data1<- log10(data1+10)  
data1<-scale(data1)
data1<- log10(data1+10)  
data1<-scale(data1)
data1<- log10(data1+10)  
data1<-scale(data1)
#data1<- log10(data1+10)  
#data1<-scale(data1)
#data1<- log10(data1+10)  
#data1<-scale(data1)
#data1<- log10(data1+10)  
#data1<-scale(data1)
#data1<- log10(data1+10)  
#data1<-scale(data1)
#data1<- log10(data1+10)  
#data1<-scale(data1)
#data1<- log10(data1+10)  
#data1<-scale(data1)

p<-pheatmap(data1, cluster_rows = F, cluster_cols = F,     #Rows are (not) clustered, columns are (not) clustered
            show_rownames = F,       #Line names are not displayed
            show_colnames = F,      #The display specifies that angle_row="15", row names are rotated by 15 degrees, and column names are similar
            #annotation_col = annotation_c,  #To annotate a column is to group the columns
            scale= "none",             #The data were normalized by row
            border=F,               #Does not display a small lattice border
            color = colorRampPalette(c("#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#98CADD","#EAEFF6","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C"))(100))
#color = colorRampPalette(c("#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#547297","#EAEFF6","white","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C","#DC143C"))(100))
p 
dev.off()
dev.new()


