##############################################################################
# Script information                                                      
# Title: Continuous variables treatment
# Date: 2024-4-1
##############################################################################

##### 
# Input: Radiomics, metabolomics, and microbiomics data
iris1<-read.table("num.txt" ,header = T,sep = "\t")
rownames(iris1) <- iris1[,1]
iris1<-as.matrix(iris1)
# Convert continuous variables into numerical values
iris1_vec=as.numeric(as.vector(iris1))
iris1_vec<-scale(iris1_vec)
# Sort the variances and select the top 20% of variables
iris1_vec_variances <- apply(iris1_vec, 2, var)
top_20_percent_vars <- names(sort(iris1_vec_variances, decreasing = TRUE)[1:floor(length(iris1_vec_variances) * 0.2)])
iris1_1<- iris1[, top_20_percent_vars]
