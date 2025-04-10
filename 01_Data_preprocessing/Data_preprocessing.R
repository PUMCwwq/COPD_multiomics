##############################################################################
# Script information                                                      
# Title: Data_preprocessing
# Date: 2024-4-1
##############################################################################

##### 1.KNN method for imputing missing values.
install.packages("xts_0.14.0.tar.gz",repos = NULL)
install.packages("quantmod_0.4.26.tar.gz",repos = NULL)
install.packages("TTR_0.24.4.tar.gz",repos = NULL)
install.packages("DMwR_0.4.1.tar.gz",repos = NULL)
library(xts)
library(TTR)
library(quantmod)
library(DMwR)

#Import data sets
data <- read.delim("CT_159.txt",quote = "", 
                   row.names = 1, 
                   stringsAsFactors = FALSE)
filled_dataset <- knnImputation(data, k = 10)

print(filled_dataset)
write.csv(filled_dataset,file="CT_knn.csv")


##### 2.Continuous variables were scaled collectively.
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


##### 3.Categorical variables - SNPs associated with COPD or lung function were filtered.
# Input: SNPs associated with COPD or lung function
# Clumping
# Snp&p 2 columns, SNP format x:position:ref:alt
df<-read.table("clump/METAsig.txt" ,header = T,sep = "\t")
for (i in c(1:22, 'X')){
  tmp<-df[df$CHR==i,]
  bfile = paste0("1000G_hg38/ALL.chr",i)
  bim<-read.table(paste0("1000G_hg38/ALL.chr",i,'.bim'))
  plink_bin = "plink" 
  # The sites to clump are written to the following file
  write.table(tmp,paste0('clump/','toclump.',i,'.tsv'),
              sep='\t',quote=F,col.names=T,row.names=F)
  fun2 <- paste0(
    shQuote(plink_bin),
    " --noweb",
    " --bfile ", shQuote(bfile),
    " --num_threads 8",
    paste0(' --clump clump/','toclump.',i,'.tsv'),#The site to clump
    " --clump-kb 500", 
    " --clump-p1 1",
    " --clump-r2 0.2",
    paste0(" --out clump/","clumped",i)#clump Results
  )
  system(fun2)
}
# The results of clump on 22 chromosomes were read
f=list.files('clump/')
f<-f[grep('.clumped',f)]
clumped<-read.table(paste0('clump/',f[1]),
                    header = T)
for(i in f){
  file=read.table(paste0('/clump/',i),
                  header = T)
  clumped<-rbind(clumped,file)
}
write.csv(clumped, file="clumped.csv")
