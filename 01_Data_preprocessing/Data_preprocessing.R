##############################################################################
# Script information                                                      
# Title: Data_preprocessing
# Date: 2025-3-1
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


##### 2.Continuous variables were normalized collectively.




##### 3.Categorical variables - SNPs associated with COPD or lung function were filtered.
# SNPs associated with COPD or lung function
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
