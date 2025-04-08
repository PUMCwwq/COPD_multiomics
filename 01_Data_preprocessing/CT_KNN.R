########## COPD multiomics ##########
##### CT KNN
rm(list=ls())
################################# The KNN method was used to impute missing CT features
BiocManager::install("TTR")
install.packages("/NASdata_2/wuwq/COPD/packages/xts_0.14.0.tar.gz",repos = NULL)
install.packages("/NASdata_2/wuwq/COPD/packages/quantmod_0.4.26.tar.gz",repos = NULL)
install.packages("/NASdata_2/wuwq/COPD/packages/TTR_0.24.4.tar.gz",repos = NULL)
install.packages("/NASdata_2/wuwq/COPD/packages/DMwR_0.4.1.tar.gz",repos = NULL)
library(xts)
library(TTR)
library(quantmod)
library(DMwR)
#####Import data sets
data <- read.delim("/NASdata_2/wuwq/COPD/OMIC_ALL/CT_159.txt",quote = "", 
                   row.names = 1, 
                   stringsAsFactors = FALSE)

filled_dataset <- knnImputation(data, k = 10)

print(filled_dataset)
write.csv(filled_dataset,file="/NASdata_2/wuwq/COPD/OMIC_ALL/CT_knn.csv")

