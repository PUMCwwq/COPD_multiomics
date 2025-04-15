##############################################################################
# Script information                                                      
# Title: imputing_missing_values
# Date: 2024-4-1
##############################################################################

##### KNN method for imputing missing values.
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

