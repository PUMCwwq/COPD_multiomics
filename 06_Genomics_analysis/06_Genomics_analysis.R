#########################################################   SNP
####################################### Chi-square test
library(knitr)
setwd("/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/SNP")
df = read.table('again_risk_allele_recessive.txt', header = TRUE)
kable(df)

out <- data.frame()
for (i in 1:nrow(df)){
  t <- chisq.test(matrix(as.vector(t(df[i, 4:9])), nrow=2, ncol=3),correct = TRUE) 
  #t <- fisher.test(matrix(as.vector(t(df[i, 1:4])), nrow=2, ncol=2))
  d <- df[i, ]
  d$p.value <- t$p.value
  d$OR <- t$estimate[[1]]
  d$OR.lower95 <- t$conf.int[1]
  d$OR.upper95 <- t$conf.int[2]
  out <- rbind(out, d)
}

kable(out)
write.csv(out,file="out_again_risk_allele_recessive.csv")


######################################### PRS
rm(list=ls())
setwd("/NASdata_2/wuwq/COPD/OMIC_ALL_SUM/SNP")
library(tidyverse)
library(readxl)
library(ggsci)
library(viridis)
library(RColorBrewer)
df_data <- read.delim("violin_PRS.txt",header = T,row.names = 1)
df_data$CLUSTER_MATCH<-as.character(df_data$CLUSTER_MATCH)
###
theme_clean <- function(){
  theme_bw() %+replace%    
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank()
    )
}


########Draw a box plot
p1<-ggplot(df_data,aes(x=CLUSTER_MATCH, y=df_data$PRS.0.001, color=CLUSTER_MATCH, fill=CLUSTER_MATCH))+
  geom_boxplot()
p2 <- p1 +
  geom_jitter(width = 0.2, size = 0.9)+
  scale_color_viridis(end = 0.8, discrete = T)
p3 <- p2 +
  scale_fill_viridis(end = 0.8, discrete = T, alpha = 0.5)
p3
p3 + 
  scale_y_continuous(limits = c(-40,-30))+    # Agatston (0,1850), calcification (0,610), aera (14,1200)
  guides(color="none")+
  theme_clean()



