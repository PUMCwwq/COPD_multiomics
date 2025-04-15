##############################################################################
# Script information                                                      
# Title: Continuous variables treatment
# Date: 2024-4-1
##############################################################################

#####clumping
rm(list=ls())
#########################################
### Snp&p 2 columns, SNP format x:position:ref:alt
df<-read.table("/NASdata_2/wuwq/COPD/OMIC_ALL/clump/METAsig.txt" ,header = T,sep = "\t")
for (i in c(1:22, 'X')){
  tmp<-df[df$CHR==i,]
  bfile = paste0("/NASdata/wangjx/data_resource/1000G_hg38/ALL.chr",i)
  bim<-read.table(paste0("/NASdata/wangjx/data_resource/1000G_hg38/ALL.chr",i,'.bim'))
  plink_bin = "/NASdata/wangjx/tools/plink_1.9/plink" #/NASdata_2/wuwq/COPD/OMIC_ALL/plink-1.07-x86_64/plink #/NASdata/wangjx/tools/plink_1.9/plink
  # The sites to clump are written to the following file
  write.table(tmp,paste0('/NASdata_2/wuwq/COPD/OMIC_ALL/clump/','toclump.',i,'.tsv'),
              sep='\t',quote=F,col.names=T,row.names=F)
  fun2 <- paste0(
    shQuote(plink_bin),
    " --noweb",
    " --bfile ", shQuote(bfile),
    " --num_threads 8",
    paste0(' --clump /NASdata_2/wuwq/COPD/OMIC_ALL/clump/','toclump.',i,'.tsv'),#The site to clump
    " --clump-kb 500", 
    " --clump-p1 1",
    " --clump-r2 0.2",
    paste0(" --out /NASdata_2/wuwq/COPD/OMIC_ALL/clump/","clumped",i)#clump Results
  )
  system(fun2)
}
##The results of clump on 22 chromosomes were read
f=list.files('/NASdata_2/wuwq/COPD/OMIC_ALL/clump/')
f<-f[grep('.clumped',f)]
clumped<-read.table(paste0('/NASdata_2/wuwq/COPD/OMIC_ALL/clump/',f[1]),
                    header = T)
for(i in f){
  file=read.table(paste0('/NASdata_2/wuwq/COPD/OMIC_ALL/clump/',i),
                  header = T)
  clumped<-rbind(clumped,file)
}
write.csv(clumped, file="/NASdata_2/wuwq/COPD/OMIC_ALL/clump/clumped.csv")
