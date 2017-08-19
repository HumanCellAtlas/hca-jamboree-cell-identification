# Task 3- identifying dead/dying cells
# Will Townes

library(Matrix)
library(scater)
library(biomaRt)
library(tidyverse)

a<-read.table("jamboree/kolodziejczyk_2015/counts/kolodziejczk_annotations.tsv",header=TRUE)
pd <- new("AnnotatedDataFrame", data = a)
rownames(pd)<-a$cell_name
m<-read.table("jamboree/kolodziejczyk_2015/counts/kolodziejczk_counttable_es.tsv")
m<-Matrix(as.matrix(m),sparse=TRUE)
Z<- m>0
zrs<- Matrix::rowSums(Z)
m<- m[zrs>5,]
#sc<-newSCESet(countData = m, phenoData = pd)
gd<-data.frame(ensembl_id=rownames(m))
mart<-useMart("ENSEMBL_MART_ENSEMBL","mmusculus_gene_ensembl")
bm<-getBM(attributes= c("ensembl_gene_id","hgnc_symbol","chromosome_name","gene_biotype"),mart=mart)
mtg<-bm$ensembl_gene_id[bm$chromosome_name=="MT"]
rbp<-read.table("http://www.genenames.org/cgi-bin/genefamilies/set/1054/download/branch",header=TRUE,sep="\t")
colnames(rbp)[colnames(rbp)=="Approved.Symbol"]<-"hgnc_symbol"
rbp<-plyr::join(rbp,bm,by="hgnc_symbol")
bt<-unique(bm$gene_biotype)
c("Mt_rRNA","Mt_tRNA","rRNA","")




