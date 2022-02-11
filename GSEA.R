library(clusterProfiler)
genlist<-read.table("d:/20211220_CJJRNAseq/counts/plasma_mature_genset.Gsea.1640160833967/ranked_gene_list_P_versus_N_1640160833967.tsv",sep = "\t",header = TRUE)
gene_list<-genlist[1:55291,c(1,3)]
geneList = gene_list[,2]
names(geneList) = as.character(gene_list[,1])
geneList = sort(geneList, decreasing = TRUE)
genset<-read.table("c:/Users/xjmik/Desktop/geneset/GSE4142.top.table.tsv",sep = "\t",header = TRUE,row.names = 1)
library(tidyr)
genset_new<-separate(genset,Gene.symbol,sep = "/",into = c("a","b"))
library(dplyr)
genset_new1<-filter(genset_new,a !="")
genset_new2<-genset_new1[!duplicated(genset_new1$a),]
colnames(genset_new2)<-c("adj","p","logfc","gene","blank")
genset_new3<-filter(genset_new2,adj < 0.05 & logfc >0)
genset_new4<-genset_new3[order(genset_new3$logfc,decreasing = TRUE),]
genset_new5<-toupper(genset_new4$gene)
genset_new6<-as.data.frame(genset_new5)
geneanno<-data.frame(rep("A",length(genset_new5)))
genset_new7<-cbind(geneanno,genset_new6)
colnames(genset_new7)<-c("term","gene")
C<-data.frame()

for (i in length(rownames(genset_new7)):1) {
  genset_new8<-genset_new7[1:i,]
  B<-GSEA(geneList = geneList,TERM2GENE = genset_new8,maxGSSize = 5000)
  if (length(rownames(B@result)) != 0) {
    if(B@result$p.adjust < 0.05){
      D<-data.frame(length(rownames(genset_new8)),B@result$NES,B@result$p.adjust)
      colnames(D)<-c("Length","NES","padj")
      C<-rbind(C,D)
    }
  }
}


