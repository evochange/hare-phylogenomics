library(dplyr)
library(scales)
library(biomaRt) 

# Read the three tables with fd distributions.
AB_files<-c("L_othus_vs_L_americanus_fdscan_50k5k.out","L_townsendii_vs_Lamericanus_fdscan_50k5k.out"
,"L_timidus_vs_Lamericanus_fdscan_50k5k.out")
AB_tables = lapply(AB_files, read.csv)
head(AB_tables[[2]])

#intables=AB_files
# convert all fd values to 0 where D is negative:
for (x in 1:length(AB_tables)){
  AB_tables[[x]]$fd = ifelse(AB_tables[[x]]$D < 0, 0, AB_tables[[x]]$fd)
}
# Find the quantiles for each scan.
alpha=0.995
numberrow=as.numeric(length(AB_tables))
output<-matrix(ncol=2,nrow=numberrow)
for(i in 1:length(AB_tables)){
  table_i<-AB_tables[[i]]
  output[i,1]=AB_files[i]
  output[i,2]=as.numeric(quantile(table_i$fd,alpha,na.rm=T))
}
print(data.frame(output))


sp1<-AB_tables[[1]]
sp2<-AB_tables[[2]]
sp3<-AB_tables[[3]]

paste("chr_",sp1$scaffold,sep="")->sp1$scaffold
paste("chr_",sp2$scaffold,sep="")->sp2$scaffold
paste("chr_",sp3$scaffold,sep="")->sp3$scaffold

colnames(sp1)<-c("scaffold","start","end","mid","sites",
                 "sitesUsed","ABBA","BABA","D_sp1","fd_sp1",
                 "fdM_sp1")
colnames(sp2)<-c("scaffold","start","end","mid","sites",
                 "sitesUsed","ABBA","BABA","D_sp2","fd_sp2",
                 "fdM_sp2")
colnames(sp3)<-c("scaffold","start","end","mid","sites",
                 "sitesUsed","ABBA","BABA","D_sp3","fd_sp3",
                 "fdM_sp3")

sp1 %>%
  mutate(significance_sp1 = ifelse(fd_sp1 >=as.numeric(output[1,2]),"Sig","NotSig"))->sp1
sp2 %>%
  mutate(significance_sp2 = ifelse(fd_sp2 >=as.numeric(output[2,2]),"Sig","NotSig"))->sp2
sp3 %>%
  mutate(significance_sp3 = ifelse(fd_sp3 >=as.numeric(output[3,2]),"Sig","NotSig"))->sp3

merge(sp1,sp2,by=c("scaffold","start","end"))->fd_temp1
merge(sp3,fd_temp1,by=c("scaffold","start","end"))->fd_temp2

colnames(fd_temp2)

# Find windows that are only significant in one of the tests.
private_sig1<-subset(fd_temp2,fd_temp2$significance_sp1=="Sig" & fd_temp2$significance_sp2 =="NotSig" & fd_temp2$significance_sp3 =="NotSig")
private_sig2<-subset(fd_temp2,fd_temp2$significance_sp1=="NotSig" & fd_temp2$significance_sp2 =="Sig" & fd_temp2$significance_sp3 =="NotSig")
private_sig3<-subset(fd_temp2,fd_temp2$significance_sp1=="NotSig" & fd_temp2$significance_sp2 =="NotSig" & fd_temp2$significance_sp3 =="Sig")


# Find windows that are significant in all tests.
collective<-subset(fd_temp2,fd_temp2$significance_sp1=="Sig" & fd_temp2$significance_sp2 =="Sig" & fd_temp2$significance_sp3 =="Sig")

# Let's output these windows.
rows_to_include<-rownames(subset(fd_temp2,fd_temp2$significance_sp3 =="Sig" | fd_temp2$significance_sp2 =="Sig" | fd_temp2$significance_sp1 =="Sig"))
as.numeric(rows_to_include)->rows_to_include
fd_temp2[rows_to_include,]->sigwindows

write.table(sigwindows,file="significant_windows_in_sp1_sp2_sp3.txt",col.names=T,row.names=F,quote=F,sep="\t")

#To generate the plot in the paper (Figure 5):

data.frame(count(private_sig1,scaffold))->counts_sig1
data.frame(count(private_sig2,scaffold))->counts_sig2
data.frame(count(private_sig3,scaffold))->counts_sig3
data.frame(count(collective,scaffold))->counts_col

counts_mainchr_sig1<-counts_sig1[- grep("GL", counts_sig1$scaffold),]
counts_mainchr_sig2<-counts_sig2[- grep("GL", counts_sig2$scaffold),]
counts_mainchr_sig3<-counts_sig3[- grep("GL", counts_sig3$scaffold),]
counts_mainchr_col<-counts_col[- grep("GL", counts_col$scaffold),]

chrorder<-paste("chr_",seq(1:22),sep="")
chrorder<-c(chrorder,"chr_X")

factor(counts_mainchr_col$scaffold, chrorder, ordered=TRUE)->counts_mainchr_col$scaffold
counts_mainchr_col[do.call(order, counts_mainchr_col[, c("scaffold","n")]), ]->counts_mainchr_col

factor(counts_mainchr_sig1$scaffold, chrorder, ordered=TRUE)->counts_mainchr_sig1$scaffold
counts_mainchr_sig1[do.call(order, counts_mainchr_sig1[, c("scaffold","n")]), ]->counts_mainchr_sig1
factor(counts_mainchr_sig2$scaffold, chrorder, ordered=TRUE)->counts_mainchr_sig2$scaffold
counts_mainchr_sig2[do.call(order, counts_mainchr_sig2[, c("scaffold","n")]), ]->counts_mainchr_sig2
factor(counts_mainchr_sig3$scaffold, chrorder, ordered=TRUE)->counts_mainchr_sig3$scaffold
counts_mainchr_sig3[do.call(order, counts_mainchr_sig3[, c("scaffold","n")]), ]->counts_mainchr_sig3


plot.default(0,axes=F,col="white",
             ylab="number of outlier windows",xlab="chromosome",ylim=c(0,80),xlim=c(1,23))
axis(1,at=seq(1,23),labels = c(1:22,"X"))
axis(2)
abline(v = seq(1,23),col="gray75")
points(counts_mainchr_sig1$scaffold,counts_mainchr_sig1$n,pch=19,col=alpha(colour="blue",alpha=0.5),cex=1)
points(counts_mainchr_sig2$scaffold,counts_mainchr_sig2$n,pch=19,col=alpha(colour="black",alpha=0.5),cex=1)
points(counts_mainchr_sig3$scaffold,counts_mainchr_sig3$n,pch=19,col=alpha(colour="purple",alpha=0.5),cex=1)
points(counts_mainchr_col$scaffold,counts_mainchr_col$n,pch=19,col=alpha(colour="orange",alpha=0.5),cex=2)

# Also, this bit will annotate whatever data.frame you would like using biomaRt.
# recent ensembl:
ensembl = useEnsembl(biomart="ensembl",dataset="ocuniculus_gene_ensembl")

significance=99.5

paste(private_sig2$scaffold,private_sig2$start,private_sig2$end,sep=":")->sig_ensembl
gsub("chr_",replacement = "",x = sig_ensembl)->sig_ensembl
#print(sig_ensembl_i)
sig_query_i<-getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position','external_gene_name','go_id','name_1006','definition_1006','namespace_1003'),
                   filters=c('chromosomal_region'),
                   values=sig_ensembl,
                   mart=ensembl)

name<-paste("outfile_name.txt",sep="")
write.table(as.data.frame(sig_query_i),file=name,col.names = T,row.names = F,quote=F,sep="\t")
