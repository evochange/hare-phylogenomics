# Script to calculate fb_C from a list of f_hom results
library(data.table)
library(reshape2)
library(gplots)
#let's calculate fb_C:

f_G_table<-data.table(read.table("fG_results_species.txt",sep="\t",header=F,col.names = c("b","a","A","B","C","P3a","P3b","f_hom","n_blocks","n_positions","f_G_sd","Z_score")))
f_G_table_comp<-data.table(read.table("fb_C_comparisons_f_G.txt",sep="\t",header=F,col.names = c("b","a","B","A","C","P3a","P3b")))

# I rename the table, because in a previous version I calculted fbC with f_hom instead of f_G.
# Mantaining the names to avoid mistakes.
f_G_table->f_hom_table
f_hom_table_temp<-f_hom_table
f_hom_table$P3a<-NULL
f_hom_table$P3b<-NULL
f_hom_table$n_blocks<-NULL
f_hom_table$n_positions<-NULL

# When fG values are negative, they should be set to zero because they have no meaning.
f_hom_table[f_hom_table$f_hom < 0]$f_hom <-0
f_hom_table[f_hom_table$Z_score < 0]$Z_score <-0

fb_Cs<-unique(f_hom_table[,c("b","C")])
fb_Cs_results=as.data.frame(matrix(ncol=5,nrow=nrow(fb_Cs)))

for(i in 1:nrow(fb_Cs)){
  temp<-f_hom_table[b==fb_Cs[i,]$b & C==fb_Cs[i,]$C]
  As_i<-unique(temp$A)
  #print(As_i)
  fb_AC_i<-as.data.frame(matrix(ncol=8,nrow=length(As_i)))
  for(j in 1:length(As_i)){
    n<-temp[A==As_i[j],lapply(.SD,function(x) which.min(x)),.SDcols="f_hom"]$f_hom
    fb_AC_i[j,]<-temp[A==As_i[j]][n]
  }
  #print(fb_AC_i)
  fb_C_i<-median(fb_AC_i$V6)
  
  fb_C_Zscore_i<-median(fb_AC_i$V8)
  p_value_i<-2*pnorm(-abs(fb_C_Zscore_i))
  
  cat(paste(as.character(fb_Cs[i,]$b),as.character(fb_Cs[i,]$C),fb_C_i,fb_C_Zscore_i,"\n",sep=" "))
  fb_Cs_results[i,]$V1<-as.character(fb_Cs[i,]$b)
  fb_Cs_results[i,]$V2<-as.character(fb_Cs[i,]$C)
  fb_Cs_results[i,]$V3<-fb_C_i
  fb_Cs_results[i,]$V4<-fb_C_Zscore_i
  fb_Cs_results[i,]$V5<-p_value_i
  #<-matrix(as.character(fb_Cs[i,]$b),as.character(fb_Cs[i,]$C),fb_C_i,fb_C_Zscore_i,ncol=4)
}
colnames(fb_Cs_results)<-c("b","C","fbC","Zscore","pvalue")

acast(fb_Cs_results,b~C,value.var="fbC")->fbC_matrix
acast(fb_Cs_results,b~C,value.var="Zscore")->fbC_Zscore_matrix
acast(fb_Cs_results,b~C,value.var="pvalue")->fbC_pvalue_matrix

matrix(rep(NA,15),nrow=1,ncol=15)->L_callotis_matrix
colnames(L_callotis_matrix)<-colnames(fbC_matrix)
rownames(L_callotis_matrix)<-c("L_callotis")
rbind(L_callotis_matrix,fbC_matrix)->fbC_matrix
rbind(L_callotis_matrix,fbC_pvalue_matrix)->fbC_pvalue_matrix

# Roworder for analysis by species:
roworder<-c("L_callotis","L_californicus","n27","L_americanus","n26","n25","L_europaeus","n24","n23","L_habessinicus","L_starki","n20","L_fagani","L_capensis","n16","n6","L_castroviejoi","L_corsicanus","n15","L_mandschuricus","n14","L_granatensis","n13","L_townsendii","n12","L_othus","L_timidus")
colorder=c("L_timidus","L_othus","L_townsendii","L_granatensis","L_mandschuricus","L_corsicanus","L_castroviejoi","L_capensis","L_fagani","L_starki","L_habessinicus","L_europaeus","L_americanus","L_californicus","L_callotis")

fbC_matrix[,colorder]->fbC_matrix
fbC_matrix[roworder,]->fbC_matrix
fbC_pvalue_matrix[,colorder]->fbC_pvalue_matrix
fbC_pvalue_matrix[roworder,]->fbC_pvalue_matrix

bfr_cor<-0.05/910 # by species
fbC_pvalue_matrix[fbC_pvalue_matrix>bfr_cor]<-NA
round(fbC_pvalue_matrix,1)->fbC_pvalue_matrix
fbC_pvalue_matrix[fbC_pvalue_matrix==0]<-'*'

#### plot heatmap ####

matrix.min  <- min(fbC_matrix[ fbC_matrix!=0 ], na.rm=TRUE)
matrix.max <- max( fbC_matrix, na.rm = TRUE )
pairs.breaks <- c(0, seq( matrix.min, matrix.max, length.out=50) )
color.palette= colorRampPalette(c("white","#FFBA3C"))(length(pairs.breaks) - 2)
color.palette= colorRampPalette(c("white",'#FF7F2A'))(length(pairs.breaks))

#color.palette= colorRampPalette(c("snow","tomato4"))(length(pairs.breaks) - 2)

#lmat_temp=rbind(c(4,3,0),c(2,1,0),c(0,0,0))
dev.off()
lmat_temp=rbind(c(3,4,0),c(2,1,0),c(0,0,0))
lhei_temp=c(2,6,1)
lwid_temp=c(0.3,3,0.4)
heatmap.2(fbC_matrix,Colv=NA,Rowv=NA,col=c(color.palette),
          trace="none",dendrogram="none",
          density.info="none",
          lmat=lmat_temp,lhei=lhei_temp,lwi=lwid_temp,
          sepcolor="black",colsep=1:ncol(fbC_matrix),rowsep=1:nrow(fbC_matrix),
          sepwidth=c(0.005,0.005),
          cellnote=fbC_pvalue_matrix,notecol="black",notecex=2,
          labCol = F,labRow = F,
          na.color="#66666650",
          key=F)
