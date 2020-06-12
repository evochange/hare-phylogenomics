# code testing f_G
source("genomics_general/jackknife.R")

library(stringi)
library(data.table)

f.stat <- function(p1, p2, p3,p3a, p3b) {
  ABBA_numerator <- (1 - p1) * p2 * p3
  BABA_numerator <- p1 * (1 - p2) * p3
  
  ABBA_denominator <- (1 - p1) * p3b * p3a
  BABA_denominator <- p1 * (1 - p3b) * p3a
  
  (sum(ABBA_numerator) - sum(BABA_numerator)) /
    (sum(ABBA_denominator) - sum(BABA_denominator))
}

# get command line arguments:
args <- commandArgs(trailingOnly = TRUE)
allele_freqs_file_name <- args[1]
output_file_name <- args[2]
b<- args[3]
a<- args[4]
P2 <- args[5]
P1 <- args[6]
P3 <- args[7]
P3a<- args[8]
P3b<- args[9]


P1 <- "L_granatensis"
P2 <- "L_townsendii"
P3 <- "L_callotis"
P3a<- "LCL_1956_A"
P3b<- "LCL_1956_B"

cat(paste("Calculating Dstat for P1:",P1,"P2:",P2,"P3:",P3,"and branch b:",b,"and branch a:",a,"\n",sep=" "))

freq_table = read.table(allele_freqs_file_name, header=T, as.is=T)

freq_table = read.table("all_frequencies_filtered.freq.tsv", header=T, as.is=T)

#Open the output file:
output_file<- file(output_file_name,"w")

#remove NAs
subset_freq_table<-na.omit(freq_table[,c("scaffold","position",P1,P2,P3,P3a,P3b)])
rownames(subset_freq_table) <- NULL

f_G<-f.stat(subset_freq_table[,P1],subset_freq_table[,P2],subset_freq_table[,P3],subset_freq_table[,P3a],subset_freq_table[,P3b])

cat(paste("f_G:",f_G,"\n",sep=" "))

#Hypothesis 1 for jackknife:
## compute jackknife:
cat(paste("Computing genome blocks...\n"))

block_indices <- get_block_indices(block_size=1000000,
                                   positions=subset_freq_table$position,
                                   chromosomes=subset_freq_table$scaffold)

# remove blocks with missing data:
block_indices[lapply(block_indices,length)>0]->block_indices

n_blocks <- length(block_indices)

print(paste("Genome divided into", n_blocks, "blocks."))

#jackknife# 
cat(paste("Computing Z-score...\n"))

P3a_sub<-stri_sub(P3a,-2,-1)
P3b_sub<-stri_sub(P3b,-2,-1)

# We have two jacknaffing strategies. The first statement will only work for the L.callotis, because it is the only
# pseudo-haploid data we are dealing with. It will randomly shuffle the L.callotis genotypes for each site
# for each jackkniffing round. This aims to randomize any phasing biases created when we split 
# each diploid chromosome into two haploid samples.
if(P3a_sub=="_A" & P3b_sub=="_B"){
  print(P3a_sub)
  print(P3b_sub)

  overall_mean <- f_G
  values<-c()
  args<-list(subset_freq_table[,P1],subset_freq_table[,P2],subset_freq_table[,P3],subset_freq_table[,P3a],subset_freq_table[,P3b])
  for(i in 1:n_blocks){
    jackknife_list<-lapply(args,function(a) a[-block_indices[[i]]])
    jackknife_table<-cbind(jackknife_list[[1]],jackknife_list[[2]],jackknife_list[[3]],jackknife_list[[4]],jackknife_list[[5]])
    colnames(jackknife_table)<-c("V1","V2","V3","V4","V5")
   for(i in 1:nrow(jackknife_table)){
      jackknife_table[i,c(4,5)]<-sample(jackknife_table[i,c(4,5)],size=2,replace=FALSE)
    }
    jackknife_list[[4]]<-jackknife_table[,4]
    jackknife_list[[5]]<-jackknife_table[,5]
    values<-c(values,overall_mean*n_blocks - do.call(f.stat,jackknife_list)*(n_blocks-1))
  }
  f_G_sd <- sd(values)
}else{
  f_G_sd <- get_jackknife_sd(block_indices=block_indices,
                             FUN=f.stat,
                             subset_freq_table[,P1],subset_freq_table[,P2],subset_freq_table[,P3],subset_freq_table[,P3a],subset_freq_table[,P3b])
}

cat(paste("f_G sd:",f_G_sd,"\n",sep=" "))

D_err <- f_G_sd/sqrt(n_blocks)
D_Z <- f_G / D_err

cat(paste("D_err:",D_err,"\n",sep=" "))
cat(paste("Z-score:",D_Z,"\n",sep=" "))

##### Output ######
header<-c("b","a","P1","P2","P3","P3a","P3b","f_G","n_blocks","n_positions","f_G_sd","Z_score")
output_df<-as.data.frame(matrix(ncol=length(header)),nrow=1)
colnames(output_df)<-header

output_df$b<-b
output_df$a<-a
output_df$P1<-P1
output_df$P2<-P2
output_df$P3<-P3
output_df$P3a<-P3a
output_df$P3b<-P3b
output_df$f_G<-f_G
output_df$n_blocks<-n_blocks
output_df$n_positions<-nrow(subset_freq_table)
output_df$f_G_sd<-f_G_sd
output_df$Z_score<-D_Z

write.table(output_df,file=output_file,col.names = T,sep="\t",row.names=F,quote=F)

# Close the output file.
close(output_file)

# Feedback.
cat(paste("\nWrote results to file ", output_file_name, ".\n\n", sep=""))
