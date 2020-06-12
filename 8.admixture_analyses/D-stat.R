
#setwd("~/Documents/OneDrive - Universidade do Porto/BIOLOGIA/Doutoramento/Chapter 1 Phylogenomics/2-Pipeline/Dstats")
library(data.table)
source("/scratch/1/users/m/mf239628e/f_hom/test_GG_fd/genomics_general-master/jackknife.R")
#source("genomics_general-master/jackknife.R")

D.stat <- function(p1, p2, p3) {
  ABBA <- (1 - p1) * p2 * p3
  BABA <- p1 * (1 - p2) * p3
  (sum(ABBA) - sum(BABA)) / (sum(ABBA) + sum(BABA))
}

# get command line arguments:
args <- commandArgs(trailingOnly = TRUE)
allele_freqs_file_name <- args[1]
output_file_name <- args[2]
spc_p1 <- args[3]  # A
spc_p2 <- args[4]  # B
spc_p3 <- args[5]  # C

cat(paste("Calculating Dstat for P1:",spc_p1,"P2:",spc_p2,"P3:",spc_p3,"\n",sep=" "))

freq_table = read.table(allele_freqs_file_name, header=T, as.is=T)

# Open the output file.
output_file <- file(output_file_name, "w")

#remove NAs
subset_freq_table<-na.omit(freq_table[,c("scaffold","position",spc_p1,spc_p2,spc_p3)])
rownames(subset_freq_table) <- NULL

cat(paste("Computing unfiltered D-stat...\n"))
## get normal D-stat value with all positions
D_stat<-D.stat(subset_freq_table[,spc_p1], subset_freq_table[,spc_p2], subset_freq_table[,spc_p3])

cat(paste("D_stat:",D_stat,"\n",sep=" "))

## compute jackknife:
cat(paste("Computing genome blocks for D...\n"))


block_indices <- get_block_indices(block_size=1000000,
                                   positions=subset_freq_table$position,
                                   chromosomes=subset_freq_table$scaffold)

# remove blocks with missing data:
block_indices[lapply(block_indices,length)>0]->block_indices

n_blocks <- length(block_indices)
cat(paste("Genome divided into", n_blocks, "blocks.\n"))

cat(paste("Computing Z-score for D...\n"))
Dstat_sd <- get_jackknife_sd(block_indices=block_indices,FUN=D.stat,subset_freq_table[,spc_p1],subset_freq_table[,spc_p2],subset_freq_table[,spc_p3])

cat(paste("Dstat sd:",Dstat_sd,"\n",sep=" "))

D_err <- Dstat_sd/sqrt(n_blocks)
D_Z <- D_stat / D_err

cat(paste("D_err:",D_err,"\n",sep=" "))
cat(paste("Z-score:",D_Z,"\n",sep=" "))

##### Output ######
header<-c("P1","P2","P3","D_stat","n_blocks","n_positions","Z_score","D_stat_f","n_blocks_f","n_positions_f","Z_score_f")
output_df<-as.data.frame(matrix(ncol=length(header)),nrow=1)
colnames(output_df)<-header

output_df$P1<-spc_p1
output_df$P2<-spc_p2
output_df$P3<-spc_p3
output_df$D_stat<-D_stat
output_df$n_blocks<-n_blocks
output_df$n_positions<-nrow(subset_freq_table)
output_df$Z_score<-D_Z

write.table(output_df,file=output_file,col.names = T,sep="\t",row.names=F,quote=F)

# Close the output file.
close(output_file)

# Feedback.
cat(paste("\nWrote results to file ", output_file_name, ".\n\n", sep=""))


