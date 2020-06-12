#### get DMIN ####
library(data.table)

# Results by species:
D_stat_300<-data.table(read.table(file="D_stat_by_species_results_FINAL.txt",header=F,as.is=T))

# Add column names.
colnames(D_stat_300)<-c("P1","P2","P3","D_stat","n_blocks","n_positions","Z_score")

# Read the combinations file.
combinations<-read.table(file="species_combinations_Dmin.txt",header=T,as.is=T)

# Create the empty dataframe where we will write Dmin.
Dmin=as.data.frame(matrix(ncol=ncol(D_stat_300),nrow=nrow(D_stat_300)/3))
colnames(Dmin)<-c("P1","P2","P3","D_stat","n_blocks","n_positions","Z_score")

# Calculate Dmin
indexes=seq(1,nrow(combinations),by=3)

for(r in 1:length(indexes)){
  i=indexes[[r]]
  A=combinations[i,1]
  B=combinations[i,2]
  C=combinations[i,3]
  
  
  result_1=D_stat_300[P1==A & P2==B & P3==C]
  result_2=D_stat_300[P1==A & P2==C & P3==B]
  result_3=D_stat_300[P1==B & P2==C & P3==A]
  
  temp=rbind(result_1,result_2,result_3)
  
  D_stat1=abs(result_1$D_stat)
  D_stat2=abs(result_2$D_stat)
  D_stat3=abs(result_3$D_stat)
  
  n=which.min(c(D_stat1,D_stat2,D_stat3))
  
  Dmin[r,]=temp[n,]

}


# Calculate the p-value for each Dmin value.
Dmin$p_value<-2*pnorm(-abs(Dmin$Z_score))
# We set the threshold for the total amount of D-stat values calculated.
threshold<-0.05/nrow(D_stat_300)
# This gives 401 out of 455 significant Dmin values.
nrow(Dmin[Dmin$p_value<threshold,]) 
# Retrieve significant values.
Dmin[Dmin$p_value<threshold,]->Dmin_significant

write.table(Dmin_significant,"Dmin_by_species_sig.txt",col.names=T,row.names=F,quote=F)