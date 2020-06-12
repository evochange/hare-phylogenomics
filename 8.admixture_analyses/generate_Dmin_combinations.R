library(treeman)
library(gtools)

# Read your tree file
tree<-readTree(file="astral_species_all_gene_trees.tree")
# Get the tips from your species tree
tree['tips']->species
# Remove outgroups, because we don't need them for the combinations.
species[c(-1,-17)]->species
# Make all possible combinations of 3 with our list of species.
combinations<-combinations(length(species),3,species)
colnames(combinations)<-c("A","B","C")
# Now, let's generate the ABBA-BABA conformations for each combination of species.
all_combinations<-as.data.frame(matrix(ncol=3,nrow=1))
colnames(all_combinations)<-c("A","B","C")

combinations<-comb
for(i in 1:nrow(combinations)){
  A=combinations[i,1]
  B=combinations[i,2]
  C=combinations[i,3]
  
  c1=as.data.frame(matrix(c(A,B,C),ncol=3,nrow=1))
  c2=as.data.frame(matrix(c(A,C,B),ncol=3,nrow=1))
  c3=as.data.frame(matrix(c(B,C,A),ncol=3,nrow=1))
  
  colnames(c1)<-c("A","B","C")
  colnames(c2)<-c("A","B","C")
  colnames(c3)<-c("A","B","C")
  
  all_combinations<-rbind(all_combinations,c1)
  all_combinations<-rbind(all_combinations,c2)
  all_combinations<-rbind(all_combinations,c3)
  
}

# Write the output.
write.table(all_combinations[-1,],file="species_combinations_Dmin.txt",sep="\t",quote=F,col.names = T,row.names = F)
