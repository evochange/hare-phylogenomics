#script to calculate the f_hom combinations from a given species tree


library(treeman)
library(ape)
library(purrr)

'%ni%' <- Negate('%in%')

#tree<-randTree(9)
#summary(tree)
#tree_test<-tree
# to plot tree convert to phylo
tree2<-as(tree,"phylo")
plot(tree2)

# Astral_tree:
# It's important that your tree is rooted! Reroot it, for example, in FigTree and export
# the newick tree
astral_tree<-read.tree(file="astral_species_TCs5_gene_trees_rerooted.tree")

# Check if tree is correctly rooted as you want
# 
plot(as(astral_tree,"phylo"),show.node.label = T)
# or
plot(astral_tree,show.node.label = T,)

# Drop your roots:
drop.tip(astral_tree,c("B_idahoensis","O_cuniculus"))->astral_tree

# Convert to TreeMan
tree=as(astral_tree,"TreeMan")
df=as.data.frame(matrix(ncol=5,nrow=1))
colnames(df)<-c("branch_b","branch_a","B_descendants","A_descendants","C")

# Loop through the tree.
for(i in tree['all']){
  tips=tree['tips']
  kids=getNdKids(tree,i)
  
  sister=getNdsSstr(tree,i)
  kids_sister=getNdKids(tree,sister)
  cousins=append(kids_sister,kids)
  
  # find the root:
  if(identical(length(tips),length(intersect(kids,tips)))==TRUE){
    cat(paste(i,"is root","\n",sep=" "))
  }
  # find the two first nodes/tips that do not allow any Cs
  else if (length(tips)==length(cousins) | length(tips)-1==length(cousins)){
    print(sister)
    print(kids_sister)
    print(cousins)
    cat(paste(i,"does not have Cs to test against\n",sep=" "))
  }
  else if (i %in% tree['tips']){
    B=i
    #cat(paste(B,"is a tip","\n", sep=" "))
    A=getNdsSstr(tree,i)
    if (A %in% tree['tips']){
      C=tips[ tips != B & tips != A]
      temp_df=expand.grid(B,A,B,A,C)
      colnames(temp_df)<-colnames(df)
      df=rbind(df,temp_df)
    }
    else{
      kidsA=getNdKids(tree,A)
      C=tips[ tips != B & tips %ni% kidsA]
      #print(expand.grid(B,A,B,kidsA,C))
      temp_df=expand.grid(B,A,B,kidsA,C)
      colnames(temp_df)<-colnames(df)
      df=rbind(df,temp_df)
    }

  }
  else {
    b=i
    #cat(paste(b,"is an inner node and branch b\n",sep=" "))
    kidsB=getNdsKids(tree,i)
    #cat(paste(kidsB,"are descendants of branch b\n",sep=" "))
    a=getNdsSstr(tree,i)
   # cat(paste(a,"is the a branch\n",sep=" "))
    kidsA=getNdsKids(tree,a)
    #cat(paste(kidsA,"are descendants of branch a\n",sep=" "))
    if(is_empty(kidsA[[1]])==TRUE){
      #print("sister branch a is a tip")
      C=tips[ tips != a & tips %ni% kidsB[[1]]]
      temp_df=expand.grid(b,a,kidsB[[1]],a,C)
      colnames(temp_df)<-colnames(df)
      df=rbind(df,temp_df)
    }else {
      #print("sister branch a is node")
      C=tips[ tips %ni% kidsA[[1]] & tips %ni% kidsB[[1]] ]
      temp_df=expand.grid(b,a,kidsB[[1]],kidsA[[1]],C)
      colnames(temp_df)<-colnames(df)
      df=rbind(df,temp_df)
    }
  }
}

#remove the first line
df<-df[-1,]
write.table(df,file="fd_C_comparisonst_f_G.txt",sep="\t",quote = F,col.names = TRUE,row.names = FALSE)

