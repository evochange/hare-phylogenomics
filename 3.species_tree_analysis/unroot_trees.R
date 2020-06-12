#!/usr/bin/env Rscript

library(ape)

#### TO UNROOT THE SPECIES TREE ####
trees_list<-list.files(pattern="*bestTree*")

for (i in 1:length(trees_list)){
  tr_i<-read.tree(trees_list[i])
  tr_unroot_i<-unroot(tr_i)
  outtree2=paste(trees_list[i],"_unrooted.tree",sep="")
  write.tree(tr_unroot_i,file=outtree2)
}

#### TO UNROOT THE BOOTSTRAP TREES ####
#bs_trees_list<-list.files(pattern="*bootstrap*")
#	for (i in 1:length(bs_trees_list)){
#	bs_tr_i<-read.tree(bs_trees_list[i])
#	bs_tree_unrooted_i<-unroot(bs_tr_i)
#	bs_out_tree=paste(bs_trees_list[i],"_unrooted.tree",sep="")
#	write.tree(bs_tree_unrooted_i,file=bs_out_tree)
#}
