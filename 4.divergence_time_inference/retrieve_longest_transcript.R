# This script will load packages biomaRt and data.ble.
# biomaRt will retrieve ENSEMBL information for the coordinates of your targets, such
# as gene ID, transcript ID, transcript length, etc. The objective is to retrieve the
# longest transcript if, for a given target region, you have more than one transcript.
# Author: Mafalda S. Ferreira.

library(biomaRt)
library(data.table)

# Load the ensembl data base
ensembl = useEnsembl(biomart="ensembl",dataset="ocuniculus_gene_ensembl")

# Load the file that has the coordinates of your targets.
regions<-read.table("141216_Lepus_Ex1_MJ_EZ_HX1_capture_targets.regions",col.names = "Regions")


# Create a query object with the information you will want to retrive for each one
# of your regions. Don't forget "transcript_length".
query<-getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','transcript_length','chromosome_name','external_gene_name','genomic_coding_start','genomic_coding_end'),
             filters=c('chromosomal_region'),
             values=regions,
             mart=ensembl)

# Create a backup, because sometimes I make mistakes.
query_backup<-query

# How many genes do we have in our capture?
genes<-unique(query$ensembl_gene_id)
length(genes) # 18798

# Convert the query into a data.table to process it.
query<-data.table(query)
filtered_query<-as.data.frame(matrix(ncol=ncol(query),nrow=1))
colnames(filtered_query)<-colnames(query)

# This bit of code will, for each gene id, obtain the id of the longest transcript available.
for(i in 1:length(genes)){
  gene<-genes[i]
  temp<-query[ensembl_gene_id==gene]
  
  # Find the transcript with maximum length:
  n<-temp[,lapply(.SD,function(x) which.min(x)),.SDcols="transcript_length"]$transcript_length
  
  transcript_id<-temp[n,]$ensembl_transcript_id
  # Sse data for that transcript to fill in the new filtered_query
  output<-query[ensembl_transcript_id==transcript_id]
  
  filtered_query<-rbind(filtered_query,output)
  
}

# Remove NAs from genomic_coding_start and genomic_coding_end  because it's stupid...
na.omit(filtered_query)->filtered_query

# Write this table:
write.table(filtered_query,file="Lepus_Ex1_MJ_EZ_HX1_capture_targets.longest_transcript.regions",col.names = T,row.names = F,quote=F, sep="\t")
# Write just the list of transcript ensembl IDS:
write.table(unique(filtered_query$ensembl_transcript_id),file="Lepus_Ex1_MJ_EZ_HX1_capture_targets.longest_transcript.list",col.names = F,row.names = F,quote = F)
