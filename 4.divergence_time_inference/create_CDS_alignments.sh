#!/bin/bash
# this script extracts CDS alignments from a gff3 file and the consensus chromosome
# for each individual
# the analysis is done by transcript and chromosome
# HASHTAG OUT reverse_concatenate_reverse_CDS.py WHEN USING THE FORWARD TRANSCRIPTS

transcript=$1
chromosome=$2

echo "generating alignments for $transcript of chromosome $chromosome"

cd $transcript

# Call consensus chromosome fasta sequences for each sample.
parallel --colsep="\t" 'ln -s /path/to/consensus_per_individual/{1}/{1}_{2}.fasta.fa ./' :::: ../samples.txt ::: $chromosome
# Unzip files, because getfasta does not handle gziped fastas.
for i in $(ls *.consensus.fa.gz); do bgzip -d $i; done

# Retrieve outgroup rabbit sequence.
faidx -x ~/reference/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa $chromosome
# Name it as the rest of the samples.
mv $chromosome.fa OryCun_$chromosome.consensus.fa

# Retrieve CDS sequence.
parallel -j 4 --colsep="\t" 'bedtools getfasta -fi {1}_{2}.fasta.fa -bed ../{3}.gff3 -fo {1}_{3}.CDS.fa -s' :::: ../samples.txt ::: $chromosome ::: $transcript

######## PART TO HASHTAG OUT ########
# Concatenate all the CDS regions in each fasta file.
# If analyzing the REVERSE STRAND transcripts, hashtag this out:
# parallel -j 4 'reverse_concatenate_reverse_CDS.py {1}_{2}.CDS.fa {1}_{2}.CDS.cat.fa' :::: ../samples.txt ::: $transcript
# If analyzing the FORWARD STRAND, hastag this out:
# for i in $(ls *CDS.fa); do echo -e ">"${i%_*}"\n"$(grep -v '^>' $i | paste -s -d "") > ${i/CDS.fa/CDS.cat.fa} ; done

# Cat fasta:
cat *.cat.fa > "$transcript".CDS_algn.fa
# Alter Ns to ?, as required by PAML/MCMCtree
sed -i.tmp '/^>/! s/N/?/g' "$transcript".CDS_algn.fa
# Convert from fasta alignment to phylip
AMAS.py convert -i "$transcript".CDS_algn.fa -f fasta -u phylip -d dna

# Clean data:
rm *.CDS.fa
rm *.CDS.cat.fa
rm *.fasta.fa.fai
rm *.fasta.fa
