#!/bin/bash
# This script loops through chromosome files to generate whole-genome alignments
# per chromosome (972) for the samples in file samples.txt
# it converts the fasta alignment into a geno file, which is the input
# for Genomics Generals

#for chromosome as $1
chromosome=$1
parallel 'echo "creating Genomics Generals alignments for {}"' ::: $chromosome
#counting time for each chromosome
start=`date +%s`

cd chr_"$chromosome"

#retrieve the consensus genome wide chromosome sequence for each individual (fasta.fa file):
parallel --colsep '\t' 'ln -s /path/to/consensus_per_individual/{1}/{1}_{2}.fasta.fa ./' :::: ../samples.txt ::: $chromosome
#retrieve the Oryctolagus sequence:
faidx -x ~/reference/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa $chromosome
mv $chromosome.fa OryCun_$chromosome.fasta.fa
# change the headers
for i in $(ls *.fasta.fa); do header=${i/_"$chromosome".fasta.fa/}; o=${i/.fasta.fa/.new.fa}; sed "s/$chromosome/$header/" "$i" > $o; done
rm *.fasta.fa 
# convert fasta into linear format (fasta2geno.pl needs linear format:
for i in $(ls *.new.fa); do o=${i/.fa/_l.fa}; fasta_formatter -w 0 -i $i -o $o; done

# concatenate all sequences in one alignment:
cat *.new_l.fa > chr_"$chromosome"_full_align.fa

# remove extra files:

rm *.new.fa
rm *.new_l.fa

# Convert the alignment to geno type file:
/usr/bin/python ~/my_programs/genomics_general/seqToGeno.py -s chr_"$chromosome"_full_align.fa -g chr_"$chromosome"_full_align.geno -f fasta --mode samples -C $chromosome

# I decided to apply a 1/3 filter:
/usr/bin/python ~/my_programs/genomics_general/filterGenotypes.py --infile chr_"$chromosome"_full_align.geno --minCalls 11 -t 8 -o chr_"$chromosome"_filtered.geno



