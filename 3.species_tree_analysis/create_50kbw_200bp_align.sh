#!/bin/bash
###### Written by Mafalda S. Ferreira
###### this scripts creates 50 kb window alignments from chromosome alignments with my 
###### individuals in "sample.txt". it will remove regions of the genome out of the targets 
###### + 200bp used for exome capture


# This is the chromosome we are processing now:
$chromosome=$1

# Enter the chromosome folder
cd chr_"$chromosome"_200

# Retrieve the consensus genome wide chromosome sequence for each individual (fasta.fa file):
parallel --colsep '\t' 'ln -s /path/to/consensus/sequences/{1}/{1}_{2}.fasta.fa ./' :::: ../samples.txt ::: $chromosome

# Retrieve Oryctolagus sequence.
faidx -x ~/reference/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa $chromosome
# Name it as the rest of the files.
mv $chromosome.fa OryCun_$chromosome.fasta.fa

# Mask regions of chromosome not contained inside target files
## 1. Bed file with coordinates per chromosome
grep -w "$chromosome" ~/reference/Oryctolagus.chromsizes > chr_"$chromosome".chromsize
## 2. Extend targets and create complement bed file to mask
bedtools slop -i ~/bed_file/"$chromosome".bed -g chr_"$chromosome".chromsize -b 200 > chr_"$chromosome"_200bp.bed
# 3. Merge overlapping features and bookended features (finishing and starting at the same position)
bedtools merge -i chr_"$chromosome"_200bp.bed > chr_"$chromosome"_200bp_merged.bed
bedtools complement -i chr_"$chromosome"_200bp_merged.bed -g chr_"$chromosome".chromsize > chr_"$chromosome"_target_complement_200bp.bed

# Mask regions outside extended targets:
parallel 'bedtools maskfasta -fi {1}_{2}.fasta.fa -bed chr_{2}_target_complement_200bp.bed -fo {1}.masked.fa -mc ?' :::: ../samples.txt ::: $chromosome

# Concatenate everything into one sequence
for i in $(ls *.masked.fa); do header=${i/.masked.fa/}; o=${i/.masked.fa/.masked.new.fa}; sed "s/$chromosome/$header/" "$i" > $o; done

# Remove excess files.
rm *.masked.fa

# Create whole chromosome alignment.
cat *.masked.new.fa > chr_"$chromosome"_full_algn.fa

# Remove excess files.
rm *.masked.new.fa

# Split all chromosome alignments in windows of 50kb. 
## The output fasta files will have names with coordinates. We can use the coordinates to filter alignments outside extended targets later on.
~/my_programs/phast/bin/msa_split chr_"$chromosome"_full_algn.fa --windows 50000,0 --in-format FASTA --out-format FASTA --out-root chr_"$chromosome"
# Zip the chromosome alignment because it occupies a lot of space!
find . -name "*_full_algn.fa" | xargs pigz 

# Grab the bed file that has the target coordinates for this chromosome
ln -s ~/bed_file/"$chromosome".bed ./

# Grab the bed file that has the coordinates of 50kb windows across a chromosome
ln -s /scratch/3/mafalda/50kbwindows_ref/"$chromosome"_output.bed ./

# The intersection will tell you what windows you have targets
# Is there any intersection between the extended targets and the 50kb windows across the genome?
bedtools intersect -a "$chromosome"_output.bed -b chr_"$chromosome"_200bp_merged.bed -wa | uniq > chr_"$chromosome"_overlap.bed

# The intervals in the fasta names will be 1 based and not 0 based
# This awk will add 1 to all the intervals in the bed file, so we can match the names of the fasta files and the bed intervals 
awk -v var="$chromosome" 'BEGIN {FS="\t"} {OFS="-"} {print $2+1, $3}' chr_"$chromosome"_overlap.bed > chr_"$chromosome"_windows_to_keep.txt
# Do a reverse grep to find the windows we should remove that do not overlap with any target
ls chr_"$chromosome".*-*.fa | grep -f chr_"$chromosome"_windows_to_keep.txt -w -v > chr_"$chromosome"_windows_to_remove.txt
# Filter windows
parallel 'rm {}' :::: chr_"$chromosome"_windows_to_remove.txt

# Convert your fasta alignments to single lines.
for i in $(ls chr_*.*-*.fa); do o=${i/.fa/_l.fasta}; seq2line.py $i > $o; done
# Remove question marks (missing data) from your alignments.
for i in $(ls *_l.fasta); do o=${i/_l.fasta/_new.fasta}; sed 's/?//g' $i > $o ; done
# Rm excess files
rm *chr_*.*-*.fa 
rm *_l.fasta
rm *.fasta.fa

# Use TriSeq to filter for coverage (33% = 9 individuals)
for i in $(ls *_new.fasta); do o=${i/.fasta/_filter}; TriSeq -in $i -of fasta -o $o --upper-case --missing-filter 0 9; done

# Rm excess files
rm *chr_*.*-*.fa 
rm *_l.fasta
rm *.fasta.fa
rm *_new.fasta