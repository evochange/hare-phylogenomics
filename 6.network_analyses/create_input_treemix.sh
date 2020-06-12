#!/bin/bash
# This list of commands will generate SNP alignments that will be used as input for treemix.
# The pipeline is similar to create_svdquartets_input.sh except we do not allow missing data 
# in our SNPs! So, filter_missing_vcf.py is set to 0. And then we proceed on extending targets 
# and selecting SNPs every 10 kb. We also use L. callotis as outgroup, so we do not call the 
# rabbit reference.

#for chromosome as $1
chromosome=$1
parallel 'echo "creating svdquartets input for {}"' ::: $chromosome

cd chr_"$chromosome"

############## pipeline with Lepus only (LCL will be outgroup) ##################

# Symbolic link all the necessary fasta files inside each chromosome file
parallel --colsep '\t' 'ln -s /path/to/consensus_per_individual/{1}/{1}_{2}.fasta.fa ./' :::: ../samples.txt ::: $chromosome
# Fix names of chromosomes inside individual files:
# >10 to >Name_of_sample
for i in $(ls *.fasta.fa); do header=${i/.fasta.fa/}; o=${i/.fasta.fa/.new.fa}; sed "s/$chromosome/$header/" "$i" > $o; done
rm *.fasta.fa
# Construct alignment
cat *.new.fa > chr_"$chromosome"_full_algn.fa
# Remove extra outgroups from the alignments
# Call snps to a vcf file
snp-sites-2.3.3 -v -o chr_"$chromosome"_full_algn.vcf chr_"$chromosome"_full_algn.fa
# Filter positions that have too much missing data:
filter_missing_vcf.py chr_"$chromosome"_full_algn.vcf 0 chr_"$chromosome"_full_algn_filtered.vcf

# If the file has more than the 4 header lines:

num=$(wc -l < chr_"$chromosome"_full_algn_filtered.vcf)

if [[ "$num" -gt 4 ]]; then
	
	echo "HURRAY, you have SNPS!"
	# Convert the vcf file to bed file
	sed -e 's/chr//' chr_"$chromosome"_full_algn_filtered.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' > chr_"$chromosome"_full_algn_snps_filtered.bed
	# Change the CHROM column in Bed from 1 to chromosome name
	awk -v var="$chromosome" 'BEGIN {FS="\t"} {OFS="\t"} {print var, $2, $3, $4, $5}' chr_"$chromosome"_full_algn_snps_filtered.bed > chr_"$chromosome"_full_snps.bed

	#### EXTEND TARGETS ######
	ln -s /home/mafalda_ferreira/reference/Oryctolagus.chromsizes ./
	grep -w "$chromosome" Oryctolagus.chromsizes > chr_"$chromosome".chromsize
	# Extend targets with bedtools slop.
	bedtools slop -i ~/bed_file/"$chromosome".bed -g chr_"$chromosome".chromsize -b 200 > chr_"$chromosome"_200bp.bed
	# This will merge overlapping features and bookended features (finishing and starting at the same position)
	bedtools merge -i chr_"$chromosome"_200bp.bed > chr_"$chromosome"_200bp_merged.bed
	# Intersect my bed_file for the targets with the positions of snps for chromosome 10
	bedtools intersect -a chr_"$chromosome"_200bp_merged.bed -b chr_"$chromosome"_full_snps.bed > chr_"$chromosome"_target_snps.bed
	
	# Extract snps that dist 10kb from each other in bed file
	filter_bed.py chr_"$chromosome"_target_snps.bed chr_"$chromosome"_target_snps_w10kb.bed
	# Change back all the >Name_of_sample to >chromosome name in all fastas (reverse of this command:)
	for i in $(ls *.new.fa); do header=${i/.new.fa/}; o=${i/.new.fa/.new.2.fa}; sed "s/$header/$chromosome/" "$i" > $o; done
	# Apply getfasta from bedtools
	ls *.new.2.fa | parallel "bedtools getfasta -fi {} -bed chr_'$chromosome'_target_snps_w10kb.bed -fo {}_snps" 
	
	# Concatenate all positions into one sequence
	ls *.new.2.fa_snps | parallel 'concatenate_snps.py {} {}_seq'
	# Produce a new alignment for chromosome 10
	cat *.fa_snps_seq > chr_"$chromosome"_final_snps_LepusOnly.fa
	# Remove extra files
	rm *.new.2.fa
	rm *.new.2.fa_snps
	rm *.new.2.fa_snps_seq
	rm *.new.fa
	rm *.fai
	
else
	echo "VCF file of $chromosome is empty. No SNPs remain."
fi
