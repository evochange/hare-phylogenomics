#!/bin/bash
# I wanted to call SNPs for my ingroups and then call the same positions for Ory and Brachylagus
# filter_missing_vcf.py is set to 10 (10/32=0.31)
# and then we proceed on extending targets and selecting SNPs
# every 10 kb.

parallel 'echo "creating svdquartets input for {}"' ::: $chromosome

chromosome=$1

# Enter the folder.
cd chr_"$chromosome"_200

############## pipeline with Lepus + call Ory variance ##################

# symbolic link all the necessary fasta files inside each chromosome file
parallel -j 1 --colsep '\t' 'ln -s /path/to/consensus/sequences/{1}/{1}_{2}.fasta.fa ./' :::: ../samples.txt ::: $chromosome
# fix names of chromosomes inside individual files:
# >10 to >Name_of_sample
for i in $(ls *.fasta.fa); do header=${i/.fasta.fa/}; o=${i/.fasta.fa/.new.fa}; sed "s/$chromosome/$header/" "$i" > $o; done
rm *.fasta.fa
# construct alignment
cat *.new.fa > chr_"$chromosome"_full_algn.fa
# remove extra outgroups from the alignments
# call snps to a vcf file
snp-sites-2.3.3 -v -o chr_"$chromosome"_full_algn.vcf chr_"$chromosome"_full_algn.fa
# filter positions that have too much missing data:
filter_missing_vcf.py chr_"$chromosome"_full_algn.vcf 10 chr_"$chromosome"_full_algn_filtered.vcf

# if the file has more than the 4 header lines:

num=$(wc -l < chr_"$chromosome"_full_algn_filtered.vcf)

if [[ "$num" -gt 4 ]]; then
	
	echo "HURRAY, you have SNPS!"
	# convert the vcf file to bed file
	sed -e 's/chr//' chr_"$chromosome"_full_algn_filtered.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' > chr_"$chromosome"_full_algn_snps_filtered.bed
	# change the CHROM column in Bed from 1 to chromosome name
	awk -v var="$chromosome" 'BEGIN {FS="\t"} {OFS="\t"} {print var, $2, $3, $4, $5}' chr_"$chromosome"_full_algn_snps_filtered.bed > chr_"$chromosome"_full_snps.bed

	#### EXTEND TARGETS ######q
	grep -w "$chromosome" /home/mafalda_ferreira/reference/Oryctolagus.chromsizes > chr_"$chromosome".chromsize
	# extend targets 
	bedtools slop -i ~/bed_file/"$chromosome".bed -g chr_"$chromosome".chromsize -b 200 > chr_"$chromosome"_200bp.bed
	# this will merge overlapping features and bookended features (finishing and starting at the same position)
	bedtools merge -i chr_"$chromosome"_200bp.bed > chr_"$chromosome"_200bp_merged.bed
	# intersect my bed_file for the targets with the positions of snps for chromosome 10
	bedtools intersect -a chr_"$chromosome"_200bp_merged.bed -b chr_"$chromosome"_full_snps.bed > chr_"$chromosome"_target_snps.bed
	
	# intersect my bed_file for the targets with the positions of snps for chromosome 10
	bedtools intersect -a ~/bed_file/"$chromosome".bed -b chr_"$chromosome"_full_snps.bed > chr_"$chromosome"_target_snps.bed
	#extract snps that dist 10kb from each other in bed file
	filter_bed.py chr_"$chromosome"_target_snps.bed chr_"$chromosome"_target_snps_w10kb.bed
	# change back all the >Name_of_sample to >chromosome name in all fastas (reverse of this command:)
	for i in $(ls *.new.fa); do header=${i/.new.fa/}; o=${i/.new.fa/.new.2.fa}; sed "s/$header/$chromosome/" "$i" > $o; done
	# apply getfasta from bedtools
	ls *.new.2.fa | parallel "bedtools getfasta -fi {} -bed chr_'$chromosome'_target_snps_w10kb.bed -fo {}_snps" 
	
	## grab Oryc reference: 
	faidx -x ~/reference/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa $chromosome
	mv $chromosome.fa OryCun_$chromosome.new.2.fa
	bedtools getfasta -fi OryCun_$chromosome.new.2.fa -bed chr_"$chromosome"_target_snps_w10kb.bed -fo OryCun_$chromosome.new.2.fa_snps
	
	# grab Brachylagus idahoensis reference:
	ln -s /path/to/consensus/sequences/BID_RC_444/BID_RC_444_"$chromosome".fasta.fa ./
	mv BID_RC_444_"$chromosome".fasta.fa BID_RC_444_"$chromosome".new.2.fa
	bedtools getfasta -fi BID_RC_444_"$chromosome".new.2.fa -bed chr_"$chromosome"_target_snps_w10kb.bed -fo BID_RC_444_"$chromosome".new.2.fa_snps
	
	#concatenate all positions into one sequence
	ls *.new.2.fa_snps | parallel 'concatenate_snps.py {} {}_seq'
	#produce a new alignment for chromosome 10
	cat *.fa_snps_seq > chr_"$chromosome"_final_snps_OryBID_outgroup.fa
	##remove extra files
	rm *.new.2.fa
	rm *.new.2.fa_snps
	rm *.new.2.fa_snps_seq
	rm *.new.fa
	rm *.fai
	
else
	echo "VCF file of $chromosome is empty. No SNPs remain."
fi