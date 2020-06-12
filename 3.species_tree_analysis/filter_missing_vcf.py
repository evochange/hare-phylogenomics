#!/usr/bin/python
# a python to filter the output of snps to sites (yes, a vcf file!)
# 21 March 2018
# MafaldaSF

import sys
from collections import Counter

input=open(sys.argv[1],'r')
min_no_individuals=sys.argv[2]
output=open(sys.argv[3],'w')

#this function will index the position of the missing data in the ALT column
def symbol_of_missing(vcf):
	
	symbol_dict={}
	
	for line in vcf:
		if line.startswith("#"):
			pass
		else:
			line=line.strip().split()
			#access the ALT column
			POS=line[1]
			ALT=line[4]
			#split ALT by column:
			alleles=ALT.split(',')
			
			#indexing a line with no missing data will result
			#in an Error, so we try the code, and if there is no
			# *, we atribute NA to the symbol
			try:
				symbol=int(alleles.index("*"))+1
				symbol=str(symbol)
				symbol_dict.setdefault(POS,symbol)
			except ValueError:			
				symbol_dict.setdefault(POS,'NA')
			
	return symbol_dict

#retrieve the symbol for missing data for each position			
dict=symbol_of_missing(input)

#open input file because why not (it doesn't work if I don't!!)
input_vcf=open(sys.argv[1],'r')

for line in input_vcf:
	if line.startswith("#"):
		output.write(line)
	else:
		input_line=line.strip().split('\t')
		POS=input_line[1]

		if dict[POS]=="NA": # there is no missing data, so write line to output:
			output.write(line)
		else:
			genotypes=[]
		
			for field in input_line[9:len(input_line)]:
				genotypes.append(field)
	
			count_of_missing=Counter(genotypes)
			
			symbol=dict[POS]
			
			if int(count_of_missing[symbol])>int(min_no_individuals):
				pass
			else:
				output.write(line)
		
