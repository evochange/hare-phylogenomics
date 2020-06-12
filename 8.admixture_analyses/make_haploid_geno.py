#!/usr/bin/python

import pandas as pd
import sys
import random

# define input and output
input=sys.argv[1]
output=sys.argv[2]

def IUPAC(letter):
	dict={}
	
	dict.setdefault("R",["A","G"])
	dict.setdefault("Y",["C","T"])
	dict.setdefault("S",["G","C"])
	dict.setdefault("W",["A","T"])
	dict.setdefault("K",["G","T"])
	dict.setdefault("M",["A","C"])
	
	alleles=dict[letter]
	
	return alleles	

print "Reading input file %s..." % (input)
#read input file into dataframe
dataset=pd.read_csv(input,delimiter="\t",dtype={'#CHROM':object},compression="gzip")
#retrieve column names
names=list(dataset)[2:len(list(dataset))]
#retrieve length of alignment
b=len(dataset)


new_dataset=pd.DataFrame()
new_dataset=new_dataset.append(dataset.iloc[:,0:2])

for i in names:
	if i=="OryCun":
		new_dataset["OryCun"]=dataset.loc[:,"OryCun"]
	else:
		sequence=''.join(dataset.loc[:,i])
		names1=''.join(i+"_A")
		names2=''.join(i+"_B")
	
		names1_sequence=[]
		names2_sequence=[]
	
		for letter in sequence:
			if letter in 'RYSWKM': #if its a het
				print i
				alleles=IUPAC(letter)
				random.shuffle(alleles)
				names1_sequence.append(alleles[0])
				names2_sequence.append(alleles[1])
			else:
				#print letter
				names1_sequence.append(letter)
				names2_sequence.append(letter)
	
		new_dataset[names1]=names1_sequence
		new_dataset[names2]=names2_sequence
	
new_dataset.to_csv(output,sep="\t",index=False,compression="gzip")
			
	
			




