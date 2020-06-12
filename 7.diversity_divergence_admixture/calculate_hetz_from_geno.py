#!/usr/bin/python
# count number of heterozygout sites in a geno file, by individual.
# Mafalda S. Ferreira, Missoula, MT, Dec 4 2018

import sys
import pandas as pd
from collections import Counter
from decimal import *

getcontext().prec= 6

input=sys.argv[1]
output=sys.argv[2]

def sequence_info(sequence,name):
	info=[]
	t_length=len(sequence)
	cnt=Counter(sequence)
	#print t_length
	#for letter in 'ACTGRYSWKMN':
	
	hetz=cnt['R']+cnt['Y']+cnt['S']+cnt['W']+cnt['K']+cnt['M']
	missing=cnt['N']
	
	p_hetz=Decimal(hetz)/Decimal(t_length)
	p_missing=Decimal(missing)/Decimal(t_length)
	
	info.extend([name,str(t_length),str(hetz),str(Decimal(hetz)/Decimal(t_length)),str(missing),str(Decimal(missing)/Decimal(t_length))])
	
	return info

print "Reading input file %s..." % (input)
dataset=pd.read_csv(input,delimiter="\t",dtype={'#CHROM':object})
names=list(dataset)[2:len(list(dataset))]

outpufile=open(output,"w")
outpufile.write('name\ttotal length\tno het sites\tprop of het sites\tno of missing sites\tprop missing sites\n')

print "Retreaving information from input file..."
for i in range(2,len(list(dataset)),1):
	n=i-2
	name=names[n]
	sequence=''.join(dataset.iloc[0:,i].tolist())
	info=sequence_info(sequence,name)
	print "Writing info of sample %s to output" % (name)
	outpufile.write('\t'.join(info))
	outpufile.write('\n')
	
print "Done!"
outpufile.close()