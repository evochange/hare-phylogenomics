#!/usr/bin/python
# count number of singletons and number of shared hetz amoung individuals in a geno file
# the output is a file with five columns:
# column 1 = no. private hets
# column 2 = no. non private hets
# column 3 = no. private hets /total no. hets
# column 4 = no. non private hets / total no. hets
# column 5 = total no. hets
# Mafalda S. Ferreira, Missoula, MT, Dec 6 2018

import sys
import pandas as pd
from collections import Counter
from decimal import *

#change the number of decimals to present to six
getcontext().prec= 6

#define input and output
input=sys.argv[1]
popfile=sys.argv[2]
output=sys.argv[3]

# loop through each column (individual) a determine how many het sites are shared 
# with other individuals and how many are singletons
def make_sp_dict(poptable):
	poptable=open(poptable,"r")
	mydict={}
	
	for line in poptable:
		line=line.strip().split()
		indv=line[1]
		species=line[0]
		
		mydict.setdefault(species,[]).append(indv)
		
	return mydict

def make_indv_dict(poptable):
	poptable=open(poptable,"r")
	mydict={}
	
	for line in poptable:
		line=line.strip().split()
		indv=line[1]
		species=line[0]
		
		mydict.setdefault(indv,species)
		
	return mydict
	
	
def what_are_singletons(slice,name,df,b):
	info=[]
	singletons=0
	nonsingletons=0
	
	length=len(slice)
	
	for i in range(length):
		letter=slice.loc[i]
		if letter in 'RYSWKM': #if its het
			seq=''.join(df.iloc[i,2:].tolist()) #check the position
			cnt=Counter(seq)
			if cnt[letter]>1: #if there are other hets
				nonsingletons+=1
			else:
				singletons+=1
				
	info.extend([name,str(singletons),str(nonsingletons),str(Decimal(singletons)/Decimal(b)),str(Decimal(nonsingletons)/Decimal(b)),str(b)])
	return info
	
print "Reading input file %s..." % (input)
#read input file into dataframe
dataset=pd.read_csv(input,delimiter="\t",dtype={'#CHROM':object})
#retrieve column names
names=list(dataset)[2:len(list(dataset))]
#retrieve length of alignment
b=len(dataset)

indvdict=make_indv_dict(popfile)
speciesdict=make_sp_dict(popfile)


print "\nRemoving non-het sites from dataframe..."
#restrict alignment to just het sites:
hetz_dataframe=pd.DataFrame(columns=list(dataset))

for i in range(len(dataset)):
	sequence=''.join(dataset.iloc[i,2:].tolist())
	cnt=Counter(sequence)
	hetz=cnt['R']+cnt['Y']+cnt['S']+cnt['W']+cnt['K']+cnt['M']
	if hetz!=0:
		hetz_dataframe=hetz_dataframe.append(dataset.iloc[i,:], ignore_index=True)

print "\nRetreaving heterozigosity information..."
# do the math:
outputfile=open(output,"w")


for h in range(2,len(list(hetz_dataframe)),1):
	n=h-2
	name=names[n]
	#print name
	sp=indvdict.get(name)
	indvs=speciesdict.get(sp)
	#print name
	#print sp
	#print indvs
	
	if len(indvs)==1:
		#print name
		slice=hetz_dataframe.loc[:,name]
		info=what_are_singletons(slice,name,hetz_dataframe,b)
		print "\nWriting info of sample %s to output" % (name)
		outputfile.write('\t'.join(info))
		outputfile.write('\n')
	else:
		namestokeep=[]
		namestokeep.append("#CHROM")
		namestokeep.append("POS")
		for item in names:
			if item in indvs:
				pass
			else:
				namestokeep.append(item)
		namestokeep.append(name)
		hetz_dataframe_filt=hetz_dataframe[namestokeep]
		slice=hetz_dataframe_filt.loc[:,name]
		info=what_are_singletons(slice,name,hetz_dataframe_filt,b)
		print "\nWriting info of sample %s to output" % (name)
		outputfile.write('\t'.join(info))
		outputfile.write('\n')		

outputfile.close()
print "\nDone"
	
