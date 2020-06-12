#!/usr/bin/python
#A script that will read in a fasta alignment and concatenate the sequences inside.

import sys
import re

input=sys.argv[1]
output=sys.argv[2]


infasta=open(input,"r")
outfasta=open(output,"w")

outfasta.write("".join(">"+re.sub(".agp.fa","",sys.argv[1])+"\n"))

new_fasta=[]

for line in infasta:
	line=line.strip()
		
	if line.startswith(">"):
		pass
	else:
		new_fasta.append(line)
		
		
sequence=''.join(new_fasta)
outfasta.write(sequence+'\n')
