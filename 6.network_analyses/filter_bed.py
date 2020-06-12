#!/usr/bin/python
# choose snps from bed file that dist X bps from each other

import sys

snps_file=sys.argv[1]
output_file=sys.argv[2]

X=0

output=open(output_file,"w")

for inline in open(snps_file,"r"):
	line=inline.strip().split()
	if X==0:
		X=int(line[1])
		output.write(inline)
	else:
		if int(line[1])-X >= 10000:
			output.write(inline)
			X=int(line[1])
			
		
	