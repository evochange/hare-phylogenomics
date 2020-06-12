#!/usr/bin/python

#converting a fasta alignment into two sequences per individual to split
#heterozygout alleles
#3 Abril 2018, Mafalda S. Ferreira
#usage: fake_phase_fasta_algn.py fasta_algn.fa output_fasta_algn.fa

import sys
import Bio
from Bio import AlignIO


### read inputs ###
fasta_align=sys.argv[1]
output_fasta=sys.argv[2]

### make a dictionary with the iupac codes for my HETS ####
dictionary_iupac={}
dictionary_iupac.setdefault("R",["A","G"])
dictionary_iupac.setdefault("Y",["C","T"])
dictionary_iupac.setdefault("S",["G","C"])
dictionary_iupac.setdefault("W",["A","T"])
dictionary_iupac.setdefault("K",["G","T"])
dictionary_iupac.setdefault("M",["A","C"])

### read the alignment ###
align=AlignIO.read(fasta_align,"fasta")
length_align=align.get_alignment_length()
number_inds=len(align)

### open output ###
new_fasta=open(output_fasta,"w")

### convert my hets and write them to different sequences ###
for record in align:
	seq_a = list()
	seq_b = list()
	for w in xrange(length_align):
		base=record.seq[w]
		if base in ['B', 'D', 'H', 'V']:
			print "Adding N to site {0} - you have a triallelic SNP in {1}".format(w,record.id)
			seq_a.append('N')
			seq_b.append('N')
		elif base not in dictionary_iupac:
			seq_a.append(base)
			seq_b.append(base)
		else:
			nuc=dictionary_iupac.get(base)
			seq_a.append(nuc[0])
			seq_b.append(nuc[1])
	
	new_header_a=">{0}_a\n".format(record.id)
	new_header_b=">{0}_b\n".format(record.id)
	new_fasta.write("".join(new_header_a))
	new_fasta.write("".join(seq_a))
	new_fasta.write("\n")
	new_fasta.write("".join(new_header_b))
	new_fasta.write("".join(seq_b))
	new_fasta.write("\n")
	
new_fasta.close()
	

