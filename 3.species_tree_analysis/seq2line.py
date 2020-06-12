#!/usr/bin/python
# Converts multiple line seq from fasta files to single line seqs
# Copyright (c) 2011, Jose Melo Ferreira (CIBIO, UP), 30NOV2011
# USAGE: python seq2line.py fasta_file > new_fasta

from Bio import SeqIO
import sys

try: #If no or insufficient args are passed prints USAGE message and exits

#Opens and parses the fasta file using Biopython
	handle_seqs = open(sys.argv[1], 'rU')
	parsed_fasta = SeqIO.parse(handle_seqs, 'fasta')

except:
	print "Converts multiple line seqs from fasta files to single line seqs"
	print "USAGE: python seq2line.py fasta_file > new_fasta"
	raise SystemExit

#Prints fasta with sequences as a single line
for seq_record in parsed_fasta:	
	print '>'+seq_record.id
	print seq_record.seq

sys.exit()
handle_seqs.close()

