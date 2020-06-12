#!/usr/bin/perl

use strict; # help 

## my input parameters 
my $vcf = $ARGV[0];
my $chrlen = $ARGV[1];


## define nucleotide codes
my %nucleotides = (
	"AA" => "A",
	"TT" => "T",
	"GG" => "G",
	"CC" => "C",
	"AG" => "R",
	"CT" => "Y",
	"GC" => "S",
	"AT" => "W",
	"GT" => "K",
	"AC" => "M",
	"GA" => "R",
	"TC" => "Y",
	"CG" => "S",
	"TA" => "W",
	"TG" => "K",
	"CA" => "M",
	".." => "N"
);

## define raw fasta sequence
my @fasta = ("N") x $chrlen;

my $chr;
## read VCF and replace fasta sequence
open (vcf, '<', $vcf) or die "can't open $vcf $!\n";
while (my $line = <vcf>) {
    chomp($line);
    # ignone comment lines
    if ($line =~ /^\#\#/) { next; }
	#
	else {
		my @cols = split("\t", $line);
		$chr = $cols[0];
		my $pos = $cols[1];
		my $gen = $cols[2];
		my $rep = $nucleotides{$gen};
		my $seqpos = $pos - 1; 
		$fasta[$seqpos] = $rep;
	}
}
	
print ">", $chr,"\n";
print join("", @fasta),"\n";

