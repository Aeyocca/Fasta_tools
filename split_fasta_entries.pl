#! /usr/bin/perl -w

# split_fasta_entries.pl
# take fasta file and split large fasta entries into smaller fragments size selected on command line, if not divisible by specified number, last entry will be the truncated one
# Alan E. Yocca
# 05-04-18

use strict;
use Getopt::Std;

my $usage = "\n$0 -f <input fasta> -k <# of nucleotides to split by> -o <output fasta> \n\n";

our ($opt_f, $opt_k, $opt_o);
getopts('f:k:o:') or die "$usage";

if ( (!(defined $opt_f)) || (!(defined $opt_k)) || (!(defined $opt_o)) ) {
  print "$usage";
  exit;
}

open (my $fasta_fh, '<', $opt_f) || die "Cannot open the fasta file: $opt_f\n\n";
open (my $out_fh, '>', $opt_o) || die "Cannot open the output file: $opt_o\n\n";

my $header;

while (my $line = <$fasta_fh>) {
	chomp $line;
	my @trans = split(">",$line);
	if ($trans[1]) {
		$header = $line;
	}
	else {
		$line =~ s/.{$opt_k}\K(?=.)/ /sg;
#		$line =~ s/\s$//;
#		print $out_fh "$line\n";
		my @splitup = split(/\s/,$line);
		for (my $i=0; $i < @splitup; $i++) {
			my $split_header = $header . "_" . $i;
			print $out_fh "$split_header\n";
			print $out_fh "$splitup[$i]\n";
		}
	}
}

close $out_fh;
close $fasta_fh;

exit;