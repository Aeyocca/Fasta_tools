#! /usr/bin/perl -w

# subset_fa.pl
# take list of genes and pull out of fa
# Alan E. Yocca
# 2-19-18

use strict;
use Getopt::Std;

my $usage = "\n$0 -i <input gene list> -f <fa file> -o <output> \n\n";

our ($opt_i, $opt_f, $opt_o);
getopts('i:f:o:') or die "$usage";

if ( (!(defined $opt_i)) || (!(defined $opt_f)) || (!(defined $opt_o)) ) {
  print "$usage";
  exit;
}

open (my $gene_list_fh, '<', $opt_i) || die "Cannot open the gene list: $opt_i\n\n";
open (my $fasta_fh, '<', $opt_f) || die "Cannot open the fasta file: $opt_f\n\n";
open (my $out_fh, '>', $opt_o) || die "Cannot open output: $opt_o\n\n";

my %genes;

while (my $line = <$gene_list_fh>) {
	chomp $line;
	$genes{$line} = 1;
}

my $loop = 0;
my $count = 0;

while (my $line = <$fasta_fh>) {
	chomp $line;
	my @trans = split(">",$line);
	if ($trans[1]) {
		my @line = split(" ",$trans[1]);
		if ($genes{$line[0]}) {
			if ($count == 0) {
				print $out_fh "$line\n";
				$loop = 1;
				$count = $count + 1;
			}
			else {
				print $out_fh "\n$line\n";
				$loop = 1;
				$count = $count + 1;
			}
		}
		else {
			$loop = 0;
		}
	}
	else {
		if ($loop) {
			print $out_fh "$line";
		}
	}
}

print "genes match line: $count\n";

close $gene_list_fh;
close $fasta_fh;
close $out_fh;

exit;