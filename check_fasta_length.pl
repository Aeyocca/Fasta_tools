#! /usr/bin/perl -w

# check_fasta_length.pl
# 03-08-18
# Alan E. Yocca

use strict;
use Getopt::Std;

my $usage = "\n$0 -c <coge fasta > -l <mine fasta (one line)> \n\n";

our ($opt_c, $opt_l);
getopts('c:l:') or die "$usage";

if ( (!(defined $opt_c)) || (!(defined $opt_l)) ) {
  print "$usage";
  exit;
}

if (!(-e $opt_c)) {
  die "\nThe alpha dup file, $opt_c, does not exist.\n\n";
}
open (my $coge_fh, '<', $opt_c) || die "Cannot open the input file $opt_c\n\n";
open (my $mine_fh, '<', $opt_l) || die "Cannot open $opt_l\n\n";

#>1||11864||12940||AT1G01030.1||-1||CDS||306206338||4

my %coge_length;
my %mine_length;
my $gene_coge;

while (my $line = <$coge_fh>) {
	chomp $line;
	if ($line =~ m/^>/){
		my @info = split(/\|\|/, $line);
		$gene_coge = $info[3];
	}
	else {
		my @sequence = split('', $line);
		my $length = @sequence;
		$coge_length{$gene_coge} = $length;
	}
}

my $gene_mine;

while (my $line = <$mine_fh>) {
	chomp $line;
	if ($line =~ m/^>/){
		my @info = split(/\|\|/, $line);
		$gene_mine = $info[3];
	}
	else {
		my @sequence = split('', $line);
		my $length = @sequence;
		$mine_length{$gene_mine} = $length;
	}
}

my $mismatch_gene = 0;
my $diff_length = 0;

foreach my $key (keys %mine_length) {
	if ($coge_length{$key}) {
		if ($coge_length{$key} == $mine_length{$key}) {
			next;
		}
		else {
			$diff_length = $diff_length + 1;
			#print "CULPRIT: $key\nCoge length:$coge_length{$key}\nMine length:$mine_length{$key}\n";
		}
	}
	else {
		$mismatch_gene = $mismatch_gene + 1;
		print "$key\n";
		last;
	}
}
print "cant find gene: $mismatch_gene\nGenes with different length: $diff_length\n";


close $coge_fh;
close $mine_fh;

exit;