#! /usr/bin/perl -w

# check concordance.pl
# Alan E. Yocca
# 02-28-18

use strict;
use Getopt::Std;

my $usage = "\n$0 -c <coge fasta headers> -l <my version fasta headers> -o <output check> \n\n";

our ($opt_l, $opt_c, $opt_o);
getopts('l:c:o:') or die "$usage";

if ( (!(defined $opt_l)) || (!(defined $opt_c)) || (!(defined $opt_o)) ) {
  print "$usage";
  exit;
}

open (my $coge_fh, '<', $opt_c) || die "Cannot open coge $opt_c\n\n";
open (my $mine_fh, '<', $opt_l) || die "Cannot open mine: $opt_l\n\n";
open (my $out_fh, '>', $opt_o) || die "Cannot open out: $opt_o\n\n";

# Chr2    AT2G20100.2     10841   10841   Chr5    AT5G57440       72244   72244   0.019
# a16911_1        1||25069727||25095526||AT1G67120.1||-1||CDS||306229466||7593||99.85     7593    7593    b37385_1        1||24616810||24642619||CDS:AT1G67120-Can_0-mGene.NGS.1||-1||CDS||1719345002||15142||99.85       15142   15142   0.0


my %coge_gn_col0;
my $no_split = 0;


while (my $line = <$coge_fh>) {
	chomp $line;
	my @gene_col0 = split(/\|\|/, $line);;
	if ($gene_col0[3] && $gene_col0[7]) {
		$coge_gn_col0{$gene_col0[3]} = $gene_col0[7];
	}
	else {
		$no_split = $no_split + 1;
	}
}

print "can't split coge: $no_split\n";

my %mine_gn_col0;

while (my $line = <$mine_fh>) {
	chomp $line;
	my @gene_col0 = split(/\|\|/, $line);;
	if ($gene_col0[3] && $gene_col0[7]) {
		$mine_gn_col0{$gene_col0[3]} = $gene_col0[7];
	}
}

my $match = 0;
my $gene_diff = 0;

foreach my $key (sort { $mine_gn_col0{$a} <=> $mine_gn_col0{$b} } keys %mine_gn_col0) {
	if ($coge_gn_col0{$key}) {
		if ($coge_gn_col0{$key} == $mine_gn_col0{$key}) {
			$match = $match + 1;
		}
		else {
#			print "Col0:\nMatched $match until exit\nCULPRIT == $key\nCoGe number: $coge_gn_col0{$key}\nMy number: $mine_gn_col0{$key}\n\n";
#			last;
			print $out_fh "$key\nCoGe number: $coge_gn_col0{$key}\nMy number: $mine_gn_col0{$key}\n\n";
		}
	}
	else {
		$gene_diff = $gene_diff + 1;
	}
}

my $match_acc = 0;
my $gene_diff_acc = 0;

close $mine_fh;
close $coge_fh;

exit;