#! /usr/bin/perl -w

# gen_fasta_aey_acc_test.pl
# take GFFs and genomes from accession and create fasta with CDS with headers formatted for CoGe DAGChainer pipe
# Alan E. Yocca
# 03-09-18

use strict;
use Getopt::Std;

my $usage = "\n$0 -g <acc_gff> -c <acc Genome> -o <fasta_output> \n\n";

our ($opt_g, $opt_c, $opt_o);
getopts('g:c:o:') or die "$usage";

if ( (!(defined $opt_g)) || (!(defined $opt_c)) || (!(defined $opt_o)) ) {
  print "$usage";
  exit;
}

if (-e $opt_o ) {
print "File: $opt_o or exist, is it okay to overwrite it?\n"; 
my $answer = <STDIN>;
	if ($answer =~ /^y(?:es)?$/i) {
		print "Excellent!\n";
	}
	else {
		die "fine, I will not overwrite your files\n";
	}
}

open (my $acc_gff_fh, '<', $opt_g) || die "Cannot open the acc GFF file $opt_g\n\n";
open (my $out_fh, '>', $opt_o) || die "Cannot open output fasta: $opt_o\n\n";

#loop through GFF, get lowest / highest start / end coordinate from each gene model based on CDS feature

my @acc_gff;

while (my $line = <$acc_gff_fh>) {
	chomp $line;
	my @info = split("\t",$line);
	if ($info[2]) {
		if ($info[0] ne "C" && $info[0] ne "M" && $info[2] eq "CDS") {
			#add $acc_gff[$f][9] containing the decimal number for a further sort
			my @first = split(";", $info[8]);
			my @gene = split("=", $first[0]);
			my @decimal = split(/\./, $gene[1]);
#			my @coge_id = split("coge_fid=",$info[8]);
			$info[8] = $gene[1];
			push (@info, $decimal[2]);
			push (@acc_gff, \@info);
		}
#		elsif ($info[0] ne "C" && $info[0] ne "M" && $info[2] =~ m/[rt]RNA/) {
			#add $acc_gff[$f][9] containing the decimal number for a further sort
#			my @first = split("Alias", $info[8]);
#			my @second = split(",", $first[1]);
#			my @gene = split(";", $second[1]);
#			my @decimal = split(/\./, $gene[0]);
#			my @coge_id = split("coge_fid=",$info[8]);
#			push (@info, $decimal[1], $coge_id[1]);
#			push (@acc_gff, \@info);
#		}
		else {
			next;
		}
	}
}

my @sorted_acc_gff = sort {
	$a->[0] cmp $b->[0] ||
	$a->[3] <=> $b->[3] ||
	$a->[9] <=> $b->[9];
} @acc_gff;

my %seen_it_acc;
my %acc_high;
my %acc_low;
my %go_acc;
my $count = 1;
	
#fill variables
#change this loop to only assign gene order, sorted acc gff has all the info we need to do it by cds
#also rewrite gene name so can sort by it
#add in optimization for start / end so can add that metainfo to fasta header
	
for (my $f = 0 ; $f < @sorted_acc_gff ; $f++) {
	my $gene;
	if ($sorted_acc_gff[$f][2] eq "CDS") {
#		my @first = split(",", $sorted_acc_gff[$f][8]);
#		my @second = split("Alias", $first[0]);
#		@gene = split("=", $second[1]);
		$gene = $sorted_acc_gff[$f][8];
#		$sorted_acc_gff[$f][8] = $gene[1];
	}
#	else {
#		my @first = split("Alias", $sorted_acc_gff[$f][8]);
#		my @second = split(",", $first[1]);
#		@gene = split(";", $second[1]);
#		$gene[1] = $gene[0];
#		$sorted_acc_gff[$f][8] = $gene[1];
#	}
	if ($seen_it_acc{$gene}) {
			if ($sorted_acc_gff[$f][3] <= $acc_low{$gene}) {
				$acc_low{$gene} = $sorted_acc_gff[$f][3];
			}
			if ($sorted_acc_gff[$f][4] >= $acc_high{$gene}) {
				$acc_high{$gene} = $sorted_acc_gff[$f][4];
			}
	}
	else {
#		my $f_count;
#		if ($count) {
#			$seen_it_acc{$gene} = 1;
#			$f_count = 3*$count + 1;
#		}
#		else {
#			$f_count = 1;
#		}
		$seen_it_acc{$gene} = 1;
		$go_acc{$gene} = $count;
		$count = $count + 1;
		$acc_low{$gene} = $sorted_acc_gff[$f][3];
		$acc_high{$gene} = $sorted_acc_gff[$f][4];
	}	
}

my @sorted_acc_for_samtools = sort {
	$a->[8] cmp $b->[8] ||
	$a->[3] <=> $b->[3];
} @sorted_acc_gff;

foreach my $key (keys %seen_it_acc) {
	delete $seen_it_acc{$key};
}

my $first_in_file = 1;

for (my $f = 0 ; $f < @sorted_acc_for_samtools ; $f++) {
	my $gene;
	if ($sorted_acc_for_samtools[$f][2] eq "CDS") {
		$gene = $sorted_acc_for_samtools[$f][8];
	}
#	else {
#		$gene = $sorted_acc_for_samtools[$f][8];
#	}
	my $strand = 1;
	if ($sorted_acc_for_samtools[$f][6] eq "-") {
		$strand = -1;
	}
#	my $chromosome = "Chr" . $sorted_acc_for_samtools[$f][0];
	(my $chrom_num = $sorted_acc_for_samtools[$f][0]) =~ s/Chr//;
	if ($seen_it_acc{$gene}){
		open (SAMTOOLSA,"/opt/software/SAMTools/1.2--GCC-4.4.5/bin/samtools faidx /mnt/home/yoccaala/04_Edger/athal_cns/$opt_c $sorted_acc_for_samtools[$f][0]:$sorted_acc_for_samtools[$f][3]-$sorted_acc_for_samtools[$f][4] |");
		while (my $samtool = <SAMTOOLSA>) {
			chomp $samtool;
			if ($samtool =~ m/^>/){
				next;
			}
			else {
				print $out_fh "$samtool";
			}
		}
		close SAMTOOLSA;
		next;
	}
	else {
		$seen_it_acc{$gene} = 1;
		open (SAMTOOLSA,"/opt/software/SAMTools/1.2--GCC-4.4.5/bin/samtools faidx /mnt/home/yoccaala/04_Edger/athal_cns/$opt_c $sorted_acc_for_samtools[$f][0]:$sorted_acc_for_samtools[$f][3]-$sorted_acc_for_samtools[$f][4] |");
#		print $out_fh "$sorted_acc_for_samtools[$f][0]:$sorted_acc_for_samtools[$f][3]-$sorted_acc_for_samtools[$f][4]\n";
				while (my $samtool = <SAMTOOLSA>) {
			chomp $samtool;
			if ($samtool =~ m/^>/) {
				$samtool =~ s/$sorted_acc_for_samtools[$f][0]:$sorted_acc_for_samtools[$f][3]-$sorted_acc_for_samtools[$f][4]/$chrom_num||$acc_low{$gene}||$acc_high{$gene}||$gene||$strand||$sorted_acc_for_samtools[$f][2]||69||$go_acc{$gene}/;
				if ($first_in_file) {
					print $out_fh "$samtool\n";
					$first_in_file = 0;
				}
				else {
					print $out_fh "\n$samtool\n";
					next;
				}
			}
			else {
				print $out_fh "$samtool";
			}
		}
		close SAMTOOLSA;
	}
}

close $out_fh;
close $acc_gff_fh;

exit;