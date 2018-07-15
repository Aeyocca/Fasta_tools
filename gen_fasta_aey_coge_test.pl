#! /usr/bin/perl -w

# gen_fasta_aey_coge_test.pl
# take GFFs and genomes and make fastas, this version is local
# Alan E. Yocca

use strict;
use Getopt::Std;

my $usage = "\n$0 -c <col0_gff> -a <genome_a> -o <fasta_output> \n\n";

our ($opt_c, $opt_g, $opt_a, $opt_b, $opt_f, $opt_o);
getopts('c:a:o:') or die "$usage";

if ( (!(defined $opt_c)) || (!(defined $opt_a)) || (!(defined $opt_o)) ) {
  print "$usage";
  exit;
}

#if (-e $opt_o || -e $opt_f) {
#print "File: $opt_o or exist, is it okay to overwrite them?\n"; #$opt_f exist, is it okay to overwrite them?\n";
#my $answer = <STDIN>;
#	if ($answer =~ /^y(?:es)?$/i) {
#		print "Excellent\n";
#	}
#	else {
#		die "fine, I will not overwrite your files\n";
#	}
#}

open (my $col_gff_fh, '<', $opt_c) || die "Cannot open the col0 GFF file $opt_c\n\n";
#open (my $acc_gff_fh, '<', $opt_g) || die "Cannot open the accession GFF file $opt_g\n\n";
#open (my $ga_fh, '<', $opt_a) || die "Cannot open genome 1: $opt_a\n\n";
#open (my $gb_fh, '<', $opt_b) || die "Cannot open genome 2: $opt_b\n\n";
open (my $temp_fh, '>', "temp.txt") || die "Cannot open temp\n\n";
open (my $outa_fh, '>', $opt_o) || die "Cannot open output fasta 1: $opt_o\n\n";
#open (my $out_fh, '>', $opt_o) || die "Cannot open output fasta 2: $opt_o\n\n";

#loop through GFF, get lowest / highest start / end coordinate from each gene model based on CDS feature

my @col0_gff;

while (my $line = <$col_gff_fh>) {
	chomp $line;
	my @info = split("\t",$line);
	if ($info[2]) {
		if ($info[0] ne "C" && $info[0] ne "M" && $info[2] eq "CDS") {
			#if ($info[6] eq "-") {
			#	my $temp = $info[3];
			#	$info[3] = $info[4];
			#	$info[4] = $temp;
			#}
			#add $col0_gff[$f][9] containing the decimal number for a further sort
			my @first = split(",", $info[8]);
			my @second = split("Alias", $first[0]);
			my @gene = split("=", $second[1]);
			my @decimal = split(/\./, $gene[1]);
			my @coge_id = split("coge_fid=",$info[8]);
			push (@info, $decimal[1], $coge_id[1]);
			push (@col0_gff, \@info);
		}
		elsif ($info[2] eq "CNS") {
			my @decimal = ("nonsense", 1);
			my @coge_id = ("nonsense", 69);
			push (@info, $decimal[1], $coge_id[1]);
			push (@col0_gff, \@info);
		}
		elsif ($info[0] ne "C" && $info[0] ne "M" && $info[2] =~ m/[rt]RNA/) {
		#	if ($info[6] eq "-") {
		#		my $temp = $info[3];
		#		$info[3] = $info[4];
		#		$info[4] = $temp;
		#	}
			#add $col0_gff[$f][9] containing the decimal number for a further sort
			my @first = split("Alias", $info[8]);
			my @second = split(",", $first[1]);
			my @gene = split(";", $second[1]);
			my @decimal = split(/\./, $gene[0]);
			my @coge_id = split("coge_fid=",$info[8]);
			push (@info, $decimal[1], $coge_id[1]);
			push (@col0_gff, \@info);
		}
		else {
			next;
		}
	}
}

my @sorted_col0_gff = sort {
	$a->[0] cmp $b->[0] ||
	$a->[3] <=> $b->[3] ||
	$a->[9] <=> $b->[9];
} @col0_gff;

my %seen_it_col0;
my %col0_high;
my %col0_low;
my %col0_chrom;
my %go_col0;
my $count = 1;
my %hash_for_array;
my %HoA;
	
#fill variables
#change this loop to only assign gene order, sorted col0 gff has all the info we need to do it by cds
#also rewrite gene name so can sort by it
#add in optimization for start / end so can add that metainfo to fasta header
	
for (my $f = 0 ; $f < @sorted_col0_gff ; $f++) {
	my @gene;
	if ($sorted_col0_gff[$f][2] eq "CDS") {
		my @first = split(",", $sorted_col0_gff[$f][8]);
		my @second = split("Alias", $first[0]);
		@gene = split("=", $second[1]);
		$sorted_col0_gff[$f][8] = $gene[1];
	}
	elsif ($sorted_col0_gff[$f][2] eq "CNS") {
		my @first = split(";", $sorted_col0_gff[$f][8]);
		my @second = split("Alias", $first[3]);
		@gene = split("=", $second[1]);
		$sorted_col0_gff[$f][8] = $gene[1];
	}
	else {
		my @first = split("Alias", $sorted_col0_gff[$f][8]);
		my @second = split(",", $first[1]);
		@gene = split(";", $second[1]);
		$gene[1] = $gene[0];
		$sorted_col0_gff[$f][8] = $gene[1];
	}
	if ($seen_it_col0{$gene[1]}) {
			if ($sorted_col0_gff[$f][3] <= $col0_low{$gene[1]}) {
				$col0_low{$gene[1]} = $sorted_col0_gff[$f][3];
			}
			if ($sorted_col0_gff[$f][4] >= $col0_high{$gene[1]}) {
				$col0_high{$gene[1]} = $sorted_col0_gff[$f][4];
			}
	}
	else {
		if ($sorted_col0_gff[$f][2] eq "CNS") {
			$seen_it_col0{$gene[1]} = 1;
#			$go_col0{$gene[1]} = $count; #don't need to call on gene order here, going to look up gene order of associated gene later
#			$count = $count + 1; #leave out to keep GO to the gene it is closest to
			$col0_low{$gene[1]} = $sorted_col0_gff[$f][3];
			$col0_high{$gene[1]} = $sorted_col0_gff[$f][4];
		}
		else {
			$seen_it_col0{$gene[1]} = 1;
			$go_col0{$gene[1]} = $count;
			$count = $count + 1;
			$col0_low{$gene[1]} = $sorted_col0_gff[$f][3];
			$col0_high{$gene[1]} = $sorted_col0_gff[$f][4];
		}
	}	
}

my $undef_cds = 1;
my $undef_rna = 1;

my @sorted_col0_for_samtools = sort {
	$a->[8] cmp $b->[8] ||
	$a->[3] <=> $b->[3];
} @sorted_col0_gff;

foreach my $key (keys %seen_it_col0) {
	delete $seen_it_col0{$key};
}

my $first_in_file = 1;
my $sam_tool_count = 0;

#for (my $f = 0 ; $f < @sorted_col0_for_samtools ; $f++) {
#	my $gene;
#	$gene = $sorted_col0_for_samtools[$f][8];
#	my $strand = 1;
#	if ($sorted_col0_for_samtools[$f][6] eq "-") {
#		$strand = -1;
#	}
#	my $chromosome = "Chr" . $sorted_col0_for_samtools[$f][0];
#	if ($gene) {
#		for (my $g = 0 ; $g < @{$sorted_col0_for_samtools[$f]} ; $g++) {
#			print $outa_fh "$sorted_col0_for_samtools[$f][$g]\t";
#		}
#		print $outa_fh "\n";
#	}
#	else {
#		for (my $g = 0 ; $g < @{$sorted_col0_for_samtools[$f]} ; $g++) {
#			print $outa_fh "$sorted_col0_for_samtools[$f][$g]\t";
#		}
#		print $outa_fh "\n";
#	}
	#	print $outa_fh "$sorted_col0_for_samtools[$f][0]||$col0_low{$gene}||$col0_high{$gene}||$gene||$strand||$sorted_col0_for_samtools[$f][2]||$sorted_col0_for_samtools[$f][10]||$go_col0{$gene}/n";

#}


#*******find lowest gene order of associated isoforms**********
#my $gene_test = "AT1G01020";
#my @isoforms = grep {/$gene_test/} keys %go_col0;
#my $lowest;
#MIN: foreach my $isoform (@isoforms) {
#	if ($lowest) { 
#		if ($go_col0{$isoform} <= $go_col0{$lowest}) {
#			$lowest = $go_col0{$isoform};
#		}
#		else {
#			next MIN;
#		}
#	}
#	else {
#		$lowest = $go_col0{$isoform};
#	}
#}
#print "$lowest\n";

my $no_gene = 0;
	
for (my $f = 0 ; $f < @sorted_col0_for_samtools ; $f++) {
	my $gene;
	my $gene_go;
	if ($sorted_col0_for_samtools[$f][2] eq "CNS") {
		$gene = $sorted_col0_for_samtools[$f][8];
#		if ($gene) {
		#	print "$sorted_col0_for_samtools[$f][2]\t$sorted_col0_for_samtools[$f][8]\t$gene\n";
#		}
#		else {
#			print "$sorted_col0_for_samtools[$f][8]\n";
#			$no_gene = $no_gene + 1;	
#		}
		my ($number, $easier_temp_than_work_around) = split("_",$gene);
		my $gene_base = $easier_temp_than_work_around;
		my @isoforms = grep {/$gene_base/} keys %go_col0;
#		if ($isoforms[0]) {
#			next;
#		}
#		else {
#			print $outa_fh "$gene\n";
#			$no_gene = $no_gene + 1;
#		}
MIN: 	foreach my $isoform (@isoforms) {
			if ($gene_go) { 
				if ($go_col0{$isoform} <= $go_col0{$gene_go}) {
					$gene_go = $isoform;
				}
				else {
					next MIN;
				}
			}
			else {
				$gene_go = $isoform;
			}
		}
	}
	else {
#		my @first = split("Alias", $sorted_col0_for_samtools[$f][8]);
#		my @second = split(",", $first[1]);
#		@gene = split(";", $second[1]);
#		$gene[1] = $gene[0];
		$gene = $sorted_col0_for_samtools[$f][8];
		my $int = $gene;
		$gene_go = $int;		
	}
	my $strand = 1;
	if ($sorted_col0_for_samtools[$f][6] eq "-") {
		$strand = -1;
	}
	my $chromosome = "Chr" . $sorted_col0_for_samtools[$f][0];	
	if ($seen_it_col0{$gene}){
		open (SAMTOOLSA,"/opt/software/SAMTools/1.2--GCC-4.4.5/bin/samtools faidx /mnt/home/yoccaala/04_Edger/athal_cns/$opt_a $chromosome:$sorted_col0_for_samtools[$f][3]-$sorted_col0_for_samtools[$f][4] |");
		while (my $samtool = <SAMTOOLSA>) {
			chomp $samtool;
			if ($samtool =~ m/^>/){
				next;
			}
			else {
				print $outa_fh "$samtool";
			}
		}
		close SAMTOOLSA;
	}
#	if ($sorted_col0_for_samtools[$f][0] && $sorted_col0_for_samtools[$f][2] ) {
#		print $outa_fh "All Defined\n";
#	}
#	else {
#		print $outa_fh "something undefined\n";
#	}
#print $outa_fh "$sorted_col0_for_samtools[$f][0]||$col0_low{$gene}||$col0_high{$gene}||$gene||$strand||$sorted_col0_for_samtools[$f][2]||$sorted_col0_for_samtools[$f][10]||$go_col0{$gene}\n"
	else {
		$seen_it_col0{$gene} = 1;
		open (SAMTOOLSA,"/opt/software/SAMTools/1.2--GCC-4.4.5/bin/samtools faidx /mnt/home/yoccaala/04_Edger/athal_cns/$opt_a $chromosome:$sorted_col0_for_samtools[$f][3]-$sorted_col0_for_samtools[$f][4] |");
		while (my $samtool = <SAMTOOLSA>) {
			chomp $samtool;
			#$sam_tool_count = $sam_tool_count + 1;
			#print "$samtool\t$sam_tool_count\n";
			if ($samtool =~ m/^>/) {
				$samtool =~ s/$chromosome:$sorted_col0_for_samtools[$f][3]-$sorted_col0_for_samtools[$f][4]/$sorted_col0_for_samtools[$f][0]||$col0_low{$gene}||$col0_high{$gene}||$gene||$strand||$sorted_col0_for_samtools[$f][2]||$sorted_col0_for_samtools[$f][10]||$go_col0{$gene_go}/;
				if ($first_in_file) {
					print $outa_fh "$samtool\n";
					$first_in_file = 0;
				}
				else {
					print $outa_fh "\n$samtool\n";
				}
			}
			else {
				print $outa_fh "$samtool";
			}
		}
		close SAMTOOLSA;
	}
}

#print "tried to remove repeats: $sam_tool_count\n";

#close $temp_fh;

#open (my $tempr_fh, '<', "temp.txt") || die "Cannot open temp a second time\n\n";

#while (my $line = <$tempr_fh>) {
#	chomp $line;
#	if ($line =~m/^>/) {
#		my @info = split(/\|\|/,$line);
#		$seen_it_col0{$info[3]} = 1;
#	}
#	else {
#	
#	}
#}
#print "$no_gene\n";

close $outa_fh;
close $col_gff_fh;

exit;
