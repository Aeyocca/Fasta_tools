#! /usr/bin/perl -w

# gen_fasta_aey_hpcc.pl
# take GFFs and genomes and make fastas, this version is local
# Alan E. Yocca

use strict;
use Getopt::Std;

my $usage = "\n$0 -c <col0_gff> -g <acc_gff> -a <genome_a> -b <genome_b> -f <fasta_output_a> -o <fasta_output_b> \n\n";

our ($opt_c, $opt_g, $opt_a, $opt_b, $opt_f, $opt_o);
getopts('c:g:a:b:f:o:') or die "$usage";

if ( (!(defined $opt_c)) || (!(defined $opt_g)) || (!(defined $opt_a)) || (!(defined $opt_b)) || (!(defined $opt_f)) || (!(defined $opt_o)) ) {
  print "$usage";
  exit;
}

open (my $col_gff_fh, '<', $opt_c) || die "Cannot open the col0 GFF file $opt_c\n\n";
open (my $acc_gff_fh, '<', $opt_g) || die "Cannot open the accession GFF file $opt_g\n\n";
open (my $ga_fh, '<', $opt_a) || die "Cannot open genome 1: $opt_a\n\n";
open (my $gb_fh, '<', $opt_b) || die "Cannot open genome 2: $opt_b\n\n";
open (my $outa_fh, '>', $opt_f) || die "Cannot open output fasta 1: $opt_f\n\n";
open (my $outb_fh, '>', $opt_o) || die "Cannot open output fasta 2: $opt_o\n\n";

#loop through GFF, get lowest / highest start / end coordinate from each gene model based on CDS feature

my %col0_model;
my %col0_up;
my %col0_down;

while (my $line = <$col_gff_fh>) {
	chomp $line;
	my @info = split("\t",$line); #chromosome source feature 5' 3' score +/- frame attribute
	my @first = split(",", $info[8]);
	my @gene = split("=", $first[0]);
	if ($info[2] eq "CDS") {
		if ($col0_model{$gene[1]}) {
			if ($info[3] <= $col0_down{$gene[1]}) {
				$col0_down{$gene[1]} = $info[3];
			}
			if ($info[3] >= $col0_up{$gene[1]}) {
				$col0_up{$gene[1]} = $info[3];
			}
			if ($info[4] <= $col0_down{$gene[1]}) {
				$col0_down{$gene[1]} = $info[4];
			}
			if ($info[4] >= $col0_up{$gene[1]}) {
				$col0_up{$gene[1]} = $info[4];
			}
		}
		elsif ($info[2] && $gene[1]) {
			$col0_model{$gene[1]} = 1;
			$col0_down{$gene[1]} = $info[3];
			$col0_up{$gene[1]} = $info[4];
		}
	}
	else {
		next;
	}
}

my %acc_model;
my %acc_up;
my %acc_down;

#ID=CDS:AT1G01040-Can_0-mGene.NGS.1;Parent=Transcript:AT1G01040-Can_0-mGene.NGS.1

while (my $line = <$acc_gff_fh>) {
	chomp $line;
	my @info = split("\t",$line); #chromosome source feature 5' 3' score +/- frame attribute
	my @first = split(";", $info[8]);
	my @gene = split(":", $first[0]);
	if ($info[2] eq "CDS") {
		if ($acc_model{$gene[1]}) {
			if ($info[3] <= $acc_down{$gene[1]}) {
				$acc_down{$gene[1]} = $info[3];
			}
			if ($info[3] >= $acc_up{$gene[1]}) {
				$acc_up{$gene[1]} = $info[3];
			}
			if ($info[4] <= $acc_down{$gene[1]}) {
				$acc_down{$gene[1]} = $info[4];
			}
			if ($info[4] >= $acc_up{$gene[1]}) {
				$acc_up{$gene[1]} = $info[4];
			}
		}
		elsif ($info[2] && $gene[1]) {
			$acc_model{$gene[1]} = 1;
			$acc_down{$gene[1]} = $info[3];
			$acc_up{$gene[1]} = $info[4];
		}
	}
	else {
		next;
	}
}

foreach my $key (keys %col0_model) {
	if ($col0_up{$key} >= $col0_down{$key}) {
		open (SAMTOOLSA,"/opt/software/SAMTools/1.2--GCC-4.4.5/bin/samtools faidx /mnt/home/yoccaala/04_Edger/athal_cns/$opt_a $key:$col0_down{$key}-$col0_up{$key} |");
		print "$col0_up{$key}\t$col0_down{$key}\t$key\n";
		while (my $samtool = <SAMTOOLSA>) {
			#my $info = $a_id;
			#$samtool =~ s/Chr$a_chrom/$info/;
			print $outa_fh "$samtool";
		}
		close SAMTOOLSA;
	}
	if ($col0_down{$key} >= $col0_up{$key}) {
		open (SAMTOOLSA,"/opt/software/SAMTools/1.2--GCC-4.4.5/bin/samtools faidx /mnt/home/yoccaala/04_Edger/athal_cns/$opt_a $key:$col0_up{$key}-$col0_down{$key} |");
		print "$col0_up{$key}\t$col0_down{$key}\t$key\n";
		while (my $samtool = <SAMTOOLSA>) {
			#my $info = $a_id;
			#$samtool =~ s/Chr$a_chrom/$info/;
			print $outa_fh "$samtool";
		}
		close SAMTOOLSA;
	}
}

foreach my $key (keys %acc_model) {
	if ($acc_up{$key} >= $acc_down{$key}) {
		open (SAMTOOLSA,"/opt/software/SAMTools/1.2--GCC-4.4.5/bin/samtools faidx /mnt/home/yoccaala/04_Edger/athal_cns/$opt_a $key:$acc_down{$key}-$acc_up{$key} |");
		print "$acc_up{$key}\t$acc_down{$key}\t$key\n";
		while (my $samtool = <SAMTOOLSA>) {
			#my $info = $a_id;
			#$samtool =~ s/Chr$a_chrom/$info/;
			print $outb_fh "$samtool";
		}
		close SAMTOOLSA;
	}
	if ($acc_down{$key} >= $acc_up{$key}) {
		open (SAMTOOLSA,"/opt/software/SAMTools/1.2--GCC-4.4.5/bin/samtools faidx /mnt/home/yoccaala/04_Edger/athal_cns/$opt_a $key:$acc_up{$key}-$acc_down{$key} |");
		print "$acc_up{$key}\t$acc_down{$key}\t$key\n";
		while (my $samtool = <SAMTOOLSA>) {
			#my $info = $a_id;
			#$samtool =~ s/Chr$a_chrom/$info/;
			print $outb_fh "$samtool";
		}
		close SAMTOOLSA;
	}
}

close $outa_fh;
close $outb_fh;
close $ga_fh;
close $gb_fh;
close $col_gff_fh;
close $acc_gff_fh;

exit;
