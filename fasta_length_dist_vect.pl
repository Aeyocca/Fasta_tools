#! /usr/bin/perl -w

# fasta_length_dist_vect.pl
# take fasta file, output a comma separated vector (list??, its just one line) of the lengths of each fasta entry
# think of adding an ignore N switch at some point
# WILL NOT DO ANY CHECKS, simply get number of characters in fasta entry
# should make it work on wrapped fasta
# 06-28-18
# Alan E. Yocca

use strict;
use Getopt::Std;

my $usage = "\n$0 -f <input fasta> -o <output> \n\n";

our ($opt_f, $opt_o);
getopts('f:o:') or die "$usage";

if ( (!(defined $opt_f)) || (!(defined $opt_o)) ) {
  print "$usage";
  exit;
}

open (my $fasta_fh, '<', $opt_f) || die "Cannot open the fasta file: $opt_f\n\n";
open (my $out_fh, '>', $opt_o) || die "Cannot open the output file: $opt_o\n\n";

#my $unambig = 0;
#my $ns = 0;
my $fasta_entries = 0; #just a check to see if working properly
my $length = 0;

while (my $line = <$fasta_fh>) {
	chomp $line;
	if ($line =~ m/^>/) {
		if ($fasta_entries != 0) { #skips first fasta entry, so have no length
			print $out_fh "$length,";
		}
		$fasta_entries = $fasta_entries + 1;
		$length = 0;
		next;
	}
	else {
		#these all from other script used as template for this one
#		my $loop_unambig = $line =~ tr/actgACTG//;
#		my $loop_ns = $line =~ tr/nN//;
#		$unambig = $unambig + $loop_unambig;
#		$ns = $ns + $loop_ns;
		$length = $length + length($line);
	}
}

#get the last one
print $out_fh "$length";

#stdout summary
print "fasta entries:\t$fasta_entries\n";

close $fasta_fh;
close $out_fh;

exit;
