#! /usr/bin/perl -w
# fasta_one_line.pl
# convert wrapped fasta file to single line per entry
# 06-01-18
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

my $header;
my $sequence;

while (my $line = <$fasta_fh>) {
        chomp $line;
        if ($line =~ m/^>/) {
        		if ($header){ #ran into a new fasta header so can write out the entry
	        		print $out_fh "$header\n$sequence\n";
    			}
                $header = $line;
                $sequence = ""; #reset sequence
                next;
        }
        else {
			$sequence = $sequence . $line
		}
}
#get that last sequence to print out
print $out_fh "$header\n$sequence\n";

close $fasta_fh;
close $out_fh;

exit;