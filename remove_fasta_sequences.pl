#! /usr/bin/perl -w

# remove_fasta_sequences.pl
# take list of fasta headers and remove them and their sequence
# Alan E. Yocca
# 03-19-18

use strict;
use Getopt::Std;

my $usage = "\n$0 -i <fasta> -r <list of headers to be removed> -o <output> \n\n";

our ($opt_i, $opt_r, $opt_o);
getopts('i:r:o:') or die "$usage";

if ( (!(defined $opt_i)) || (!(defined $opt_r)) || (!(defined $opt_o)) ) {
  print "$usage";
  exit;
}

#if (-e $opt_o ) {
#print "File: $opt_o or exist, is it okay to overwrite it?\n"; 
#my $answer = <STDIN>;
#	if ($answer =~ /^y(?:es)?$/i) {
#		print "Excellent!\n";
#	}
#	else {
#		die "fine, I will not overwrite your files\n";
#	}
#}

open (my $fasta_fh, '<', $opt_i) || die "Cannot open the file $opt_i\n\n";
open (my $remove_fh, '<', $opt_r) || die "Cannot open the to be removed list $opt_r\n\n";
open (my $out_fh, '>', $opt_o) || die "Cannot open output: $opt_o\n\n";

my %remove;

while (my $line = <$remove_fh>){
	chomp $line;
	$remove{$line} = 1;
}

my $loop_logic = 0;

FASTA: while (my $line = <$fasta_fh>) {
	chomp $line;
	if ($loop_logic) {
		$loop_logic = 0;
		next;
	}
	else {
		my @header = split(">", $line);
		if ($header[1] && $remove{$header[1]}) {
			$loop_logic = 1;
			next FASTA;
		}
		else {
			print $out_fh "$line\n";
		}
	}
}

close $fasta_fh;
close $remove_fh;
close $out_fh;

exit;