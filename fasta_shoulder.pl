#! /usr/bin/perl -w
# fasta_shoulder.pl 
# collect shoulder fasta sequences by length or percent, specify minimum cutoffs and things
# Alan E. Yocca
# 05-29-18

use strict;
use warnings;
use Getopt::Long;

my $usage = "./fasta_shoulder.pl\n-f <input fasta>\n-o <output fasta>\n--min <minimum length of sequences to output, default: 30>\n--length <if specified, this length of each end of the input fasta sequences will be kept.>\n--percent <percent of each end of the sequence to output.>\n\n";

my $min_length = 30; #don't think I want sequences shorter than 30 bp because assembly will be bad, right?
my $length;
my $percent;
my $opt_f;
my $opt_o;

GetOptions (	'min=i' => \$min_length, 		
		'length=i' => \$length,
		'percent=i' => \$percent,
		'f=s' => \$opt_f,
		'o=s' => \$opt_o
) or die "$usage\n";

if ($length && $percent) {
	print "only specify length OR percent:\n$usage";
	exit;
}

if (!defined $length && !defined $percent) {
	print "$usage";
	exit;
}

#$length = 100; #set default length to 100 bp on either side
#$percent = 10; #set default percent to 10% on either side

open (my $fasta_fh, '<', $opt_f) || die "Cannot open the fasta file: $opt_f\n\n";
open (my $out_fh, '>', $opt_o) || die "Cannot open the output file: $opt_o\n\n";

my $header;
#my @ambiguous_sequence_distribution;
my $too_short = 0;
my $written = 0;

while (my $line = <$fasta_fh>) {
	chomp $line;
	if ($line =~ m/^>/) {
		$header = $line;
		next;
	}
	else {
		#could if else length percent here
		if ($length) {
			my $left_shoulder = substr($line, 0, $length);
			my $right_shoulder = substr($line, -$length);
			my $ls_length = length $left_shoulder;
			my $rs_length = length $right_shoulder;
			if ($ls_length >= $min_length && $rs_length >= $min_length) {
				my $header_left = $header . "_1"; #I think 1 and 2 will be easier to parse later than left and right
				my $header_right = $header . "_2";
				print $out_fh "$header_left\n$left_shoulder\n$header_right\n$right_shoulder\n";
				$written = $written + 1;
			}
			else {
				$too_short = $too_short + 1;
			}
		} elsif ($percent) {
			my $sequence_length = length $line;
			my $percent_decimal = $percent / 100;
			my $shoulder_length = $sequence_length*$percent_decimal;
			my $shoulder_length_rounded = sprintf "%.0f",$shoulder_length;
			my $left_shoulder = substr($line, 0, $shoulder_length_rounded);
			my $right_shoulder = substr($line, -$shoulder_length_rounded);
			my $ls_length = length $left_shoulder;
			my $rs_length = length $right_shoulder;
			if ($ls_length >= $min_length && $rs_length >= $min_length) {
				my $header_left = $header . "_1"; #I think 1 and 2 will be easier to parse later than left and right
				my $header_right = $header . "_2";
				print $out_fh "$header_left\n$left_shoulder\n$header_right\n$right_shoulder\n";
				$written = $written + 1;
			}
			else {
				$too_short = $too_short + 1;
			}
		}
	}
}

#should I do some summary stats? This should do it for me, need to write usage statement and stuff
print "filtered out for being too short:\t$too_short\n";
print "pairs written out to $opt_o:\t$written\n";

close $fasta_fh;
close $out_fh;




