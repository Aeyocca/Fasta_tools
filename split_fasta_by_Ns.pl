#! /usr/bin/perl -w

# split_fasta_by_Ns.pl
# getting really lazy with my names.. take in fasta split everything by ambiguous characters (either N or n) and output fasta of contiguous unambiguous sequences
# Alan E. Yocca
# 05-09-18
# 07-03-18
# added option for minimum length to output


use strict;
use warnings;
use Getopt::Long;

my $usage = "\n$0 -f <input fasta> -o <output> --min <sequences shorter than this length will be dropped [30]>\n\n";

my $min_length = 30;

GetOptions (	'min=i' => \$min_length, 		
		'f=s' => \$opt_f,
		'o=s' => \$opt_o
) or die "$usage\n";

open (my $fasta_fh, '<', $opt_f) || die "Cannot open the fasta file: $opt_f\n\n";
open (my $out_fh, '>', $opt_o) || die "Cannot open the output file: $opt_o\n\n";

my $header;
my $test = 0;
my @ambiguous_sequence_distribution;
my $dropped = 0;
my $length_dropped = 0;

while (my $line = <$fasta_fh>) {
	chomp $line;
	if ($line =~ m/^>/) {
		$header = $line;
		next;
	}
	else {
		my @sequence = split(/[Nn]/,$line);
		#then loop through array to produce statistics, thinking you can load all the lengths into an array to summarize
		my $ambiguous_sequence_length = 1;
		my $sub_header_count = 0;
		my $sequence_length = @sequence;
		my $sequence_max_index = $sequence_length - 1;
		for (my $i=0; $i < @sequence; $i++) {
			if ($sequence[$i]) {
				#found some unambiguous sequence!
				if ($i != 0) {
					if ($i == 1 && $line =~ m/^[Nn]/) {
						#this means there was just a single N at the beginning followed by actual sequence
						#reset ambiguous_sequence_length to 1
						$ambiguous_sequence_length = 1;
						push(@ambiguous_sequence_distribution, $ambiguous_sequence_length);
					}
					else {
						push(@ambiguous_sequence_distribution, $ambiguous_sequence_length);
						$ambiguous_sequence_length = 1;
					}
				}
				if (length $sequence[$i] >= $min_length) {
					$sub_header_count = $sub_header_count + 1;
					my $sub_header = $header . "_" . $sub_header_count;
					print $out_fh "$sub_header\n$sequence[$i]\n";
				}
				else {
					$dropped = $dropped + 1:
					$length_dropped = $length_dropped + length($sequence[$i]);
				}
				if ($i == $sequence_max_index && $line =~ m/[nN]$/) {
					#this means we only have a single trailing N or n at the end of the sequence so we need to count it
					push(@ambiguous_sequence_distribution, $ambiguous_sequence_length);
				}
			}
			else {
				if ($i == 0){
					#don't add another one
				}
				else {
					$ambiguous_sequence_length = $ambiguous_sequence_length + 1;
				}
			}
		}
	}
}

my $number_of_splits = @ambiguous_sequence_distribution;
my $asd_max_index = $number_of_splits - 1;

my $sum = 0;
for (my $i = 0; $i < @ambiguous_sequence_distribution; $i++) {
	$sum = $sum + $ambiguous_sequence_distribution[$i];
}

my $avg = $sum / $number_of_splits;

my @sorted_asd = sort {$a <=> $b } @ambiguous_sequence_distribution;
my $largest = $sorted_asd[$asd_max_index];
my $median;
if (0 == $number_of_splits % 2) {
	my $med_index = $number_of_splits / 2;
	my $med_index_minus_one = $med_index - 1;
	$median = ($sorted_asd[$med_index_minus_one] + $sorted_asd[$med_index]) / 2; #take average of the two numbers around the median	
}
else {
	my $middle = $number_of_splits / 2;
	my $med_index = $middle - 0.5;
	$median = $sorted_asd[$med_index];
}

my @asd_not_one;
my $singletons = 0;
my $sum_not_one = 0;
my $splits_not_one = 0;

#for this next section tagged a bunch of variables n1 for not one so the length of stretches of N that aren't singletons. probably could have just overwritten the variables used earlier but memory not a concern here and not sure if something would go wrong if I tried to overwrite, this one a little less ambiguous also... man I seem to love that word today
for (@ambiguous_sequence_distribution) {
	if (/^1$/) {
		$singletons = $singletons + 1;
		next;
	}
	else {
		$sum_not_one = $sum_not_one + $_;
		$splits_not_one = $splits_not_one + 1;
		push(@asd_not_one, $_);
	}
}

my $avg_not_1 = $sum_not_one / $splits_not_one;

my @sorted_asd_n1 = sort {$a <=> $b } @asd_not_one;

my $median_n1;
if (0 == $splits_not_one % 2) {
	my $med_index_n1 = $splits_not_one / 2;
	my $med_index_minus_one_n1 = $med_index_n1 - 1;
	$median_n1 = ($sorted_asd_n1[$med_index_minus_one_n1] + $sorted_asd_n1[$med_index_n1]) / 2; #take average of the two numbers around the median	
	
}
else {
	my $middle_n1 = $splits_not_one / 2;
	my $med_index_n1 = $middle_n1 - 0.5;
	$median_n1 = $sorted_asd_n1[$med_index_n1];
}

print "\n";
print "Number of gaps split: $number_of_splits\n";
print "Average size of gaps: $avg\n";
print "Average size of gaps >1: $avg_not_1\n";
print "Singletones = $singletons\n";
print "median size of gaps: $median\n";
print "median size of gaps >1: $median_n1\n";
print "Largest gap: $largest\n";
print "Segments dropped for falling below --min ($min_length): $dropped\n";
print "Total length of segments dropped: $length_dropped\n";

close $fasta_fh;

exit;
