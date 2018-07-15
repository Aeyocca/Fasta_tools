#! /usr/bin/perl -w

# filter_fasta_length.pl
# filter fasta length either set max length or min length or both
#***************currently only supports one line fasta entries (eg, not wrapped sequences), think about changing later***************
# Alan E. Yocca
# 05-12-18

use strict;
use Getopt::Long;

my $usage = "\n$0 -f <input fasta> -o <output> --min <set minimum length (inclusive)> --max <set max length (inclusive)> \n\n";

my ($opt_f, $opt_o, $min, $max);

GetOptions ("f=s" => \$opt_f,
			"o=s" => \$opt_o,
			"min=i" => \$min,
			"max=i" => \$max
);

if ( (!(defined $opt_f)) || (!(defined $opt_o)) ) {
  print "$usage";
  exit;
}

if ( (!(defined $max)) && (!(defined $min)) ) {
  print "$usage";
  exit;
}

open (my $fasta_fh, '<', $opt_f) || die "Cannot open the fasta file: $opt_f\n\n";
open (my $out_fh, '>', $opt_o) || die "Cannot open the output file: $opt_o\n\n";

my $written = 0;
my $max_filtered_out = 0;
my $min_filtered_out = 0;
my $header;

FASTA_ENTRY: while (my $line = <$fasta_fh>) {
	chomp $line;
	if ($line =~ m/^>/) {
		$header = $line;
		next;
	}
	else {
		my $fasta_length = length($line);
		if ($max && $fasta_length > $max) {
			$max_filtered_out = $max_filtered_out + 1;
			next FASTA_ENTRY;
		}
		if ($min && $fasta_length < $min) {
			$min_filtered_out = $min_filtered_out + 1;
			next FASTA_ENTRY;
		}
		else {
			print $out_fh "$header\n$line\n";
			$written = $written + 1;
		}
	}
}

if ($min) {
	print "Sequences filtered out < $min bases:\t$min_filtered_out\n";
}

if ($max) {
	print "Sequences filtered out > $max bases:\t$max_filtered_out\n";
}

print "Sequences written to $opt_o:\t$written\n";

close $fasta_fh;
close $out_fh;

exit;