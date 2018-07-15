#! /usr/bin/perl -w

# fasta_rev_comp.pl
# reverse complement your fasta sequences with some added functionality
# Alan E. Yocca
# 07-11-18
# should work fine on wrapped fastas, but will return unwrapped

use strict;
use warnings;
use Getopt::Long;

#my $usage = "\n$0 -f <input fasta> -o <output> --list <list of fasta entries to rev complement. one per line please> --entries <comma-separated list of fasta entries, handy if only have like two so don't want to make list file> --all <reverse complement everything> --rev_only <reverse only, do not complement, handy for dot plot stuffs> --drop <keep only sequences to be reverse complemented\n\n";

my $opt_f;
my $opt_o;
my $list;
my $entries;
my $all = 0;
my $rev_only = 0;
my $drop = 0;

GetOptions (	'list=s' => \$list, 
		'entries=s' => \$entries,		
		'all=i' => \$all,
		'rev_comp=i' => \$rev_only,
		'drop=i' => \$drop,
		'f=s' => \$opt_f,
		'o=s' => \$opt_o
) or die usage();

if (!defined $opt_o || !defined $opt_f) {
	usage();
	die;
}

if (!defined $list && !defined $entries && !defined $all) {
	print "Please specify either a file with a list of fasta headers to reverse complement, or specify manually on command line, or explicitly tell me you want all of them switched. See usage:\n";
	usage();
	die;
}

open (my $fasta_fh, '<', $opt_f) || die "Cannot open the fasta file: $opt_f\n\n";
open (my $out_fh, '>', $opt_o) || die "Cannot open the output file: $opt_o\n\n";

my $header;
my $test = 0;
my $fasta_sequence;
my %list_hash;
my %entry_hash;
my $skip_not_spec = 0;
my $dropped = 0;

if ($list) {
	print "reading in $list\n";
	open (my $list_fh, '<', $list) || die "Cannot open the list file: $list\n\n";
	while (my $list_line = <$list_fh>) {
		chomp $list_line;
		if ($list_line =~ /^>/) {
			$list_hash{$list_line} = 1;
		}
		else {
			my $fasta_line = ">" . $list_line;
			$list_hash{$fasta_line} = 1;
		}
	}
	close $list_fh;
}
else {
	#nothing!
}

if (defined $entries) {
	my @entries = split(/,/,$entries);
	for (my $i = 0; $i < @entries; $i++) {
		if ($i =~ /^>/) {
			$entry_hash{$entries[$i]} = 1;
		}
		else {
			my $fasta_i = ">" . $entries[$i];
			$entry_hash{$fasta_i} = 1;
		}
	}
}
else {
	#nothing!
}

while (my $line = <$fasta_fh>) {
	chomp $line;
	if ($line =~ m/^>/) {
#		print "line -> $line\n";
#		print "header -> $header\n";
		if (not defined $fasta_sequence) { # first fasta sequence or empty one somewhere
#			print "should print round 1\n";
			$header = $line;
			next;
		}
		else { # loaded in full sequence
			if ($list_hash{$header} || $entry_hash{$header} || $all) {
				my $rev = reverse($fasta_sequence);
				if (not defined $rev_only) {
					$rev =~ tr/ATCGatcgNn/TAGCtagcnN/;
				} else {
					#do not complement
				}		
				print $out_fh "$header\n$rev\n";
			}
			else {
				$skip_not_spec = $skip_not_spec + 1;
				if ($drop) {
					$dropped = $dropped + 1;
				}
				else {
					print $out_fh "$header\n$fasta_sequence\n";
				}
			}
		}
		$header = $line;
#		print "defining header as $line\n";
		$fasta_sequence = "";
		next;
	}
	else {
		$fasta_sequence = $fasta_sequence . $line;		
	}
}

# get last sequence
if ($list_hash{$header} || $entry_hash{$header} || $all) {
	my $rev = reverse($fasta_sequence);
	if (not defined $rev_only) {
		$rev =~ tr/ATCGatcgNn/TAGCtagcnN/;
	} else {
		#do not complement
	}		
	print $out_fh "$header\n$rev\n";
}
else {
	$skip_not_spec = $skip_not_spec + 1;
	if ($drop) {
		$dropped = $dropped + 1;
	}
	else {
		print $out_fh "$header\n$fasta_sequence\n";
	}
}

print "skipped because not specified: $skip_not_spec\n";
print "dropped from $opt_o: $dropped\n";

close $fasta_fh;
close $out_fh;

sub usage {
	print "\n$0 -f <input fasta>\n";
	print "-o <output>\n";
	print "--list <list of fasta entries to rev complement. one per line please>\n";
	print "--entries <comma-separated list of fasta entries, handy if only have like two so don't want to make list file>\n";
	print "--all <reverse complement everything>\n";
	print "--rev_only <reverse only, do not complement, handy for dot plot stuffs>\n";
	print "--drop <keep only sequences to be reverse complemented>\n";
	print "\n";
}

exit;
