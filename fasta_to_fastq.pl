#! /usr/bin/perl -w
#Copyright (c) 2010 LUQMAN HAKIM BIN ABDUL HADI (csilhah@nus.edu.sg)
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files 
#(the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, 
#merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
#OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
#LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
#IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#06-01-18
#Alan E. Yocca edit to make base quality 30 instead of 40, not sure how that will affect my MaSuRCA assembly by o whale
#added rant:
#Also added getopt so it takes -f and -o flag, works better for my pipe because I run through it many times and if it fails it will still create the output file if i > pipe into it then when I fix the pipe and run again it it will fail because empty file still exists and cannot overwrite,,, I guess I could just specify the overwrite output.... no thats too easy


use strict;
use Getopt::Std;

my $usage = "\n$0 -f <input fasta> -o <output fastq> \nASSIGNS BASE QUALITY OF 30 IN PHRED33!!!\n\n";

our ($opt_f, $opt_o);
getopts('f:o:') or die "$usage";

if ( (!(defined $opt_f)) || (!(defined $opt_o)) ) {
  print "$usage";
  exit;
}

open (my $fasta_fh, '<', $opt_f) || die "Cannot open the fasta file: $opt_f\n\n";
open (my $out_fh, '>', $opt_o) || die "Cannot open the output file: $opt_o\n\n";

#my $file = $ARGV[0];
#open FILE, $file;

my ($header, $sequence, $sequence_length, $sequence_quality);
while(my $line = <$fasta_fh>) {
        chomp $line;
        if ($line =~ /^>(.+)/) {
                if(length $header) {
                        print $out_fh "\@".$header."\n";
                        print $out_fh "$sequence\n";
                        print $out_fh "+\n";
                        print $out_fh "$sequence_quality\n";
                }
                $header = $1;
		$sequence = "";
		$sequence_length = "";
		$sequence_quality = "";
        }
	else { 
		$sequence .= $line;
		$sequence_length = length($line); 
		for(my $i=0; $i<$sequence_length; $i++) {$sequence_quality .= "?"} 
	}
}
close $fasta_fh;
print $out_fh "\@".$header."\n";
print $out_fh $sequence."\n";
print $out_fh "+"."\n";
print $out_fh $sequence_quality."\n";

close $out_fh;
exit;
