# Fasta_tools
A bunch of perl scripts I wrote to do various things with fasta files

Heads up! just putting these up here after I got them to work for my purposes.
Mostly for documentation of things I have done.
Some scripts may depend on specific header formats.
Hell, some might not even do what they specify in cases that I did not encounter.
Also, will probably not handle wrapped fasta files unless specified (can use fa_to_one_line.pl to convert it, see below)
I do not claim responsibility for any bit of code.
I simply stand on the shoulders of giants before me, heavily reliant on others

But I do think several can be useful to others
For all scripts, executing them without any arguments will print the usage statement to STDOUT
It may also be handy to open the script and read some of the documentation at the top

To note a few:

# fasta_one_line.pl
convert wrapped fasta file to single line per entry
Since some of these scripts can't handle wrapped fastas (unless specified), this will be a good first step

# fasta_length_dist_vect.pl
take fasta file, output a comma separated vector (list??, its just one line) of the lengths of each fasta entry
think of adding an ignore N switch at some point
WILL NOT DO ANY CHECKS, simply get number of characters in fasta entry
should make it work on wrapped fasta
Useful for making a histogram of the lengths of your fasta entries

# fasta_to_fastq.pl
Convert fasta to fastq, assigning a single specified Phred33 score to every single base-pair

To run these scripts:
check you have Perl v5.10.1 (or anything compatible, it was written and tested in this version)

Either:
from the command line, in the directory where the perl script is located;
./fasta_to_fastq.pl

OR

perl ./fasta_to_fastq.pl

To get the first method to work, you may have to change the file permissions to make it executable

Any questions just email me:
aeyap42@gmail.com
