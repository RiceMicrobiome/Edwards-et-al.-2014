#!/usr/bin/perl
#fq_group_extract.pl
##################################
# Written by Joe Edwards
# Sundaresan Lab
# Department of Plant Biology
# University of California - Davis
# 2013
##################################

# This script is intended to take the group file from the demultiplex.pl script output
# and use it to pull out the specific lines of another FQ sequencing file.  This script
# can only be run on one FQ file at a time.  For our purposes we have been using MiSeq
# runs with paired end reads and have to run this script over each file.

use strict; use warnings;
use Getopt::Std;

#==============================================================================
## Get the options from command line
my %opts;
getopts('g:q:', \%opts);


#==============================================================================
## Usage statement
my $usage = "
fq_group_extract.pl
Last Updated: 

Written by Joe Edwards
PhD Student
Sundaresan Lab
Department of Plant Biology
University of California - Davis

Usage
fq_group_extract.pl -g groups.group -q fastq.fq

Options
";

#==============================================================================
open GROUP, "$opts{g}" or die "Cannot open group file $opts{g}\n";
open FQ, "$opts{q}" or die "Cannot open fastq file $opts{q}\n";

#==============================================================================
my $found = 0;																			# Start a counter to monitor progress
while (my $line = <GROUP>) {														# Open the group file
	chomp $line;																		# Remove the line break from the incoming line
	my ($id, $group, $bc) = split(/\t/, $line);								# Split the line into header, sample id, and barcode sequence and save into scalars 

	# Move through fastq file until correct header is found
	while (my $header = <FQ>) {													# Start the FQ reading loop				
		chomp $header;																	# Remove LB from the incoming line
		my $seq = <FQ>; chomp $seq;												# Save the next line as the sequence and remove LB
		my $skip = <FQ>;																# Skip the next line
		my $qual = <FQ>; chomp $qual;												# Save the next line as the quality and remove LB
		if ($header =~ /$id/) {														# If the header matches the header on the group file proceed
			print "$header\n$seq\n+\n$qual\n";									# And print this crap
			$found++;																	# Add to the amount of sequences foundfound
			last;																			# Exit the loop so it goes back to the correct place when resuming going through the FQ file
		} else {
			next;																			# If it's not the right header, proceed to the next line in the FQ file
		}
	}
warn "$found sequences found so far\n" if ($found % 10000 == 0);		# Update the user in real time of how many are found so they know the progress
}


