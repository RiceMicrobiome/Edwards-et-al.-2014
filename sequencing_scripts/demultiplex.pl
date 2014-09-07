#!/usr/bin/perl 
# demultiplex.pl
  
#################################
# Written by Joe Edwards
# Sundaresan Lab
# Department of Plant Biology
# University of California - Davis
# 2013
##################################

# This script takes the fused barcode reads from a miseq run and
# uses a qiime formatted mapping file to attribute the barcodes
# to specific libraries. The results of this script will be fed into
# the script fq_group_extract.pl
# WARNING: This script does EXACT matching of the barcodes.

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('b:m:', \%opts); # b for barcode fastq and m for mapping file


open MAP, "$opts{m}" or die "Cannot open the mapping file\n"; 	# Shut down and warn if no mapping file
open BAR, "$opts{b}" or die "Cannot open the barcode file\n"; 	# Shut down and warn if no barcode file

# Open the mapping file and create a hash that hash with barcode as key and the sample name as the value
my %bcs = ();																	
while (my $line = <MAP>) {													# Start while loop
	chomp $line;																# Get rid of line break at the end of the line
	next if $line =~ /^#/;													# Skip the line if it's commented out
	my @fields = split(/\t/, $line);										# Split the line by tabs and save into array
	my $bc = $fields[1];														# Second member of the array is the barcode
	my $samp = $fields[0];													# First member of the array is the sample name
	$bcs{$bc} = $samp;														# Load everything into the hash
}

my $found = 0;																	# Create a counter for the amount of barcodes found

while (my $header = <BAR>) {     										# Start looping through the barcode FQ file
	chomp $header;																# Get rid of the line breaks
	my $seq = <BAR>; chomp $seq;											# Save the sequence and remove LB
	my $skip = <BAR>;															# Skip the next line
	my $qual = <BAR>; chomp $qual;										# Save the quality into a scalar and remove LB
  
	if (exists $bcs{$seq}) {												# See if the barcode exists in the previous hash
		print "$header\t$bcs{$seq}\t$seq\n";							# If it does, print out the header and the sample id and the barcode sequence
		$found++;																# Add to the amount of barcodes found
	} else {
		next;																		# If not found, go to the next line
	}
	warn "$found bcs found\n" if ($found % 10000 == 0);			# Update in real time how many barcodes are found
}
