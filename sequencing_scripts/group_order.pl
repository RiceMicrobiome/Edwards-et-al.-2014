#!/usr/bin/perl
#group_order.pl
##################################
# Written by Joe Edwards
# Sundaresan Lab
# Department of Plant Biology
# University of California - Davis
# 2012
##################################

# This script reorders a group file based on the order of the headers in a fastq file

use strict; use warnings;
use Getopt::Std;

#==============================================================================
## Get the options from command line
my %opts;
getopts('g:q:', \%opts);


#==============================================================================
## Usage statement
my $usage = "
group_order.pl
Last Updated: 

Written by Joe Edwards
PhD Student
Sundaresan Lab
Department of Plant Biology
University of California - Davis

sage
group_order.pl -g group.file -q fastq.file

Options
";

#==============================================================================
# Open the two files to work on
open FQ, "$opts{q}" or die "Cannot open fastq file $opts{q}\n";
open GROUP, "$opts{g}" or die "Cannot open group file $opts{g}\n";

#==============================================================================
my %group_hash = ();													# Initialize hash to save the contents of the group file into
while (my $line = <GROUP>) {										# Start to loop through the group file
	chomp $line;														# Take the line break off the incoming line
	my ($id, $group, $barcode) = split(/\t/, $line);		# Split the line in the group file into id, group, and barcode
	$group_hash{$id} = "$group\t$barcode";						# ...and load it into the hash
}
warn "Finished constructing hash\n";
#==============================================================================
# Order the group file based on the order of the fastq file
while (my $header = <FQ>) {										# Start looping through the FASTQ file
	chomp $header;														# Remove the LB
	my $skip = <FQ>; $skip = <FQ>; $skip = <FQ>;				# Skip everything but the header
	print "$header\t$group_hash{$header}\n";					# Print to a new file
}



