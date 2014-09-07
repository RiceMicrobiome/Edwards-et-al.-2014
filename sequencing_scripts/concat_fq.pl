#!/usr/bin/perl
#concat_fq.pl
##################################
# Written by Joe Edwards
# Sundaresan Lab
# Department of Plant Biology
# University of California - Davis
# 2012
##################################

# This script concatenates the reads of two fastq files together.  I use this
# to concatenate the barcodes together from miseq runs.  The results of this 
# can be fed into the demultiplex.pl script.

use strict; use warnings;
use Getopt::Std;

#==============================================================================
## Get the options from command line
my %opts;
getopts('f:r:', \%opts);


#==============================================================================
## Usage statement
my $usage = "
concat_fq.pl
Last Updated: 

Written by Joe Edwards
PhD Student
Sundaresan Lab
Department of Plant Biology
University of California - Davis

Usage
concat_fq.pl -f fwd.fq -r rvs.fq

Options
	-f:	forward read
	-r:	reverse read
";

#==============================================================================
## Ensure kosherness on the commmand line
die "Please include forward read $usage\n" unless defined $opts{f};
die "Please include reverse read $usage\n" unless defined $opts{r};

#==============================================================================
## Open the fq files and get down to business
open F, "$opts{f}" or die "Cannot open file $opts{f}";
open R, "$opts{r}" or die "Cannot open file $opts{r}";

while (my $fhead = <F>) {								
	my $rhead = <R>;
	my $fseq = <F>;
	my $rseq = <R>;
	my $skip = <F>; $skip = <R>;
	my $fqual = <F>;
	my $rqual = <R>;
	
	chomp $fhead;
	chomp $rhead;
	chomp $fseq;
	chomp $rseq;
	chomp $fqual;
	chomp $rqual;

	my ($seqname) = $fhead =~ /@(.*)\s/;
	my $seq = $fseq . revcomp($rseq);
	my $qual = $fqual . reverse($rqual);

	print "\@$seqname\n$seq\n+\n$qual\n";
	
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Subroutines Follow
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub revcomp {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr/ATCGatcg/TAGCtagc/;
	return $seq;
#	return (reverse ($_[0] =~ tr/ATCGatcg/TAGCtagc/));
}



















