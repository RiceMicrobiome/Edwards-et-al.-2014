#!/usr/bin/perl
#fa_to_fq.pl
##################################
# Written by Joe Edwards
# Sundaresan Lab
# Department of Plant Biology
# University of California - Davis
# 2012
##################################
use lib "/home/ledwards/Perl_Scripts/";
use strict; use warnings;
use Getopt::Std;

#==============================================================================
## Get the options from command line
my %opts;
getopts('q:', \%opts);


#==============================================================================
## Usage statement
my $usage = "
fa_to_fq.pl
Last Updated: 

Written by Joe Edwards
PhD Student
Sundaresan Lab
Department of Plant Biology
University of California - Davis

Usage
fa_to_fq.pl -q [quality] file.fa

Options
	-q:	Specify the quality you want to use
";

#==============================================================================
die "Please specify a quality\n$usage" unless defined $opts{q};
die "$usage" unless @ARGV == 1;

#==============================================================================
open FA, "$ARGV[0]" or die "Cannot open file $ARGV[0]\n";
while (my $header = <FA>) {
	chomp $header;
	$header =~ tr/>/@/;
	my $seq = <FA>; chomp $seq;
	my $qual = "";
	$qual .= $opts{q} for (1..length $seq);
	print "$header\n$seq\n+\n$qual\n";
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Subroutines Follow
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++










