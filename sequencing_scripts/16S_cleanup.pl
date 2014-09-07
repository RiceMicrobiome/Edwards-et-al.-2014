#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

## Get the options from the command line
my %opts;
getopts('f:m:vu', \%opts);

## Set up usage statement
my $usage = "Usage:\n16S_cleanup.pl -f input.fa -m 275 > cleaned.fasta
	-f\tinput fasta file
	-m\tmaximum sequence length
	-v\tverbose
	-u\tprint usage statment and die\n";

## Print the usage statement if the user wants it
die "$usage" if $opts{u};

## Open the fasta file
open FA, "$opts{f}" or die "[ERROR] Could not open fast file\n$usage";


my $processed = 0;
my $removed_for_length = 0;
my $ambigs = 0;
my $kept = 0;
my $max = 275;
if ($opts{m}) {
	$max = $opts{m};
	print "Using $max bp as the maximum length for sequences.";
} else {
	print "Defaulting to 275 bp as the maximum length for sequences.";
}

## Loop through fasta file
while (my $header = <FA>) {
	chomp $header;
	my $seq = <FA>; chomp $seq;
	if ($processed % 10000 == 0) {
		warn "$processed total\t$removed_for_length removed due to length\t$ambigs removed due to ambiguities\n" if $opts{v};
	}
	if (length $seq > $opts{m}) {
		$processed++;
		$removed_for_length++;
		next;
	}
	## Look for ambiguous base calls in reads
	$seq =~ s/N$//;
	if ($seq =~ /N/) {
		$processed++;
		$ambigs++;
		next;
	}

	## Print it if it's all okay
	$processed++;
	print "$header\n$seq\n";
}

warn "Processed $processed total sequences.
Removed $ambigs sequences due to ambiguities.
Removed $removed_for_length due to being too long.";
