#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

my %opts;
getopts('f:m:', \%opts);

open FA, "$opts{f}" or die;

my $processed = 0;
my $removed_for_length = 0;
my $ambigs = 0;
my $kept = 0;
while (my $header = <FA>) {
	chomp $header;
	my $seq = <FA>; chomp $seq;
	warn "$processed total\t$removed_for_length removed due to length\t$ambigs removed due to ambiguities\n" if $processed % 10000 == 0;
	if (length $seq > $opts{m}) {
		$processed++;
		$removed_for_length++;
		next;
	}
	$seq =~ s/N$//;
	if ($seq =~ /N/) {
		$processed++;
		$ambigs++;
		next;
	}
	$processed++;
	print "$header\n$seq\n";
}
