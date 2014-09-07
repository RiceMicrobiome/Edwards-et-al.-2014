#!/usr/bin/perl
#qiime_taxer.pl
##################################
# Written by Joe Edwards
# Sundaresan Lab
# Department of Plant Biology
# University of California - Davis
# 2012
##################################

use strict;
use warnings;

my $file = $ARGV[0];

open TAX, "$file" or die "Cannot open the file\n";

print "OTU\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
while (my $line = <TAX>) {
        chomp $line;
        my($otu, $taxonomy, $number) = (split(/\t/, $line));
        $taxonomy =~ s/ /_/g;
        $taxonomy =~ s/.__//g;
        print "$otu\t";
        my @tax = split(/;/, $taxonomy);
        for(0..6) {
                if(exists $tax[$_]) {
                        if ($tax[$_] =~ /\w/) {
                                print "$tax[$_]\t";
                        } else {
                                print "unclassified\t";
                        }
                } else {
                        print "unclassified\t";
                }
        }
        print "\n";
}
