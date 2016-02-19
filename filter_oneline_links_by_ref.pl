#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Basename;

my %options;
getopts("i:r:", \%options);
my $infile = $options{i} or &usage;
my $refId = $options{r} or &usage;
sub usage {die "USAGE: " . basename($0) . " [-i Circos one-line links infile to process] [-r ref id to filter by].\n";}

my @line;
open (IN, $infile) or die "ERROR: could not open $infile.\n";
while (@line = split(/\s/, <IN>)) {

	if ($line[0] =~ /^$refId$/i || $line[3] =~ /^$refId$/i) {
		print join(" ", @line), "\n";		
		}	

	}	
close (IN) or die "ERROR: could not close $infile.\n";

