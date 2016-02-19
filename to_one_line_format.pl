#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Basename;

my %options;
getopts("i:", \%options);
my $infile = $options{i} or &usage;
sub usage {die "USAGE: " . basename($0) . " [-i two-line circos linkage infile to convert].\n";}

my @line1;
my @line2;
open (IN, $infile) or die "ERROR: could not open $infile.\n";
while (@line1 = split(/\s/, <IN>)) {
        @line2 = split(/\s/, <IN>);

	# sanity check
	die if $line1[0] ne $line2[0];

	printf(
		"%s %d %d %s %d %d %s\n",
		$line1[1],
		$line1[2],
		$line1[3],
		$line2[1],
                $line2[2],
                $line2[3],
		$line1[4]
		);

        }
close (IN) or die "ERROR: could not close $infile.\n";

