#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Basename;

my %options;
getopts("i:m:g:r:", \%options);
my $infile 	= $options{i} or &usage;
my $min		= $options{m} or &usage;
my $gap		= $options{g} or &usage;
my $roundup	= $options{r} || "false";
sub usage {die "USAGE: " . basename($0) . " [-i mpileup infile] [-m min mean coverage of region] [-g max gap between blocks] [-r if true, roundup (default, false)]\n";}

my $start = undef;
my $end = undef;
my $currentChr = undef;
my $covSum = 0;
open (IN, $infile) or die "ERROR: Could not open $infile.\n";

while (<IN>) {
	chomp;
	/^(\S+)\t(\d+)\t(\w)\t(\d+)/ or die "ERROR: regex failed with line '$_'.\n";
	my $chr = $1;
	my $pos = $2;
	my $cov = $4;	

	if (!defined $start) {
		$start = $pos-1;
		$end = $pos;
		$currentChr = $chr;
		}

	elsif ($currentChr ne $chr) {

		process($covSum, $roundup, $currentChr, $start, $end);
		
		$start = $pos-1;
		$currentChr = $chr;
		$covSum = 0;
		}

	elsif ($end+$gap < $pos) {

        	process($covSum, $roundup, $currentChr, $start, $end);        

		$start = $pos-1;		
		$covSum = 0;
		}

	

	$end = $pos;
	$covSum += $cov;


	}
 
process($covSum, $roundup, $currentChr, $start, $end) unless !defined $start;

close (IN) or die "ERROR: could not close $infile.\n";

###############################################################################

sub process {

	my $covSum 	= shift;
	my $roundup 	= shift;
	my $currentChr	= shift;	
	my $start 	= shift;
	my $end		= shift;

	my $rawMeanCov = $covSum / ($end-$start);

	my $meanCov = $rawMeanCov;
	if ($roundup =~ /^t(rue)$/i) {
        	$meanCov = sprintf("%.0f", $rawMeanCov);
        	}
	
	printf("%s\t%d\t%d\t%.3f\n", $currentChr, $start, $end, $rawMeanCov) if $meanCov >= $min;

	} # End of method.

###############################################################################
