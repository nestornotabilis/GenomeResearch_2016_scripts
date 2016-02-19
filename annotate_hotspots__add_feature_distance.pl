#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Basename;

my %options;
getopts("i:f:l:", \%options);
my $infile 	= $options{i} or &usage;
my $bedfile 	= $options{f} or &usage;
my $label	= $options{l} or &usage;
sub usage {die "USAGE:" . basename($0) . " [-i hotspot table infile] [-f feature file bed format] [-l label].\n";}

#-----------------------------------------------------------------------------
# Process each hotspot in table file.
open (IN, $infile) or die "ERROR: could not open $infile.\n";
while (my $line = <IN>) {
	chomp($line);
	if ($line =~ /^Location/) {
		print "$line\t$label\n";
		next;
		} 

	my (
		$location, 
		$meanCov, 
		$range, 	
		$noReads, 
		$freqPlus, 
		$freqMinus, 
		$orientationFreq, 
		$meanMAPQ, 
		$meanSoftclipped, 
		$percentErrors, 
		$mateLocation, 
		$mateRange
		) = split(/\t/, $line);

	my $shortestDistance = process($location, $bedfile);

	if (!defined $shortestDistance) {
		$shortestDistance = "-";
		}
	elsif ($shortestDistance == -1) {
		$shortestDistance = "overlap";
		}


	print "$line\t$shortestDistance\n";

	}
close (IN) or die "ERROR: could not close $infile.\n";

###############################################################################

sub process {

	# Location of hotspot within chromosome.	
	my $location = shift;
	# Bed filee containing features to measure distance from.
	my $infile = shift;

	# Deconstruct location
	$location =~ /^(\S+):(\d+)-(\d+)$/ or die "ERROR: regex failed with '$location'.\n";
        my $chr         = $1;
        my $start       = $2;
        my $end         = $3;
	
	# Stores location of closest feature 
	my $shortestDistance;

	# Work through each feature in turn.
	open (FEATURE, $infile) or die "ERROR: could not open $infile.\n";
	
	# For current feature need to determine its proximity to hotspot
	# and whether that location is the closest found so far.
	while (<FEATURE>) {

		chomp;
		# Get feature position.
                my ($featureChr, $featureStart, $featureEnd) = split(/\t/, $_);

		# Convert from 0 (BED format) to 1 indexing.
                $featureStart++;
               
		 # Ignore if feature falls on a different chromosome.
                next if $featureChr ne $chr;

		# Determine whether feature is before, after or overlapping
     		# hotspot and process accordingly.
		if ($featureEnd < $start) {
			$shortestDistance = processDownstreamFeature($shortestDistance, $featureEnd, $start);
			}
		elsif ($featureStart > $end) {
			$shortestDistance = processUpstreamFeature($shortestDistance, $featureStart, $end);
			}
		else {
			# Assume overlapping.
			$shortestDistance = -1;
			}	

		}
	close (FEATURE) or die "ERROR: could not close $infile.\n";

	return ($shortestDistance);

	} # End of method.

###############################################################################

sub processDownstreamFeature {
	
	my $shortestDistance 	= shift;
	my $featureEnd	 	= shift;
	my $start 		= shift;	

	my $score = ($start - $featureEnd) - 1;

	if (!defined $shortestDistance) {
		return $score;
		}

	elsif ($score < $shortestDistance) {
		return $score;
		}
	else {
		return $shortestDistance;
		}

	} # End of method.

###############################################################################

sub processUpstreamFeature {
	
	my $shortestDistance 	= shift;
	my $featureStart	= shift;
	my $end			= shift;	

	my $score = ($featureStart  - $end) - 1;
	
	if (!defined $shortestDistance) {
                return $score;
                }

        elsif ($score < $shortestDistance) {
                return $score;
                }
        else {  
                return $shortestDistance;
                }

	} # End of method.

###############################################################################
1;
