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
	# Bed file containing features to measure distance from.
	my $infile = shift;

	# Deconstruct location
	$location =~ /^(\S+):(\d+)-(\d+)$/ or die "ERROR: regex failed with '$location'.\n";
        my $chr         = $1;
        my $start       = $2;
        my $end         = $3;
	
	# Stores location of closest feature 
	my %shortest;

	# Work through each feature in turn.
	open (FEATURE, $infile) or die "ERROR: could not open $infile.\n";
	
	# For current feature need to determine its proximity to hotspot
	# and whether that location is the closest found so far.
	while (<FEATURE>) {

		chomp;
		# Get feature position.
                my ($featureChr, $featureStart, $featureEnd, $label) = split(/\t/, $_);

		# Convert from 0 (BED format) to 1 indexing.
                $featureStart++;
               
		 # Ignore if feature falls on a different chromosome.
                next if $featureChr ne $chr;

		# Determine whether feature is before, after or overlapping
     		# hotspot and process accordingly.
		if ($featureEnd < $start) {
			
			my $score = ($start - $featureEnd) - 1;

			processFeature(
				\%shortest,
				$score,
				$label,
				"downstream"
				);

			}
		elsif ($featureStart > $end) {
				
			my $score = ($featureStart  - $end) - 1;
				
			processFeature(
				\%shortest, 
				$score, 
				$label,
				"upstream"
				);

			}
		else {
			# Assume overlapping - shouldn't happen in this context 
			# so throw error exception

			# Might overlap centromere - this is ok
			if ($label eq "CENTROMERE") {
				warn "WARNING: '$location' is overlapping centromere '$featureChr:$featureStart-$featureEnd, $label'.\n";
				$shortest{downstream}->{distance} = -1;
				$shortest{downstream}->{label} = $label;
				$shortest{upstream}->{distance} = -1;
                                $shortest{upstream}->{label} = $label;	
				}
			# Shouldn't overlap telomere?
			else {
				die "ERROR: unexpected overlap with $label. '$location'\n";
				}

			}	

		}
	close (FEATURE) or die "ERROR: could not close $infile.\n";

	if (!exists $shortest{downstream}) {
		return;
		}

	# If downstream is centromere, return upstream (obviously telomere) distance.
	elsif (exists $shortest{downstream}->{label} && $shortest{downstream}->{label} eq "CENTROMERE") {
		return $shortest{upstream}->{distance};
		}
	# if upstream is centromere, return downstream (obviously telomere) distance.
	elsif (exists $shortest{upstream}->{label} && $shortest{upstream}->{label} eq "CENTROMERE") {
		return $shortest{downstream}->{distance};
		}
	elsif ($shortest{downstream}->{distance} == -1) {
		return -1;	
		}
	else {
		warn "WARNING: neither upstream nor downstream is centromere. '$location'.\n";
		
		if ($shortest{downstream}->{distance} < $shortest{upstream}->{distance}) {
			return $shortest{downstream}->{distance};
			} 
		else {
                        return $shortest{upstream}->{distance};
                	}
		}

	} # End of method.

###############################################################################

sub processFeature {
	
	my $shortest	 	= shift;
	my $score		= shift;	
	my $label		= shift;
	my $key			= shift;

	if (!defined $shortest->{$key}->{distance}) {
		$shortest->{$key}->{distance} = $score;
		$shortest->{$key}->{label} = $label;
		}

	elsif ($score < $shortest->{$key}->{distance}) {
		$shortest->{$key}->{distance} = $score;
		$shortest->{$key}->{label} = $label;
		}
	else {
		# do nothing.
		}

	} # End of method.

###############################################################################
1;
