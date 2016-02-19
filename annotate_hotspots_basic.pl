#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Basename;

my %options;
getopts("f:b:r:", \%options);
my $bedfile = $options{f} or &usage;
my $bamfile = $options{b} or &usage;
my $reffile = $options{r} or &usage;
sub usage {die "USAGE:" . basename($0) . " [-f hotspot feature bedfile] [-b bamfile] [-r reference].\n";}

printf(
	"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
#                "chr",
#                "start",
#                "end",
		"Location",
                "mean cov",
                "range",
                "no reads",
                "freq plus",
                "freq minus",
		"orientation freq",
                "mean MAPQ",
		"mean softclipped",
		"% errors",
		"mate location",
		"mate range"
                );


#-----------------------------------------------------------------------------
# Process each hotspot in BED file.
open (BED, $bedfile) or die "ERROR: could not open $bedfile.\n";
while (my $bedLine = <BED>) {
	chomp($bedLine);
	if ($bedLine =~ /^ID/) {next;} 
	
	my ($chr, $start, $end, $meanCoverage) = split(/\t/, $bedLine);

	# Only consider standard reference sequences and stong sequences.
        next if $chr !~ /^chr/ && $chr !~ /_stong$/;
	
	my $size = $end - $start;
	$start = $start + 1;
	
	my $isPlus = 0;
	my $isMinus = 0;
	my @mapqs;
	my @softClipLengths;
	my $n = 0;
	my %mateMappings;
	# Analyse reads associated with current hotspot
	open (PIPE, "samtools view $bamfile $chr:$start-$end |") or die;
	while (my $bamLine = <PIPE>) {
		my ($qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual) = split(/\t/, $bamLine);

		# Note where mate maps to.
		push (@{$mateMappings{$rnext}}, $pnext);	
		
		# Count number of reads associated with hotspot 
		$n++;

		# Determine orientation
		if (($flag & 16) == 16) {
			$isMinus++;
			}	
		else {
			$isPlus++;
			}	

		# Collect mapq scores for later to calculate average.
		push(@mapqs, $mapq);

		# Determine whether read is soft clipped and by how much - collect lengths for later to calculate average
		my $length = 0;
		# For now just consider 3' end
		if ($cigar =~ /(\d+)S$/) {
			$length = $1;
			}

		push(@softClipLengths, $length);

		}
	close (PIPE);

	# Estimate number of errors.
	my $totalBases = 0;
	my $totalErrors = 0;
	open (PIPE2, "samtools view -b $bamfile $chr:$start-$end | samtools mpileup -A -f $reffile - 2> /dev/null |") or die;
	while (<PIPE2>) {
		my ($chr, $pos, $ref, $cov, $string) = split(/\t/, $_);
		$totalBases += $cov;
		$totalErrors += countErrors($string, $cov);	
		}
	close (PIPE2) or die;


	printf(
		"%s:%d-%d\t%.2f\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n", 
		$chr,
		$start,
		$end,
		$meanCoverage,
		$size,
		$n,
		$isPlus / $n,
		$isMinus / $n,
		&orientationFreq($isPlus, $isMinus),
		&mean(\@mapqs),
		&mean(\@softClipLengths),
		$totalErrors/$totalBases * 100,
		&processMateMappings(\%mateMappings)
		);


	}
close (BED) or die "ERROR: could not close $bedfile.\n";
#------------------------------------------------------------------------------

###############################################################################

sub processMateMappings {
	
	my $mateMappings = shift;

	my @buffer;
	foreach my $chr (sort keys %$mateMappings) {
				
		my $positions = $mateMappings->{$chr};
		my @sorted = sort {$a <=> $b} @$positions;
	
		my $start = $sorted[0];
		my $end = $sorted[@sorted-1];

		push(@buffer, "$chr:$start-$end\t" . ($end-$start+1) . "");		

		}

	return join("\t", @buffer);

	} # End of method.

###############################################################################

sub orientationFreq {
	
	my ($largest, $smallest) = sort {$b <=> $a} @_;
	my $n = $largest + $smallest;
	return $largest/ $n;
	
	} # End of method.

###############################################################################

sub countErrors {
	
	my $string = shift;
	my $cov = shift;

	chomp($string);

	my @buffer = split(//, $string);
	
	my $n = 0;
	my $error = 0;
	foreach (@buffer) {
		$n++ if /[\.\,ACGTN]/i;
		$error ++ if /[ACGTN]/i;
		}

	if ($cov != $n) {
		#warn "WARNING: coverage $cov != $n: $string\n";
		}

	return $error;

	} # End of method.

###############################################################################

sub mean {
	
	my $data = shift;
	
	my $n = 0;
	my $sum = 0;
	foreach (@$data) {	
		$sum += $_;
		$n++;
		}

	return $sum / $n;

	} # End of method.

###############################################################################
1;
