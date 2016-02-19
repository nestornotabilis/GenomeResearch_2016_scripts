#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Basename;

my %options;
getopts("b:i:",\%options);
my $bamfile = $options{b} or &usage;
my $infile = $options{i} or &usage;
sub usage {die "USAGE: " . basename($0) . " [-b bam file] [-i hotspot table infile].\n";}

# Process each hotspot listed in table infile:
my $i = 0;
open (IN, $infile) or die "ERROR: could not open $infile.\n";
while (my $line = <IN>) {
        chomp($line);
        if ($line =~ /^Location/) {
                print "$line\tconsensus\n";
                next;
                }

	# Just interested in location.
        my ($location) = split(/\t/, $line);

	# Construct consensus
        my $consensus = processBAM($bamfile, $location, $i++);

	# Add to table
        print "$line\t$consensus\n";

        }
close (IN) or die "ERROR: could not close $infile.\n";

###############################################################################

sub processBAM {
	
	my $bamfile 	= shift;
	my $pos 	= shift;
	my $i 		= shift;

	my $TMP_FAS_FILE 	= "TMP.$i.fas";
        my $TMP_ALN_FILE 	= "TMP.$i.aln.fas";
        my $TMP_CONSENSUS_FILE	= "TMP.$i.consensus.fas";

	# Create TMP file for storing reads
	open (OUT, ">$TMP_FAS_FILE") or die "ERROR: could not open $TMP_FAS_FILE";

	# Extract reads from BAM
	my $readCount = 0;
	open (BAM, "samtools view $bamfile $pos |") or die "ERROR: could not open pipe.\n";
	while (<BAM>) {
		chomp;
		my @cells = split(/\t/, $_);
		my $id = $cells[0];
		my $read = $cells[9];
			
		print OUT ">$id\n";
		print OUT "$read\n";
		$readCount++;
		}
	close (BAM) or die "ERROR: could not close pipe.\n";
	close (OUT) or die "ERROR: could not close $TMP_FAS_FILE.\n"; 


	# Align extracted reads and create a TMP alignment file
	system("muscle -in $TMP_FAS_FILE -out $TMP_ALN_FILE -quiet");   

	# Now create TMP consensus file from alignment.
	# If no more than 5,000 reads (avoiding segmentation fault)
	if ($readCount > 1 && $readCount <= 5000) {
		system("consambig -sequence $TMP_ALN_FILE -outseq $TMP_CONSENSUS_FILE -name $pos\_consensus 2> /dev/null");
		}
	elsif ($readCount > 5000) {
		# Avoiding segmentation fault 
		# by ignoring very deep coverage
		#----------------------------
		unlink($TMP_FAS_FILE);
        	unlink($TMP_ALN_FILE);
        	unlink($TMP_CONSENSUS_FILE);
		warn "WARNING: Hotspot too deep at $readCount x; ignoring.\n";
		return "TOO DEEP TO CALCULATE";
		#----------------------------
		}
	# If only one read
	else {
		`cp $TMP_ALN_FILE $TMP_CONSENSUS_FILE`;
		}

	# Finally extract consensus from TMP consensus file.
	my $consensus;
	open (CONSENSUS, $TMP_CONSENSUS_FILE) or die "ERROR: could not open $TMP_CONSENSUS_FILE.\n";
	$/ = ">";
	<CONSENSUS>;
	while (<CONSENSUS>) {
        	chomp;
        	my ($header, @lines) = split(/\n/, $_);
        	$consensus = join("", @lines);
       		}
	$/ = "\n";
	close (CONSENSUS) or die "ERROR: could not close $TMP_CONSENSUS_FILE.\n";

	# Tidy up.
	unlink($TMP_FAS_FILE);
	unlink($TMP_ALN_FILE);
	unlink($TMP_CONSENSUS_FILE);

	# return consensus.
	return uc($consensus);

	} # End of method.

###############################################################################
1;
