#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Basename;

my %options;
getopts("i:", \%options);
my $infile = $options{i} or &usage;
sub usage {die "USAGE: " . basename($0) . " [-i Scenario A, B or C Bam infile (only properly paired reads should be present)].\n";}

# Store PE-read information
my %reads;
open (PIPE, "java -jar /share/apps/WGP-Toolkit/ReadPositionsFromBAM.jar -b $infile |") or die "ERROR: could not open PIPE for $infile.\n";
while (<PIPE>) {
        chomp;
        my @line 	= split(/\t/, $_);
	my $id 		= $line[0];
        my $flag 	= $line[1];
        my $chr 	= $line[2];
        my $start 	= $line[3];
        my $end 	= $line[4];

	# First in pair
	if (($flag & 64) == 64) {
		$reads{$id}->{'first'}->{'ref'} = $chr;
		$reads{$id}->{'first'}->{'end'} = $end;	
		$reads{$id}->{'first'}->{'flag'} = $flag; 		
                }

        # Second in pair
        elsif (($flag & 128) == 128) {
		$reads{$id}->{'second'}->{'ref'} = $chr;
                $reads{$id}->{'second'}->{'end'} = $end;
		$reads{$id}->{'second'}->{'flag'} = $flag;
                }
        else {
                die "ERROR: unexpected flag = $flag\n";
                }

        }
close (PIPE) or die "ERROR: could not close pipe.\n";




my %data;
foreach my $id (sort keys %reads) {

	my $ref1 = $reads{$id}->{first}->{ref};
	my $pos1 = $reads{$id}->{first}->{end};
	my $ref2 = $reads{$id}->{second}->{ref};
	my $pos2 = $reads{$id}->{second}->{end};

	process($id, $ref1, $pos1, $ref2, $pos2, \%data);

	} # End of while loop.



foreach my $id (sort keys %data) {

	if (!defined $data{$id}->{first}->{pos}) {
		warn "WARNING: No r1 read for '$id'\n";
		next;
		}

	 if (!defined $data{$id}->{second}->{pos}) {
                warn "WARNING: No r2 read for '$id'\n";
                next;
                }

	printf(
		"%s %s %d %d z=0,color=%s\n",
		$id,		
		$data{$id}->{first}->{ref},
		$data{$id}->{first}->{pos},
		$data{$id}->{first}->{pos},
		lc($data{$id}->{first}->{ref}) # lowercase to conform to circos requiements.
		);

	printf(
                "%s %s %d %d z=0,color=%s\n",
                $id,
                $data{$id}->{second}->{ref},
                $data{$id}->{second}->{pos},
                $data{$id}->{second}->{pos},
                lc($data{$id}->{first}->{ref})  # lowercase to conform to circos requiements.
                );

	}


###############################################################################

sub process {

	my $id = shift;
	my $ref1 = shift;
	my $pos1 = shift;
	my $ref2 = shift;
	my $pos2 = shift;
	my $data = shift;
	
	if (exists $data->{$id}->{first}) {
		
		die "ERROR: $data->{$id}->{first}->{ref} ne $ref1 or $data->{$id}->{first}->{pos} != $pos1.\n" if $data->{$id}->{first}->{ref} ne $ref1 || $data->{$id}->{first}->{pos} != $pos1;

		} 
	else {
		$data->{$id}->{first}->{ref} = $ref1;
		$data->{$id}->{first}->{pos} = $pos1;
		} 


	if (exists $data->{second}) {
                
                die "ERROR: $data->{$id}->{second}->{ref} ne $ref2 or $data->{$id}->{second}->{pos} != $pos2.\n" if $data->{$id}->{second}->{ref} ne $ref2 || $data->{$id}->{second}->{pos} != $pos2;

                } 
        else {
                $data->{$id}->{second}->{ref} = $ref2;
                $data->{$id}->{second}->{pos} = $pos2;
                } 


	} # End of method.

###############################################################################
1;
