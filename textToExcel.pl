#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use Excel::Writer::XLSX;

my %options;
getopts("i:o:", \%options);
my $infile = $options{i} or &usage;
my $outfile = $options{o} or &usage;
sub usage {die "USAGE: " . basename($0) . " [-i text infile (.txt)] [-o Excel outfile (.xlsx)].\n";}


# Create a new Excel workbook
my $workbook = Excel::Writer::XLSX->new($outfile);

$workbook->set_tempdir( '/tmp' );

# Add a worksheet
my $worksheet = $workbook->add_worksheet();

#  Header format
my $headerFormat = $workbook->add_format();
$headerFormat->set_bold();
$headerFormat->set_border(1);

#  Cell format
my $cellFormat = $workbook->add_format();
$cellFormat->set_border(1);

# Process infile
open (IN, $infile) or die "ERROR: Could not open '$infile'.\n";

my $rowNumber = 0;
while (my $row = <IN>) {
	
	my @cells = split(/\t/, $row);

	my $columnNumber = 0;
	foreach my $cell (@cells) {
		if ($rowNumber == 0) {
			$worksheet->write($rowNumber, $columnNumber, $cell, $headerFormat);
			}
		else {
			$worksheet->write($rowNumber, $columnNumber, $cell, $cellFormat);
			}
		$columnNumber++;
		}
	$rowNumber++;
	}


close(IN) or die "ERROR: could not close $infile.\n";

