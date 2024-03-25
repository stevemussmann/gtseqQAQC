#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;
use Data::Dumper;

# kill program and print help if no command line arguments were given
if( scalar( @ARGV ) == 0 ){
	&help;
	die "Exiting program because no command line options were used.\n\n";
}

# take command line arguments
my %opts;
getopts( 'hf:o:q:', \%opts );

# if -h flag is used, or if no command line arguments were specified, kill program and print help
if( $opts{h} ){
	&help;
	die "Exiting program because help flag was used.\n\n";
}

# parse the command line
my( $orig, $out, $qaqc ) = &parsecom( \%opts );

my @qaqcArr;
my %qaqcHash;
my @origArr;
my %origHash;

# read input files
&filetoarray( $qaqc, \@qaqcArr );
&filetoarray( $orig, \@origArr );

# pull headers off input files
my @qaqcHead = split( /,/, shift( @qaqcArr ) );
my @origHead = split( /,/, shift( @origArr ) );

&getNames( \@qaqcArr, \%qaqcHash );
&getNames( \@origArr, \%origHash );

# capture missing values
my %qaqcMiss;
my %origMiss;
my %bothMiss;
my %origPercMiss;
my %qaqcPercMiss;

# capture percentage values
my @rawArr;
my @adjArr;
my @origPercMissArr;
my @qaqcPercMissArr;
my @locRawArr;
my @locAdjArr;

# capture individuals not found in original file
my @missingInds;

# capture counts of successfully genotyped loci
my %origLociTotal; #total individuals attempted to genotype
my %qaqcLociTotal;
my %origLociCount; #total individuals missing genotypes
my %qaqcLociCount;
my %locusMismatch;
my %locusCompare;
my %locusCompareMissOrig;
my %locusCompareMissQAQC;
my %locusCompareMissBoth;

# populate locus hashes
foreach my $q( sort keys %qaqcHash ){
	foreach my $l( sort keys %{$qaqcHash{$q}} ){
		$origLociTotal{$l}=0;
		$qaqcLociTotal{$l}=0;
		$origLociCount{$l}=0;
		$qaqcLociCount{$l}=0;
		$locusMismatch{$l}=0;
		$locusCompare{$l}=0;
		$locusCompareMissOrig{$l}=0;
		$locusCompareMissQAQC{$l}=0;
		$locusCompareMissBoth{$l}=0;
	}
}

my $indout = "$out.qaqc.individuals.txt";
open( OUT, '>', $indout ) or die "Can't open $indout: $!\n\n";

# write new header to print to file
print OUT "Individual\tMismatch_Loci\tOrig_Missing\tQAQC_Missing\tBoth_Missing\tAdjusted_Loci\tRaw_Loci\tAdjusted_%_Missing\tRaw_%_Missing\n";

# parse data
foreach my $q( sort keys %qaqcHash ){
	if( exists( $origHash{$q} ) ){
		#print "exists\n";
		my $mismatch=0; # initialize new mismatch value 
		my $count=0; # initialize count of loci genotyped in both 
		$qaqcMiss{$q} = 0;
		$origMiss{$q} = 0;
		$bothMiss{$q} = 0;
		$origPercMiss{$q} = 0;
		$qaqcPercMiss{$q} = 0;
		foreach my $l( sort keys %{$qaqcHash{$q}} ){
			# locus data
			if( exists( $origHash{$q}{$l} ) ){
				$origLociTotal{$l}++;
				if( $origHash{$q}{$l} eq "00" ){
					$origLociCount{$l}++;
				}
			}
			if( exists( $qaqcHash{$q}{$l} ) ){
				$qaqcLociTotal{$l}++;
				if( $qaqcHash{$q}{$l} eq "00" ){
					$qaqcLociCount{$l}++;
				}
			}

			# individual data
			if( exists( $origHash{$q}{$l} ) ){
				#print $l, "\n";
				if( $origHash{$q}{$l} ne "00" and $qaqcHash{$q}{$l} ne "00" ){
					$locusCompare{$l}++;
					$count++; #add to number of loci genotyped in both original and qaqc data
					if( $origHash{$q}{$l} ne $qaqcHash{$q}{$l} ){
						$mismatch++;
						$locusMismatch{$l}++;
					}
				}else{
					if( $origHash{$q}{$l} eq "00" and $qaqcHash{$q}{$l} eq "00" ){
						$bothMiss{$q}++;
						$locusCompareMissBoth{$l}++;
					}elsif( $origHash{$q}{$l} eq "00" ){
						$origMiss{$q}++;
						$locusCompareMissOrig{$l}++;
					}elsif( $qaqcHash{$q}{$l} eq "00" ){
						$qaqcMiss{$q}++;
						$locusCompareMissQAQC{$l}++;
					}
				}
			}
		}
		## total loci
		my $total = $count + $bothMiss{$q} + $origMiss{$q} + $qaqcMiss{$q};

		## calculate mismatch percents
		# raw mismatch percent
		my $raw = sprintf( "%.2f", (($mismatch/$total)*100) );
		push( @rawArr, $raw );
		$raw = join( "", $raw, "%" );

		# adjusted mismatch percent
		my $adj = "NA";
		if( $count != 0 ){
			$adj = sprintf( "%.2f", (($mismatch/$count)*100) );
			push( @adjArr, $adj );
			$adj = join( "", $adj, "%");
		}
		
		## missing data
		# calculate original file missing data percents
		my $origMissData = $origMiss{$q} + $bothMiss{$q};
		$origPercMiss{$q} = sprintf("%.2f", (($origMissData / $total)*100) );
		push( @origPercMissArr, $origPercMiss{$q} );

		# calculate qaqc file missing data percents
		my $qaqcMissData = $qaqcMiss{$q} + $bothMiss{$q};
		$qaqcPercMiss{$q} = sprintf("%.2f", (($qaqcMissData / $total)*100) );
		push( @qaqcPercMissArr, $qaqcPercMiss{$q} );


		print OUT $q, "\t", $mismatch, "\t", $origMiss{$q}, "\t", $qaqcMiss{$q}, "\t", $bothMiss{$q}, "\t", $count, "\t", $total, "\t", $adj, "\t", $raw, "\n";
	}else{
		push( @missingInds, $q );
		print OUT $q, "\t", "NA", "\t", "NA", "\t", "NA", "\t", "NA", "\t", "NA", "\t", "NA", "\t", "NA", "\t", "NA", "\n";
	}
}

close OUT;

# print out summaries for individual loci
my $locout = "$out.qaqc.loci.txt";
open( OUT, '>', $locout ) or die "Can't open $locout: $!\n\n";
print OUT "Locus\tMismatch_Individuals\tOrig_Missing\tQAQC_Missing\tBoth_Missing\tAdjusted_Individuals\tRaw_Individuals\tAdjusted_%_Missing\tRaw_%_Missing\n";
foreach my $l( sort keys %qaqcLociTotal ){
	my $locAdj;
	if( $locusCompare{$l} == 0 ){
		$locAdj = "NA";
	}else{
		$locAdj = sprintf("%.2f", ($locusMismatch{$l}/$locusCompare{$l})*100);
		push( @locAdjArr, $locAdj );
		$locAdj = join( "", $locAdj, "%" );
	}
	
	my $locusTot = $locusCompare{$l}+$locusCompareMissOrig{$l}+$locusCompareMissQAQC{$l}+$locusCompareMissBoth{$l};
	my $locTot;
	if( $locusTot == 0 ){
		$locTot = "NA";
	}else{
		$locTot = sprintf("%.2f", ($locusMismatch{$l}/$locusTot)*100);
		push( @locRawArr, $locTot );
		$locTot = join( "", $locTot, "%" );
	}
	print OUT $l, "\t", $locusMismatch{$l}, "\t", $locusCompareMissOrig{$l}, "\t", $locusCompareMissQAQC{$l}, "\t", $locusCompareMissBoth{$l}, "\t", $locusCompare{$l}, "\t", $locusTot,"\t", $locAdj, "\t", $locTot, "\n";
}

close OUT;

# calculate individual based stats
my $rawMean = &calcMean( \@rawArr );
my $adjMean = &calcMean( \@adjArr );
my $origMissMean = &calcMean( \@origPercMissArr );
my $qaqcMissMean = &calcMean( \@qaqcPercMissArr );

my $rawStdev = &calcStdev( \@rawArr, $rawMean );
my $adjStdev = &calcStdev( \@adjArr, $adjMean );
my $origMissStdev = &calcStdev( \@origPercMissArr, $origMissMean );
my $qaqcMissStdev = &calcStdev( \@qaqcPercMissArr, $qaqcMissMean );

my $rawMedian = &calcMedian( \@rawArr );
my $adjMedian = &calcMedian( \@adjArr );
my $origMissMedian = &calcMedian( \@origPercMissArr );
my $qaqcMissMedian = &calcMedian( \@qaqcPercMissArr );

my $rawMin = &calcMin( \@rawArr );
my $adjMin = &calcMin( \@adjArr );
my $origMissMin = &calcMin( \@origPercMissArr );
my $qaqcMissMin = &calcMin( \@qaqcPercMissArr );

my $rawMax = &calcMax( \@rawArr );
my $adjMax = &calcMax( \@adjArr );
my $origMissMax = &calcMax( \@origPercMissArr );
my $qaqcMissMax = &calcMax( \@qaqcPercMissArr );

# calculate locus based stats
my $locRawMean = &calcMean( \@locRawArr );
my $locAdjMean = &calcMean( \@locAdjArr );

my $locRawStdev = &calcStdev( \@locRawArr, $locRawMean );
my $locAdjStdev = &calcStdev( \@locAdjArr, $locAdjMean );

my $locRawMedian = &calcMedian( \@locRawArr );
my $locAdjMedian = &calcMedian( \@locAdjArr );

my $locRawMin = &calcMin( \@locRawArr );
my $locAdjMin = &calcMin( \@locAdjArr );

my $locRawMax = &calcMax( \@locRawArr );
my $locAdjMax = &calcMax( \@locAdjArr );

print("\nMISSING INDIVIDUAL SUMMARY:\n");
# print missing individuals
if( scalar(@missingInds>0) ){
	print("WARNING - The following individuals from the QAQC file ($qaqc) were not found among the original genotypes file ($orig):\n");
	foreach my $ind( @missingInds ){
		print $ind, "\n";
	}
}else{
	print("All individuals from the QAQC file ($qaqc) were found in the original genotype file ($orig).\n")
}

print("\nMISSING DATA SUMMARIES:\n");
# print missing data stats for original and qaqc genotypes
print("Missing data stats for QAQC'd individuals from original genotype file ($orig):\n");
$origMissMean = sprintf("%.2f", $origMissMean);
$origMissStdev = sprintf("%.2f", $origMissStdev);
print "Mean = ", $origMissMean, "%\n";
print "StDev = ", $origMissStdev, "%\n";
print "Median = ", $origMissMedian, "%\n";
print "Min = ", $origMissMin, "%\n";
print "Max = ", $origMissMax, "%\n";

print("\nMissing data stats for QAQC genotype file ($qaqc):\n");
$qaqcMissMean = sprintf("%.2f", $qaqcMissMean);
$qaqcMissStdev = sprintf("%.2f", $qaqcMissStdev);
print "Mean = ", $qaqcMissMean, "%\n";
print "StDev = ", $qaqcMissStdev, "%\n";
print "Median = ", $qaqcMissMedian, "%\n";
print "Min = ", $qaqcMissMin, "%\n";
print "Max = ", $qaqcMissMax, "%\n";

print("\nMISMATCHED GENOTYPES SUMMARIES (Individuals):\n");
# print stats for raw mismatching genotypes in individuals
print("Stats for raw percentage of mismatched genotypes (individuals):\n");
$rawMean = sprintf("%.2f", $rawMean);
$rawStdev = sprintf("%.2f", $rawStdev);
$rawMedian = sprintf("%.2f", $rawMedian);
$rawMin = sprintf("%.2f", $rawMin);
$rawMax = sprintf("%.2f", $rawMax);
print "Mean = ", $rawMean, "%\n";
print "StDev = ", $rawStdev, "%\n";
print "Median = ", $rawMedian, "%\n";
print "Min = ", $rawMin, "%\n";
print "Max = ", $rawMax, "%\n";

# print stats for adjusted mismatching genotypes in individuals
print("\nStats for adjusted percentage of mismatched genotypes (individuals):\n");
$adjMean = sprintf("%.2f", $adjMean);
$adjStdev = sprintf("%.2f", $adjStdev);
$adjMedian = sprintf("%.2f", $adjMedian);
$adjMin = sprintf("%.2f", $adjMin);
$adjMax = sprintf("%.2f", $adjMax);
print "Mean = ", $adjMean, "%\n";
print "StDev = ", $adjStdev, "%\n";
print "Median = ", $adjMedian, "%\n";
print "Min = ", $adjMin, "%\n";
print "Max = ", $adjMax, "%\n";

print("\nMISMATCHED GENOTYPES SUMMARIES (Loci):\n");
# print stats for raw mismatching genotypes in Loci
print("Stats for raw percentage of mismatched genotypes (loci):\n");
$locRawMean = sprintf("%.2f", $locRawMean);
$locRawStdev = sprintf("%.2f", $locRawStdev);
$locRawMedian = sprintf("%.2f", $locRawMedian);
$locRawMin = sprintf("%.2f", $locRawMin);
$locRawMax = sprintf("%.2f", $locRawMax);
print "Mean = ", $locRawMean, "%\n";
print "StDev = ", $locRawStdev, "%\n";
print "Median = ", $locRawMedian, "%\n";
print "Min = ", $locRawMin, "%\n";
print "Max = ", $locRawMax, "%\n";

# print stats for adjusted mismatching genotypes in loci
print("\nStats for adjusted percentage of mismatched genotypes (loci):\n");
$locAdjMean = sprintf("%.2f", $locAdjMean);
$locAdjStdev = sprintf("%.2f", $locAdjStdev);
$locAdjMedian = sprintf("%.2f", $locAdjMedian);
$locAdjMin = sprintf("%.2f", $locAdjMin);
$locAdjMax = sprintf("%.2f", $locAdjMax);
print "Mean = ", $locAdjMean, "%\n";
print "StDev = ", $locAdjStdev, "%\n";
print "Median = ", $locAdjMedian, "%\n";
print "Min = ", $locAdjMin, "%\n";
print "Max = ", $locAdjMax, "%\n";


print("\nIndividual-specific data printed to file: $indout.\n\n");
print("Locus-specific data printed to file: $locout.\n\n");


#print Dumper( \%qaqcHash );
#print Dumper( \%origHash );

exit;
#####################################################################################################
############################################ Subroutines ############################################
#####################################################################################################

# subroutine to print help
sub help{

  print "\ngtseqQAQC.pl is a perl script developed by Steve Mussmann\n\n";
  print "To report bugs send an email to mussmann\@uark.edu\n";
  print "When submitting bugs please include all input files, options used for the program, and all error messages that were printed to the screen\n\n";
  print "Program Options:\n";
  print "\t\t[ -h | -f | -o | -q ]\n\n";
  print "\t-h:\tUse this flag to display this help message.\n";
  print "\t\tThe program will die after the help message is displayed.\n\n";
  print "\t-f:\tSpecify your .csv file of original genotypes output from Nate Campbell's gtseq pipeline.\n\n";
  print "\t-o:\tSpecify the output file prefix.\n\n";
  print "\t-q:\tSpecify the .csv file of QAQC genotypes output from Nate Campbell's gtseq pipeline.\n\n";

}

#####################################################################################################
# subroutine to parse the command line options

sub parsecom{

  my( $params ) =  @_;
  my %opts = %$params;

  # set default values for command line arguments
  my $str = $opts{f} || die "No input file of original genotypes specified.\n\n"; #used to specify original .csv input file name.
  my $out = $opts{o} || "gtseqQAQCoutput"  ; #used to specify output file prefix.
  my $map = $opts{q} || die "No input file of qaqc genotypes specified.\n\n"; #used to specify qaqc .csv input file name

  return( $str, $out, $map );

}

#####################################################################################################
# subroutine to put file into an array

sub filetoarray{

  my( $infile, $array ) = @_;


  # open the input file
  open( FILE, $infile ) or die "Can't open $infile: $!\n\n";

  # loop through input file, pushing lines onto array
  while( my $line = <FILE> ){
    chomp( $line );
    next if($line =~ /^\s*$/);
    push( @$array, $line );
  }
  
  # close input file
  close FILE;

}

#####################################################################################################
# subroutine to get sample names from data lines in array

sub getNames{

	my( $array, $hash ) = @_;

	foreach my $line( @$array ){
		my @temp = split( /,/, $line );
		my @name = split( /_/, $temp[0] );
		my $samp = $name[0];
		for( my $i=6; $i<@temp; $i++ ){
			$$hash{$samp}{$origHead[$i]} = $temp[$i];
		}
	}

}

#####################################################################################################
# subroutine to calculate mean of an array

sub calcMean{

	my( $array ) = @_;

	my $total = 0;
	my $count = scalar( @$array );

	foreach my $item( @$array ){
		$total+=$item;
	}
	my $mean = $total/$count;

	return $mean;

}

#####################################################################################################
# subroutine to find median of an array

sub calcMedian{

	my( $array ) = @_;

	my @sorted = sort{ $a <=> $b } @$array;

	my $len = scalar(@sorted);
	my $i = (($len-1)/2);

	my $median = 0;

	if( ($len % 2) == 0){
		$median = $sorted[$i];
	}else{
		$median = (($sorted[$i] + $sorted[$i+1])/2);
	}

	return $median;

}

#####################################################################################################
# subroutine to find standard deviation of an array

sub calcStdev{

	my( $array, $mean ) = @_;

	my @devArr;

	my $stdev = 0;

	foreach my $val( @$array ){
		my $dev = ($val - $mean)**2;
		push( @devArr, $dev );
	}

	my $sum = 0;
	foreach my $val( @devArr ){
		$sum+=$val;
	}

	if( (scalar(@$array)-1)>0 ){
		my $temp = $sum/(scalar(@$array)-1);
		$stdev = sqrt($temp);
	}

	return $stdev;

}

#####################################################################################################
# subroutine to find minimum value of an array

sub calcMin{

	my( $array ) = @_;

	my @sorted = sort{ $a <=> $b } @$array;

	return shift(@sorted);

}

#####################################################################################################
# subroutine to find maximum value of an array

sub calcMax{

	my( $array ) = @_;

	my @sorted = sort{ $a <=> $b } @$array;

	return pop(@sorted);

}

#####################################################################################################
