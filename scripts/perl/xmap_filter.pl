#!/usr/bin/env perl
# Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences
# MOMS is licensed under the Mulan PSL v1.
# You can use this software according to the terms and conditions of the Mulan PSL v1.
# You may obtain a copy of Mulan PSL v1 at:
#    http://license.coscl.org.cn/MulanPSL
# THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND, EITHER EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT, MERCHANTABILITY OR FIT FOR A PARTICULAR
# PURPOSE.
# See the Mulan PSL v1 for more details.
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
use List::Util qw[min max];

BEGIN{
	my $progpath = abs_path(dirname($0));
	unshift @INC, $progpath;
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

use BNG::Utility;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for filtering alignments in a XMAP
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-x, --xmap <str>    The XMAP file containing the alignment information (REQUIRED)
	-c, --conf <float>  The mininum allowed value of confidence (default: 11)
	-d, --delta <float> The maximum allowed difference between the confidence with the maximum (default: 2.0)
	-r, --ratio <float> The minimum allowed ratio in terms of QryLen / max(QryLen) for the alignments with the same QryContigID (default: 0.9)
	-h, --help          Help

Example:
	$program -x alignment.xmap > alignment_filtered.xmap
USAGE

use vars qw($opt_x $opt_c $opt_d $opt_r $opt_h);
GetOptions( "x|xmap=s" => \$opt_x,
			"c|conf=f" => \$opt_c,
			"d|delta=f" => \$opt_d,
			"r|ratio=f" => \$opt_r,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

die("**ERROR: -x option must be specified\n") unless(defined $opt_x);

my $confi = (defined $opt_c) ? $opt_c : 11;
my $delta = (defined $opt_d) ? $opt_d : 2.0;
my $ratio = (defined $opt_r) ? $opt_r : 0.9;

my $xmap = readXMap($opt_x);
my $QryContigIDs = $xmap->{hits}->{QryContigID};
my $Confidences = $xmap->{hits}->{Confidence};
my $QryLens = $xmap->{hits}->{QryLen};
die("**ERROR: irregular content in $opt_x\n") unless($xmap->{totalHits} == scalar(@$QryContigIDs));

# sort the alignments by QryContigID
my @indexes = sort{ $QryContigIDs->[$a] <=> $QryContigIDs->[$b] } 0..$#{$QryContigIDs};

my ($i, $j, $startI, $ii, $jj);
my ($oldQryID, $qryID);
my ($maxConfidence, $confidence);
my ($maxQryLen, $QryLen);

# print out the headers
for($i=0; $i<scalar(@{$xmap->{headers}}); $i++){
	print "$xmap->{headers}->[$i]\n";
}

# processing the alignments
my @fields =@{ $xmap->{dataName} };
my $nField = scalar(@fields);
my $k;
our $no = 0;
my $confi_threshold;
for($oldQryID=-1, $i=0; $i<scalar(@$QryContigIDs); $i++){
	$j = $indexes[$i];
	$qryID = $QryContigIDs->[$j];
	$confidence = $Confidences->[$j];
	$QryLen = $QryLens->[$j];
	if($qryID != $oldQryID){
		if($oldQryID > 0){
			for($ii=$startI; $ii<$i; $ii++){
				$jj = $indexes[$ii];
				$confi_threshold = max($maxConfidence - $delta, $confi);
				next if($Confidences->[$jj] <= $confi_threshold);
				next if($QryLens->[$jj] < $maxQryLen * $ratio);
				&printRow($xmap, $jj);
			}
		}
		$maxConfidence = $confidence;
		$maxQryLen = $QryLen;
		$startI = $i; # set the start index
		$oldQryID = $qryID;
	}
	else{
		if($confidence > $maxConfidence){
			$maxConfidence = $confidence;
		}
		if($QryLen > $maxQryLen){
			$maxQryLen = $QryLen;
		}
	}
}

if($oldQryID > 0){
	for($ii=$startI; $ii<$i; $ii++){
		$jj = $indexes[$ii];
		next if($Confidences->[$jj] <= $maxConfidence - $delta);
		next if($QryLens->[$jj] < $maxQryLen * $ratio);
		&printRow($xmap, $jj);
	}
}

exit 0;

sub printRow
{
	my ($xmap, $i) = @_;
	our $no;
	my $fields =$xmap->{dataName};
	my $nField = scalar(@{$fields});
	for($k=0; $k<$nField; $k++){
		if($k == 0){
			print ++$no;
		} else {
			print $xmap->{hits}->{$fields[$k]}->[$i];
		}
		if($k < $nField - 1){
			print "\t";
		} else {
			print "\n";
		}
	}
}
