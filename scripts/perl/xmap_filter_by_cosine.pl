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
$program: A perl script for filtering alignments in a XMAP by cosine
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-x, --xmap <str>   The XMAP file containing the alignment information (REQUIRED)
	-r, --ref <str>    The reference CMAP file (REQUIRED)
	-q, --qry <str>    The query CMAP file (REQUIRED)
	-h, --help         Help

Example:
	$program -x aligned.xmap -r multicolor.cmap -q unused_contigs_multicolor.cmap > filtered.xmap
USAGE

use vars qw($opt_x $opt_r $opt_q $opt_h);
GetOptions( "x|xmap=s" => \$opt_x,
			"r|ref=s" => \$opt_r,
			"q|qry=s" => \$opt_q,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

die("**ERROR: -x option must be specified.\n") unless(defined $opt_x);
die("**ERROR: -r option must be specified.\n") unless(defined $opt_r);
die("**ERROR: -q option must be specified.\n") unless(defined $opt_q);

my $xmap = readXMap($opt_x);
my ($refCmap, $refNum, $refLength) = readCMap($opt_r);
my ($qryCmap, $qryNum, $qryLength) = readCMap($opt_q);

my $QryContigIDs = $xmap->{hits}->{QryContigID};
my $RefContigIDs = $xmap->{hits}->{RefContigID};
my $QryStartPos = $xmap->{hits}->{QryStartPos};
my $QryEndPos = $xmap->{hits}->{QryEndPos};
my $RefStartPos = $xmap->{hits}->{RefStartPos};
my $RefEndPos = $xmap->{hits}->{RefEndPos};

# print out the headers
for(my $i=0; $i<scalar(@{$xmap->{headers}}); $i++){
	print "$xmap->{headers}->[$i]\n";
}
my $data_names = $xmap->{"dataName"};
my $numc = scalar(@$data_names);
my ($labels1, $labels2);
my ($channels1, $channels2, $union_channels);
my ($vec1, $vec2, $vec3, $vec4);
my $cos;
my $max = max(scalar(keys %{$refCmap->{"channels"}}), scalar(keys %{$qryCmap->{"channels"}}));

my $k = 0;
for(my $i=0; $i<$xmap->{totalHits}; $i++){
	($labels1, $channels1) = retrieveLabels($qryCmap, $QryContigIDs->[$i], $QryStartPos->[$i], $QryEndPos->[$i]);
	($labels2, $channels2) = retrieveLabels($refCmap, $RefContigIDs->[$i], $RefStartPos->[$i], $RefEndPos->[$i]);
	$union_channels = union($channels1, $channels2);
	($vec1, $vec3) = calculateEigenVectors($labels1, $max);
	($vec2, $vec4) = calculateEigenVectors($labels2, $max);
	if(!compatible($vec1, $vec2)){
		next;
	}
	$cos = cosine($vec3, $vec4);
	if($cos < 0.1){
		next;
	}
	$k++;
	print "$k\t";
	for (my $j=1; $j < $numc; $j++) {
		print $xmap->{hits}->{$data_names->[$j]}->[$i];
		if ($j < $numc - 1) {
			print "\t";
		}else{
			print "\n";
		}
	}
}

exit 0;

sub retrieveLabels
{
	my ($cmap, $id, $start, $end) = @_;
	my @labels = ();
	my %channels = ();
	my $ctg = $cmap->{contigs}->{$id};
	my $numsites = $ctg->{NumSites};
	if($start > $end){
		my $tmp = $start; $start = $end; $end = $tmp;
	}
	my $i;
	my $positions = $ctg->{Position};
	my $labelchannels = $ctg->{LabelChannel};
	for($i=0; $i<$numsites; $i++){
		my $pos = $positions->[$i];
		if($pos >= $start && $pos <= $end){
			my $ichannel = $labelchannels->[$i];
			push(@labels, $ichannel);
			$channels{$ichannel} = 1;
		}
	}
	return (\@labels, \%channels);
}

sub union
{
	my ($channels1, $channels2) = @_;
	my %union = ();
	foreach (keys %$channels1, keys %$channels2){
		$union{$_} = 1;
	}
	return \%union;
}

sub calculateEigenVectors
{
	my ($labels, $max) = @_;
	my @vector = (0) x $max;
	my @vector2 = (0) x (($max * ($max+1)) / 2);
	my ($idx, $idx1, $idx2);
	for(my $i=0; $i<=$#$labels; $i++){
		$idx = $labels->[$i];
		$vector[$idx-1]++;
		if($i < $#$labels){
			if($labels->[$i] >= $labels->[$i+1]){
				$idx1 = $labels->[$i];
				$idx2 = $labels->[$i+1];
			}
			else{
				$idx1 = $labels->[$i+1];
				$idx2 = $labels->[$i];
			}
			$vector2[($idx1-1)*$idx1/2+$idx2-1]++;
		}
	}

	return (\@vector, \@vector2);
}

sub compatible
{
	my ($vec1, $vec2) = @_;
	if(scalar(@$vec1) != scalar(@$vec2)){
		return 0;
	}
	for(my $i=0; $i<scalar(@$vec1); $i++){
		if( ($vec1->[$i] == 0) ^ ($vec2->[$i] == 0) ){
			return 0;
		}
	}
	return 1;
}

sub cosine
{
	my ($vec1, $vec2) = @_;
	die("vector sizes do not match\n") if(scalar(@$vec1) != scalar(@$vec2));
	return dotProduct($vec1, $vec2) / sqrt(dotProduct($vec1, $vec1) * dotProduct($vec2, $vec2));
}

sub dotProduct
{
	my ($vec1, $vec2) = @_;
	die("vector sizes do not match\n") if(scalar(@$vec1) != scalar(@$vec2));
	my $product = 0;
	for(my $i=0; $i<scalar(@$vec1); $i++){
		$product += $vec1->[$i] * $vec2->[$i];
	}
	return $product;
}
