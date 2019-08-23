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
use XML::Simple;

BEGIN{
	my $progpath = abs_path(dirname($0));
	unshift @INC, $progpath;
	select(STDERR);	$| = 1;
	select(STDOUT);	$| = 1;
}

use BNG::Utility;

my $program=basename($0);
my $usage = << "USAGE";
$program: A perl script for classifying the alignments between CMAPs
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-i, --xmap <str>      The XMAP file containing the alignment information (REQUIRED)
	-r, --reference <str> The reference CMAP file for alignment (REQUIRED)
	-q, --query <str>     The query CMAP file for alignment (REQUIRED)
	-r0, --ref0  <str>    The original reference CMAP file for alignment (REQUIRED)
	-q0, --qry0  <str>    The original query CMAP file for alignment (REQUIRED)
	-o, --output <str>    The output path (REQUIRED)
	-h, --help            Help

Options for invoking \$bionano/binary/RefAligner:
	-x, --xml <str>     Parameter file in XML format (default: \$bionano/xml/conflictsArguments.xml)

Example:
	$program -i align.xmap -r align_r.cmap -q align_q.cmap -r0 ngs_assembly.cmap -q0 bng_assembly.cmap -o classfied
USAGE

use vars qw($opt_i $opt_r $opt_q $opt_r0 $opt_q0 $opt_o $opt_x $opt_h);
GetOptions( "i|xmap=s" => \$opt_i,
			"r|reference=s" => \$opt_r,
			"q|query=s" => \$opt_q,
			"r0|ref0=s" => \$opt_r0,
			"q0|qry0=s" => \$opt_q0,
			"o|output=s" => \$opt_o,
			"x|xml=s" => \$opt_x,
			"h|help"  => \$opt_h);

die($usage) if($opt_h);

my $bionano = abs_path(dirname($0) . "/../bionano");
my $optxml = (defined $opt_x) ? $opt_x : "$bionano/xml/conflictsArguments.xml";
die("**ERROR: parameter file \"$optxml\" does not exist") unless(-f $optxml);

# check the required arguments
die("**ERROR: -i option must be specified\n") if(!defined $opt_i);
die("**ERROR: -r option must be specified\n") if(!defined $opt_r);
die("**ERROR: -q option must be specified\n") if(!defined $opt_q);
die("**ERROR: -r0 option must be specified\n") if(!defined $opt_r0);
die("**ERROR: -q0 option must be specified\n") if(!defined $opt_q0);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);

my ($xmap_in, $ref_cmap_in, $qry_cmap_in, $ref0_cmap_in, $qry0_cmap_in) = ($opt_i, $opt_r, $opt_q, $opt_r0, $opt_q0);

# check the input files
die("**ERROR: XMAP file \"$xmap_in\" does not exist") unless(-f $xmap_in);
die("**ERROR: reference cmap \"$ref_cmap_in\" does not exist") unless(-f $ref_cmap_in);
die("**ERROR: query cmap \"$qry_cmap_in\" does not exist") unless(-f $qry_cmap_in);
die("**ERROR: original reference cmap \"$ref0_cmap_in\" does not exist") unless(-f $ref0_cmap_in);
die("**ERROR: original query cmap \"$qry0_cmap_in\" does not exist") unless(-f $qry0_cmap_in);

# set the output directory and prefix
my ($outdir, $outpre);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = $opt_o;
	$outpre = "classified";
}
else{
	$outdir = dirname($opt_o);
	$outpre = basename($opt_o);
}

# read in the parameter configuration
my $XML = new XML::Simple(KeyAttr=>[]);
my $configs = $XML->XMLin($optxml);
my $paraDict = parseConfig($configs, "conflicts_identify");
my ($T_cutoff, $max_overhang) = (-log($paraDict->{'T_cutoff'}{val})/log(10), $paraDict->{'max_overhang'}{val});

my $xmap = readXMap($xmap_in);
my ($ref_cmap) = readCMap($ref_cmap_in);
my ($qry_cmap) = readCMap($qry_cmap_in);
my ($orig_ref_cmap) = readCMap($ref0_cmap_in);
my ($orig_qry_cmap) = readCMap($qry0_cmap_in);

my ($stickyXMap, $filtered_ref_cmap, $filtered_qry_cmap, $breakpoints) = &createStickXMap($xmap, $ref_cmap, $qry_cmap, $orig_ref_cmap, $orig_qry_cmap, $T_cutoff, $max_overhang);

# output the result
&outputResults($outdir, $outpre, $stickyXMap, $filtered_ref_cmap, $filtered_qry_cmap, $breakpoints);

exit;

sub outputResults
{
	my ($outdir, $outpre, $xmap, $filtered_ref, $filtered_qry, $breakpoints) = @_;

	mkdir $outdir if (! -e $outdir);

	writeXMapFile($xmap, "$outdir/$outpre.xmap");
	writeCMapFile($filtered_ref, "$outdir/${outpre}_r.cmap");
	writeCMapFile($filtered_qry, "$outdir/${outpre}_q.cmap");

	&outputBreakPoints($breakpoints, "$outdir/${outpre}_conflicts.txt");
}

sub createStickXMap
{
	my ($xmap, $ngs_cmap, $bn_cmap, $orig_ngs_cmap, $orig_bn_cmap, $T_cutoff, $max_overhang) = @_;

	my $stickyXMap = {};
	my @breakpoints = ();

	my @refLeftLabCnt = ();
	my @refRightLabCnt = ();
	my @qryLeftLabCnt = ();
	my @qryRightLabCnt = ();

	my ($numL, $numR);
	my $shiftbp = 20;
	for(my $i=0; $i<$xmap->{totalHits}; $i++){
		($numL, $numR) = countOverhangLabels($ngs_cmap, $xmap->{hits}->{RefContigID}->[$i],
					 $xmap->{hits}->{RefStartPos}->[$i] - $shiftbp, $xmap->{hits}->{RefEndPos}->[$i] + $shiftbp);
		push(@refLeftLabCnt, {numLabel => $numL, bkpt => $xmap->{hits}->{RefStartPos}->[$i]});
		push(@refRightLabCnt, {numLabel => $numR, bkpt => $xmap->{hits}->{RefEndPos}->[$i]});

		if( $xmap->{hits}->{Orientation}->[$i] eq "+"){
			($numL, $numR) = countOverhangLabels($bn_cmap, $xmap->{hits}->{QryContigID}->[$i],
					$xmap->{hits}->{QryStartPos}->[$i] - $shiftbp, $xmap->{hits}->{QryEndPos}->[$i] + $shiftbp);
		}
		else{
			($numR, $numL) = countOverhangLabels($bn_cmap, $xmap->{hits}->{QryContigID}->[$i],
					$xmap->{hits}->{QryEndPos}->[$i] - $shiftbp, $xmap->{hits}->{QryStartPos}->[$i] + $shiftbp);
		}
		push(@qryLeftLabCnt, {numLabel => $numL, bkpt => $xmap->{hits}->{QryStartPos}->[$i]});
		push(@qryRightLabCnt, {numLabel => $numR, bkpt => $xmap->{hits}->{QryEndPos}->[$i]});
	}

	$stickyXMap->{headers} = $xmap->{headers};
	$stickyXMap->{dataName} = $xmap->{dataName};
	$stickyXMap->{dataType} = $xmap->{dataType};
	$stickyXMap->{version} = $xmap->{version};
	my @xmap_data_name = @{$xmap->{dataName}};
	my $numc = scalar(@xmap_data_name);
	for(my $i=0; $i<$numc; $i++){
		$stickyXMap->{hits}->{$xmap_data_name[$i]} = [];
	}
	my $stickyNGSContigs = {};
	my $stickyBNContigs = {};
	my $tg = 0;
	for(my $i=0; $i<$xmap->{totalHits}; $i++){
		if( ($xmap->{hits}->{Confidence}->[$i] >= $T_cutoff) && (
			($refLeftLabCnt[$i]{numLabel} > $max_overhang && $qryLeftLabCnt[$i]{numLabel} > $max_overhang) ||
			($refRightLabCnt[$i]{numLabel} > $max_overhang && $qryRightLabCnt[$i]{numLabel} > $max_overhang) ) ){
			#significant overhang
			$tg++;
			for(my $m=0; $m<$numc; $m++){
				my $a = $stickyXMap->{hits}->{$xmap_data_name[$m]};
				my $b = $xmap->{hits}->{$xmap_data_name[$m]};
				my $c = $b->[$i];
				push(@$a, $c);
				$stickyXMap->{hits}->{$xmap_data_name[$m]}=$a;
			}
			$stickyNGSContigs->{ $xmap->{hits}->{RefContigID}->[$i] } = 1;
			$stickyBNContigs->{ $xmap->{hits}->{QryContigID}->[$i] } = 1;

			my ($leftRefBkpt, $leftQryBkpt) = ($refLeftLabCnt[$i]{numLabel} > $max_overhang && $qryLeftLabCnt[$i]{numLabel} > $max_overhang) ? ($refLeftLabCnt[$i]{bkpt}, $qryLeftLabCnt[$i]{bkpt}) : (-1, -1);
			my ($rightRefBkpt, $rightQryBkpt) = ($refRightLabCnt[$i]{numLabel} > $max_overhang && $qryRightLabCnt[$i]{numLabel} > $max_overhang) ? ($refRightLabCnt[$i]{bkpt}, $qryRightLabCnt[$i]{bkpt}) : (-1, -1);
			push(@breakpoints, {xmapId => $xmap->{hits}->{XmapEntryID}->[$i], refId => $xmap->{hits}->{RefContigID}->[$i], qryId => $xmap->{hits}->{QryContigID}->[$i], leftRefBkpt => $leftRefBkpt, leftQryBkpt => $leftQryBkpt, rightRefBkpt => $rightRefBkpt, rightQryBkpt => $rightQryBkpt, alignOrientation => $xmap->{hits}->{Orientation}->[$i]});
		}
	}
	$stickyXMap->{totalHits} = $tg;

	my $filtered_ngs_cmap = getSubsetCMap($orig_ngs_cmap, $stickyNGSContigs);
	my $filtered_bionano_cmap = getSubsetCMap($orig_bn_cmap, $stickyBNContigs);

	return ($stickyXMap, $filtered_ngs_cmap, $filtered_bionano_cmap, \@breakpoints);
}

sub outputBreakPoints
{
	my ($breakpoints, $filename) = @_;
	open(OUT, ">$filename") or die("Can not open \"$filename\" for writing.\n");

	my $shiftbp = 30;
	# shift certain amount of bp to avoid losing of the last/first label of the alignment
	# this is okay because there must be at least $max_overhang unaligned labels beyond the aligned region
	print OUT "# xMapId\trefQry\trefId\tleftRefBkpt\trightRefBkpt\talignmentOrientation\trefQry\tqryId\tleftQryBkpt\trightQryBkpt\talignmentOrientation\n";
	for(my $i=0; $i<scalar(@{$breakpoints}); $i++){
		my $bref = $breakpoints->[$i];
		my $orientation = $bref->{alignOrientation};
		my ($leftRefBkpt, $rightRefBkpt) = ($bref->{leftRefBkpt}, $bref->{rightRefBkpt});
		$leftRefBkpt -= $shiftbp unless($leftRefBkpt == -1);
		$rightRefBkpt += $shiftbp unless($rightRefBkpt == -1);
		my ($leftQryBkpt, $rightQryBkpt) = ($bref->{leftQryBkpt}, $bref->{rightQryBkpt});
		$leftQryBkpt = (($orientation eq "+") ? ($leftQryBkpt - $shiftbp) : ($leftQryBkpt + $shiftbp)) unless($leftQryBkpt == -1);
		$rightQryBkpt = (($orientation eq "+") ? ($rightQryBkpt + $shiftbp) : ($rightQryBkpt - $shiftbp)) unless($rightQryBkpt == -1);
		print OUT join("\t", $bref->{xmapId}, "ref", $bref->{refId}, $leftRefBkpt, $rightRefBkpt, $orientation, "qry", $bref->{qryId}, $leftQryBkpt, $rightQryBkpt, $orientation) . "\n";
	}

	close OUT;
}
