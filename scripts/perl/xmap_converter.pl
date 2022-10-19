#!/usr/bin/env perl
# Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences
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

BEGIN{
	my $progpath = abs_path(dirname($0));
	unshift @INC, $progpath;
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

use BNG::Utility;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for transform coordinates in a XMAP file and output a new XMAP file
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-i, --input <str>   The input XMAP file (Required)
	-o, --output <str>  The output XMAP file (Required)
	-t, --trans <str>   The coordinate transformation file for the Query contigs of the XMAP file (Required)
	-c, --cmap <str>    The CMAP file of the Query before coordinate transformation (Required)
	-s, --switch        To switch the pairs in alignment strings (default: no)
	-h, --help          Help

Exmaple:
	$program -i input.xmap -t coord_translation.txt -c ori_qry.cmap -o output.xmap
USAGE

use vars qw($opt_i $opt_o $opt_t $opt_c $opt_s $opt_h);
GetOptions( "i|input=s" => \$opt_i,
			"o|output=s" => \$opt_o,
			"t|trans=s" => \$opt_t,
			"c|cmap=s" => \$opt_c,
			"s|switch" => \$opt_s,
			"h|help" => \$opt_h);

die($usage) if($opt_h);
die("**ERROR: -i option must be specified\n") unless(defined $opt_i);
die("**ERROR: -o option must be specified\n") unless(defined $opt_o);
die("**ERROR: -t option must be specified\n") unless(defined $opt_t);
die("**ERROR: -c option must be specified\n") unless(defined $opt_c);

my ($xmap_file, $coord_translation, $cmap_file, $new_xmap_file) = ($opt_i, $opt_t, $opt_c, $opt_o);
my $xmap = readXMap($xmap_file);
die("**ERROR: some of the required field is missing in \"$xmap_file\"\n") unless(&checkXmapFields($xmap->{dataName}));
my ($cmap) = readCMap($cmap_file);
$xmap = &transformCoords($xmap, $coord_translation, $cmap, defined($opt_s));
writeXMapFile($xmap, $new_xmap_file);

exit 0;

sub transformCoords
{
	my ($xmap, $coord_file, $cmap, $switch) = @_;
	$switch = 0 unless(defined $switch);
	open(IN, "<$coord_file") || die("ERROR: cannot open $coord_file for reading\n");
	my $line;
	my @tokens;
	my %map = ();
	my %length = ();
	my $leading_nsites;
	my @data_name = @{ $cmap->{dataName} };
	while($line = <IN>){
		chomp($line);
		next if($line !~ /^\d/); # skip header
		@tokens = split(/\t/, $line);
		next if(scalar(@tokens) < 6);
		my $ctg = $cmap->{contigs}->{$tokens[0]};
		my $numsites = $ctg->{$data_name[2]};
		for ($leading_nsites=0; $leading_nsites <= $numsites; $leading_nsites++) {
			my $position = $ctg->{$data_name[5]}->[$leading_nsites];
			if($tokens[1] <= $position){
				last;
			}
		}
		$map{$tokens[3]} = {'id' => $tokens[0], 'offset' => $tokens[1], 'len' => ($tokens[5] - $tokens[4] + 1), 'leading_sites' => $leading_nsites};
		if(!defined $length{$tokens[0]} or ($tokens[2] + 1) > $length{$tokens[0]}){
			$length{$tokens[0]} = ($tokens[2] + 1);
		}
	}
	close IN;
	my $QryContigIDs = $xmap->{hits}->{QryContigID};
	my $QryStartPoses = $xmap->{hits}->{QryStartPos};
	my $QryEndPoses = $xmap->{hits}->{QryEndPos};
	my $QryLens = $xmap->{hits}->{QryLen};
	my $nHits = $xmap->{totalHits};
	my $Alignments = $xmap->{hits}->{Alignment};
	my ($qryId, $startPos, $endPos, $alignment);
	for(my $i=0; $i<$nHits; $i++){
		$qryId = $QryContigIDs->[$i];
		my $item = $map{$qryId};
		if(!defined $item){
			warn "TransformCoords: old ID for $qryId is undefined\n";
			next;
		}
		($startPos, $endPos) = ($QryStartPoses->[$i], $QryEndPoses->[$i]);
		if( ($startPos > $item->{len}) || ($endPos > $item->{len}) ){
			warn "TransformCoords: Out of bounds detected for QryID $qryId in the input XMAP\n";
			next;
		}
		my $oldId = $item->{id};
		($QryContigIDs->[$i], $QryStartPoses->[$i], $QryEndPoses->[$i], $QryLens->[$i]) = ($oldId, $startPos + $item->{offset}, $endPos + $item->{offset}, $length{$oldId});
		$alignment = &modifyAlignment($Alignments->[$i], $switch, $item->{leading_sites});
		$Alignments->[$i] = $alignment;
	}

	return $xmap;
}

sub modifyAlignment
{
	my ($alignment, $switch, $leading_sites) = @_;
	$alignment =~ s/^\(//;
	$alignment =~ s/\)$//;
	my @pairs = ();
	foreach (split(/\)\(/, $alignment)){
		push(@pairs, [split(/,/, $_)]);
	}
	foreach my $pair (@pairs){
		if($switch){
			my $tmp = $pair->[0];
			$pair->[0] = $pair->[1];
			$pair->[1] = $tmp;
		}
		$pair->[0] += $leading_sites;
	}
	@pairs = sort { $a->[1] <=> $b->[1] } @pairs;
	my @aligns = ();
	foreach my $pair (@pairs){
		push(@aligns, join(",", @$pair));
	}
	$alignment = "(" . join(")(", @aligns) . ")";
	return $alignment;
}

sub checkXmapFields
{
	my ($data_name) = @_;
	my %hasField = ();
	for(my $i=0; $i<scalar(@{$data_name}); $i++){
		$hasField{$data_name->[$i]} = 1;
	}
	my @required_common = ("QryContigID", "QryLen");
	foreach(@required_common){
		if(!defined $hasField{$_}){
			return 0;
		}
	}
	my @required_features = ("RefContigID", "RefLen", "Confidence", "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", "Orientation");
	foreach(@required_features){
		if(!defined $hasField{$_}){
			return 0;
		}
	}
	return 1;
}
