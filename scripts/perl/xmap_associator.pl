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
$program: A perl script for associating contigs in two XMAP files, resulting in filtered XMAP files and a TSV file
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-f, --first <str>   The first XMAP file containing the alignment information (REQUIRED)
	-s, --second <str>  The second XMAP file containing the alignment information (REQUIRED)
	-t1 <str>           The coordinate transformation file for the Query contigs of the first XMAP file
	-t2 <str>           The coordinate transformation file for the Query contigs of the second XMAP file
	-m, --min <float>   The minimum overlap threshold in bp for combining two alignments (default: 50000)
	-o, --output <str>  The output path (REQUIRED)
	-h, --help          Help

Example:
	$program -f alignment1.xmap -s alignment2.xmap -o common
USAGE

use vars qw($opt_f $opt_s $opt_t1 $opt_t2 $opt_m $opt_o $opt_h);
GetOptions( "f|first=s" => \$opt_f,
			"s|second=s" => \$opt_s,
			"t1=s" => \$opt_t1,
			"t2=s" => \$opt_t2,
			"m|min=f" => \$opt_m,
			"o|output=s" => \$opt_o,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

die("**ERROR: -f option must be specified\n") unless(defined $opt_f);
die("**ERROR: -s option must be specified\n") unless(defined $opt_s);
die("**ERROR: the minumum overlap threshold must be a positive\n") if((defined $opt_m) and ($opt_m <= 0));
die("**ERROR: -o option must be specified\n") unless(defined $opt_o);

my $min_overlap = ((defined $opt_m) ? $opt_m : 50000);
my ($outdir, $outpre);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = ($opt_o eq "") ? "/" : $opt_o;
	$outpre = basename($opt_f); $outpre =~ s/\.xmap$//;
	$outpre .= "_" . basename($opt_s); $outpre =~ s/\.xmap$//;
}
else{
	$outdir = dirname($opt_o);
	$outpre = basename($opt_o);
}

my ($cmd, $retCode);
$cmd = "mkdir -p $outdir";
$retCode = system($cmd);
die("**ERROR: Can not create directory \"$outdir\"\n") if($retCode != 0);

our @required_common = ("QryContigID", "QryLen");
our @required_features = ("RefContigID", "RefLen", "Confidence", "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", "Orientation");

my $xmap1 = readXMap($opt_f);
die("**ERROR: some of the required field is missing in \"$opt_f\"\n") unless(&checkFields($xmap1->{dataName}));
my $xmap2 = readXMap($opt_s);
die("**ERROR: some of the required field is missing in \"$opt_s\"\n") unless(&checkFields($xmap2->{dataName}));

my ($filename1, $filename2);
if(defined $opt_t1){
	$xmap1 = &transformCoords($xmap1, $opt_t1);
}
if(defined $opt_t2){
	$xmap2 = &transformCoords($xmap2, $opt_t2);
}

# sort the alignments by QryContigID, QryStartPos, and QryEndPos
my ($firstEntry1, $indexes1) = &sortXMap($xmap1);
my ($firstEntry2, $indexes2) = &sortXMap($xmap2);

$filename1 = &makeFileName($outdir, "xmap", "sorted", basename($opt_f));
$filename2 = &makeFileName($outdir, "xmap", "sorted", basename($opt_s));
&writeXMap($xmap1, $indexes1, $filename1);
&writeXMap($xmap2, $indexes2, $filename2);

# TSV file containing association information between two CMAPs with different channels
my $outfile = "$outdir/$outpre.tsv";
open(OUT, ">$outfile") or die("Can not open \"$outfile\" for writing\n");

# print out the header
print OUT join("\t", ("ID", "QryContigID", "QryLen", "RefContigID.x", "RefLen.x", "RefContigID.y", "RefLen.y", "Orientation.x", "Orientation.y", "Overlap", "Confidence.x", "Confidence.y", "QryStartPos.x", "QryEndPos.x", "RefStartPos.x", "RefEndPos.x", "QryStartPos.y", "QryEndPos.y", "RefStartPos.y", "RefEndPos.y", "OriContigID.x", "offset.x", "len.x", "OriContigID.y", "offset.y", "len.y")) . "\n";

# processing the alignments
my $iret;
my $id = 0;
my ($overlap, $idx1, $idx2);
my ($entry1, $entry2) = ($firstEntry1, $firstEntry2);
my ($hits1, $hits2) = ($xmap1->{hits}, $xmap2->{hits});
my (%keep1, %keep2) = ((), ());
while( (defined $entry1) and (defined $entry2) ){
	$iret = &compare($entry1, $entry2);
	if($iret < 0){
		$entry1 = $entry1->{next};
		next;
	}
	if($iret > 0){
		$entry2 = $entry2->{next};
		next;
	}
	# match
	$overlap = &getOverlap($entry1, $entry2);
	if($overlap >= $min_overlap){
		$idx1 = $entry1->{idx};
		$idx2 = $entry2->{idx};
		$keep1{$idx1} = 1;
		$keep2{$idx2} = 1;
		$id++;
		print OUT "$id\t$entry1->{contigID}\t$hits1->{QryLen}->[$idx1]\t$hits1->{RefContigID}->[$idx1]\t$hits1->{RefLen}->[$idx1]\t$hits2->{RefContigID}->[$idx2]\t$hits2->{RefLen}->[$idx2]\t$hits1->{Orientation}->[$idx1]\t$hits2->{Orientation}->[$idx2]\t$overlap\t$hits1->{Confidence}->[$idx1]\t$hits2->{Confidence}->[$idx2]\t$hits1->{QryStartPos}->[$idx1]\t$hits1->{QryEndPos}->[$idx1]\t$hits1->{RefStartPos}->[$idx1]\t$hits1->{RefEndPos}->[$idx1]\t$hits2->{QryStartPos}->[$idx2]\t$hits2->{QryEndPos}->[$idx2]\t$hits2->{RefStartPos}->[$idx2]\t$hits2->{RefEndPos}->[$idx2]\t$hits1->{OriContigID}->[$idx1]\t$hits1->{map}->[$idx1]->{offset}\t$hits1->{map}->[$idx1]->{len}\t$hits2->{OriContigID}->[$idx2]\t$hits2->{map}->[$idx2]->{offset}\t$hits2->{map}->[$idx2]->{len}\n";
	}
	($entry1, $entry2) = &getNextPair($entry1, $entry2);
}

close OUT;

$filename1 = &makeFileName($outdir, "xmap", "filtered", basename($filename1));
$filename2 = &makeFileName($outdir, "xmap", "filtered", basename($filename2));

my @filtered1 = grep {$keep1{$_}} @{$indexes1};
my @filtered2 = grep {$keep2{$_}} @{$indexes2};
&writeXMap($xmap1, \@filtered1, $filename1);
&writeXMap($xmap2, \@filtered2, $filename2);

exit 0;

sub transformCoords
{
	my ($xmap, $coord_file) = @_;
	open(IN, "<$coord_file") || die("ERROR: cannot open $coord_file for reading\n");
	my $line;
	my @tokens;
	my %map = ();
	my %length = ();
	while($line = <IN>){
		chomp($line);
		next if($line !~ /^\d/); # skip header
		@tokens = split(/\t/, $line);
		next if(scalar(@tokens) < 6);
		$map{$tokens[3]} = {'id' => $tokens[0], 'offset' => $tokens[1], 'len' => $tokens[5]};
		if(!defined $length{$tokens[0]} or $tokens[2] + 1 > $length{$tokens[0]}){
			$length{$tokens[0]} = $tokens[2] + 1;
		}
	}
	close IN;
	my $QryContigIDs = $xmap->{hits}->{QryContigID};
	my $QryStartPoses = $xmap->{hits}->{QryStartPos};
	my $QryEndPoses = $xmap->{hits}->{QryEndPos};
	my $QryLens = $xmap->{hits}->{QryLen};
	$xmap->{hits}->{OriContigID} = [];
	$xmap->{hits}->{map} = [];
	my $nHits = $xmap->{totalHits};
	my ($qryId, $startPos, $endPos);
	for(my $i=0; $i<$nHits; $i++){
		$qryId = $QryContigIDs->[$i];
		$xmap->{hits}->{OriContigID}->[$i] = $qryId;
		my $item = $map{$qryId};
		if(!defined $item){
			warn "TransformCoords: old ID for $qryId is undefined\n";
			next;
		}
		$xmap->{hits}->{map}->[$i] = $item;
		($startPos, $endPos) = ($QryStartPoses->[$i], $QryEndPoses->[$i]);
		if( ($startPos > $item->{len}) || ($endPos > $item->{len}) ){
			warn "TransformCoords: Out of bounds detected for QryID $qryId in the input XMAP\n";
			next;
		}
		my $oldId = $item->{id};
		($QryContigIDs->[$i], $QryStartPoses->[$i], $QryEndPoses->[$i], $QryLens->[$i]) = ($oldId, $startPos + $item->{offset}, $endPos + $item->{offset}, $length{$oldId});
	}

	return $xmap;
}

sub checkFields
{
	my ($data_name) = @_;
	my %hasField = ();
	for(my $i=0; $i<scalar(@{$data_name}); $i++){
		$hasField{$data_name->[$i]} = 1;
	}
	our (@required_common, @required_features);
	foreach(@required_common){
		if(!defined $hasField{$_}){
			return 0;
		}
	}
	foreach(@required_features){
		if(!defined $hasField{$_}){
			return 0;
		}
	}
	return 1;
}

sub sortXMap
{
	my ($xmap) = @_;
	my $QryContigIDs = $xmap->{hits}->{QryContigID};
	my $QryStartPoses = $xmap->{hits}->{QryStartPos};
	my $QryEndPoses = $xmap->{hits}->{QryEndPos};
	my @entries = ();
	my ($pos1, $pos2);
	my $i;
	my $nHits = $xmap->{totalHits};
	for($i=0; $i<$nHits; $i++){
		($pos1, $pos2) = ($QryStartPoses->[$i], $QryEndPoses->[$i]);
		if($pos1 < $pos2){
			push @entries, {
				'contigID' => $QryContigIDs->[$i],
				'startPos' => $pos1,
				'endPos' => $pos2
			};
		} else{
			push @entries, {
				'contigID' => $QryContigIDs->[$i],
				'startPos' => $pos2,
				'endPos' => $pos1
			};
		}
	}
	my @indexes = sort {
		$entries[$a]->{contigID} <=> $entries[$b]->{contigID} or
		$entries[$a]->{startPos} <=> $entries[$b]->{startPos} or
		$entries[$b]->{endPos} <=> $entries[$a]->{endPos}
	} 0..$#entries;

	my ($firstEntry, $entry);
	for($i=$nHits-1; $i>=0; $i--){
		$entry = \%{$entries[$indexes[$i]]};
		$entry->{idx} = $indexes[$i];
		$entry->{next} = $firstEntry;
		$firstEntry = $entry;
	}
	return  ($firstEntry, \@indexes);
}

sub compare
{
	my ($entry1, $entry2) = @_;
	my ($id1, $id2) = ($entry1->{contigID}, $entry2->{contigID});
	if($id1 < $id2){
		return -1;
	}
	if($id1 > $id2){
		return 1;
	}
	# $id1 == $id2
	return 0;
}

sub getOverlap
{
	my ($entry1, $entry2) = @_;
	if($entry1->{contigID} ne $entry2->{contigID}){
		return 0;
	}
	my ($start1, $start2) = ($entry1->{startPos}, $entry2->{startPos});
	my ($end1, $end2) = ($entry1->{endPos}, $entry2->{endPos});

	my $overlap = min($end1, $end2) - max($start1, $start2) + 1;

	return (($overlap > 0) ? $overlap : 0);
}

sub getNextPair
{
	my ($entry1, $entry2) = @_;
# (defined $entry1) and (defined $entry2)
	my ($entry3, $entry4) = ($entry1->{next}, $entry2->{next});
	if(!defined $entry3 or !defined $entry4){
		return ($entry3, $entry4);
	}
# (defined $entry3) and (defined $entry4)
	my $overlap1 = &getOverlap($entry3, $entry2);
	my $overlap2 = &getOverlap($entry1, $entry4);
	return  ($overlap1 == 0 and $overlap2 == 0) ? ($entry3, $entry4) :
			($overlap1 < $overlap2) ? ($entry1, $entry4) :
			($entry3, $entry2);
}

sub makeFileName
{
	my ($outdir, $ext, $suffix, $filename) = @_;
	$filename =~ s/\.$ext$//;
	my $newname = "$outdir/$filename-$suffix\.$ext";
	return $newname;
}

sub writeXMap
{
	my ($xmap, $indexes, $filename) = @_;
	open(XMF, ">$filename") or die("Can not open \"$filename\" for writing\n");

	my $i;
	my $header = $xmap->{headers};
	for($i=0; $i<scalar(@{$header}); $i++){
		print XMF "$header->[$i]\n";
	}
	my $j;
	my @data_name = @{ $xmap->{dataName} };
	my $nFields = scalar(@data_name);
	for($i=0; $i<scalar(@{$indexes}); $i++){
		for($j=0; $j<$nFields; $j++){
			if($j == 0){
				print XMF int($i+1);
			} else {
				print XMF $xmap->{hits}->{$data_name[$j]}->[$indexes->[$i]];
			}
			if($j < $nFields - 1){
				print XMF "\t";
			} else {
				print XMF "\n";
			}
		}
	}

	close XMF;
}
