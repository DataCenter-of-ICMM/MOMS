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
use Math::Round;

BEGIN{
	my $progpath = abs_path(dirname($0));
	unshift @INC, $progpath;
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

use BNG::Utility;

my %enzymes = (
	"BspQI" => "GCTCTTC",
	"BssSI" => "CACGAG",
	"BbvCI" => "CCTCAGC",
	"BsmI"  => "GAATGC",
	"BsrDI" => "GCAATG",
	"bseCI" => "ATCGAT",
	"BamHI" => "GGATCC",
	"DLE1" => "CTTAAG"
	
	# You can add more enzymes here ...
	
);

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for converting alignment format from Bionano XMAP to OpGen XML
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-i, --input <str>   The input XMAP file containing the alignment information (Required)
	-r, --ref <str>     The reference CMAP file (Required)
	-q, --qry <str>     The query CMAP file (Required)
	-s, --switch        Switch the query and reference of the input XMAP file (default: no)
	-o, --output <str>  The output OpGen XML file (Required)
	-h, --help          Help

Example:
	$program -i bspqi.xmap -r EXP_REFEINEFINAL_BSPQI.cmap -q contigs_BSPQI.cmap -o bspqi.xml
USAGE

use vars qw($opt_i $opt_r $opt_q $opt_s $opt_o $opt_h);
GetOptions( "i|input=s" => \$opt_i,
			"r|ref=s" => \$opt_r,
			"q|qry=s" => \$opt_q,
			"s|switch" => \$opt_s,
			"o|output=s" => \$opt_o,
			"h|help");

die($usage) if($opt_h);
die("**ERROR: -i option must be specified\n") unless(defined $opt_i);
die("**ERROR: -r option must be specified\n") unless(defined $opt_r);
die("**ERROR: -q option must be specified\n") unless(defined $opt_q);
die("**ERROR: -o option must be specified\n") unless(defined $opt_o);

my ($xmap_file, $ref_file, $qry_file, $out_xml) = ($opt_i, $opt_r, $opt_q, $opt_o);
my $xmap = readXMap($xmap_file);
my ($cmap_ref) = readCMap($ref_file);
my ($cmap_qry) = readCMap($qry_file);

die("**ERROR: the Nickase Recognition Sites are incompatible for the input CMAP files\n") if( &compareChannels($cmap_ref->{channels}, $cmap_qry->{channels}) != 0 );

my $experimentInfo = &createExperiment(abs_path(dirname($0)) . "/$program", basename($xmap_file), basename($out_xml), $xmap->{totalHits});
my $enzyme_str = &getEnzymeStr($cmap_ref->{channels}, \%enzymes);
my $ref_maps = &createMaps($cmap_ref, "consensus", $enzyme_str);
my $qry_maps = &createMaps($cmap_qry, "opmap", $enzyme_str);

my $alignments = &createAlignments($xmap, $cmap_ref, $cmap_qry, $opt_s);
&writeOpGenXml($out_xml, $experimentInfo, $ref_maps, $qry_maps, $alignments);

exit 0;

sub compareChannels
{
	my ($channels1, $channels2) = @_;
	if(scalar(values %{$channels1}) ne scalar(values %{$channels2})){
		return (scalar(values %{$channels1}) < scalar(values %{$channels2})) ? -1 : 1;
	}
	my @motifs1 = sort {$a cmp $b } values %{$channels1};
	my @motifs2 = sort {$a cmp $b } values %{$channels2};
	for(my $i=0; $i<scalar(@motifs1); $i++){
		if($motifs1[$i] ne $motifs2[$i]){
			return ($motifs1[$i] < $motifs2[$i]) ? -1 : 1;
		}
	}
	return 0;
}

sub getEnzymeStr
{
	my ($channels, $enzymes) = @_;
	my %usedMotif = map{uc($_)=>1} values %{$channels};
	my @enzyme_names = grep {$usedMotif{$enzymes{$_}}} (keys %$enzymes);
	return join(",", @enzyme_names);
}

sub createExperiment
{
	my ($creator, $xmap_file, $xml_file, $nhits) = @_;
	my %expInfo = ();
	$expInfo{name} = "experiment";
	$expInfo{value} = [];
	push(@{$expInfo{value}}, {"name" => "creator", "value" => $creator});
	push(@{$expInfo{value}}, {"name" => "start_time", "value" => time()});
	push(@{$expInfo{value}}, {"name" => "soma_version", "value" => 0.8});
	push(@{$expInfo{value}}, {"name" => "uuid", "value" => "00000000-0000-0000-0000-000000000000"});
	push(@{$expInfo{value}}, {"name" => "genomes_mapset_file", "value" => $xmap_file});
	push(@{$expInfo{value}}, {"name" => "genomes_parameter_file", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "mapset_range", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "genspect_file", "value" => $xml_file});
	push(@{$expInfo{value}}, {"name" => "number_alignments", "value" => $nhits});
	push(@{$expInfo{value}}, {"name" => "minimum_score", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "minimum_aligned_chunks", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "mum_flanking_chunks", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "p_value_cutoff_probability", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "p_value_lambda", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "p_value_K", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "number_standard_error", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "size_mismatch_fraction", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "cut_penalty", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "expected_cut_mismatch", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "maximum_mismatched_cuts", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "gap_maximum_distance", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "gap_open_penalty", "value" => ""});
	push(@{$expInfo{value}}, {"name" => "gap_extend_penalty", "value" => ""});

	return \%expInfo;
}

sub createMaps
{
	my ($cmap, $type, $enzyme_str) = @_;
	my @maps = ();
	my @cmapIds = sort {$a <=> $b} keys %{ $cmap->{contigs} };
	for (my $i=0; $i < scalar(@cmapIds); $i++) {
		my $ctg = $cmap->{contigs}->{$cmapIds[$i]};
		next if($ctg->{NumSites} <= 2);
		my $mapInfo = createMap($ctg, $cmapIds[$i], $type, $enzyme_str);
		push(@maps, $mapInfo);
	}
	return \@maps;
}

sub createMap
{
	my ($contig, $id, $type, $enzyme_str) = @_;
	my %mapInfo = ();
	$mapInfo{name} = "restriction_map";
	$mapInfo{value} = [];
	
	my $numsites = $contig->{NumSites};
	my $positions = $contig->{Position};
	my @map_block = ();
	for (my $j=1; $j < $numsites; $j++) {
		my $fragsize = sprintf("%.3f", round($positions->[$j] - $positions->[$j-1]) / 1000.0);
		push(@map_block, $fragsize);
	}
	push(@{$mapInfo{value}}, {"name" => "name", "value" => (($type eq "consensus") ? "OM_$id": "$id")});
	push(@{$mapInfo{value}}, {"name" => "type", "value" => $type});
	push(@{$mapInfo{value}}, {"name" => "circular", "value" => "false"});
	push(@{$mapInfo{value}}, {"name" => "enzymes", "value" => $enzyme_str});
	push(@{$mapInfo{value}}, {"name" => "num_frags", "value" => scalar(@map_block)});
	push(@{$mapInfo{value}}, {"name" => "map_block", "value" => join(" ", @map_block)});

	return \%mapInfo;
}

sub createAlignments
{
	my ($xmap, $cmap_ref, $cmap_qry, $switch) = @_;
	my @Alignments = ();
	my $QryIDs = $xmap->{hits}->{QryContigID};
	my $RefIDs = $xmap->{hits}->{RefContigID};
	my $Orientations = $xmap->{hits}->{Orientation};
	my $Confidences = $xmap->{hits}->{Confidence};
	my $AlignmentStrings = $xmap->{hits}->{Alignment};
	my ($refId, $qryId, $orientation, $pairs);
	for (my $i=0; $i < $xmap->{totalHits}; $i++) {
		if(defined $switch){
			($refId, $qryId) = ($QryIDs->[$i], $RefIDs->[$i]);
		}
		else{
			($refId, $qryId) = ($RefIDs->[$i], $QryIDs->[$i]);
		}
		if(!defined $cmap_ref->{contigs}->{$refId}){
			warn "$refId is not defined in reference CMAP file\n";
			next;
		}
		if(!defined $cmap_qry->{contigs}->{$qryId}){
			warn "$qryId is  not defined in query CMAP file\n";
			next;
		}
		if($cmap_qry->{contigs}->{$qryId}->{NumSites} <= 1){
			warn "Can not convert alignment for qryId $qryId\n";
			next;
		}
		$orientation = $Orientations->[$i];
		$pairs = getPairs($AlignmentStrings->[$i], $switch, $orientation);
		my $alignment = createAlignment($refId, $qryId, $orientation, $Confidences->[$i], $pairs, $cmap_qry->{contigs}->{$qryId}->{NumSites} - 1);
		push(@Alignments, $alignment);
	}

	return \@Alignments;
}

sub getPairs
{
	my ($alignment, $switch, $orientation) = @_;
	$alignment =~ s/^\(//;
	$alignment =~ s/\)$//;
	my @pairs = ();
	foreach (split(/\)\(/, $alignment)){
		push(@pairs, [split(/,/, $_)]);
	}
	if(defined $switch){
		foreach my $pair (@pairs){
			my $tmp = $pair->[0];
			$pair->[0] = $pair->[1];
			$pair->[1] = $tmp;
		}
		if($orientation eq '+'){
			@pairs = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @pairs;
		}
		else{
			@pairs = sort { $a->[0] <=> $b->[0] || $b->[1] <=> $a->[1] } @pairs;
		}
	}
	return \@pairs;
}

sub createAlignment
{
	my ($refId, $qryId, $orientation, $confidence, $pairs, $nFragments) = @_;
	my %alignInfo = ();
	my $alignedTriples = createTriples($pairs, $orientation, $nFragments);
	$orientation = ($orientation eq '+') ? 'N' : 'R';
	$alignInfo{name} = "map_alignment";
	$alignInfo{value} = [];
	push($alignInfo{value}, {"name" => "reference_map", "value" => {"name" => "name", "value" => "OM_$refId"}});
	push($alignInfo{value}, {"name" => "aligned_map",
								 "value" => [{"name" => "name", "value" => $qryId},
								 			 {"name" => "orientation", "value" => $orientation}]
							});
	push($alignInfo{value}, {"name" => "soma_score", "value" => $confidence});
	push($alignInfo{value}, {"name" => "count", "value" => scalar(@$alignedTriples)});
	push($alignInfo{value}, {"name" => "", "value" => $alignedTriples});

	return \%alignInfo;
}

sub createTriples
{
	my ($pairs, $orientation, $nFragments) = @_;
	my @triples = ();
	if($orientation eq '+'){
		my ($dr, $dq);
		my $last_i = 0;
		for(my $i=1; $i<scalar(@$pairs); $i++){
			$dr = $pairs->[$i]->[0] - $pairs->[$last_i]->[0];
			$dq = $pairs->[$i]->[1] - $pairs->[$last_i]->[1];
			if($dq >= 1){
				while($dq >= 1){
					push(@triples, [$pairs->[$i]->[1] - $dq - 1, $pairs->[$last_i]->[0] - ($dr == 0 ? 1 : 0) - 1, $pairs->[$i]->[0] - 1 - 1]);
					$dq--;
				}
				$last_i = $i;
			}
		}
	}
	else{ # '-'
		my ($dr, $dq);
		my $last_i = 0;
		for(my $i=1; $i<scalar(@$pairs); $i++){
			$dr = $pairs->[$i]->[0] - $pairs->[$last_i]->[0];
			$dq = $pairs->[$last_i]->[1] - $pairs->[$i]->[1];
			if($dq >= 1){
				for(my $tdq = 1; $tdq<=$dq; $tdq++){
					push(@triples, [$nFragments - 1 - ($pairs->[$last_i]->[1] - $tdq - 1), $pairs->[$last_i]->[0]  - ($dr == 0 ? 1 : 0) - 1, $pairs->[$i]->[0] - 1 - 1]);
				}
				$last_i = $i;
			}
		}
	}
	my @tripleElements = ();
	foreach(@triples){
		push(@tripleElements, {"name" => "f", "value" => "<i>$_->[0]</i><l>$_->[1]</l><r>$_->[2]</r>"});
	}
	return \@tripleElements;
}

sub writeOpGenXml
{
	my ($out_xml, $experimentInfo, $ref_maps, $qry_maps, $alignments) = @_;
	open(OUT, ">$out_xml") or die("Can not open $out_xml for writing\n");
	print OUT "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n\n";
	print OUT "<aligned_maps_document>\n\n";
	print OUT "<version>1.0</version>\n\n";
	&writeElement(\*OUT, $experimentInfo);
	&writeElement(\*OUT, $ref_maps);
	&writeElement(\*OUT, $qry_maps);
	&writeElement(\*OUT, $alignments);
	print OUT "</aligned_maps_document>\n";
	close OUT;
}

sub writeElement
{
	my ($fh, $element, $indent) = @_;
	$indent = "" unless(defined $indent);
	if(ref($element) eq "ARRAY"){
		foreach(@{$element}){
			writeElement($fh, $_, $indent);
		}
		return;
	}
	if(ref($element) eq "HASH"){
		if(ref($element->{value}) eq ""){
			print $fh "$indent<$element->{name}>$element->{value}</$element->{name}>\n";
		}
		else{
			if($element->{name} ne ""){
				print $fh "$indent<$element->{name}>\n";
			}
			writeElement($fh, $element->{value}, "$indent  ");
			if($element->{name} ne ""){
				print $fh "$indent</$element->{name}>\n";
			}
		}
		if($indent eq ""){
			print $fh "\n";
		}
	}
}
