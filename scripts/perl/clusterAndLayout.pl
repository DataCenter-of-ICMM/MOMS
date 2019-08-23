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
use List::MoreUtils qw(uniq);
use Math::Round;
use Cwd 'abs_path';
use Switch;
use Storable qw(dclone);
use feature qw(state);

BEGIN{
	my $progpath = abs_path(dirname($0));
	unshift(@INC, $progpath);
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

use BNG::Utility;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script which clusters the nodes and make layouts of these nodes within clusters
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-r, --reference <str> The reference CMAP file(s) (REQUIRED)
	-e, --edges <str>   The TSV file containing the edge information of nodes (REQUIRED)
	-g, --glues <str>   The TSV file containing the glue information of edges (REQRUIED)
	-o, --output <str>  The path of the output table (default: STDOUT)
	-h, --help          Help

Example:
	$program -r BSPQI_r.cmap -r BSSSI_r.cmap -e edges.tsv -g glues.tsv -o layout
USAGE

use vars qw(@opt_r $opt_e $opt_g $opt_o $opt_h);
GetOptions( "r|reference=s" => \@opt_r,
			"e|edges=s" => \$opt_e,
			"g|glues=s" => \$opt_g,
			"o|output=s" => \$opt_o,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

die("**ERROR: at least two -r options must be specified\n") unless(scalar(@opt_r)>=2);
die("**ERROR: -e option must be specified\n") unless(defined $opt_e);
die("**ERROR: -g option must be specified\n") unless(defined $opt_g);

die("Only one of '-e' and '-g' can be specified as STDIN\n") if( ($opt_e eq "-") and ($opt_g eq "-") );

my $outfile;
if(defined $opt_o){
	my ($outdir, $outpre);
	if($opt_o =~ /\/$/ or $opt_o !~ /\//){
		$opt_o =~ s/\/$//;
		$outdir = ($opt_o eq "") ? "/" : $opt_o;
		$outpre = basename($opt_e); $outpre =~ s/\.tsv$//;
	}
	else{
		$outdir = dirname($opt_o);
		$outpre = basename($opt_o);
	}

	my ($cmd, $retCode);
	$cmd = "mkdir -p $outdir";
	$retCode = system($cmd);
	die("**ERROR: Can not create directory \"$outdir\"\n") if($retCode != 0);
	$outfile = "$outdir/$outpre.tsv";
}

my @cmaps = ();
for(my $i=0; $i<scalar(@opt_r); $i++){
	($cmaps[$i+1]) = readCMap($opt_r[$i]);
}
my $edgeTable = &readTable($opt_e, ["id_str", "nglues", "overlap", "intercept", "slope", "type1", "type2", "id1", "id2", "length1", "length2", "glues"]);
my $glueTable = &readTable($opt_g, ["ID", "QryContigID", "QryLen", "type1", "RefContigID.x", "RefLen.x", "type2", "RefContigID.y", "RefLen.y", "Orientation.x", "Orientation.y", "Overlap", "Confidence.x", "Confidence.y", "QryStartPos.x", "QryEndPos.x", "RefStartPos.x", "RefEndPos.x", "QryStartPos.y", "QryEndPos.y", "RefStartPos.y", "RefEndPos.y", "OriContigID.x", "offset.x", "len.x", "OriContigID.y", "offset.y", "len.y"]);

my $graph = &buildGraph($edgeTable, $glueTable, \@cmaps);
my $cluster_pool = &traverseGraph($graph);

my ($paths, $headers) = &makeLayouts($graph->{nodes}, $cluster_pool);
&outputLayouts($paths, $headers, $outfile);

exit 0;

sub readTable
{
	my ($infile, $mustFields) = @_;
	my @table = ();
	my $row;
	my $in;
	if($infile eq "-"){
		$in = \*STDIN;
	} else {
		open(IN, "$infile") or die("Can not open \"$infile\" for reading\n");
		$in = \*IN;
	}
	my @fields = ();
	my $nfields = 0;
	my @indices = ();
	my $line;
	my @columns = ();
	my ($i, $idx);
	while($line = <$in>){
		chomp($line);
		if($nfields == 0){
			@fields = split(/\t/, $line);
			$nfields = scalar(@fields);
			@indices = &getIndices(\@fields, $mustFields);
			if(scalar(@indices) != scalar(@$mustFields)){
				last;
			}
			next;
		}
		@columns = split(/\t/, $line);
		next if(scalar(@columns) != $nfields);
		$row = undef;
		for($i=0; $i<scalar(@indices); $i++){
			$idx = $indices[$i];
			$row->{$fields[$idx]} = $columns[$idx];
		}
		push(@table, $row);
	}
	close $in;

	die("**ERROR: some of the required fields is missing in \"$infile\"\n") unless(scalar(@indices) == scalar(@$mustFields));

	return \@table;
}

sub buildGraph
{
	my ($edgeTable, $glueTable, $cmaps) = @_;
	my @nodes = ();
	my @edges = ();
	my $row;
	my ($id1, $id2, $id_str);
	my ($type1, $type2);
	my ($ori_id1, $ori_id2);
	my ($length1, $length2);
	my $maxtype = 0;
	my @glueIDs;
	my ($glue, $slope, $intercept);
	my @glues;
	my $edge;
	for(my $i=0; $i<scalar(@{$edgeTable}); $i++){
		my $row = \%{$edgeTable->[$i]};
		($id_str, $ori_id1, $ori_id2) = ($row->{"id_str"}, $row->{"id1"}, $row->{"id2"});
		($id1, $id2) = split(/_/, $id_str);
		($length1, $length2) = ($row->{"length1"}, $row->{"length2"});
		if(!defined $nodes[$id1]){
			$type1 = $row->{"type1"};
			if($type1 > $maxtype){
				$maxtype = $type1;
			}
			$nodes[$id1] = &defineNode($id1, $type1, $ori_id1, $length1);
		}
		if(!defined $nodes[$id2]){
			$type2 = $row->{"type2"};
			if($type2 > $maxtype){
				$maxtype = $type2;
			}
			$nodes[$id2] = &defineNode($id2, $type2, $ori_id2, $length2);
		}
		($slope, $intercept) = ($row->{slope}, $row->{intercept});
		@glueIDs = split(/,/, $row->{glues});
		@glues = ();
		foreach my $id (@glueIDs){
			$glue = &defineGlue($glueTable->[$id - 1]);
			next if(&isOutlierGlue($glue, $slope, $intercept));
			push(@glues, $glue);
		}
		$edge = &defineEdge($id1, $id2, $row->{slope}, $row->{intercept}, $row->{overlap}, [sort{ $a->{first}->[0] <=> $b->{first}->[0] } (@glues)]);
		if(scalar(@{$edge->{glues}}) == 0){
			warn ("Edge ($id1, $id2): no valid glue found, neglected\n");
			next;
		}
		push(@edges, $edge);
	}
	my %graph = ();
	$graph{nodes} = \@nodes;
	$graph{edges} = \@edges;
	$graph{cmaps} = $cmaps;
	$graph{nchannels} = $maxtype;

	return \%graph;
}

sub defineGlue
{
	my ($row) = @_;
	my %glue = ();
	my ($dir1, $dir2) = ($row->{"Orientation.x"}, $row->{"Orientation.y"});
	my ($qryStart1, $qryEnd1, $qryStart2, $qryEnd2) = ($row->{"QryStartPos.x"}, $row->{"QryEndPos.x"}, $row->{"QryStartPos.y"}, $row->{"QryEndPos.y"});
	my ($refStart1, $refEnd1, $refStart2, $refEnd2) = ($row->{"RefStartPos.x"}, $row->{"RefEndPos.x"}, $row->{"RefStartPos.y"}, $row->{"RefEndPos.y"});
	$glue{len} = $row->{Overlap};
	$glue{id} = $row->{ID};
	if($dir1 eq "+"){
		if($dir2 eq "+"){ # ++
			if($qryStart1 < $qryStart2){
				$glue{first} = [$refStart1 + ($qryStart2 - $qryStart1), $refStart2];
			}
			else{
				$glue{first} = [$refStart1, $refStart2 + ($qryStart1 - $qryStart2)];
			}
			if($qryEnd1 < $qryEnd2){
				$glue{second} = [$refEnd1, $refEnd2 - ($qryEnd2 - $qryEnd1)];
			}
			else{
				$glue{second} = [$refEnd1 - ($qryEnd1 - $qryEnd2), $refEnd2];
			}
		}
		else{ # +-
			if($qryStart1 < $qryEnd2){
				$glue{first} = [$refStart1 + ($qryEnd2 - $qryStart1), $refEnd2];
			}
			else{
				$glue{first} = [$refStart1, $refEnd2 - ($qryStart1 - $qryEnd2)];
			}
			if($qryEnd1 < $qryStart2){
				$glue{second} = [$refEnd1, $refStart2 + ($qryStart2 - $qryEnd1)];
			}
			else{
				$glue{second} = [$refEnd1 - ($qryEnd1 - $qryStart2), $refStart2];
			}
		}
	}
	else{
		if($dir2 eq "+"){ # -+
			if($qryEnd1 < $qryStart2){
				$glue{second} = [$refEnd1 - ($qryStart2 - $qryEnd1), $refStart2];
			}
			else{
				$glue{second} = [$refEnd1, $refStart2 + ($qryEnd1 - $qryStart2)];
			}
			if($qryStart1 < $qryEnd2){
				$glue{first} = [$refStart1, $refEnd2 - ($qryEnd2 - $qryStart1)];
			}
			else{
				$glue{first} = [$refStart1 + ($qryStart1 - $qryEnd2), $refEnd2];
			}
		}
		else{ # --
			if($qryEnd1 < $qryEnd2){
				$glue{second} = [$refEnd1 - ($qryEnd2 - $qryEnd1), $refEnd2];
			}
			else{
				$glue{second} = [$refEnd1, $refEnd2 - ($qryEnd1 - $qryEnd2)];
			}
			if($qryStart1 < $qryStart2){
				$glue{first} = [$refStart1, $refStart2 + ($qryStart2 - $qryStart1)];
			}
			else{
				$glue{first} = [$refStart1 + ($qryStart1 - $qryStart2), $refStart2];
			}
		}
	}

	return \%glue;
}

sub isOutlierGlue
{
	my ($glue, $slope, $intercept) = @_;
	my $diff;
	$diff = abs($glue->{first}->[0] * $slope + $intercept - $glue->{first}->[1]);
	if( ($diff > 50000) && ($diff / $glue->{len} > 0.1) ){
		return 1;
	}
	$diff = abs($glue->{second}->[0] * $slope + $intercept - $glue->{second}->[1]);
	return ( ($diff > 50000) && ($diff / $glue->{len} > 0.1) ) ? 1 : 0;
}

sub defineNode
{
	my ($id, $type, $ori_id, $length) = @_;
	my %node = ();
	$node{id} = $id;
	$node{type} = $type;
	$node{ori_id} = $ori_id;
	$node{length} = $length;

	return \%node;
}

sub defineEdge
{
	my ($id1, $id2, $slope, $intercept, $weight, $glues) = @_;
	my %edge = ();
	$edge{id1} = $id1;
	$edge{id2} = $id2;
	$edge{slope} = $slope;
	$edge{intercept} = $intercept;
	$edge{weight} = $weight;
	$edge{glues} = $glues;

	return \%edge;
}

sub traverseGraph
{
	my ($graph) = @_;
	my $nodes = $graph->{nodes};
	my $edges = $graph->{edges};
	my @sorted = sort { $a->{weight} <=> $b->{weight} } @{$edges};

	my $cluster_pool = {nchannels=>$graph->{nchannels}, graph=>$graph, clusters=>[], frees=>[], homo=>{}};
	my $cluster;
	my %clsid = ();

	my $edge;
	my ($id1, $id2);
	my ($type1, $type2);
	my ($cls1, $cls2);

	for(my $i=$#sorted; $i>=0; $i--){
		$edge = $sorted[$i];
		handleHomoNode($edge, $cluster_pool->{homo});
		($id1, $id2) = ($edge->{id1}, $edge->{id2});
		($type1, $type2) = ($nodes->[$id1]->{type}, $nodes->[$id2]->{type});
		($cls1, $cls2) = ($clsid{$id1}, $clsid{$id2});
#		print "edge $i: ($id1, $id2)\n"; # for DEBUG
		if(defined $cls1){
			if(defined $cls2){ # cls1 & cls2
				if($cls1 != $cls2){ # belong to different clusters
					my $cluster1 = getCluster($cluster_pool, $cls1);
					my $cluster2 = getCluster($cluster_pool, $cls2);
					$cluster = mergeClusters($cluster_pool, $cluster1, $cluster2, $edge, $type1, $type2);
					if(defined $cluster){
						setClusterIDs($cluster, \%clsid);
					}
				}
			}
			else{ # cls1 & ^cls2
				$cluster = getCluster($cluster_pool, $cls1);
				if(insertAnchors($cluster_pool, $cluster, $edge, $type1, $type2)){
					$clsid{$id2} = $cluster->{id};
				}
			}
		}
		else{
			if(defined $cls2){  # ^cls1 & cls2
				$cluster = getCluster($cluster_pool, $cls2);
				if(insertAnchors($cluster_pool, $cluster, $edge, $type1, $type2)){
					$clsid{$id1} = $cluster->{id};
				}
			}
			else{ # ^cls1 & ^cls2
				$cluster = newCluster($cluster_pool);
				insertAnchors($cluster_pool, $cluster, $edge, $type1, $type2);
				$clsid{$id1} = $clsid{$id2} = $cluster->{id};
			}
		}
#		next if(!defined $cluster);
#		checkCluster($cluster); # For DEBUG
	}
	my @uniq_ids = grep {defined $cluster_pool->{clusters}->[$_]} (uniq values %clsid);
	my @sorted_ids = sort {scalar(keys %{$cluster_pool->{clusters}->[$b]->{locations}}) <=> scalar(keys %{$cluster_pool->{clusters}->[$a]->{locations}})} (@uniq_ids);
	my @clusters = ();
	my $no = 0;
	foreach my $id (@sorted_ids){
		$cluster = $cluster_pool->{clusters}->[$id];
		$cluster->{ends} = getMarginAnchorEnds($cluster->{anchor_pool}, $cluster->{locations}, $graph->{nodes});
		$cluster = simplifyCluster($cluster);
		$cluster->{id} = $no++;
		push(@clusters, $cluster);
	}
	$cluster_pool->{clusters} = \@clusters;
	$cluster_pool->{frees} = [];

	return $cluster_pool;
}

sub checkCluster
{
	my ($cluster) = @_;

	my $anchor_pool = $cluster->{anchor_pool};
	my $locations = $cluster->{locations};
	checkAnchorPool($anchor_pool);
	checkLocations($locations, $anchor_pool);
}

sub checkAnchorPool
{
	my ($anchor_pool) = @_;
	foreach my $anchor (@{$anchor_pool->{anchors}}){
		if(defined $anchor){
			for(my $end=0; $end<=1; $end++){
				next if(!defined $anchor->{$end}->{neighbor});
				my $neighbor = $anchor->{$end}->{neighbor};
				my ($nextAnchor, $nextEnd) = getAnchorEnd($anchor_pool, $neighbor);
				unless (defined $nextAnchor->{$nextEnd}->{neighbor} and
						$nextAnchor->{$nextEnd}->{neighbor}->{id} == $anchor->{id} and
						$nextAnchor->{$nextEnd}->{neighbor}->{end} == $end){
					die("inconsistent neighbor found\n");
				}
				for(my $tp=$anchor_pool->{nchannels}; $tp>=1; $tp--){
					unless(defined $anchor->{nodeId}->[$tp] and defined $nextAnchor->{nodeId}->[$tp] and
							$anchor->{nodeId}->[$tp] == $nextAnchor->{nodeId}->[$tp]){
						next;
					}
					if( round($anchor->{$end}->{coor}->[$tp] - $anchor->{$end^1}->{coor}->[$tp]) * round($nextAnchor->{$nextEnd}->{coor}->[$tp] - $anchor->{$end}->{coor}->[$tp]) < 0 ){
						die("order inconsistent for neighbor\n");
					}
					if( round($anchor->{$end}->{coor}->[$tp] - $anchor->{$end^1}->{coor}->[$tp]) * round($nextAnchor->{$nextEnd^1}->{coor}->[$tp] - $nextAnchor->{$nextEnd}->{coor}->[$tp]) < 0 ){
						die("order inconsistent for neighbor\n");
					}
				}
			}
			my ($leftAnchor, $leftEnd) = (defined $anchor->{0}->{neighbor}) ? getAnchorEnd($anchor_pool,  $anchor->{0}->{neighbor}) : (undef, undef);
			my ($rightAnchor, $rightEnd) = (defined $anchor->{1}->{neighbor}) ? getAnchorEnd($anchor_pool,  $anchor->{1}->{neighbor}) : (undef, undef);
			if(defined $leftAnchor and defined $rightAnchor){
				my @ctype = ();
				for(my $tp=$#{$leftAnchor->{nodeId}}; $tp>=1; $tp--){
					my $nodeId = $leftAnchor->{nodeId}->[$tp];
					next unless (defined $nodeId);
					next unless (defined $rightAnchor->{nodeId}->[$tp] and $rightAnchor->{nodeId}->[$tp] == $nodeId);
					if(!defined $anchor->{nodeId}->[$tp]){
						die("missing node found\n");
					}
					if($anchor->{nodeId}->[$tp] != $nodeId){
						die("conflicting node found\n");
					}
					push(@ctype, $tp);
				}
				next if(scalar(@ctype) > 0);
			}
			my @type1 = ();
			my @type2 = ();
			for(my $tp=$#{$anchor->{nodeId}}; $tp>=1; $tp--){
				my $nodeId = $anchor->{nodeId}->[$tp];
				next unless (defined $nodeId);
				if(defined $leftAnchor and defined $leftAnchor->{nodeId}->[$tp] and $leftAnchor->{nodeId}->[$tp] == $nodeId){
					push(@type1, $tp);
				}
				elsif(defined $rightAnchor and defined $rightAnchor->{nodeId}->[$tp] and $rightAnchor->{nodeId}->[$tp] == $nodeId){
					push(@type2, $tp);
				}
			}
			if( defined $leftAnchor and scalar(@type1) == 0 or
				defined $rightAnchor and scalar(@type2) == 0 ){
				die("invalid path met\n");
			}
		}
	}
}

sub checkLocations
{
	my ($locations, $anchor_pool) = @_;
	foreach my $nodeId (keys %{$locations}){
		my $loc = $locations->{$nodeId};
		my $type = $loc->{type};
		for(my $i=0; $i<2; $i++){
			my $ends = $loc->{ends}->{$i};
			my ($anchor, $end) = getAnchorEnd($anchor_pool, $ends);
			my ($nextAnchor, $nextEnd) = getNextAnchorEnd($anchor_pool, $anchor, $end^1);
			if(defined $nextAnchor->{nodeId}->[$type] and $nextAnchor->{nodeId}->[$type] == $nodeId){
				die("incomplete path found");
			}
		}
		my ($anchor, $end) = getAnchorEnd($anchor_pool, $loc->{ends}->{0});
		while(defined $anchor){
			if(!defined $anchor->{nodeId}->[$type] or $anchor->{nodeId}->[$type] != $nodeId){
				die("inconsistent path found");
			}
			($anchor, $end) = getNextAnchorEnd($anchor_pool, $anchor, $end, $type);
		}
	}
}

### CLUSTER related subroutines ###
sub newCluster
{
	my ($cluster_pool) = @_;
	my %cluster = ();
	$cluster{anchor_pool} = {nchannels=>$cluster_pool->{nchannels}, anchors=>[], frees=>[]};
	$cluster{locations} = {};
	my $id;
	$id = scalar(@{$cluster_pool->{clusters}});
	push(@{$cluster_pool->{clusters}}, \%cluster);
	$cluster{id} = $id;
	return \%cluster;
}

sub getCluster
{
	my ($cluster_pool, $clsid) = @_;
	return ($clsid >= 0 and $clsid < scalar(@{$cluster_pool->{clusters}})) ? $cluster_pool->{clusters}->[$clsid] : undef;
}

sub freeCluster
{
	my ($cluster_pool, $cluster) = @_;
	my $id = $cluster->{id};
	return if(!defined $cluster_pool->{clusters}->[$id]);
	$cluster_pool->{clusters}->[$id] = undef;
	push(@{$cluster_pool->{frees}}, $id);
}

sub setClusterIDs
{
	my ($cluster, $clsIDs) = @_;

	my $id = $cluster->{id};
	foreach my $nodeId (keys %{$cluster->{locations}}){
		$clsIDs->{$nodeId} = $id;
	}
}

sub simplifyCluster
{
	my ($cluster) = @_;

	my $anchor_pool = $cluster->{anchor_pool};
	clearSynonymousPairs($anchor_pool);
	my ($tp, $gap_len, $tp_cnt);
	my ($anchor, $end) = (getAnchor($anchor_pool, $cluster->{ends}->{0}->{id}), $cluster->{ends}->{0}->{end});
	while(defined $anchor){
		my ($nextAnchor, $nextEnd) = getNextAnchorEnd($anchor_pool, $anchor, $end);
		if(defined $nextAnchor){
			($gap_len, $tp_cnt) = (0, 0);
			for($tp=$anchor_pool->{nchannels}; $tp>=1; $tp--){
				if(defined $anchor->{nodeId}->[$tp]){
					if(defined $nextAnchor->{nodeId}->[$tp]){
						if($anchor->{nodeId}->[$tp] != $nextAnchor->{nodeId}->[$tp]){
							last;
						}
						$gap_len += abs($anchor->{$end^1}->{coor}->[$tp] - $nextAnchor->{$nextEnd}->{coor}->[$tp]);
						$tp_cnt++;
					}
					else{
						last;
					}
				}
				else{
					if(defined $nextAnchor->{nodeId}->[$tp]){
						last;
					}
				}
			}
			if( ($tp == 0) and ($tp_cnt > 0) ){
				$gap_len = sprintf("%.1f", ($gap_len / $tp_cnt));
				%{$anchor->{$end^1}} = %{dclone($nextAnchor->{$nextEnd^1})};
				$anchor->{len} = sprintf("%.1f", ($anchor->{len} + $gap_len + $nextAnchor->{len}));
				if(defined $nextAnchor->{$nextEnd^1}->{neighbor}){
					my ($tmpAnchor, $tmpEnd) = getAnchorEnd($anchor_pool, $nextAnchor->{$nextEnd^1}->{neighbor});
					$tmpAnchor->{$tmpEnd}->{neighbor} = {id=>$anchor->{id}, end=>$end^1};
				}
				setSynonymousPair($anchor_pool, $nextAnchor, $nextEnd^1, $anchor, $end^1);
				freeAnchor($anchor_pool, $nextAnchor);
				next;
			}
		}
		($anchor, $end) = ($nextAnchor, $nextEnd);
	}
	updateLocations($cluster->{locations}, $anchor_pool->{synonym});
	my $syn = $anchor_pool->{synonym}->{$cluster->{ends}->{1}->{id}}->{$cluster->{ends}->{1}->{end}};
	($cluster->{ends}->{1}->{id}, $cluster->{ends}->{1}->{end}) = ($syn->{id}, $syn->{end}) if(defined $syn);

	return $cluster;
}

### ANCHOR related subroutines ###
sub newAnchor
{
	my ($anchor_pool) = @_;
	my %anchor = ();
	if(!defined $anchor_pool){
		$anchor{0} = {coor=>[], neighbor=>undef};
		$anchor{1} = {coor=>[], neighbor=>undef};
		$anchor{nodeId} = [];
		return \%anchor;
	}
	$anchor{0} = {coor=>[(undef) x ($anchor_pool->{nchannels} + 1)], neighbor=>undef};
	$anchor{1} = {coor=>[(undef) x ($anchor_pool->{nchannels} + 1)], neighbor=>undef};
	$anchor{nodeId} = [(undef) x ($anchor_pool->{nchannels} + 1)];
	my $id;
	if(scalar(@{$anchor_pool->{frees}}) > 0){
		$id = shift(@{$anchor_pool->{frees}});
		$anchor_pool->{anchors}->[$id] = \%anchor;
	}
	else{
		$id = scalar(@{$anchor_pool->{anchors}});
		push(@{$anchor_pool->{anchors}}, \%anchor);
	}
	$anchor{id} = $id;

	return \%anchor;
}

sub freeAnchor
{
	my ($anchor_pool, $anchor) = @_;
	my $id = $anchor->{id};
	return if(!defined $anchor_pool->{anchors}->[$id]);
	$anchor_pool->{anchors}->[$id] = undef;
	push(@{$anchor_pool->{frees}}, $id);
}

sub freeAnchors
{
	my ($anchor_pool, $anchor, $end, $type) = @_;
	my $nextAnchor;
	while(defined $anchor){
		($nextAnchor, $end) = getNextAnchorEnd($anchor_pool, $anchor, $end, $type);
		freeAnchor($anchor_pool, $anchor);
		$anchor = $nextAnchor;
	}
}

sub freeAnchorsFromOne
{
	my ($anchor_pool, $anchor) = @_;
	my ($neighbor_anchor, $neighbor_end) = getNextAnchorEnd($anchor_pool, $anchor, 1);
	freeAnchors($anchor_pool, $anchor, 0);
	freeAnchors($anchor_pool, $neighbor_anchor, $neighbor_end);
}

sub getNewAnchorId
{
	my ($anchor_pool) = @_;
	return (scalar(@{$anchor_pool->{frees}}) > 0) ? @{$anchor_pool->{frees}}[0] : scalar(@{$anchor_pool->{anchors}});
}

sub getAnchorCount
{
	my ($anchor_pool) = @_;
	return (scalar(@{$anchor_pool->{anchors}}) - scalar(@{$anchor_pool->{frees}}));
}

sub clearSynonymousPairs
{
	my ($anchor_pool) = @_;
	$anchor_pool->{synonym} = {};
}

sub setNodeIdProtected
{
	my ($anchor_pool, $nodeId) = @_;
	$anchor_pool->{proteected}->{$nodeId} = 1;
}

sub setSynonymousPair
{
	my ($anchor_pool, $anchor, $end, $newAnchor, $newEnd) = @_;
	
	if(!defined $anchor_pool->{synonym}->{$anchor->{id}}->{$end}){
		$anchor_pool->{synonym}->{$anchor->{id}}->{$end} = {id=>$newAnchor->{id}, end=>$newEnd};
	}
}

sub calculateGap
{
	my ($anchor_pool, $anchor, $end, $nextAnchor, $nextEnd) = @_;

	my ($gap_len, $tp_cnt) = (0, 0);
	for(my $tp=$anchor_pool->{nchannels}; $tp>=1; $tp--){
		if(defined $anchor->{nodeId}->[$tp] and defined $nextAnchor->{nodeId}->[$tp]){
			if($anchor->{nodeId}->[$tp] == $nextAnchor->{nodeId}->[$tp]){
				$gap_len += abs($anchor->{$end^1}->{coor}->[$tp] - $nextAnchor->{$nextEnd}->{coor}->[$tp]);
				$tp_cnt++;
			}
		}
	}
	die("Can not calculate gap size\n") if($tp_cnt == 0);
	$gap_len = sprintf("%.1f", ($gap_len / $tp_cnt));
	return $gap_len;
}

sub getAnchor
{
	my ($anchor_pool, $anchor_id) = @_;
	return (defined $anchor_pool->{anchors}) ?  $anchor_pool->{anchors}->[$anchor_id] : undef;
}

sub getAnchorEnd
{
	my ($anchor_pool, $id_end) = @_;
	return (getAnchor($anchor_pool, $id_end->{id}), $id_end->{end});
}

sub getNextAnchorEnd
{
	my ($anchor_pool, $anchor, $end, $type) = @_;
	my $neighbor = $anchor->{$end^1}->{neighbor};
	if(!defined $neighbor){
		return (undef, undef);
	}
	my ($nextAnchor, $nextEnd) = getAnchorEnd($anchor_pool, $neighbor);
	if(defined $type and defined $anchor->{nodeId}->[$type]){
		if( !defined $nextAnchor->{nodeId}->[$type] or $nextAnchor->{nodeId}->[$type] != $anchor->{nodeId}->[$type]){
			return (undef, undef);
		}
	}
	if($nextAnchor->{$nextEnd}->{neighbor}->{id} != $anchor->{id} or $nextAnchor->{$nextEnd}->{neighbor}->{end} != ($end^1)){
		die("Invalid structure found\n"); # Debug
	}
	return ($nextAnchor, $nextEnd);
}

sub getMarginAnchorEnds
{
	my ($anchor_pool, $locations, $nodes) = @_;

	my ($firstAnchor, $firstEnd);
	my ($lastAnchor, $lastEnd);
	my ($anchor, $end);
	my ($nodeId, $loc);
# get one margin end
	while (($nodeId, $loc) = each (%{$locations}) ){
		for my $locEnd (values %{$loc->{ends}}){
			($anchor, $end) = getAnchorEnd($anchor_pool, $locEnd);
			if(!defined $anchor->{$end}->{neighbor}){
				($firstAnchor, $firstEnd) = ($anchor, $end);
				last;
			}
		}
	}
	my %order = ();
	my $no = 0;
# get the other margin end and mark the order of all the ends
	($anchor, $end) = ($firstAnchor, $firstEnd);
	while(defined $anchor){
		$order{$anchor->{id}} = ($end ? [$no+2, $no+1] : [$no+1, $no+2]);
		$no+=2;
		if(!defined $anchor->{$end^1}->{neighbor}){
			($lastAnchor, $lastEnd) = ($anchor, $end^1);
			last;
		}
		($anchor, $end) = getAnchorEnd($anchor_pool, $anchor->{$end^1}->{neighbor});
	}
# get the weights of two opposite directions according to current ordering
	my @weights = ( (0) x 2 );
	my $type;
	while (($nodeId, $loc) = each (%{$locations}) ){
		my $node = $nodes->[$nodeId];
		$type = $loc->{type};
		my @ends = values %{$loc->{ends}};
		my $orient = (($order{$ends[0]->{id}}->[$ends[0]->{end}] - $order{$ends[1]->{id}}->[$ends[1]->{end}]) * ($anchor_pool->{anchors}->[$ends[0]->{id}]->{$ends[0]->{end}}->{coor}->[$type] - $anchor_pool->{anchors}->[$ends[1]->{id}]->{$ends[1]->{end}}->{coor}->[$type]) > 0) ? 0 : 1;
		$weights[$orient] += $node->{length};
	}
# if the negative direction is dominant, then change the ordering
	return ($weights[0] >= $weights[1]) ? {0=>{id=>$firstAnchor->{id}, end=>$firstEnd}, 1=>{id=>$lastAnchor->{id}, end=>$lastEnd}} :
			{0=>{id=>$lastAnchor->{id}, end=>$lastEnd}, 1=>{id=>$firstAnchor->{id}, end=>$firstEnd}};
}

sub copyAnchor
{
	my ($anchor_pool, $anchor, $end) = @_;

	my $newAnchor = newAnchor($anchor_pool);
	$newAnchor->{len} = $anchor->{len};
	@{$newAnchor->{nodeId}} = @{dclone($anchor->{nodeId})};
	@{$newAnchor->{0}->{coor}} = @{dclone($anchor->{$end}->{coor})};
	@{$newAnchor->{1}->{coor}} = @{dclone($anchor->{$end^1}->{coor})};

	return $newAnchor;
}

sub collectPositions
{
	my ($ctg, $start, $end) = @_;

	my $positions = $ctg->{Position};
	my $idx1 = locate($positions, $start);
	my $idx2 = locate($positions, $end);
	if($start <= $end){
		return [@{$positions}[$idx1..($idx2-1)]];
	}

	my @reversed = ();
	my $length = $ctg->{ContigLength};
	for(my $i=($idx1-1); $i>$idx2; $i--){
		push(@reversed, $length - $positions->[$i]);
	}

	return \@reversed;
}

sub getFirstAnchorEnd
{
	my ($anchor_pool, $ends, $type) = @_;
	my ($firstAnchor, $firstEnd) = getAnchorEnd($anchor_pool, $ends);
	my $nodeId = $firstAnchor->{nodeId}->[$type];
	$firstEnd ^= 1;
	my ($nextAnchor, $nextEnd) = getNextAnchorEnd($anchor_pool, $firstAnchor, $firstEnd);
	while(defined $nextAnchor){
		last if(!defined $nextAnchor->{nodeId}->[$type] or $nextAnchor->{nodeId}->[$type] != $nodeId);
		($firstAnchor, $firstEnd) = ($nextAnchor, $nextEnd);
		($nextAnchor, $nextEnd) = getNextAnchorEnd($anchor_pool, $firstAnchor, $firstEnd);
	}
	$firstEnd ^= 1;
	return ($firstAnchor, $firstEnd);
}

sub unifyIDs4HomologousNeighbors
{
	my ($anchor_pool, $firstAnchor, $firstEnd, $type, $graph) = @_;
	my $nodes = $graph->{nodes};
	my $contigs = $graph->{cmaps}->[$type]->{contigs};

	my $nodeId = $firstAnchor->{nodeId}->[$type];
	my $ctg = $contigs->{$nodes->[$nodeId]->{ori_id}};
	my ($anchor, $end) = ($firstAnchor, $firstEnd^1);
	my ($bstart, $bend) = ($anchor->{$end}->{coor}->[$type], $anchor->{$end^1}->{coor}->[$type]);
	my $dir = ($bend - $bstart) > 0 ? 1 : -1;
	my $offset = $bend;
	
	my $length_left = (($dir > 0) ? ($nodes->[$nodeId]->{length} - $bend) : $bend);
	my $mercy = 50000;
	my ($xstart, $xend, $ystart, $yend);
	my ($nextAnchor, $nextEnd) = getNextAnchorEnd($anchor_pool, $anchor, $end);
	while(defined $nextAnchor){
		my $gap = calculateGap($anchor_pool, $anchor, $end, $nextAnchor, $nextEnd);
		my $nodeId2 = $nextAnchor->{nodeId}->[$type];
		if(!defined $nodeId2){
			last if($gap + $nextAnchor->{len} > $length_left);
			$ystart = $offset + $gap * $dir;
			$yend = $ystart + $nextAnchor->{len} * $dir;
			$offset = $yend;
			$length_left -= ($yend - $ystart) * $dir;
			$nextAnchor->{nodeId}->[$type] = $nodeId;
			$nextAnchor->{$nextEnd}->{coor}->[$type] = $ystart;
			$nextAnchor->{$nextEnd^1}->{coor}->[$type] = $yend;
			($anchor, $end) = ($nextAnchor, $nextEnd);
			($nextAnchor, $nextEnd) = getNextAnchorEnd($anchor_pool, $anchor, $end);
			next;
		}
		last if($nodes->[$nodeId2]->{length} > $length_left + $mercy);
		$xstart = $nextAnchor->{$nextEnd}->{coor}->[$type];
		$xend = $nextAnchor->{$nextEnd^1}->{coor}->[$type];
		$ystart = $offset + $gap * $dir;
		$yend = $ystart + $nextAnchor->{len} * $dir;
		($firstAnchor, $firstEnd) = ($anchor, $end);
		($anchor, $end) = ($nextAnchor, $nextEnd);
		($nextAnchor, $nextEnd) = getNextAnchorEnd($anchor_pool, $anchor, $end);
		while(defined $nextAnchor and defined $nextAnchor->{nodeId}->[$type] and $nextAnchor->{nodeId}->[$type] == $nodeId2){	
			$xend = $nextAnchor->{$nextEnd^1}->{coor}->[$type];
			$gap = calculateGap($anchor_pool, $anchor, $end, $nextAnchor, $nextEnd);
			$yend += ($gap + $nextAnchor->{len}) * $dir;
			($anchor, $end) = ($nextAnchor, $nextEnd);
			($nextAnchor, $nextEnd) = getNextAnchorEnd($anchor_pool, $anchor, $end);
		}
		$offset = $yend;
		$length_left -= ($yend - $ystart) * $dir;
		my $ctg2 = $contigs->{$nodes->[$nodeId2]->{ori_id}};
		my $positions1 = collectPositions($ctg, $ystart, $yend);
		my $positions2 = collectPositions($ctg2, $xstart, $xend);
		my $alignment = globalAlign($positions1, $positions2, 0.9);
		if(scalar(@{$alignment}) == 0){
			# no significant overlap found
			($anchor, $end) = ($firstAnchor, $firstEnd);
			last;
		}
		my ($slope, $intercept) = calculateLinearParams($xstart, $xend, $ystart, $yend);
		$anchor_pool->{homo}->{$nodeId2} = {id=>$nodeId, xstart=>$xstart, xend=>$xend, ystart=>$ystart, yend=>$yend, slope=>$slope, intercept=>$intercept};
		($anchor, $end) = ($firstAnchor, $firstEnd);
		($nextAnchor, $nextEnd) = getNextAnchorEnd($anchor_pool, $anchor, $end);
		while(defined $nextAnchor and defined $nextAnchor->{nodeId}->[$type] and $nextAnchor->{nodeId}->[$type] == $nodeId2){	
			$nextAnchor->{nodeId}->[$type] = $nodeId;
			for(my $i=0; $i<=1; $i++){
				$nextAnchor->{$i}->{coor}->[$type] = sprintf("%.1f", $nextAnchor->{$i}->{coor}->[$type] * $slope + $intercept);
			}
			($anchor, $end) = ($nextAnchor, $nextEnd);
			($nextAnchor, $nextEnd) = getNextAnchorEnd($anchor_pool, $anchor, $end);
		}
	}

	return ($anchor, $end^1);
}

# this subroutine is to fill the missing coordinates between surrounding anchors
sub fillMissingCoorsExtended
{
	my ($middleAnchor, $leftAnchor, $leftEnd, $rightAnchor, $rightEnd) = @_;

	my @btypes = ();
	my @etypes = ();
	for(my $tp=$#{$leftAnchor->{nodeId}}; $tp>=1; $tp--){
		my $nodeId = $leftAnchor->{nodeId}->[$tp];
		next unless (defined $nodeId);
		if(defined $rightAnchor->{nodeId}->[$tp] and $rightAnchor->{nodeId}->[$tp] == $nodeId){
			if(defined $middleAnchor->{nodeId}->[$tp]){
				if($middleAnchor->{nodeId}->[$tp] != $nodeId){
					return 0; # conflict found
				}
				# no need to fill
			}
			else{
				push(@etypes, $tp);
			}
			next;
		}
		if(defined $middleAnchor->{nodeId}->[$tp] and $middleAnchor->{nodeId}->[$tp] == $nodeId){
			push(@btypes, $tp);
		}
	}

	if(scalar(@etypes) > 0){
		if(scalar(@btypes) == 0){
			die("estimateInterval: invalid path met\n");
		}
		my ($xstart, $xend) = ($leftAnchor->{$leftEnd^1}->{coor}->[$btypes[0]], $leftAnchor->{$leftEnd}->{coor}->[$btypes[0]]);
		foreach my $tp (@etypes){
			$middleAnchor->{nodeId}->[$tp] = $leftAnchor->{nodeId}->[$tp];
			my ($ystart, $yend) = ($leftAnchor->{$leftEnd^1}->{coor}->[$tp], $leftAnchor->{$leftEnd}->{coor}->[$tp]);
			my ($slope, $intercept) = calculateLinearParams($xstart, $xend, $ystart, $yend);
			for(my $i=0; $i<2; $i++){
				$middleAnchor->{$i}->{coor}->[$tp] = sprintf("%.1f", $middleAnchor->{$i}->{coor}->[$btypes[0]] * $slope  + $intercept);
			}
		}
	}

	return 1;
}

sub fillMissingCoors
{
	my ($middleAnchor, $leftAnchor, $leftEnd, $rightAnchor, $rightEnd, $type) = @_;

	my ($flag, $alpha, $beta) = estimateIntervalParams($middleAnchor, $leftAnchor, $leftEnd, $rightAnchor, $rightEnd, $type);
# middle anchor and at least one of surrounding anchors have "$type" coordinates, but can not share other types of coordinates at this point
	for(my $tp=$#{$leftAnchor->{nodeId}}; $tp>=1; $tp--){
		next if($tp == $type);
		next unless (defined $leftAnchor->{nodeId}->[$tp] and defined $rightAnchor->{nodeId}->[$tp]);
		next unless ($leftAnchor->{nodeId}->[$tp] == $rightAnchor->{nodeId}->[$tp]);
		if(defined $middleAnchor->{nodeId}->[$tp]){
			return 0; # conflict found
		}
		$middleAnchor->{nodeId}->[$tp] = $leftAnchor->{nodeId}->[$tp];
		if($flag == 3){
			$middleAnchor->{0}->{coor}->[$tp] = $leftAnchor->{$leftEnd}->{coor}->[$tp];
			$middleAnchor->{0}->{coor}->[$tp] += $alpha * ($rightAnchor->{$rightEnd}->{coor}->[$tp] - $leftAnchor->{$leftEnd}->{coor}->[$tp]);
			$middleAnchor->{1}->{coor}->[$tp] = $leftAnchor->{$leftEnd}->{coor}->[$tp];
			$middleAnchor->{1}->{coor}->[$tp] += $beta * ($rightAnchor->{$rightEnd}->{coor}->[$tp] - $leftAnchor->{$leftEnd}->{coor}->[$tp]);
		}
		elsif($flag == 1){
			$middleAnchor->{0}->{coor}->[$tp] = $leftAnchor->{$leftEnd^1}->{coor}->[$tp];
			$middleAnchor->{0}->{coor}->[$tp] += $alpha * ($leftAnchor->{$leftEnd}->{coor}->[$tp] - $leftAnchor->{$leftEnd^1}->{coor}->[$tp]);
			$middleAnchor->{1}->{coor}->[$tp] = $leftAnchor->{$leftEnd^1}->{coor}->[$tp];
			$middleAnchor->{1}->{coor}->[$tp] += $beta * ($leftAnchor->{$leftEnd}->{coor}->[$tp] - $leftAnchor->{$leftEnd^1}->{coor}->[$tp]);
		}
		else{ # 2
			$middleAnchor->{0}->{coor}->[$tp] = $rightAnchor->{$rightEnd}->{coor}->[$tp];
			$middleAnchor->{0}->{coor}->[$tp] += $alpha * ($rightAnchor->{$rightEnd^1}->{coor}->[$tp] - $rightAnchor->{$rightEnd}->{coor}->[$tp]);
			$middleAnchor->{1}->{coor}->[$tp] = $rightAnchor->{$rightEnd}->{coor}->[$tp];
			$middleAnchor->{1}->{coor}->[$tp] += $beta * ($rightAnchor->{$rightEnd^1}->{coor}->[$tp] - $rightAnchor->{$rightEnd}->{coor}->[$tp]);
		}
	}
	return 1;
}

sub splitAnchor
{
	my ($anchor_pool, $anchor, $end, $alpha) = @_;
	my ($firstAnchor, $secondAnchor) = (newAnchor($anchor_pool), newAnchor());
	for(my $tp=$#{$anchor->{nodeId}}; $tp>=1; $tp--){
		next unless (defined $anchor->{nodeId}->[$tp]);
		$firstAnchor->{nodeId}->[$tp] = $anchor->{nodeId}->[$tp];
		$firstAnchor->{0}->{coor}->[$tp] = $firstAnchor->{1}->{coor}->[$tp] = $anchor->{$end}->{coor}->[$tp];
		$firstAnchor->{1}->{coor}->[$tp] += sprintf("%.1f", $alpha * ($anchor->{$end^1}->{coor}->[$tp] - $anchor->{$end}->{coor}->[$tp]));
		$secondAnchor->{nodeId}->[$tp] = $anchor->{nodeId}->[$tp];
		$secondAnchor->{$end}->{coor}->[$tp] = $firstAnchor->{1}->{coor}->[$tp];
		$secondAnchor->{$end^1}->{coor}->[$tp] = $anchor->{$end^1}->{coor}->[$tp];
	}
	$firstAnchor->{len} = sprintf("%.1f", $anchor->{len} * $alpha);
	$secondAnchor->{id} = $anchor->{id};
	$secondAnchor->{len} = sprintf("%.1f", $anchor->{len} - $firstAnchor->{len});
	%{$secondAnchor->{0}->{neighbor}} = %{dclone($anchor->{0}->{neighbor})} if(defined $anchor->{0}->{neighbor});
	%{$secondAnchor->{1}->{neighbor}} = %{dclone($anchor->{1}->{neighbor})} if(defined $anchor->{1}->{neighbor});

	return ($firstAnchor, $secondAnchor);
}

sub splitAnchorExtended
{
	my ($anchor_pool, $anchor, $end, $alpha, $rev) = @_;
	my ($firstAnchor, $secondAnchor) = (newAnchor($anchor_pool), newAnchor());
	for(my $tp=$#{$anchor->{nodeId}}; $tp>=1; $tp--){
		next unless (defined $anchor->{nodeId}->[$tp]);
		$firstAnchor->{nodeId}->[$tp] = $anchor->{nodeId}->[$tp];
		$firstAnchor->{$rev}->{coor}->[$tp] = $firstAnchor->{$rev^1}->{coor}->[$tp] = $anchor->{$end}->{coor}->[$tp];
		$firstAnchor->{$rev^1}->{coor}->[$tp] += sprintf("%.1f", $alpha * ($anchor->{$end^1}->{coor}->[$tp] - $anchor->{$end}->{coor}->[$tp]));
		$secondAnchor->{nodeId}->[$tp] = $anchor->{nodeId}->[$tp];
		$secondAnchor->{$end}->{coor}->[$tp] = $firstAnchor->{$rev^1}->{coor}->[$tp];
		$secondAnchor->{$end^1}->{coor}->[$tp] = $anchor->{$end^1}->{coor}->[$tp];
	}
	$firstAnchor->{len} = sprintf("%.1f", $anchor->{len} * $alpha);
	$secondAnchor->{id} = $anchor->{id};
	$secondAnchor->{len} = sprintf("%.1f", $anchor->{len} - $firstAnchor->{len});
	%{$secondAnchor->{0}->{neighbor}} = %{dclone($anchor->{0}->{neighbor})} if(defined $anchor->{0}->{neighbor});
	%{$secondAnchor->{1}->{neighbor}} = %{dclone($anchor->{1}->{neighbor})} if(defined $anchor->{1}->{neighbor});

	return ($firstAnchor, $secondAnchor);
}

sub joinAnchors
{
	my ($anchor_pool, $anchor1, $end1, $anchor2, $end2, $type) = @_;
	my $bSuccess = 1;
	my $newAnchor = newAnchor($anchor_pool);
	for(my $tp=$anchor_pool->{nchannels}; $tp>0; $tp--){
		if(defined $anchor1->{nodeId}->[$tp]){
			if( ($tp != $type) and defined $anchor2->{nodeId}->[$tp]){
				$bSuccess = 0; # conflict found
				last;
			}
			$newAnchor->{nodeId}->[$tp] = $anchor1->{nodeId}->[$tp];
			$newAnchor->{0}->{coor}->[$tp] = $anchor1->{$end1}->{coor}->[$tp];
			$newAnchor->{1}->{coor}->[$tp] = $anchor1->{$end1^1}->{coor}->[$tp];
			next;
		}
		if(defined $anchor2->{nodeId}->[$tp]){
			$newAnchor->{nodeId}->[$tp] = $anchor2->{nodeId}->[$tp];
			$newAnchor->{0}->{coor}->[$tp] = $anchor2->{$end2}->{coor}->[$tp];
			$newAnchor->{1}->{coor}->[$tp] = $anchor2->{$end2^1}->{coor}->[$tp];
		}
	}
	if(!$bSuccess){
		freeAnchor($anchor_pool, $newAnchor);
		return undef;
	}
	$newAnchor->{len} = sprintf("%.1f", ($anchor1->{len} + $anchor2->{len}) / 2);

	return $newAnchor;
}

sub joinAnchorsExtended
{
	my ($anchor_pool, $anchor1, $end1, $anchor2, $end2, $rev) = @_;
	my $bSuccess = 1;
	my $newAnchor = newAnchor($anchor_pool);
	for(my $tp=$anchor_pool->{nchannels}; $tp>0; $tp--){
		if(defined $anchor1->{nodeId}->[$tp]){
			$newAnchor->{nodeId}->[$tp] = $anchor1->{nodeId}->[$tp];
			if(defined $anchor2->{nodeId}->[$tp]){
				if($anchor2->{nodeId}->[$tp] != $anchor1->{nodeId}->[$tp]){
					$bSuccess = 0; # conflict found
					last;
				}
				$newAnchor->{$rev}->{coor}->[$tp] = sprintf("%.1f", ($anchor1->{$end1}->{coor}->[$tp] + $anchor2->{$end2}->{coor}->[$tp]) / 2);
				$newAnchor->{$rev^1}->{coor}->[$tp] = sprintf("%.1f", ($anchor1->{$end1^1}->{coor}->[$tp] + $anchor2->{$end2^1}->{coor}->[$tp])/ 2);
			}
			else{
				$newAnchor->{$rev}->{coor}->[$tp] = $anchor1->{$end1}->{coor}->[$tp];
				$newAnchor->{$rev^1}->{coor}->[$tp] = $anchor1->{$end1^1}->{coor}->[$tp];
			}
			next;
		}
		if(defined $anchor2->{nodeId}->[$tp]){
			$newAnchor->{nodeId}->[$tp] = $anchor2->{nodeId}->[$tp];
			$newAnchor->{$rev}->{coor}->[$tp] = $anchor2->{$end2}->{coor}->[$tp];
			$newAnchor->{$rev^1}->{coor}->[$tp] = $anchor2->{$end2^1}->{coor}->[$tp];
		}
	}
	if(!$bSuccess){
		freeAnchor($anchor_pool, $newAnchor);
		return undef;
	}
	$newAnchor->{len} = sprintf("%.1f", ($anchor1->{len} + $anchor2->{len}) / 2);

	return $newAnchor;
}

sub getAnchorInterval
{
	my ($anchor, $startEnd, $type) = @_;
	return [sprintf("%.1f", $anchor->{$startEnd}->{coor}->[$type]), sprintf("%.1f", $anchor->{$startEnd^1}->{coor}->[$type])];
}

sub estimateIntervalParams
{
	my ($middleAnchor, $leftAnchor, $leftEnd, $rightAnchor, $rightEnd, $type) = @_;
	my $nodeId = $middleAnchor->{nodeId}->[$type];
	my $it = getAnchorInterval($middleAnchor, 0, $type);

# middle anchor and at least one of the surrounding anchors have the "$type" coordinates, the corresponding situation is marked as $flag
	my $flag = 0;
	$flag |= 0x1 if(defined $leftAnchor->{nodeId}->[$type] and $leftAnchor->{nodeId}->[$type] == $nodeId);
	$flag |= 0x2 if(defined $rightAnchor->{nodeId}->[$type] and $rightAnchor->{nodeId}->[$type] == $nodeId);
	my ($leftCoor, $rightCoor, $interval);
	if($flag == 3){ # both
		$leftCoor = $leftAnchor->{$leftEnd}->{coor}->[$type];
		$rightCoor = $rightAnchor->{$rightEnd}->{coor}->[$type];
	}
	elsif($flag == 1){ # left only
		$leftCoor = $leftAnchor->{$leftEnd^1}->{coor}->[$type];
		$rightCoor = $leftAnchor->{$leftEnd}->{coor}->[$type];
	}
	elsif($flag == 2){ # right only
		$leftCoor = $rightAnchor->{$rightEnd}->{coor}->[$type];
		$rightCoor = $rightAnchor->{$rightEnd^1}->{coor}->[$type];
	}
	else{ # 0
		# definitely can not happen
		die("estimateInterval: fatal error occurs\n");
	}

	$interval = ($rightCoor - $leftCoor);
	if($interval == 0){
		# usually can't happen, because the algorithm merges the overlapped parts in case flag==3;
		#  in other cases each anchor can't have size of 0
		if($flag == 3){
			die("estimateInterval: no space for insertion\n");
		}
		else{ # 1 | 2
			die("estimateInterval: invalid size-zero anchor found\n");
		}
	}
	if( ($it->[1] - $it->[0]) * $interval < 0 ){
		# usually can't happen, because the algorithm should have adjusted the directions beforehand
		die("estimateInterval: conflicting directions\n");
	}
	my $alpha = ($it->[0] - $leftCoor) / $interval;
	my $beta = ($it->[1] - $leftCoor) / $interval;

	return ($flag, $alpha, $beta);
}

sub checkEndConflict
{
	my ($anchor1, $end1, $anchor2, $end2, $type, $dir) = @_;
	my $it1 = getAnchorInterval($anchor1, $end1, $type);
	my $it2 = getAnchorInterval($anchor2, $end2, $type);
	my $res = ($it1->[0] - $it2->[0]) * $dir;
	if($res < 0){
		return (defined $anchor2->{$end2}->{neighbor}) ? 1 : 0;
	}
	if($res > 0){
		return (defined $anchor1->{$end1}->{neighbor}) ? 1 : 0;
	}
	return (defined $anchor1->{$end1}->{neighbor} and defined $anchor2->{$end2}->{neighbor}) ? 1 : 0;
}

sub setupPaths
{
	my ($anchor_pool, $glues, $nodeId1, $nodeId2, $type1, $type2) = @_;
	my ($path1, $path2);
	my @ids = ();
	my ($anchor, $preAnchor);
	for(my $i=0; $i<scalar(@$glues); $i++){
		my $glue = $glues->[$i];
		$anchor = newAnchor($anchor_pool);
		$anchor->{len} = $glue->{len};
		$anchor->{nodeId}->[$type1] = $nodeId1;
		$anchor->{0}->{coor}->[$type1] = $glue->{first}[0];
		$anchor->{1}->{coor}->[$type1] = $glue->{second}[0];
		$anchor->{nodeId}->[$type2] = $nodeId2;
		$anchor->{0}->{coor}->[$type2] = $glue->{first}[1];
		$anchor->{1}->{coor}->[$type2] = $glue->{second}[1];
		while( defined $preAnchor and
			$anchor->{0}->{coor}->[$type2] <= $preAnchor->{0}->{coor}->[$type2] and
			$preAnchor->{1}->{coor}->[$type2] <= $anchor->{1}->{coor}->[$type2]){
			# previous anchor is covered by anchor
			my $tmpAnchor = $preAnchor;
			($preAnchor) = getNextAnchorEnd($anchor_pool, $tmpAnchor, 1);
			pop(@ids);
			freeAnchor($anchor_pool, $tmpAnchor);
		}
		if(defined $preAnchor){
			# assert $preAnchor->{0}->{coor}->[$type1] <= $anchor->{0}->{coor}->[$type1]
			if( $anchor->{1}->{coor}->[$type1] <= $preAnchor->{1}->{coor}->[$type1] or
				$preAnchor->{0}->{coor}->[$type2] <= $anchor->{0}->{coor}->[$type2] and
				$anchor->{1}->{coor}->[$type2] <= $preAnchor->{1}->{coor}->[$type2]){
				# anchor is covered by previous anchor
				freeAnchor($anchor_pool, $anchor);
				next;
			}
			foreach my $tp ($type1, $type2){ # to garantee non-overlap between adjacent anchors
				if( ($anchor->{0}->{coor}->[$tp] - $anchor->{1}->{coor}->[$tp]) * ($preAnchor->{1}->{coor}->[$tp] - $anchor->{0}->{coor}->[$tp]) < 0 ){
					my $diff = sprintf("%.1f", ($preAnchor->{1}->{coor}->[$tp] - $anchor->{0}->{coor}->[$tp]) / 2);
					$anchor->{0}->{coor}->[$tp] += $diff;
					$preAnchor->{1}->{coor}->[$tp] = $anchor->{0}->{coor}->[$tp];
				}
			}
			$preAnchor->{1}->{neighbor} = {id=>$anchor->{id}, end=>0};
			$anchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
		}
		push(@ids, $anchor->{id});
		$preAnchor = $anchor;
	}
	my $firstAnchor = getAnchor($anchor_pool, $ids[0]);
	my $lastAnchor = getAnchor($anchor_pool, $ids[-1]);
	if($firstAnchor->{0}->{coor}->[$type1] < $lastAnchor->{1}->{coor}->[$type1]){
		$path1 = {type=>$type1, ends=>{0=>{id=>$ids[0], end=>0}, 1=>{id=>$ids[-1], end=>1}}};
	}
	else{
		$path1 = {type=>$type1, ends=>{0=>{id=>$ids[-1], end=>1}, 1=>{id=>$ids[0], end=>0}}};
	}
	if($firstAnchor->{0}->{coor}->[$type2] < $lastAnchor->{1}->{coor}->[$type2]){
		$path2 = {type=>$type2, ends=>{0=>{id=>$ids[0], end=>0}, 1=>{id=>$ids[-1], end=>1}}};
	}
	else{
		$path2 = {type=>$type2, ends=>{0=>{id=>$ids[-1], end=>1}, 1=>{id=>$ids[0], end=>0}}};
	}
	return ($path1, $path2);
}

sub addLocations
{
	my ($cluster, $locations, $protected_id) = @_;
	my ($nodeId, $loc);
	my $newLocations = $cluster->{locations};
	while (($nodeId, $loc) = each (%{$locations})){
		next if($nodeId == $protected_id);
		%{$newLocations->{$nodeId}} = %{dclone($loc)};
	}
}

sub updateLocations
{
	my ($locations, $synonyms, $protected_id) = @_;
	my ($id, $end, $newId, $newEnd);
	my ($nodeId, $loc);
	while (($nodeId, $loc) = each (%{$locations}) ){
		next if(defined $protected_id and $nodeId == $protected_id);
		for (my $i=0; $i<2; $i++){
			($id, $end) = ($loc->{ends}->{$i}->{id}, $loc->{ends}->{$i}->{end});
			if(!defined $synonyms->{$id}){
				next;
			}
			my $syn = $synonyms->{$id}->{$end};
			next if(!defined $syn);
			($loc->{ends}->{$i}->{id}, $loc->{ends}->{$i}->{end}) = ($syn->{id}, $syn->{end});
		}
	}
}

sub deleteLocations
{
	my ($locations, $obsolete) = @_;
	foreach (@{$obsolete}){
		delete($locations->{$_});
	}
}

sub collectHomoPairs
{
	my ($g_homo, $homo) = @_;
	return if(!defined $homo);
	my ($nodeId, $hm);
	while (($nodeId, $hm) = each (%{$homo})){
		$g_homo->{$nodeId} = $hm;
	}
}

sub handleHomoNode
{
	my ($edge, $homo) = @_;
	my ($id1, $id2) = ($edge->{id1}, $edge->{id2});
	my ($hm1, $hm2) = ($homo->{$id1}, $homo->{$id2});
	return if(!defined $hm1 and !defined $hm2);
	$edge->{id1} = (defined $hm1) ? $hm1->{id} : $id1;
	$edge->{id2} = (defined $hm2) ? $hm2->{id} : $id2;
	my $glues = [@{dclone($edge->{glues})}];
	for(my $i=0; $i<scalar(@$glues); $i++){
		my $glue = $glues->[$i];
		if(defined $hm1){
			$glue->{first}[0] = $hm1->{slope} * $glue->{first}[0] + $hm1->{intercept};
			$glue->{second}[0] = $hm1->{slope} * $glue->{second}[0] + $hm1->{intercept};
		}
		if(defined $hm2){
			$glue->{first}[1] = $hm2->{slope} * $glue->{first}[1] + $hm2->{intercept};
			$glue->{second}[1]= $hm2->{slope} * $glue->{second}[1] + $hm2->{intercept};
		}
	}
	$edge->{glues} = $glues;
}

sub expandPath
{
	my ($anchor_pool, $path, $type, $graph) = @_;

	my ($firstAnchor, $firstEnd) = getFirstAnchorEnd($anchor_pool, $path->{ends}->{0}, $type);
	my ($lastAnchor, $lastEnd) = getFirstAnchorEnd($anchor_pool, $path->{ends}->{1}, $type);

	($firstAnchor, $firstEnd) = unifyIDs4HomologousNeighbors($anchor_pool, $firstAnchor, $firstEnd, $type, $graph);
	($lastAnchor, $lastEnd) = unifyIDs4HomologousNeighbors($anchor_pool, $lastAnchor, $lastEnd, $type, $graph);

	($path->{ends}->{0}->{id}, $path->{ends}->{0}->{end}) = ($firstAnchor->{id}, $firstEnd);
	($path->{ends}->{1}->{id}, $path->{ends}->{1}->{end}) = ($lastAnchor->{id}, $lastEnd);
	
	return (scalar(keys %{$anchor_pool->{mono}}) > 0) ? 1 : 0;
}

sub updatePathEnds
{
	my ($anchor_pool, $firstAnchorId, $firstEnd, $lastAnchorId, $lastEnd, $type) = @_;
	my ($firstAnchor, $lastAnchor);
	($firstAnchor, $firstEnd) = getFirstAnchorEnd($anchor_pool, {id=>$firstAnchorId, end=>$firstEnd}, $type);
	($lastAnchor, $lastEnd) = getFirstAnchorEnd($anchor_pool, {id=>$lastAnchorId, end=>$lastEnd}, $type);

	return {"0" => {id=>$firstAnchor->{id}, end=>$firstEnd}, "1" => {id=>$lastAnchor->{id}, end=>$lastEnd}};
}

sub getExtendedIntervals
{
	my ($preAnchor, $preEnd, $anchor1, $end1, $anchor2, $end2) = @_;
	my ($it1, $it2);
	my @tp1 = ();
	my @tp2 = ();
	for(my $tp=$#{$preAnchor->{nodeId}}; $tp>=1; $tp--){
		my $nodeId = $preAnchor->{nodeId}->[$tp];
		next if(!defined $nodeId);
		if(defined $anchor1->{nodeId}->[$tp] and $anchor1->{nodeId}->[$tp] == $nodeId){
			push(@tp1, $tp);
		}
		if(defined $anchor2->{nodeId}->[$tp] and $anchor2->{nodeId}->[$tp] == $nodeId){
			push(@tp2, $tp);
		}
	}
	if(scalar(@tp1) == 0 or scalar(@tp2) == 0){
		die("Can't calculate intervals\n");
	}
	my ($xstart, $xend) = ($preAnchor->{$preEnd^1}->{coor}->[$tp2[0]], $preAnchor->{$preEnd}->{coor}->[$tp2[0]]);
	my ($ystart, $yend) = ($preAnchor->{$preEnd^1}->{coor}->[$tp1[0]], $preAnchor->{$preEnd}->{coor}->[$tp1[0]]);
	my ($slope, $intercept) = calculateLinearParams($xstart, $xend, $ystart, $yend);
	$it1 = [sprintf("%.1f", abs($anchor1->{$end1}->{coor}->[$tp1[0]] - $yend)), sprintf("%.1f", abs($anchor1->{$end1^1}->{coor}->[$tp1[0]] - $yend))];
	$it2 = [sprintf("%.1f", abs($anchor2->{$end2}->{coor}->[$tp2[0]] * $slope + $intercept - $yend)), sprintf("%.1f", abs($anchor2->{$end2^1}->{coor}->[$tp2[0]] * $slope + $intercept - $yend))];
	return ($it1, $it2);
}

sub combineSurroundingAnchors
{
	my ($anchor_pool, $anchor1, $end1, $aux_anchor_pool, $anchor2, $end2, $lastNewAnchorId, $rev) = @_;

	my $bSuccess = 1;
	my $preAnchor = getAnchor($anchor_pool, $lastNewAnchorId);
	my ($newEnd, $newMateEnd) = ($rev, $rev^1);
	my $newAnchor;

	my ($it1, $it2);
	my $alpha;
	my $res = 0;
	while(defined $anchor1 and defined $anchor2){
		($it1, $it2) = getExtendedIntervals($preAnchor, $newMateEnd, $anchor1, $end1, $anchor2, $end2);
		$res = compareInterval($it1, $it2, 1);
		if($res < 0){
			$newAnchor = copyAnchor($anchor_pool, $anchor1, $end1^$rev);
			if(!fillMissingCoorsExtended($newAnchor, $preAnchor, $newMateEnd, $anchor2, $end2)){
				return 0;
			}
			$preAnchor->{$newMateEnd}->{neighbor} = {id=>$newAnchor->{id}, end=>$newEnd};
			$newAnchor->{$newEnd}->{neighbor} = {id=>$preAnchor->{id}, end=>$newMateEnd};
			$preAnchor = $newAnchor;
			setSynonymousPair($anchor_pool, $anchor1, $end1, $preAnchor, $newEnd);
			setSynonymousPair($anchor_pool, $anchor1, $end1^1, $preAnchor, $newMateEnd);
			($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $anchor1, $end1);
			next;
		}
		if($res > 0){
			$newAnchor = copyAnchor($anchor_pool, $anchor2, $end2^$rev);
			if(!fillMissingCoorsExtended($newAnchor, $preAnchor, $newMateEnd, $anchor1, $end1)){
				return 0;
			}
			$preAnchor->{$newMateEnd}->{neighbor} = {id=>$newAnchor->{id}, end=>$newEnd};
			$newAnchor->{$newEnd}->{neighbor} = {id=>$preAnchor->{id}, end=>$newMateEnd};
			$preAnchor = $newAnchor;
			setSynonymousPair($aux_anchor_pool, $anchor2, $end2, $preAnchor, $newEnd);
			setSynonymousPair($aux_anchor_pool, $anchor2, $end2^1, $preAnchor, $newMateEnd);
			($anchor2, $end2) = getNextAnchorEnd($aux_anchor_pool, $anchor2, $end2);
			next;
		}
		if(abs($it1->[0] - $it2->[0]) >= 1){
			if($it1->[0] < $it2->[0]){
				$alpha = ($it2->[0] - $it1->[0]) / ($it1->[1] - $it1->[0]);
				($newAnchor, $anchor1) = splitAnchorExtended($anchor_pool, $anchor1, $end1, $alpha, $rev);
				if(!fillMissingCoorsExtended($newAnchor, $preAnchor, $newMateEnd, $anchor2, $end2)){
					return 0;
				}
				$preAnchor->{$newMateEnd}->{neighbor} = {id=>$newAnchor->{id}, end=>$newEnd};
				$newAnchor->{$newEnd}->{neighbor} = {id=>$preAnchor->{id}, end=>$newMateEnd};
				$preAnchor = $newAnchor;
				setSynonymousPair($anchor_pool,$anchor1, $end1, $preAnchor, $newEnd);
			}
			else{
				$alpha = ($it1->[0] - $it2->[0]) / ($it2->[1] - $it2->[0]);
				($newAnchor, $anchor2) = splitAnchorExtended($anchor_pool, $anchor2, $end2, $alpha, $rev);
				if(!fillMissingCoorsExtended($newAnchor, $preAnchor, $newMateEnd, $anchor1, $end1)){
					return 0;
				}
				$preAnchor->{$newMateEnd}->{neighbor} = {id=>$newAnchor->{id}, end=>$newEnd};
				$newAnchor->{$newEnd}->{neighbor} = {id=>$preAnchor->{id}, end=>$newMateEnd};
				$preAnchor = $newAnchor;
				setSynonymousPair($aux_anchor_pool,$anchor2, $end2, $preAnchor, $newEnd);
			}
			($it1, $it2) = getExtendedIntervals($preAnchor, $newMateEnd, $anchor1, $end1, $anchor2, $end2);
		}
# assert: $it1->[0] == $it2->[0]
		if(abs($it1->[1] - $it2->[1]) < 1){ # to avoid rounding error
			$newAnchor = joinAnchorsExtended($anchor_pool, $anchor1, $end1, $anchor2, $end2, $rev);
			if(!defined $newAnchor){
				return 0;
			}
			$preAnchor->{$newMateEnd}->{neighbor} = {id=>$newAnchor->{id}, end=>$newEnd};
			$newAnchor->{$newEnd}->{neighbor} = {id=>$preAnchor->{id}, end=>$newMateEnd};
			$preAnchor = $newAnchor;
			setSynonymousPair($anchor_pool,$anchor1, $end1, $preAnchor, $newEnd);
			setSynonymousPair($anchor_pool,$anchor1, $end1^1, $preAnchor, $newMateEnd);
			setSynonymousPair($aux_anchor_pool,$anchor2, $end2, $preAnchor, $newEnd);
			setSynonymousPair($aux_anchor_pool,$anchor2, $end2^1, $preAnchor, $newMateEnd);
			($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $anchor1, $end1);
			($anchor2, $end2) = getNextAnchorEnd($aux_anchor_pool, $anchor2, $end2);
			next;
		}
		if($it1->[1] < $it2->[1]){
			$alpha = ($it1->[1] - $it2->[0]) / ($it2->[1] - $it2->[0]);
			($newAnchor, $anchor2) = splitAnchorExtended(undef, $anchor2, $end2, $alpha, $rev);
			$newAnchor = joinAnchorsExtended($anchor_pool, $anchor1, $end1, $newAnchor, $newEnd, $rev);
			if(!defined $newAnchor){
				return 0;
			}
			$preAnchor->{$newMateEnd}->{neighbor} = {id=>$newAnchor->{id}, end=>$newEnd};
			$newAnchor->{$newEnd}->{neighbor} = {id=>$preAnchor->{id}, end=>$newMateEnd};
			$preAnchor = $newAnchor;
			setSynonymousPair($anchor_pool,$anchor1, $end1, $preAnchor, $newEnd);
			setSynonymousPair($anchor_pool,$anchor1, $end1^1, $preAnchor, $newMateEnd);
			setSynonymousPair($aux_anchor_pool,$anchor2, $end2, $preAnchor, $newEnd);
			($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $anchor1, $end1);
		}
		else{
			$alpha = ($it2->[1] - $it1->[0]) / ($it1->[1] - $it1->[0]);
			($newAnchor, $anchor1) = splitAnchorExtended(undef, $anchor1, $end1, $alpha, $rev);
			$newAnchor = joinAnchorsExtended($anchor_pool, $newAnchor, $newEnd, $anchor2, $end2, $rev);
			if(!defined $newAnchor){
				return 0;
			}
			$preAnchor->{$newMateEnd}->{neighbor} = {id=>$newAnchor->{id}, end=>$newEnd};
			$newAnchor->{$newEnd}->{neighbor} = {id=>$preAnchor->{id}, end=>$newMateEnd};
			$preAnchor = $newAnchor;
			setSynonymousPair($anchor_pool,$anchor1, $end1, $preAnchor, $newEnd);
			setSynonymousPair($aux_anchor_pool,$anchor2, $end2, $preAnchor, $newEnd);
			setSynonymousPair($aux_anchor_pool,$anchor2, $end2^1, $preAnchor, $newMateEnd);
			($anchor2, $end2) = getNextAnchorEnd($aux_anchor_pool, $anchor2, $end2);
		}
	}
	while(defined $anchor1){
		$newAnchor = copyAnchor($anchor_pool, $anchor1, $end1^$rev);
		$preAnchor->{$newMateEnd}->{neighbor} = {id=>$newAnchor->{id}, end=>$newEnd};
		$newAnchor->{$newEnd}->{neighbor} = {id=>$preAnchor->{id}, end=>$newMateEnd};
		adjustCoors($preAnchor, $newMateEnd, $newAnchor, $newEnd);
		$preAnchor = $newAnchor;
		setSynonymousPair($anchor_pool, $anchor1, $end1, $preAnchor, $newEnd);
		setSynonymousPair($anchor_pool, $anchor1, $end1^1, $preAnchor, $newMateEnd);
		($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $anchor1, $end1);
	}
	while(defined $anchor2){
		$newAnchor = copyAnchor($anchor_pool, $anchor2, $end2^$rev);
		$preAnchor->{$newMateEnd}->{neighbor} = {id=>$newAnchor->{id}, end=>$newEnd};
		$newAnchor->{$newEnd}->{neighbor} = {id=>$preAnchor->{id}, end=>$newMateEnd};
		$preAnchor = $newAnchor;
		setSynonymousPair($aux_anchor_pool, $anchor2, $end2, $preAnchor, $newEnd);
		setSynonymousPair($aux_anchor_pool, $anchor2, $end2^1, $preAnchor, $newMateEnd);
		($anchor2, $end2) = getNextAnchorEnd($aux_anchor_pool, $anchor2, $end2);
	}
	return 1;
}

sub adjustCoors
{
	my ($preAnchor, $preEnd, $anchor, $end) = @_;
	for(my $tp=$#{$preAnchor->{nodeId}}; $tp>=1; $tp--){
		next unless(defined $preAnchor->{nodeId}->[$tp] and defined $anchor->{nodeId}->[$tp]);
		next unless($preAnchor->{nodeId}->[$tp] == $anchor->{nodeId}->[$tp]);
		if( ($preAnchor->{$preEnd}->{coor}->[$tp] - $preAnchor->{$preEnd^1}->{coor}->[$tp]) *
			($anchor->{$end}->{coor}->[$tp] - $preAnchor->{$preEnd}->{coor}->[$tp]) < 0 ){
			if( ($preAnchor->{$preEnd}->{coor}->[$tp] - $preAnchor->{$preEnd^1}->{coor}->[$tp]) *
				($anchor->{$end^1}->{coor}->[$tp] - $preAnchor->{$preEnd}->{coor}->[$tp]) < 0 ){
				$preAnchor->{$preEnd}->{coor}->[$tp] = $anchor->{$end}->{coor}->[$tp];
			}
			else{
				my $mean = sprintf("%.1f", ($preAnchor->{$preEnd}->{coor}->[$tp] + $anchor->{$end}->{coor}->[$tp]) / 2);
				$anchor->{$end}->{coor}->[$tp] = $preAnchor->{$preEnd}->{coor}->[$tp] = $mean;
			}
		}
	}
}

sub mergeClusters
{
	my ($cluster_pool, $cluster1, $cluster2, $edge, $type1, $type2) = @_;
	my ($nodeId1, $nodeId2) = ($edge->{id1}, $edge->{id2});
	my ($loc1, $loc2) = ($cluster1->{locations}->{$nodeId1}, $cluster2->{locations}->{$nodeId2});

	# form a new cluster with node1 and node2
	my $cluster = newCluster($cluster_pool);
	my $anchor_pool = $cluster->{anchor_pool};
	# initialize homo related fields
	$anchor_pool->{homo} = {};
	my ($path1, $path2) = setupPaths($anchor_pool, $edge->{glues}, $nodeId1, $nodeId2, $type1, $type2);
	$cluster->{locations}->{$nodeId1} = $path1;
	$cluster->{locations}->{$nodeId2} = $path2;

	my $anchor_pool1 = $cluster1->{anchor_pool};
	$anchor_pool1->{homo} = {};
	if(expandPath($anchor_pool1, $loc1, $type1, $cluster_pool->{graph})){
		collectHomoPairs($cluster_pool->{homo}, $anchor_pool1->{homo});
		deleteLocations($cluster1->{locations}, [keys $anchor_pool1->{homo}]);
	}
	my $anchor_pool2 = $cluster2->{anchor_pool};
	$anchor_pool2->{homo} = {};
	if(expandPath($anchor_pool2, $loc2, $type2, $cluster_pool->{graph})){
		collectHomoPairs($cluster_pool->{homo}, $anchor_pool2->{homo});
		deleteLocations($cluster2->{locations}, [keys $anchor_pool2->{homo}]);
	}

	clearSynonymousPairs($anchor_pool);
	clearSynonymousPairs($anchor_pool1);
	if(!mergeOuterPaths($anchor_pool, $path1, $anchor_pool1, $loc1)){
		freeCluster($cluster_pool, $cluster);
		return undef;
	}
	updateLocations($cluster->{locations}, $anchor_pool->{synonym}, $anchor_pool->{protected_id});
	my %locations = %{dclone($cluster1->{locations})};
	updateLocations(\%locations, $anchor_pool1->{synonym}, $anchor_pool->{protected_id});
	addLocations($cluster, \%locations, $anchor_pool->{protected_id});

	expandPath($anchor_pool, $path1, $type1, $cluster_pool->{graph});
	expandPath($anchor_pool, $path2, $type2, $cluster_pool->{graph});
	if(scalar(keys %{$anchor_pool->{homo}}) > 0){
		deleteLocations($cluster->{locations}, [keys $anchor_pool->{homo}]);
	}

	clearSynonymousPairs($anchor_pool);
	clearSynonymousPairs($anchor_pool2);
	if(!mergeOuterPaths($anchor_pool, $path2, $anchor_pool2, $loc2)){
		freeCluster($cluster_pool, $cluster);
		return undef;
	}
	updateLocations($cluster->{locations}, $anchor_pool->{synonym}, $anchor_pool->{protected_id});
	%locations = %{dclone($cluster2->{locations})};
	updateLocations(\%locations, $anchor_pool2->{synonym}, $anchor_pool->{protected_id});
	addLocations($cluster, \%locations, $anchor_pool->{protected_id});

	freeCluster($cluster_pool, $cluster1);
	freeCluster($cluster_pool, $cluster2);

	collectHomoPairs($cluster_pool->{homo}, $anchor_pool->{homo});
	return $cluster;
}

sub insertAnchors
{
	my ($cluster_pool, $cluster, $edge, $type1, $type2) = @_;

	my $anchor_pool = $cluster->{anchor_pool};
	my ($nodeId1, $nodeId2) = ($edge->{id1}, $edge->{id2});
	my ($path1, $path2) = setupPaths($anchor_pool, $edge->{glues}, $nodeId1, $nodeId2, $type1, $type2);

	clearSynonymousPairs($anchor_pool);
	my ($loc1, $loc2) = ($cluster->{locations}->{$nodeId1}, $cluster->{locations}->{$nodeId2});
	if(defined $loc1){
		if(defined $loc2){
			# both nodes exist in the same cluster
			return 0;
		}
		# node1 exists in the cluster
		# add node2 to the cluster
		$anchor_pool->{homo} = {};
		if(expandPath($anchor_pool, $loc1, $type1, $cluster_pool->{graph})){
			collectHomoPairs($cluster_pool->{homo}, $anchor_pool->{homo});
			deleteLocations($cluster->{locations}, [keys $anchor_pool->{homo}]);
		}
		if(mergeInnerPaths($anchor_pool, $loc1, $path2)){
			$cluster->{locations}->{$nodeId2} = $path2;
			updateLocations($cluster->{locations}, $anchor_pool->{synonym}, $anchor_pool->{protected_id});
			return 1;
		}
		return 0;
	}
	if(defined $loc2){
		# node2 exists in the cluster
		# add node1 to the cluster
		$anchor_pool->{homo} = {};
		if(expandPath($anchor_pool, $loc2, $type2, $cluster_pool->{graph})){
			collectHomoPairs($cluster_pool->{homo}, $anchor_pool->{homo});
			deleteLocations($cluster->{locations}, [keys $anchor_pool->{homo}]);
		}
		if(mergeInnerPaths($anchor_pool, $loc2, $path1)){
			$cluster->{locations}->{$nodeId1} = $path1;
			updateLocations($cluster->{locations}, $anchor_pool->{synonym}, $anchor_pool->{protected_id});
			return 1;
		}
		return 0;
	}
	# form a new cluster with node1 and node2
	$cluster->{locations}->{$nodeId1} = $path1;
	$cluster->{locations}->{$nodeId2} = $path2;
	return 1;
}

sub mergeOuterPaths
{
	my ($anchor_pool, $path1, $aux_anchor_pool, $path2) = @_;
	my $type = $path1->{type}; # $path1->{type} == $path2->{type}

	my ($firstAnc1, $firstEd1) = getAnchorEnd($anchor_pool, $path1->{ends}->{0});
	my ($lastAnc1, $lastEd1) = getAnchorEnd($anchor_pool, $path1->{ends}->{1});
	my ($firstAnc2, $firstEd2) = getAnchorEnd($aux_anchor_pool, $path2->{ends}->{0});
	my ($lastAnc2, $lastEd2) = getAnchorEnd($aux_anchor_pool, $path2->{ends}->{1});

	my $nodeId = $firstAnc1->{nodeId}->[$type];
	$anchor_pool->{protected_id} = $nodeId; # whose location will be calculated separately

	my $dir = ($firstAnc1->{$firstEd1}->{coor}->[$type] < $lastAnc1->{$lastEd1}->{coor}->[$type]) ? 1 : -1;
	my $dir2 = ($firstAnc2->{$firstEd2}->{coor}->[$type] < $lastAnc2->{$lastEd2}->{coor}->[$type]) ? 1 : -1;
	if($dir * $dir2 < 0){
		my ($tmpAnc, $tmpEnd);
		($tmpAnc, $tmpEnd) = ($firstAnc2, $firstEd2);
		($firstAnc2, $firstEd2) = ($lastAnc2, $lastEd2);
		($lastAnc2, $lastEd2) = ($tmpAnc, $tmpEnd);
	}
	my $warning = 0;
	if( checkEndConflict($firstAnc1, $firstEd1, $firstAnc2, $firstEd2, $type, $dir) or
		checkEndConflict($lastAnc1, $lastEd1, $lastAnc2, $lastEd2, $type, $dir * -1) ){
		$warning = 1;
	}

	my ($neighbor_anchor1, $neighbor_end1) = getNextAnchorEnd($anchor_pool, $firstAnc1, $firstEd1^1);
	my ($neighbor_anchor2, $neighbor_end2) = getNextAnchorEnd($aux_anchor_pool, $firstAnc2, $firstEd2^1);

	my ($anchor1, $end1) = ($firstAnc1, $firstEd1);
	my ($anchor2, $end2) = ($firstAnc2, $firstEd2);
	my ($lastAnchor1, $lastEnd1);
	my ($lastAnchor2, $lastEnd2);
	my ($it1, $it2);
	my $alpha;

	my $firstAnchorId = getNewAnchorId($anchor_pool);
	my ($preAnchor, $newAnchor);
	my ($bSuccess, $res) = (1, 0);
	while( defined $anchor1 and defined $anchor2 ){
		$it1 = getAnchorInterval($anchor1, $end1, $type) if($res <= 0);
		$it2 = getAnchorInterval($anchor2, $end2, $type) if($res >= 0);
		$res = compareInterval($it1, $it2, $dir);
		if($res < 0){
			$newAnchor = copyAnchor($anchor_pool, $anchor1, $end1);
			if(defined $lastAnchor2){
				if(!fillMissingCoors($newAnchor, $preAnchor, 1, $anchor2, $end2, $type)){
					$bSuccess = 0;
					last;
				}
			}
			elsif(defined $neighbor_anchor2){
				if(!fillMissingCoors($newAnchor, $neighbor_anchor2, $neighbor_end2, $anchor2, $end2, $type)){
					$bSuccess = 0;
					last;
				}
			}
			if(defined $preAnchor){
				$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
				$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
			}
			$preAnchor = $newAnchor;
			setSynonymousPair($anchor_pool, $anchor1, $end1, $preAnchor, 0);
			setSynonymousPair($anchor_pool, $anchor1, $end1^1, $preAnchor, 1);
			($lastAnchor1, $lastEnd1) = ($anchor1, $end1);
			($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $anchor1, $end1, $type);
			next;
		}
		if($res > 0){
			$newAnchor = copyAnchor($anchor_pool, $anchor2, $end2);
			if(defined $lastAnchor1){
				if(!fillMissingCoors($newAnchor, $preAnchor, 1, $anchor1, $end1, $type)){
					$bSuccess = 0;
					last;
				}
			}
			elsif(defined $neighbor_anchor1){
				if(!fillMissingCoors($newAnchor, $neighbor_anchor1, $neighbor_end1, $anchor1, $end1, $type)){
					$bSuccess = 0;
					last;
				}
			}
			if(defined $preAnchor){
				$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
				$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
			}
			$preAnchor = $newAnchor;
			setSynonymousPair($aux_anchor_pool, $anchor2, $end2, $preAnchor, 0);
			setSynonymousPair($aux_anchor_pool, $anchor2, $end2^1, $preAnchor, 1);
			($lastAnchor2, $lastEnd2) = ($anchor2, $end2);
			($anchor2, $end2) = getNextAnchorEnd($aux_anchor_pool, $anchor2, $end2, $type);
			next;
		}
		if($it1->[0] != $it2->[0]){
			if( ($it1->[0] - $it2->[0]) * $dir < 0 ){
				$alpha = ($it2->[0] - $it1->[0]) / ($it1->[1] - $it1->[0]);
				($newAnchor, $anchor1) = splitAnchor($anchor_pool, $anchor1, $end1, $alpha);
				if(defined $lastAnchor2){
					if(!fillMissingCoors($newAnchor, $preAnchor, 1, $anchor2, $end2, $type)){
						$bSuccess = 0;
						last;
					}
				}
				elsif(defined $neighbor_anchor2){
					if(!fillMissingCoors($newAnchor, $neighbor_anchor2, $neighbor_end2, $anchor2, $end2, $type)){
						$bSuccess = 0;
						last;
					}
				}
				if(defined $preAnchor){
					$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
					$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
				}
				$preAnchor = $newAnchor;
				setSynonymousPair($anchor_pool, $anchor1, $end1, $preAnchor, 0);
				$it1 = getAnchorInterval($anchor1, $end1, $type);
			}
			else{
				$alpha = ($it1->[0] - $it2->[0]) / ($it2->[1] - $it2->[0]);
				($newAnchor, $anchor2) = splitAnchor($anchor_pool, $anchor2, $end2, $alpha);
				if(defined $lastAnchor1){
					if(!fillMissingCoors($newAnchor, $preAnchor, 1, $anchor1, $end1, $type)){
						$bSuccess = 0;
						last;
					}
				}
				elsif(defined $neighbor_anchor1){
					if(!fillMissingCoors($newAnchor, $neighbor_anchor1, $neighbor_end1, $anchor1, $end1, $type)){
						$bSuccess = 0;
						last;
					}
				}
				if(defined $preAnchor){
					$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
					$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
				}
				$preAnchor = $newAnchor;
				setSynonymousPair($aux_anchor_pool, $anchor2, $end2, $preAnchor, 0);
				$it2 = getAnchorInterval($anchor2, $end2, $type);
			}
		}
# assert: $it1->[0] == $it2->[0]
		if($it1->[1] == $it2->[1]){
			$newAnchor = joinAnchors($anchor_pool, $anchor1, $end1, $anchor2, $end2, $type);
			if(!defined $newAnchor){
				$bSuccess = 0;
				last;
			}
			if(defined $preAnchor){
				$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
				$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
			}
			$preAnchor = $newAnchor;
			setSynonymousPair($anchor_pool, $anchor1, $end1, $preAnchor, 0);
			setSynonymousPair($anchor_pool, $anchor1, $end1^1, $preAnchor, 1);
			setSynonymousPair($aux_anchor_pool, $anchor2, $end2, $preAnchor, 0);
			setSynonymousPair($aux_anchor_pool, $anchor2, $end2^1, $preAnchor, 1);
			($lastAnchor1, $lastEnd1) = ($anchor1, $end1);
			($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $anchor1, $end1, $type);
			($lastAnchor2, $lastEnd2) = ($anchor2, $end2);
			($anchor2, $end2) = getNextAnchorEnd($aux_anchor_pool, $anchor2, $end2, $type);
			next;
		}
		if( ($it1->[1] - $it2->[1]) * $dir < 0 ){
			$alpha = ($it1->[1] - $it2->[0]) / ($it2->[1] - $it2->[0]);
			($newAnchor, $anchor2) = splitAnchor(undef, $anchor2, $end2, $alpha);
			$newAnchor = joinAnchors($anchor_pool, $anchor1, $end1, $newAnchor, 0, $type);
			if(!defined $newAnchor){
				$bSuccess = 0;
				last;
			}
			if(defined $preAnchor){
				$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
				$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
			}
			$preAnchor = $newAnchor;
			setSynonymousPair($anchor_pool, $anchor1, $end1, $preAnchor, 0);
			setSynonymousPair($anchor_pool, $anchor1, $end1^1, $preAnchor, 1);
			setSynonymousPair($aux_anchor_pool, $anchor2, $end2, $preAnchor, 0);
			($lastAnchor1, $lastEnd1) = ($anchor1, $end1);
			($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $anchor1, $end1, $type);
		}
		else{
			$alpha = ($it2->[1] - $it1->[0]) / ($it1->[1] - $it1->[0]);
			($newAnchor, $anchor1) = splitAnchor(undef, $anchor1, $end1, $alpha);
			$newAnchor = joinAnchors($anchor_pool, $newAnchor, 0, $anchor2, $end2, $type);
			if(!defined $newAnchor){
				$bSuccess = 0;
				last;
			}
			if(defined $preAnchor){
				$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
				$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
			}
			$preAnchor = $newAnchor;
			setSynonymousPair($anchor_pool, $anchor1, $end1, $preAnchor, 0);
			setSynonymousPair($aux_anchor_pool, $anchor2, $end2, $preAnchor, 0);
			setSynonymousPair($aux_anchor_pool, $anchor2, $end2^1, $preAnchor, 1);
			($lastAnchor2, $lastEnd2) = ($anchor2, $end2);
			($anchor2, $end2) = getNextAnchorEnd($aux_anchor_pool, $anchor2, $end2, $type);
		}
	}

# append succeeding anchors not belonging to updated path1
	($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $lastAnchor1, $lastEnd1) if(!defined $anchor1); # reach the end
	($anchor2, $end2) = getNextAnchorEnd($aux_anchor_pool, $lastAnchor2, $lastEnd2) if(!defined $anchor2); # reach the end
	$bSuccess &&= combineSurroundingAnchors($anchor_pool, $anchor1, $end1, $aux_anchor_pool, $anchor2, $end2, $preAnchor->{id}, 0);
	($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $firstAnc1, $firstEd1^1);
	($anchor2, $end2) = getNextAnchorEnd($aux_anchor_pool, $firstAnc2, $firstEd2^1);
	$bSuccess &&= combineSurroundingAnchors($anchor_pool, $anchor1, $end1, $aux_anchor_pool, $anchor2, $end2, $firstAnchorId, 1);
	if(!$bSuccess){
# free the new anchors
		freeAnchorsFromOne($anchor_pool, $preAnchor);
		return 0;
	}

# update the path1 range
	$path1->{ends} = updatePathEnds($anchor_pool, $firstAnchorId, 0, $preAnchor->{id}, 1, $type);

# free the used anchors
	freeAnchorsFromOne($anchor_pool, $firstAnc1);
	return 1;
}

sub mergeInnerPaths
{
	my ($anchor_pool, $path1, $path2) = @_;
	my ($type, $type_check) = ($path1->{type}, $path2->{type});

	my ($firstAnc1, $firstEd1) = getAnchorEnd($anchor_pool, $path1->{ends}->{0});
	my ($lastAnc1, $lastEd1) = getAnchorEnd($anchor_pool, $path1->{ends}->{1});
	my ($firstAnc2, $firstEd2) = getAnchorEnd($anchor_pool, $path2->{ends}->{0});
	my ($lastAnc2, $lastEd2) = getAnchorEnd($anchor_pool, $path2->{ends}->{1});

	my $nodeId = $firstAnc1->{nodeId}->[$type];
	$anchor_pool->{protected_id} = $nodeId; # whose location will be calculated separately

	my $dir = ($firstAnc1->{$firstEd1}->{coor}->[$type] < $lastAnc1->{$lastEd1}->{coor}->[$type]) ? 1 : -1;
	my $dir2 = ($firstAnc2->{$firstEd2}->{coor}->[$type] < $lastAnc2->{$lastEd2}->{coor}->[$type]) ? 1 : -1;
	if($dir * $dir2 < 0){
		my ($tmpAnc, $tmpEnd);
		($tmpAnc, $tmpEnd) = ($firstAnc2, $firstEd2);
		($firstAnc2, $firstEd2) = ($lastAnc2, $lastEd2);
		($lastAnc2, $lastEd2) = ($tmpAnc, $tmpEnd);
	}

	my ($neighbor_anchor, $neighbor_end) = getNextAnchorEnd($anchor_pool, $firstAnc1, $firstEd1^1);
	my ($anchor1, $end1) = ($firstAnc1, $firstEd1);
	my ($anchor2, $end2) = ($firstAnc2, $firstEd2);
	my ($lastAnchor1, $lastEnd1);
	my ($lastAnchor2, $lastEnd2);
	my ($it1, $it2);
	my $alpha;

	my $firstAnchorId = getNewAnchorId($anchor_pool);
	my ($preAnchor, $newAnchor);
	my ($bSuccess, $res) = (1, 0);
	while( defined $anchor1 and defined $anchor2 ){
		$it1 = getAnchorInterval($anchor1, $end1, $type) if($res <= 0);
		$it2 = getAnchorInterval($anchor2, $end2, $type) if($res >= 0);
		$res = compareInterval($it1, $it2, $dir);
		if($res < 0){
			$newAnchor = copyAnchor($anchor_pool, $anchor1, $end1);
			if(defined $lastAnchor2){
				if(defined $anchor1->{nodeId}->[$type_check]){
					$bSuccess = 0;
					last;
				}
				fillMissingCoors($newAnchor, $preAnchor, 1, $anchor2, $end2, $type);
			}
			if(defined $preAnchor){
				$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
				$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
			}
			$preAnchor = $newAnchor;
			setSynonymousPair($anchor_pool, $anchor1, $end1, $preAnchor, 0);
			setSynonymousPair($anchor_pool, $anchor1, $end1^1, $preAnchor, 1);
			($lastAnchor1, $lastEnd1) = ($anchor1, $end1);
			($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $anchor1, $end1, $type);
			next;
		}
		if($res > 0){
			$newAnchor = copyAnchor($anchor_pool, $anchor2, $end2);
			if(defined $lastAnchor1){
				if(!fillMissingCoors($newAnchor, $preAnchor, 1, $anchor1, $end1, $type)){
					$bSuccess = 0;
					last;
				}
			}
			elsif(defined $neighbor_anchor){
				if(!fillMissingCoors($newAnchor, $neighbor_anchor, $neighbor_end, $anchor1, $end1, $type)){
					$bSuccess = 0;
					last;
				}
			}
			if(defined $preAnchor){
				$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
				$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
			}
			$preAnchor = $newAnchor;
			setSynonymousPair($anchor_pool, $anchor2, $end2, $preAnchor, 0);
			setSynonymousPair($anchor_pool, $anchor2, $end2^1, $preAnchor, 1);
			($lastAnchor2, $lastEnd2) = ($anchor2, $end2);
			($anchor2, $end2) = getNextAnchorEnd($anchor_pool, $anchor2, $end2, $type);
			next;
		}
		if(defined $anchor1->{nodeId}->[$type_check]){
			$bSuccess = 0; # by default it fails as conflict found
			if(!defined $anchor2->{$end2^1}->{neighbor} and defined $preAnchor){ # if it is the last anchor and there's already something merged
				# then small overlap is allowed for boundary anchor
				my $overlap_len = (($it1->[0] - $it2->[0]) * $dir < 0) ? abs($it2->[1] - $it2->[0]) : abs($it2->[1] - $it1->[0]);
				if($overlap_len < $anchor1->{len} * 0.01 and $overlap_len < $anchor2->{len} * 0.01){
					$newAnchor = copyAnchor($anchor_pool, $anchor2, $end2);
					$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
					$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
					$preAnchor = $newAnchor;
					setSynonymousPair($anchor_pool, $anchor2, $end2, $preAnchor, 0);
					setSynonymousPair($anchor_pool, $anchor2, $end2^1, $preAnchor, 1);
					($lastAnchor2, $lastEnd2) = ($anchor2, $end2);
					($anchor2, $end2) = getNextAnchorEnd($anchor_pool, $anchor2, $end2, $type);
					$bSuccess = 1;
				}
			}
			last;
		}
		if($it1->[0] != $it2->[0]){
			if( ($it1->[0] - $it2->[0]) * $dir < 0 ){
				$alpha = ($it2->[0] - $it1->[0]) / ($it1->[1] - $it1->[0]);
				($newAnchor, $anchor1) = splitAnchor($anchor_pool, $anchor1, $end1, $alpha);
				if(defined $lastAnchor2){
					fillMissingCoors($newAnchor, $preAnchor, 1, $anchor2, $end2, $type);
				}
				if(defined $preAnchor){
					$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
					$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
				}
				$preAnchor = $newAnchor;
				setSynonymousPair($anchor_pool, $anchor1, $end1, $preAnchor, 0);
				$it1 = getAnchorInterval($anchor1, $end1, $type);
			}
			else{
				$alpha = ($it1->[0] - $it2->[0]) / ($it2->[1] - $it2->[0]);
				($newAnchor, $anchor2) = splitAnchor($anchor_pool, $anchor2, $end2, $alpha);
				if(defined $lastAnchor1){
					fillMissingCoors($newAnchor, $preAnchor, 1, $anchor1, $end1, $type);
				}
				elsif(defined $neighbor_anchor){
					fillMissingCoors($newAnchor, $neighbor_anchor, $neighbor_end, $anchor1, $end1, $type);
				}
				if(defined $preAnchor){
					$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
					$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
				}
				$preAnchor = $newAnchor;
				setSynonymousPair($anchor_pool, $anchor2, $end2, $preAnchor, 0);
				$it2 = getAnchorInterval($anchor2, $end2, $type);
			}
		}
# assert: $it1->[0] == $it2->[0]
		if($it1->[1] == $it2->[1]){
			$newAnchor = joinAnchors($anchor_pool, $anchor1, $end1, $anchor2, $end2, $type);
			if(defined $preAnchor){
				$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
				$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
			}
			$preAnchor = $newAnchor;
			setSynonymousPair($anchor_pool, $anchor1, $end1, $preAnchor, 0);
			setSynonymousPair($anchor_pool, $anchor1, $end1^1, $preAnchor, 1);
			setSynonymousPair($anchor_pool, $anchor2, $end2, $preAnchor, 0);
			setSynonymousPair($anchor_pool, $anchor2, $end2^1, $preAnchor, 1);
			($lastAnchor1, $lastEnd1) = ($anchor1, $end1);
			($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $anchor1, $end1, $type);
			($lastAnchor2, $lastEnd2) = ($anchor2, $end2);
			($anchor2, $end2) = getNextAnchorEnd($anchor_pool, $anchor2, $end2, $type);
			next;
		}
		if( ($it1->[1] - $it2->[1]) * $dir < 0 ){
			$alpha = ($it1->[1] - $it2->[0]) / ($it2->[1] - $it2->[0]);
			($newAnchor, $anchor2) = splitAnchor(undef, $anchor2, $end2, $alpha);
			$newAnchor = joinAnchors($anchor_pool, $anchor1, $end1, $newAnchor, 0, $type);
			if(defined $preAnchor){
				$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
				$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
			}
			$preAnchor = $newAnchor;
			setSynonymousPair($anchor_pool, $anchor1, $end1, $preAnchor, 0);
			setSynonymousPair($anchor_pool, $anchor1, $end1^1, $preAnchor, 1);
			setSynonymousPair($anchor_pool, $anchor2, $end2, $preAnchor, 0);
			($lastAnchor1, $lastEnd1) = ($anchor1, $end1);
			($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $anchor1, $end1, $type);
		}
		else{
			$alpha = ($it2->[1] - $it1->[0]) / ($it1->[1] - $it1->[0]);
			($newAnchor, $anchor1) = splitAnchor(undef, $anchor1, $end1, $alpha);
			$newAnchor = joinAnchors($anchor_pool, $newAnchor, 0, $anchor2, $end2, $type);
			if(defined $preAnchor){
				$preAnchor->{1}->{neighbor} = {id=>$newAnchor->{id}, end=>0};
				$newAnchor->{0}->{neighbor} = {id=>$preAnchor->{id}, end=>1};
			}
			$preAnchor = $newAnchor;
			setSynonymousPair($anchor_pool, $anchor1, $end1, $preAnchor, 0);
			setSynonymousPair($anchor_pool, $anchor2, $end2, $preAnchor, 0);
			setSynonymousPair($anchor_pool, $anchor2, $end2^1, $preAnchor, 1);
			($lastAnchor2, $lastEnd2) = ($anchor2, $end2);
			($anchor2, $end2) = getNextAnchorEnd($anchor_pool, $anchor2, $end2, $type);
		}
	}

# append succeeding anchors not belonging to updated path1
	($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $lastAnchor1, $lastEnd1) if(!defined $anchor1); # reach the end
	($anchor2, $end2) = getNextAnchorEnd($anchor_pool, $lastAnchor2, $lastEnd2) if(!defined $anchor2); # reach the end
	$bSuccess &&= combineSurroundingAnchors($anchor_pool, $anchor1, $end1, $anchor_pool, $anchor2, $end2, $preAnchor->{id}, 0);
	($anchor1, $end1) = getNextAnchorEnd($anchor_pool, $firstAnc1, $firstEd1^1);
	($anchor2, $end2) = getNextAnchorEnd($anchor_pool, $firstAnc2, $firstEd2^1);
	$bSuccess &&= combineSurroundingAnchors($anchor_pool, $anchor1, $end1, $anchor_pool, $anchor2, $end2, $firstAnchorId, 1);
	if(!$bSuccess){
# free the new anchors and the temporary anchors
		freeAnchorsFromOne($anchor_pool, $preAnchor);
		freeAnchorsFromOne($anchor_pool, $firstAnc2);
		return 0;
	}
# update the path1 range
	$path1->{ends} = updatePathEnds($anchor_pool, $firstAnchorId, 0, $preAnchor->{id}, 1, $type);

# free the used anchors
	freeAnchorsFromOne($anchor_pool, $firstAnc1);
	freeAnchorsFromOne($anchor_pool, $firstAnc2);

	return 1;
}

sub compareInterval
{
	my ($it1, $it2, $dir) = @_;
	if( ($it1->[1] - $it2->[0]) * $dir <= 0 ){
		return -1;
	}
	if( ($it1->[0] - $it2->[1]) * $dir >= 0 ){
		return 1;
	}
	return 0;
}

### LAYOUT function ###
sub makeLayouts
{
	my ($nodes, $cluster_pool) = @_;
	my @layouts = ();
	foreach my $cluster (@{$cluster_pool->{clusters}}){
		my $layout = &makeLayout($nodes, $cluster);
		push(@layouts, $layout);
	}

	return (\@layouts, ["node_id", "ori_id", "type", "start", "end", "dir", "xstart", "xend", "xlength"]);
}

sub makeLayout
{
	my ($nodes, $cluster) = @_;
	my @layout = ();

	my $anchor_pool = $cluster->{anchor_pool};

# calculate the offsets of all the anchor ends
	my $offset = 0;
	my ($nextAnchor, $nextEnd);
	my ($anchor, $end) = (getAnchor($anchor_pool, $cluster->{ends}->{0}->{id}), $cluster->{ends}->{0}->{end});
	while(defined $anchor){
		$anchor->{$end}->{offset} = $offset;
		$offset += $anchor->{len};
		$anchor->{$end^1}->{offset} = $offset;
		($nextAnchor, $nextEnd) = getNextAnchorEnd($anchor_pool, $anchor, $end);
		if(defined $nextAnchor){
			$offset += calculateGap($anchor_pool, $anchor, $end, $nextAnchor, $nextEnd);
		}
		($anchor, $end) = ($nextAnchor, $nextEnd);
	}
# for each node in the cluster, tentatively calculate its layout information
	my $locations = $cluster->{locations};
	my ($nodeId, $loc);
	my ($type, $dir, $xstart, $xend, $ystart, $yend);
	while (($nodeId, $loc) = each (%{$locations}) ){
		my $node = $nodes->[$nodeId];
		$type = $node->{type};
		my ($anchor1, $end1) = getAnchorEnd($anchor_pool, $loc->{ends}->{0});
		my ($anchor2, $end2) = getAnchorEnd($anchor_pool, $loc->{ends}->{1});
		my ($offset1, $offset2) = ($anchor1->{$end1}->{offset}, $anchor2->{$end2}->{offset});
		my ($coor1, $coor2) = ($anchor1->{$end1}->{coor}->[$type], $anchor2->{$end2}->{coor}->[$type]);
		$dir = ($coor1 - $coor2) * ($offset1 - $offset2) > 0 ? "+" : "-";
		if($offset1 < $offset2){
			push(@layout, {node_id=>$nodeId, ori_id=>$node->{ori_id}, type=>$type, dir=>$dir, start=>$offset1, end=>$offset2, xstart=>$coor1, xend=>$coor2, xlength=>$node->{length}});
		}
		else{
			push(@layout, {node_id=>$nodeId, ori_id=>$node->{ori_id}, type=>$type, dir=>$dir, start=>$offset2, end=>$offset1, xstart=>$coor2, xend=>$coor1, xlength=>$node->{length}});
		}
	}
	@layout = sort{ $a->{start} <=> $b->{start} } @layout;
	return \@layout;
}

sub outputLayouts
{
	my ($paths, $headers, $outfile) = @_;
	my $out;
	if(defined $outfile){
		open(OUT, ">$outfile") or die("Can not open \"$outfile\" for writing\n");
		$out = \*OUT;
	}
	else{
		$out = \*STDOUT;
	}

	print $out "path_id\t" . join("\t", @{$headers}) . "\n";
	my ($i, $j, $k);
	my $path;
	for($i=0; $i<scalar(@{$paths}); $i++){
		$path = $paths->[$i];
		for($j=0; $j<scalar(@{$path}); $j++){
			print $out ($i+1);
			for($k=0; $k<scalar(@{$headers}); $k++){
				print $out "\t". $path->[$j]->{$headers->[$k]};
			}
			print $out "\n";
		}
	}
	close $out;
}

sub getIndices
{
	my ($fields, $mustFields) = @_;

	my %indices = ();
	for(my $i=0; $i<scalar(@$fields); $i++){
		$indices{$fields->[$i]} = $i;
	}
	my @subIndices = ();
	foreach my $field (@{$mustFields}){
		if(defined $indices{$field}){
			push(@subIndices, $indices{$field});
		}
	}

	return @subIndices;
}

sub globalAlign
{
	my ($positions1, $positions2, $threshold) = @_;
	my $alignment = [];
	my ($m, $n) = (scalar(@{$positions1}), scalar(@{$positions2}));
	if( ($m <= 3) or ($n <= 3) ){
		return [];
	}
	my ($aligned, $total) = (0, 0);
	if($m <= $n){
		my ($base, $size) = (0, ($m > 50) ? 50 : $m);
		my ($base2, $size2) = (0, ($n > 150) ? 150 : $n);
		my $left = $m;
		while($left > 0 and $base2 < $n){
			my $algn = alignLCS([@{$positions1}[$base..($base+$size-1)]], [@{$positions2}[$base2..($base2+$size2-1)]]);

			$left -= $size;
			$total += $size;
			if(scalar(@{$algn}) / $size >= $threshold){
				foreach (@{$algn}){
					push(@{$alignment}, [$_->[0]+$base, $_->[1]+$base2]);
				}
				$aligned += scalar(@{$algn});

				$base += $algn->[-1]->[0] + 1;
				$base2 += $algn->[-1]->[1] + 1;
				$size = ($left > 50) ? 50 : $left;
				$size2 = ($n - $base2 > 150) ? 150 : ($n - $base2);
			}
			else{
				$base += $size;
				$base2 = $base - 100 if($base2 < $base - 100);
				$size = ($left > 50) ? 50 : $left;
				$size2 = ($n - $base2 > $base - $base2 + 150) ? ($base - $base2 + 150) : ($n - $base2);
			}
			if( ($aligned + $left) / ($total + $left) < $threshold ){
				last;
			}
		}
	}
	else{
		my ($base, $size) = (0, ($n > 50) ? 50 : $n);
		my ($base2, $size2) = (0, ($m > 150) ? 150 : $m);
		my $left = $n;
		while($left > 0 and $base2 < $m){
			my $algn = alignLCS([@{$positions2}[$base..($base+$size-1)]], [@{$positions1}[$base2..($base2+$size2-1)]]);

			$left -= $size;
			$total += $size;
			if(scalar(@{$algn}) / $size >= $threshold){
				foreach (@{$algn}){
					push(@{$alignment}, [$_->[1]+$base2, $_->[0]+$base]);
				}
				$aligned += scalar(@{$algn});

				$base += $algn->[-1]->[0] + 1;
				$base2 += $algn->[-1]->[1] + 1;
				$size = ($left > 50) ? 50 : $left;
				$size2 = ($m - $base2 > 150) ? 150 : ($m - $base2);
			}
			else{
				$base += $size;
				$base2 = $base - 100 if($base2 < $base - 100);
				$size = ($left > 50) ? 50 : $left;
				$size2 = ($m - $base2 > $base - $base2 + 150) ? ($base - $base2 + 150) : ($m - $base2);
			}
			if( ($aligned + $left) / ($total + $left) < $threshold ){
				last;
			}
		}
	}
	return ($aligned / $total >= $threshold) ? $alignment : [];
}
