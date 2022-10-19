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
use Math::Round;
use Cwd qw(abs_path);
use Statistics::LineFit;
use Storable qw(dclone);
use List::Util qw[min max];

BEGIN{
        my $progpath = abs_path(dirname($0));
        unshift(@INC, $progpath);
        select(STDERR); $| = 1;
        select(STDOUT); $| = 1;
}

use BNG::Utility;
use Getopt::Std;
use Tie::File;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for exporting AGP file and FASTA file, given relevant NGS and multi-channel BNG information
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options];
	-i, --input <str>  The input XMAP file (REQUIRED)
	-o, --output <str> The output path containing the output directory and prefix (REQUIRED)
	-c, --cmap <str>   The input CMAP file from hybrid assembly (REQUIRED)
	-s, --seq <str>    The input FASTA file from NGS assembly (REQUIRED)
	-m, --map <str>    The file containing the mapping information from ID to sequence name (REQUIRED)
	-t, --trans <str>  The file containing the coordinate transformation information when resolving conflicts
	-g, --gap <int>    The length of gap to be inserted between overlapping NGS contigs in a scaffold (default: 13)
	-h, --help         Help
USAGE

#global varaibles

#enzyme cut site sequences
our %enzyme = (
	"BSPQI" => "GCTCTTC",
	"BSSSI" => "CACGAG",
	"BBVCI" => "CCTCAGC",
	"BSMI"  => "GAATGC",
	"BSRDI" => "GCAATG",
	"BSECI" => "ATCGAT",
	"BAMHI" => "GGATCC",
	"DLE1" => "CTTAAG"
	# You can add more enzymes here ...
	
);

use vars qw($opt_i $opt_o $opt_c $opt_s $opt_m $opt_t $opt_g $opt_h);
if(!GetOptions( "i|input=s" => \$opt_i,
	"o|output=s" => \$opt_o,
	"c|cmap=s" => \$opt_c,
	"s|seq=s" => \$opt_s,
	"m|map=s" => \$opt_m,
	"t|trans=s" => \$opt_t,
	"g|gap=i" => \$opt_g,
	"h|help" => \$opt_h)){
	print ("Please try -h for more details\n");
	exit 1;
}

die($usage) if($opt_h);

# check the required arguments
die("**ERROR: -i option must be specified\n") if(!defined $opt_i);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);
die("**ERROR: -c option must be specified\n") if(!defined $opt_c);
die("**ERROR: -s option must be specified\n") if(!defined $opt_s);
die("**ERROR: -m option must be specified\n") if(!defined $opt_m);

# get the XMAP file name
my $xmap_file = $opt_i;

# set the output directory and prefix
my ($outdir, $outpre);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = $opt_o;
	$outpre = basename($opt_i); $outpre =~ s/\.xmap$//;
}
else{
	$outdir = dirname($opt_o);
	$outpre = basename($opt_o);
}
my ($cmd, $retCode);
$cmd = "mkdir -p $outdir";
$retCode = system($cmd);
die("**ERROR: Can not create directory \"$outdir\"\n") if($retCode != 0);

# get the other arguments
my $fasta_file = $opt_s;
my $ngs_namemap_file = $opt_m;
my $cut_coord_file = $opt_t;
my $hybrid_cmap_file = $opt_c;

###################Main entry point of the export script########################
my ($hybridCmap, $numcontig, $contigLength) = readCMap($hybrid_cmap_file);
print("Number of contigs in hybrid: $numcontig\n");
my $seqTable = readFasta($ngs_namemap_file, $fasta_file);
my $contigInfo = readContigInfo($seqTable, $cut_coord_file);
my $alignTable = readAlign($xmap_file);
my $scaffoldInfo = processAlign($alignTable, $contigInfo, $seqTable);

# output Gap file and AGP files
my $gap_out_file = "$outdir/${outpre}.gap";
my $agp_out_file = "$outdir/$outpre.agp";
my $begin_end_file = "$outdir/${outpre}_trimHeadTailGap.coord";
printGapFile($scaffoldInfo, $gap_out_file);
printAGPFile($scaffoldInfo, $agp_out_file);
printAGPAuxFile($scaffoldInfo, $begin_end_file);

# output FASTA files
my $fasta_out_file = "$outdir/$outpre.fasta";
my $NCBI_fasta_out_file = "$outdir/${outpre}_NCBI.fasta";
my $channels = $hybridCmap->{channels};
my $cut_sites = [@{$channels}{sort { $a <=> $b } keys (%$channels)}];
my $max_len = 1;
for(@$cut_sites){
	my $len = length($_);
	$max_len = $len if $len > $max_len;
}
my @NCBI_cut_sites = ('N' x $max_len) x scalar(@{$cut_sites});
my $padding_gap_len = ((!defined $opt_g) ? 13 : (($opt_g < 1) ? 13 : $opt_g));
my $padding_gap = ("N" x $padding_gap_len);
printFasta($scaffoldInfo, $NCBI_fasta_out_file, \@NCBI_cut_sites, $padding_gap);

my $unused_fasta_out = "$outdir/${outpre}_NOT_SCAFFOLDED.fasta";
printUnUsedNGS($scaffoldInfo, $unused_fasta_out);

exit 0;

sub readFasta
{
	my ($namemap_file, $fasta_file) = @_;
	my %seqTable = ();
	my ($seqid, $name);
	my %rev_map = (); # temporary map for retrieving ID given NAME
	my $line;
	open(NM, "<$namemap_file") or die("Can not open \"$namemap_file\" for reading\n");
	while($line = <NM>){
		chomp($line);
		next if($line !~ /^\d/);
		my @tokens = split(/\t/, $line);
		next if(scalar(@tokens) < 3);
		$seqid = $tokens[0];
		$name = (split(/\s+/, $tokens[1]))[0];
		$seqTable{$seqid}{name} = $name;
		$seqTable{$seqid}{len} = $tokens[2];
		$rev_map{$name} = $seqid;
	}
	close NM;
	my $seq = "";
	open(FA, "<$fasta_file") or die("Can not open \"$fasta_file\" for reading\n");
	while($line = <FA>){
		chomp($line);
		if($line =~ /^>(\S+)/){
			if($seq ne ""){
				$seqid = $rev_map{$name};
				if(defined $seqid){
					$seqTable{$seqid}{seq} = $seq;
				}
				else{
					warn "Warning: $name is not recognized to be a valid NAME, check $namemap_file or $fasta_file\n";
				}
				$seq = "";
			}
			$name = $1;
			next;
		}
		$seq .= $line;
	}
	close FA;
	if($seq ne ""){
		$seqid = $rev_map{$name};
		if(defined $seqid){
			$seqTable{$seqid}{seq} = $seq;
		}
		else{
			warn "Warning: $name is not recognized to be a valid NAME, check $namemap_file or $fasta_file\n";
		}
	}
	return \%seqTable;
}

sub readContigInfo
{
	my ($seqTable, $coord_file) = @_;
	my $line;
	my %contigInfo = ();
	my %subseqs = ();
	my ($id, $name, $newId, $start, $end, $len);
	open(CR, "<$coord_file") or die ("Can not open \"$coord_file\" for reading\n");
	while($line = <CR>){
		chomp($line);
		next if($line =~ /^oldId/); # skipping header
		my @tokens = split(/\t/, $line);
		next if(scalar(@tokens) < 6);
		($id, $newId) = ($tokens[0], $tokens[3]);
		($start, $end) = ($tokens[1], $tokens[2]);
		$len = $end - $start + 1;
		$contigInfo{$newId} = {seqid=>$id, start=>$start, end=>$end, len=>$len, length=>$seqTable->{$id}->{len}, sibling=>{}};
		if(!exists $subseqs{$id}){
			$subseqs{$id} = [];
		}
		push(@{$subseqs{$id}}, $newId);
	}
	close CR;

	foreach my $group (values %subseqs){
		next unless (scalar(@$group) > 1);
		my @sorted_group = sort { $contigInfo{$a}->{start} <=> $contigInfo{$b}->{start} || $contigInfo{$b}->{end} <=> $contigInfo{$a}->{end} } @{$group};
		my ($i, $j, $id1, $id2);
		my ($start1, $end1, $start2, $end2);
		my ($offset, $overlap);
		my ($sibling1, $sibling2);
		for($i=0; $i<$#sorted_group; $i++){
			$id1 = $sorted_group[$i];
			$sibling1 = $contigInfo{$id1}->{sibling};
			($start1, $end1) = ($contigInfo{$id1}->{start}, $contigInfo{$id1}->{end});
			for($j=$i+1; $j<scalar(@sorted_group); $j++){
				$id2 = $sorted_group[$j];
				$start2 = $contigInfo{$id2}->{start};
				last unless ($start2 <= $end1);
				# found an overlap
				$end2 = $contigInfo{$id2}->{end};
				($offset, $overlap) = ($start2 - $start1, min($end1, $end2) - $start2 + 1);
				$sibling2 = $contigInfo{$id2}->{sibling};
				$sibling1->{$id2} = {offset => $offset, overlap => $overlap};
				$sibling2->{$id1} = {offset => -$offset, overlap => $overlap};
			}
		}
	}
	return \%contigInfo;
}

# read the alignments from XMAP file and only retain part of the information for further processing
sub readAlign
{
	my ($xmap_file) = @_;
	my $xmap = readXMap($xmap_file);
	my @mustFields = ("QryContigID", "RefContigID", "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", "Orientation", "Confidence", "QryLen", "RefLen", "LabelChannel");
	my $hits = $xmap->{hits};
	foreach my $field (@mustFields){
		if(!defined $hits->{$field}){
			die("**ERROR: field \"$field\" is missing in \"$xmap_file\"\n")
		}
	}
	my %alignTable = ();
	$alignTable{alignments} = [];
	for(my $i=0; $i < $xmap->{totalHits}; $i++){
		my $alignment = [];
		for(my $j=0; $j<scalar(@mustFields); $j++){
			push(@$alignment, $hits->{$mustFields[$j]}->[$i]);
		}
		push(@{$alignTable{alignments}}, $alignment);
	}
	$alignTable{nchannels} = $xmap->{nchannels};
	return \%alignTable;
}

# process the alignments and construct the scaffolds
sub processAlign
{
	my ($alignTable, $contigInfo, $seqTable) = @_;
	my %scaffoldTable = ();
	my $scaffold;
	my ($refId, $iChannel);
	my ($refStartPos, $refEndPos, $extRefStartPos, $extRefEndPos);
	my ($qryId, $seqid, $qryStartPos, $qryEndPos, $extQryStartPos, $extQryEndPos, $qryLen);
	my ($orientation, $slope, $intercept);
	my $ctgInfo;
	my $nchannels = $alignTable->{nchannels};
	# 0	QryContigID
	# 1	RefContigID
	# 2	QryStartPos
	# 3	QryEndPos
	# 4	RefStartPos
	# 5	RefEndPos
	# 6	Orientation
	# 7	Confidence
	# 8	QryLen
	# 9	RefLen
	#10	LabelChannel
	foreach my $alignment (@{$alignTable->{alignments}}){
		($refId, $iChannel) = ($alignment->[1], $alignment->[10]); # RefContigID, LabelChannel
		if(!defined $scaffoldTable{$refId}){ # initialization for scaffold 'refId'
			$scaffold = ();
			$scaffold->{regions}[0] = undef;
			for(my $type=0; $type<=$nchannels; $type++){
				$scaffold->{regions}[$type] = [];
			}
			$scaffoldTable{$refId} = $scaffold;
		}
		$scaffold = $scaffoldTable{$refId};
		$scaffold->{len} = round($alignment->[9]); # RefLen

		$orientation = $alignment->[6]; # Orientation
		($refStartPos, $refEndPos) = ($alignment->[4], $alignment->[5]); # RefStartPos, RefEndPos
		($qryStartPos, $qryEndPos, $qryLen) = (round($alignment->[2]), round($alignment->[3]), round($alignment->[8])); # QryStartPos, QryEndPos, QryLen
		($slope, $intercept) = getPerfectLinearParams($qryStartPos, $qryEndPos, $refStartPos, $refEndPos);
# adjust end positions in reference contig so that the aligned pair of fragments have the same length
		$refStartPos = round($qryStartPos * $slope + $intercept);
		$refEndPos = round($refStartPos + $slope * ($qryEndPos - $qryStartPos));

		if($orientation eq "+"){
			$extRefStartPos = round($refStartPos - ($qryStartPos - 1));
			$extRefEndPos = round($refEndPos + ($qryLen - $qryEndPos));
			$extQryStartPos = 1;
			$extQryEndPos = $qryLen;
		}
		else{
			$extRefStartPos = round($refStartPos - ($qryLen - $qryStartPos));
			$extRefEndPos = round($refEndPos + ($qryEndPos - 1));
			$extQryStartPos = $qryLen;
			$extQryEndPos = 1;
		}
		$qryId = $alignment->[0]; # QryContigID
		$ctgInfo = $contigInfo->{$qryId};
		$seqid = $ctgInfo->{seqid};
		if($seqid != $qryId and $ctgInfo->{start} != 0){
# transform to coordinates in pre-cutting contigs
			$qryStartPos += $ctgInfo->{start};
			$qryEndPos += $ctgInfo->{start};
			$extQryStartPos += $ctgInfo->{start};
			$extQryEndPos += $ctgInfo->{start};
			($slope, $intercept) = getPerfectLinearParams($qryStartPos, $qryEndPos, $refStartPos, $refEndPos); # update the parameters
		}
		push(@{$scaffold->{regions}[$iChannel]}, {extStart=>$extRefStartPos, extEnd=>$extRefEndPos, coreStart=>$refStartPos, coreEnd=>$refEndPos, qryStart=>$qryStartPos, qryEnd=>$qryEndPos, extQryStart=>$extQryStartPos, extQryEnd=>$extQryEndPos, qryid=>$qryId, id=>$seqid, len=>$qryLen, dir=>$orientation, slope=>$slope, intercept=>$intercept});
		$contigInfo->{$qryId}->{used} = 1; # used in the scaffolding
		while ( my ($sibling_id, $sibling) = each(%{$contigInfo->{$qryId}->{sibling}}) ){
			next unless ($sibling->{overlap} / $contigInfo->{$sibling_id}->{len} > 0.75); # this criterion can be further tuned
			$contigInfo->{$sibling_id}->{used} = 1;
		}
	}
#	foreach $scaffold (values %scaffoldTable){
	foreach my $refId (sort {$a <=> $b} keys %scaffoldTable){
		$scaffold = $scaffoldTable{$refId};
		for(my $type=0; $type<=$nchannels; $type++){
			next if(scalar(@{$scaffold->{regions}[$type]}) == 0);
			# sort all the regions and remove embedded fragments
			$scaffold->{regions}[$type] = removeEmbedded([sort {$a->{coreStart} <=> $b->{coreStart} || $b->{coreEnd} <=> $a->{coreEnd}} (@{$scaffold->{regions}[$type]})]);
		}
		$scaffold->{rgns} = mergeRegions(\@{$scaffold->{regions}});
	}
	my %scaffoldInfo = ();
	$scaffoldInfo{scaffolds} = \%scaffoldTable;
	$scaffoldInfo{contigs} = $contigInfo;
	$scaffoldInfo{seqs} = $seqTable;
	return \%scaffoldInfo;
}

sub removeEmbedded
{
	my ($regions) = @_;
	my ($i, $k) = (0, 0);
	while( $i < scalar(@$regions) ){
		if( ($i == 0) || $regions->[$i]->{coreEnd} > $regions->[$k-1]->{coreEnd} ){
			$regions->[$k++] = $regions->[$i];
		}
		$i++;
	}
	splice @$regions, $k;

	return $regions;
}

sub mergeRegions
{
	my ($regions) = @_;

# adjust regions in multiple channels by comparison of contigs
	adjustRegions($regions);

	my $combined_regions = [];
# process the regions in an order
	my @heap = ();
	my ($type, $j);
	my $rgns;
	for($type=0; $type<scalar(@$regions); $type++){
		$rgns = $regions->[$type];
		next if(scalar(@$rgns) == 0);
		push(@heap, {rgn=>$rgns->[0], i=>$type, j=>1}); # indices of channel and next element
	}
	minHeap(\@heap); # create a min heap
	my ($element, $rgn, $last_rgn);
	while(scalar(@heap) > 0){
		$element = $heap[0];
		$rgn = $element->{rgn};
# process the region
		if(!defined $last_rgn){
			$rgn->{setStart} = $rgn->{extStart}; # no conflict found, then extend to the entire part
			push(@{$combined_regions}, $rgn);
			$last_rgn = $rgn;
		}
		else{
			if($rgn->{extStart} > $last_rgn->{extEnd}){ # non-overlapped
				$last_rgn->{setEnd} = $last_rgn->{extEnd};
				$rgn->{setStart} = $rgn->{extStart}; # conflict found, then extend to the entire part
				push(@{$combined_regions}, $rgn);
				$last_rgn = $rgn;
			}
			else{ # overlapped
				if($last_rgn->{id} == $rgn->{id} and $last_rgn->{dir} eq $rgn->{dir}){ # overlapped in the same original contig
# they use the same coordinates, so they can be merged seamlessly
					if($rgn->{coreEnd} > $last_rgn->{coreEnd}){
						if($last_rgn->{dir} eq "+"){
							$last_rgn->{qryEnd} += $rgn->{coreEnd} - $last_rgn->{coreEnd};
						}
						else{
							$last_rgn->{qryEnd} -= $rgn->{coreEnd} - $last_rgn->{coreEnd};
						}
						$last_rgn->{coreEnd} = $rgn->{coreEnd};
					}
					if($rgn->{extEnd} > $last_rgn->{extEnd}){
						if($last_rgn->{dir} eq "+"){
							$last_rgn->{extQryEnd} += $rgn->{extEnd} - $last_rgn->{extEnd};
						}
						else{
							$last_rgn->{extQryEnd} -= $rgn->{extEnd} - $last_rgn->{extEnd};
						}
						$last_rgn->{extEnd} = $rgn->{extEnd};
					}
				}
				else{ # overlapped for different contigs; TODO: alignment
					if($rgn->{extEnd} > $last_rgn->{extEnd}){ # not embedded
						my $tmp;
						if($rgn->{coreStart} > $last_rgn->{coreEnd}){ # non-overlapped in core part
							if($last_rgn->{extEnd} < $rgn->{coreStart}){
								$tmp = $last_rgn->{extEnd};
							}
							elsif($rgn->{extStart} > $last_rgn->{coreEnd}){
								$tmp = $rgn->{extStart} - 1;
							}
							else{
								$tmp = round(($rgn->{coreStart} - 1 + $last_rgn->{coreEnd})/2); 
							}
						}
						else{ # overlapped in core part
							$tmp = $rgn->{coreStart} - 1;
						}
						$last_rgn->{setEnd} = $tmp;
						$rgn->{setStart} = $tmp + 1;
						push(@{$combined_regions}, $rgn);
						$last_rgn = $rgn;
					}
				}
			}
		}
		($type, $j) = ($element->{i}, $element->{j});
		if($j < scalar(@{$regions->[$type]})){
			$heap[0] = {rgn=>$regions->[$type]->[$j], i=>$type, j=>$j+1};
			heapify(\@heap, 0);
		}
		else{
			shift(@heap);
			minHeap(\@heap); # adjust the heap
		}
	}
	if(defined $last_rgn){
		if($last_rgn->{extEnd} < $last_rgn->{coreEnd} + 100000){
			$last_rgn->{setEnd} = $last_rgn->{extEnd}; # no conflict found, then extend to the entire part
		}
		else{
			$last_rgn->{setEnd} = $last_rgn->{coreEnd} + 100000;
		}
	}
	return $combined_regions;
}

# adjust regions within a scaffold
sub adjustRegions
{
	my ($regions) = @_;

	my %occurrences = ();
	my $occurs;
	my ($type, $j, $id);
	my $ntypes = scalar(@$regions)-1;
# collect all the occurrences, grouped by sequence contig
	for($type=0; $type<=$ntypes; $type++){
		my $rgns = $regions->[$type];
		my ($last_rgn, $rgn);
		for($j=0; $j<scalar(@$rgns); $j++){
			$rgn = $rgns->[$j];
			$id = $rgn->{id};
			if(!defined $occurrences{$id}){
				$occurrences{$id} = {};
			}
			$occurs = $occurrences{$id};
			if(!defined $occurs->{$type}){
				$occurs->{$type} =[];
			}
			push(@{$occurs->{$type}}, $rgn); # already sorted by {coreStart}, {coreEnd}
			if(defined $last_rgn){
				$last_rgn->{next} = $rgn;
				$rgn->{prev} = $last_rgn;
			}
			$last_rgn = $rgn;
		}
	}
	my @ids = ();
	my @ids2 = ();
	while( ($id, $occurs) = each (%occurrences) ){
		my $cnt = scalar(keys %{$occurs});
		if($cnt < $ntypes){
			if($cnt <= 1){
				push(@ids, $id);
			}
			else{
				push(@ids2, $id);
			}
			next;
		}
		calculateRgnOffsets($occurs); # calculate offsets for each rgn using k-merging sort
		adjustRgns4Contig($occurs); # adjust regions for each rgn
	}
	foreach $id (@ids2){
		$occurs = $occurrences{$id};
		calculateRgnOffsetsByAdjacency($occurs);
		adjustRgns4Contig($occurs, 1);

		calculateRgnOffsets($occurs); # calculate offsets for each rgn using k-merging sort
		adjustRgns4Contig($occurs); # adjust regions for each rgn
	}
	foreach $id (@ids){
		$occurs = $occurrences{$id};
		calculateRgnOffsetsByAdjacency($occurs);
		adjustRgns4Contig($occurs);
	}
}

sub calculateRgnOffsets
{
	my ($occurs) = @_;
# process the regions in an order
	my @heap = ();
	my ($type, $j, $rgns);
	while( ($type, $rgns) = each (%$occurs) ){
		push(@heap, {rgn=>$rgns->[0], i=>$type, j=>1}); # indices of channel and next element
	}
	minHeap(\@heap); # create a min heap
	my ($element, $rgn, $first_rgn, $last_rgn);
	my ($last_coor, $coor);
	my ($point, $offset);
	my %offsets = ();
	my @regions = ();
	while(scalar(@heap) > 0){
		$element = $heap[0];
		($type, $j) = ($element->{i}, $element->{j});
		$rgn = $element->{rgn};
# process the region
		if(defined $last_rgn){
			for(my $i=$#regions; $i>=0; $i--){
				my $test_rgn = $regions[$i];
				$point = commonPoint($test_rgn, $rgn);
				if(defined $point){
					$last_rgn = $test_rgn;
					last;
				}
			}
			if(!defined $point){
				$point = looselyCommonPoint($last_rgn, $rgn);
			}
			if(defined $point){
				$last_coor = round($last_rgn->{slope} * $point + $last_rgn->{intercept});
				$coor = round($rgn->{slope} * $point + $rgn->{intercept});
				$offset = $last_rgn->{offset} + ($coor - $last_coor);
			}
			else{
				$offset = $last_rgn->{offset};
			}
			$rgn->{offset} = $offset;
			if(!defined $offsets{$type}){
				$offsets{$type} = [0, 0];
			}
			$offsets{$type}->[0]++;
			$offsets{$type}->[1] += $rgn->{offset};
			push(@regions, $rgn);
		}
		else{ # the 1st one
			$rgn->{offset} = 0;
			push(@regions, $rgn);
		}
		$last_rgn = $rgn;

		if($j < scalar(@{$occurs->{$type}})){
			$heap[0] = {rgn=>$occurs->{$type}->[$j], i=>$type, j=>$j+1};
			heapify(\@heap, 0);
		}
		else{
			shift(@heap);
			minHeap(\@heap); # adjust the heap
		}
	}
	$offset = 0;
	if(scalar(@regions) > 0){
		my $values;
		while( ($type, $values) = each (%offsets) ){
			$offset += $values->[1] / $values->[0];
		}
		$offset = round($offset / scalar(@regions));
	}
	foreach $rgn (@regions){
		$rgn->{offset} -= $offset;
	}
}

sub adjustRgns4Contig
{
	my ($occurs, $bUpdate) = @_;
	my ($type, $rgns);
	while( ($type, $rgns) = each (%$occurs) ){
		foreach (@{$rgns}){
			adjustRgn($_, $bUpdate);
		}
	}
}

sub adjustRgn
{
	my ($rgn, $bUpdate) = @_;
	if(!defined $rgn->{offset}){
		return;
	}
	if($rgn->{offset} != 0){
		my $offset = $rgn->{offset};
		$rgn->{extStart} -= $offset;
		$rgn->{extEnd} -= $offset;
		$rgn->{coreStart} -= $offset;
		$rgn->{coreEnd} -= $offset;
	}
	if($bUpdate){
		# update the parameters
		($rgn->{slope}, $rgn->{intercept}) = getPerfectLinearParams($rgn->{qryStart}, $rgn->{qryEnd}, $rgn->{coreStart}, $rgn->{coreEnd});
		$rgn->{offset} = undef;
	}
}

sub calculateRgnOffsetsByAdjacency
{
	my ($occurs) = @_;
	foreach my $rgns (values %$occurs){
		foreach my $rgn (@$rgns){
			getOffsetByAdjacency($rgn);
		}
	}
}

sub getOffsetByAdjacency
{
	my ($rgn) = @_;
	if(defined $rgn->{offset}){
		return;
	}
	my ($prev, $next) = ($rgn->{prev}, $rgn->{next});
	my @tofill = ();
	while(defined $prev){
		last if(defined $prev->{offset});
		unshift(@tofill, $prev);
		$prev = $prev->{prev};
	}
	while(defined $next){
		last if(defined $next->{offset});
		push(@tofill, $next);
		$next = $next->{next};
	}
	my $offset;
	if(defined $prev){
		if(defined $next){
			$offset = round(($prev->{offset} + $next->{offset})/2);
		}
		else{
			$offset = $prev->{offset};
		}
	}
	else{
		if(defined $next){
			$offset = $next->{offset};
		}
		else{
			$offset = 0;
		}
	}
	$rgn->{offset} = $offset;
	foreach(@tofill){
		$_->{offset} = $offset;
	}
}

sub commonPoint
{
	my ($rgn, $rgn2) = @_;
	my ($start, $end);
	if($rgn->{slope} != $rgn2->{slope}){
		return undef;
	}
	my $slope = $rgn->{slope};
	if($slope > 0){
		$start = max($rgn->{extQryStart}, $rgn2->{extQryStart});
		$end = min($rgn->{extQryEnd}, $rgn2->{extQryEnd});
	}
	else{
		$start = min($rgn->{extQryStart}, $rgn2->{extQryStart});
		$end = max($rgn->{extQryEnd}, $rgn2->{extQryEnd});
	}
	return ($slope * ($end - $start) < 0) ? undef : round( ($start + $end) / 2 );
}

sub looselyCommonPoint
{
	my ($rgn, $rgn2) = @_;
	my ($start, $end);
	if($rgn->{slope} != $rgn2->{slope}){
		if($rgn2->{extStart} > $rgn->{extEnd} + 50000){ # the gap is too big for adjustment
			return undef;
		}
		# if loosely overlapped
		($start, $end) = ($rgn2->{qryStart}, $rgn2->{qryEnd});
		return round( ($start + $end) / 2 ); # return the middle point so that the direction can be ignored
	}
	my $slope = $rgn->{slope};
	if($slope > 0){
		$start = max($rgn->{extQryStart}, $rgn2->{extQryStart}) - 25000;
		$end = min($rgn->{extQryEnd}, $rgn2->{extQryEnd}) + 25000;
	}
	else{
		$start = min($rgn->{extQryStart}, $rgn2->{extQryStart}) + 25000;
		$end = max($rgn->{extQryEnd}, $rgn2->{extQryEnd}) - 25000;
	}
	return ($slope * ($end - $start) < 0) ? undef : round( ($start + $end) / 2 );
}

#### subroutines for k-merging sort
sub minHeap
{
	my ($array) = @_;
	my $len = scalar(@$array);
	for(my $i=int($len/2) - 1; $i>=0; $i--){
		heapify($array, $i);
	}
}

sub heapify
{
	my ($heap, $start) = @_;
	my $dad = $start;
	my $son = $dad * 2 + 1;
	my $end = $#{ $heap };
	while($son <= $end){
		if($son + 1 <= $end && compareElem($heap->[$son+1]->{rgn}, $heap->[$son]->{rgn}) < 0){
			$son++; # choose the smaller one
		}
		if(compareElem($heap->[$dad]->{rgn}, $heap->[$son]->{rgn}) < 0){ # dad is already the smallest one
			return;
		}
		swapElem($heap, $dad, $son);
		$dad = $son;
		$son = $dad * 2 + 1;
	}
}

sub compareElem
{
	my ($rgn1, $rgn2) = @_;

	if($rgn1->{coreStart} < $rgn2->{coreStart}){
		return -1;
	}
	if($rgn1->{coreStart} > $rgn2->{coreStart}){
		return 1;
	}
	if($rgn1->{coreEnd} > $rgn2->{coreEnd}){
		return -1;
	}
	if($rgn1->{coreEnd} < $rgn2->{coreEnd}){
		return 1;
	}
	return 0;
}

sub swapElem
{
	my ($array, $first, $second) = @_;
	my $temp = $array->[$first];
	$array->[$first] = $array->[$second];
	$array->[$second] = $temp;
}
######

sub getNameAndLen
{
	my ($name, $start, $end, $length) = @_;
	my $len = ($end >= $start) ? ($end - $start + 1) : ($start - $end + 1);
	if($len < $length){
		$name .= "_subseq_$start:$end";
	}
	return ($name, $len);
}

#print the gap informaton from alignments to a file
sub printGapFile
{
	my ($scaffoldInfo, $out_file) = @_;
	print "GAP output file $out_file\n";
	open(OUT, ">$out_file") or die("Can not open \"$out_file\" for writing\n");
	print OUT "#NGSId1\tNGSId2\tSuperScaffoldId\tXmapGapLength\tAdjustedGapLength\tNGSLength1\tNGSLength2\n";
	my $scaffoldTable = $scaffoldInfo->{scaffolds};
	my $contigInfo = $scaffoldInfo->{contigs};
	my $seqTable = $scaffoldInfo->{seqs};
	foreach my $refId (sort {$a <=> $b} keys %$scaffoldTable){
		my $regions = $scaffoldTable->{$refId}->{rgns};
		my ($curRgn, $nextRgn);
		my ($curCtg, $nextCtg);
		my ($name1, $name2, $len1, $len2);
		my ($gap, $oriGap);
		$curRgn = $regions->[0];
		$curCtg = $seqTable->{$curRgn->{id}};
		($name1, $len1) = getNameAndLen($curCtg->{name}, $curRgn->{extQryStart}, $curRgn->{extQryEnd}, $curCtg->{len});
		for(my $i=1; $i<scalar(@$regions); $i++){
			$nextRgn = $regions->[$i];
			$nextCtg = $seqTable->{$nextRgn->{id}};
			($name2, $len2) = getNameAndLen($nextCtg->{name}, $nextRgn->{extQryStart}, $nextRgn->{extQryEnd}, $nextCtg->{len});
			($gap, $oriGap) = ($nextRgn->{setStart} - $curRgn->{setEnd} - 1, $nextRgn->{coreStart} - $curRgn->{coreEnd} - 1);
			print OUT "$name1\t$name2\t$refId\t$oriGap\t$gap\t$len1\t$len2\n";
			$curRgn = $nextRgn;
			($name1, $len1) = ($name2, $len2);
		}
		# the last one
		my $lastRgn = $regions->[-1];
		my $lastCtg = $seqTable->{$lastRgn->{id}};
		($name1, $len1) = getNameAndLen($lastCtg->{name}, $lastRgn->{extQryStart}, $lastRgn->{extQryEnd}, $lastCtg->{len});
		print OUT "$name1\t$name1\t$refId\t0\t0\t$len1\t$len1\n";
	}
	close OUT;
}

sub printAGPFile
{
	my ($scaffoldInfo, $agp_out_file) = @_;
	print "AGP output file $agp_out_file\n";
	open(OUT, ">$agp_out_file") or die("Can not open \"$agp_out_file\" for writing.\n");
# print AGP header
	print OUT "##agp-version\t2.0\n".
			  "# Organism:  \n".
			  "# Platform:     \n".
			  "# Model:        \n".
			  "# Enzyme(s):    \n".
			  "# BioSample:    \n".
			  "# BioProject:   \n".
              "# Obj_Name\tObj_Start\tObj_End\tPartNum\tCompnt_Type\tCompntId_GapLength\tCompntStart_GapType\tCompntEnd_Linkage\tOrientation_LinkageEvidence\n";
# print content for scaffolded contigs
	my $scaffoldTable = $scaffoldInfo->{scaffolds};
	my $seqTable = $scaffoldInfo->{seqs};
	foreach my $refId (sort {$a <=> $b} keys %$scaffoldTable){
		my $regions = $scaffoldTable->{$refId}->{rgns};
		my ($curRgn, $nextRgn);
		my ($curCtg, $nextCtg);
		my ($name1, $name2, $len1, $len2);
		$curRgn = $regions->[0];
		$curCtg = $seqTable->{$curRgn->{id}};
		($name1, $len1) = getNameAndLen($curCtg->{name}, $curRgn->{extQryStart}, $curRgn->{extQryEnd}, $curCtg->{len});
		my $k = 0;
		print OUT "Super-Scaffold_$refId\t$curRgn->{setStart}\t$curRgn->{setEnd}\t" . ++$k . "\tW\t$name1\t1\t$len1\t$curRgn->{dir}\n";
		for(my $i=1; $i<scalar(@$regions); $i++){
			$nextRgn = $regions->[$i];
			$nextCtg = $seqTable->{$nextRgn->{id}};
			($name2, $len2) = getNameAndLen($nextCtg->{name}, $nextRgn->{extQryStart}, $nextRgn->{extQryEnd}, $nextCtg->{len});
			print OUT "Super-Scaffold_$refId\t" . ($curRgn->{setEnd} + 1) . "\t" . ($nextRgn->{setStart} - 1) . "\t" . ++$k . "\tN\t" . ($nextRgn->{setStart} - $curRgn->{setEnd} - 1) . "\tscaffold\tyes\tmap\n";
			print OUT "Super-Scaffold_$refId\t$nextRgn->{setStart}\t$nextRgn->{setEnd}\t" . ++$k . "\tW\t$name2\t1\t$len2\t$nextRgn->{dir}\n";
			$curRgn = $nextRgn;
			($name1, $len1) = ($name2, $len2);
		}
	}
# print content for singleton contigs
	my $contigInfo = $scaffoldInfo->{contigs};
	foreach my $qryId (sort {$a <=> $b} keys %$contigInfo){
		my $ctg = $contigInfo->{$qryId};
		next if($ctg->{used}); # skip those used in scaffolds
		my $seqid = $ctg->{seqid};
		my $contig = $seqTable->{$seqid};
		my ($name, $len) = getNameAndLen($contig->{name}, $ctg->{start}, $ctg->{end}, $contig->{len});
		print OUT "${name}_obj\t1\t$len\t1\tW\t$name\t1\t$len\t+\n";
	}
	close OUT;
}

sub printAGPAuxFile
{
	my ($scaffoldInfo, $begin_end_file) = @_;
	print "Auxiliary AGP output file $begin_end_file\n";
	open(OUT, ">$begin_end_file") or die("Can not open \"$begin_end_file\" for writing.\n");
# print AGP Aux file header
	print OUT "##agp-version\t2.0\n".
			  "# Organism:   \n".
			  "# Platform:     \n".
			  "# Model:        \n".
			  "# Enzyme(s):    \n".
			  "# BioSample:    \n".
			  "# BioProject:   \n".
			  "Obj_Id\tHeadTrimmedLength\tTailTrimmedLength\n";
# print content for scaffolds
	my $scaffoldTable = $scaffoldInfo->{scaffolds};
	foreach my $refId (sort {$a <=> $b} keys %$scaffoldTable){
		my $scaffold = $scaffoldTable->{$refId};
		my $regions = $scaffold->{rgns};
		print OUT "Super-Scaffold_$refId\t$regions->[0]->{setStart}\t" . ($scaffold->{len} - $regions->[-1]->{setEnd}) . "\n";
	}
# print content for singleton contigs
	my $seqTable = $scaffoldInfo->{seqs};
	my $contigInfo = $scaffoldInfo->{contigs};
	foreach my $qryId (sort {$a <=> $b} keys %$contigInfo){
		my $ctg = $contigInfo->{$qryId};
		next if($ctg->{used}); # skip those used in scaffolds
		my $seqid = $ctg->{seqid};
		my $contig = $seqTable->{$seqid};
		my ($name, $len) = getNameAndLen($contig->{name}, $ctg->{start}, $ctg->{end}, $contig->{len});
		print OUT "${name}_obj\t0\t0\n";
	}
	close OUT;
}

sub printFasta
{
	my ($scaffoldInfo, $fasta_out_file, $cut_sites, $padding_gap) = @_;
	my $scaffoldTable = $scaffoldInfo->{scaffolds};
	my $seqTable = $scaffoldInfo->{seqs};
	print "FASTA output file $fasta_out_file\n";
	open(OUT, ">$fasta_out_file") or die("Can not open \"$fasta_out_file\" for writing.\n");
	foreach my $refId (sort {$a <=> $b} keys %$scaffoldTable){
		my $regions = $scaffoldTable->{$refId}->{rgns};
		my ($curRgn, $nextRgn);
		my ($seq, $gap_size, $gap, $len);
		$curRgn = $regions->[0];
		my $k = 0;
		print OUT ">Super-Scaffold_$refId\n";
		$len = ($curRgn->{setEnd} - $curRgn->{setStart} + 1);
		if($curRgn->{dir} eq "+"){
			$seq = retrieveSequence($curRgn->{id}, $seqTable, $curRgn->{dir},
						 $curRgn->{setStart} - $curRgn->{extStart} + $curRgn->{extQryStart} - 1, $len);
		}
		else{
			$seq = retrieveSequence($curRgn->{id}, $seqTable, $curRgn->{dir},
						 $curRgn->{extEnd} - $curRgn->{setEnd} + $curRgn->{extQryEnd} - 1, $len);
		}
		print OUT "$seq";
		for(my $i=1; $i<scalar(@$regions); $i++){
			$nextRgn = $regions->[$i];
			$gap_size = ($nextRgn->{setStart} - $curRgn->{setEnd} - 1);
			if($gap_size > 0){
				$gap = ("N" x $gap_size);
				print OUT "$gap";
			}
			$curRgn = $nextRgn;
			$len = ($curRgn->{setEnd} - $curRgn->{setStart} + 1);
			if($curRgn->{dir} eq "+"){
				$seq = retrieveSequence($curRgn->{id}, $seqTable, $curRgn->{dir},
							 $curRgn->{setStart} - $curRgn->{extStart} + $curRgn->{extQryStart} - 1, $len);
			}
			else{
				$seq = retrieveSequence($curRgn->{id}, $seqTable, $curRgn->{dir},
							 $curRgn->{extEnd} - $curRgn->{setEnd} + $curRgn->{extQryEnd} - 1, $len);
			}
			print OUT "$seq";
		}
		print OUT "\n";
	}
	close OUT;
}

sub retrieveSequence
{
	my ($seqid, $seqTable, $dir, $offset, $len) = @_;
	$len = $seqTable->{$seqid}{len} if(!defined $len);
	if($len <= 0){
		return "";
	}
	$dir = "+" if(!defined $dir);
	$offset = 0 if(!defined $offset or $offset < 0);
	my $seq = substr($seqTable->{$seqid}{seq}, $offset, $len);
	if($dir ne "+"){
		$seq = reverseComplement($seq);
	}
	return $seq;
}

sub reverseComplement
{
	my ($str) = @_;
	$str =~ tr/ACGTacgt/TGCAtgca/;
	return reverse($str);
}

sub printUnUsedNGS
{
	my ($scaffoldInfo, $fasta_out_file) = @_;
	my $contigInfo = $scaffoldInfo->{contigs};
	my $seqTable = $scaffoldInfo->{seqs};
	print "NONE-SCAFFOLDED output file $fasta_out_file\n";
	open(OUT, ">$fasta_out_file") or die("Can not open \"$fasta_out_file\" for writing.\n");
	foreach my $qryId (sort {$a <=> $b} keys %$contigInfo){
		my $ctg = $contigInfo->{$qryId};
		next if($ctg->{used}); # skip those used in scaffolds
		my $seqid = $ctg->{seqid};
		my $contig = $seqTable->{$seqid};
		my ($name, $len) = getNameAndLen($contig->{name}, $ctg->{start}, $ctg->{end}, $contig->{len});
		print OUT ">$name\n";
		my $seq = retrieveSequence($seqid, $seqTable, "+", $ctg->{start}, $len);
		print OUT "$seq\n";
		while ( my ($sibling_id, $sibling) = each (%{$contigInfo->{$qryId}->{sibling}}) ){
			next unless ($sibling->{overlap} / $contigInfo->{$sibling_id}->{len} > 0.75); # this criterion can be further tuned
			$contigInfo->{$sibling_id}->{used} = 1;
		}
	}
	close OUT;
}
