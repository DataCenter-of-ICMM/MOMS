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
use Parallel::ForkManager;
use Cwd 'abs_path';
use POSIX qw/ceil/;
use Storable qw(dclone);
use List::Util qw[min max];

BEGIN{
	my $progpath = abs_path(dirname($0));
	unshift(@INC, $progpath);
	select(STDERR); $|= 1;
	select(STDOUT); $|= 1;
}

use BNG::Utility;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for filling gaps within a CMAP file
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-i, --cmap <str>    The CMAP file whose gaps are to be filled (REQUIRED)
	-o, --output <str>  The output path (STDOUT)
	-a, --add <str>     Additional CMAP(s) containing small contigs for gap filling
	--addonly           Only use additional CMAP(s)
	-t, --thread <num>  The number of threads for parallel computing (default: 4)
	-h, --help          Help

Example:
	$program -i BSPQI.hybrid.cmap -a BSPQI.hybrid-BNG-non_used.cmap -o BSPQI.filled
USAGE

use vars qw($opt_i $opt_o @opt_a $opt_ao $opt_t $opt_h);
GetOptions( "i|cmap=s" => \$opt_i,
			"o|output=s" => \$opt_o,
			"a|add=s" => \@opt_a,
			"addonly" => \$opt_ao,
			"t|thread=i" => \$opt_t,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

# check the required arguments
die("**ERROR: -i option must be specified\n") if(!defined $opt_i);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);
if(defined $opt_ao){
	die("**ERROR: -a option must be specified if --addonly has been specified\n") if(scalar(@opt_a)<=0);
}

# get the arguments
my ($outdir, $outfile);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = ($opt_o eq "") ? "/" : $opt_o;
	$outfile = basename($opt_i);
	$outfile =~ s/\.cmap$/.filled/;
}
else{
	$outdir = dirname($opt_o);
	$outfile = basename($opt_o);
}
$outfile .= ".cmap";
my $nthreads = (defined $opt_t and ($opt_t > 0)) ? $opt_t : 4;

# input the cmap file
my ($cmap) = readCMap($opt_i);
print "Read in $opt_i completed with $cmap->{nContigs} cmaps.\n";

my @cmaps = ($cmap);
for(my $i=0; $i<scalar(@opt_a); $i++){
	($cmaps[$i+1]) = readCMap($opt_a[$i]);
	print " Read in $opt_a[$i] completed with $cmaps[$i+1]->{nContigs} cmaps.\n";
}

# prepare for output
my ($cmd, $retCode);
$cmd = "mkdir -p $outdir";
$retCode = system($cmd);
die("**ERROR: can not create directory \"$outdir\"\n") if($retCode != 0);

my $mean_distance = meanLabelDistance($cmap);
my $gappedFragments = retrieveGappedFragments($cmap, $mean_distance);
my $smallContigs = gatherSmallContigs(\@cmaps, $mean_distance, (defined $opt_ao ? 1 : 0));

my $matches = matchFragments($gappedFragments, $smallContigs, $nthreads);
my $new_cmap = fillCMap(\@cmaps, $matches);

writeCMapFile($new_cmap, "$outdir/$outfile", 1);

exit(0);

sub meanLabelDistance
{
	my ($cmap) = @_;
	my ($totalLength, $num_sites) = (0, 0);
	foreach my $ctg (values %{$cmap->{contigs}}){
		$totalLength += $ctg->{ContigLength};
		$num_sites += $ctg->{NumSites};
	}
	return $totalLength / $num_sites;
}

sub retrieveGappedFragments
{
	my ($cmap, $mean_distance) = @_;
	my @fragments = ();
	my $minGapSize = sprintf("%.1f", 10 * $mean_distance);
	my $strictMinGapSize = sprintf("%.1f", 0.618 * $minGapSize);
	my @cmapIds = sort {$a <=> $b} keys %{ $cmap->{contigs} };
	my ($positions, $prePos, $pos, $gapSize);
	for (my $i=0; $i<scalar(@cmapIds); $i++) {
		my $id = $cmapIds[$i];
		my $ctg = $cmap->{contigs}->{$id};
		my $num_sites = $ctg->{NumSites};
		next if($num_sites < 51); # only check those relatively large contigs
		$positions = $ctg->{Position};
		$prePos= $positions->[1];
		my $lastJ;
		for (my $j=2; $j<$num_sites-1; $j++){
			$pos = $positions->[$j];
			$gapSize = ($pos - $prePos);
			if( $gapSize >= $strictMinGapSize){
				if( $gapSize >= $minGapSize ){ # should be large enough to be treated as a gap
					my $k = (defined $lastJ) ? (($j - 29 > $lastJ) ? ($j - 29) : $lastJ) : (($j > 29) ? $j - 29 : 0);
					if($j - $k >= 3){
						my $frag = {cmapid=>$id, offset=>$k, positions=>[], offset2=>$j, gapsize=>$gapSize};
						while($k < $j){
							push(@{$frag->{positions}}, $positions->[$k++]);
						}
						my $maxK = ($num_sites - $j) > 29 ? (29 +  $j) : $num_sites;
						my $prePos2 = $pos;
						while($k < $maxK){
							my $pos2 = $positions->[$k];
							last if($pos2 - $prePos2 >= $strictMinGapSize);
							push(@{$frag->{positions}}, $pos2);
							$prePos2 = $pos2;
							$k++;
						}
						if($k - $j >= 3 and scalar(@{$frag->{positions}}) >= 10){
							push(@fragments, $frag);
						}
					}
				}
				$lastJ = $j;
			}
			$prePos = $pos;
		}
	}
	return \@fragments;
}

sub gatherSmallContigs
{
	my ($cmaps, $mean_distance, $startK) = @_;
	my @fragments = ();
	my $minFragmentSize = sprintf("%.1f", 20 * $mean_distance);
	my $maxFragmentSize = sprintf("%.1f", 20 * 20 * $mean_distance);
	for (my $k=$startK; $k<scalar(@{$cmaps}); $k++){
		my $cmap = $cmaps->[$k];
		my @cmapIds = sort {$a <=> $b} keys %{ $cmap->{contigs} };
		for (my $i=0; $i<scalar(@cmapIds); $i++) {
			my $id = $cmapIds[$i];
			my $ctg = $cmap->{contigs}->{$id};
			my $length = $ctg->{ContigLength};
			# the contigs for gap-filling should be large enough but NOT too large
			next if($length < $minFragmentSize or $length > $maxFragmentSize);
			next if($ctg->{NumSites} < 10 or $ctg->{NumSites} > 250);
			my $ext_id = ($k == 0) ? $id : "${k}_$id";
			my $frag = {cmapid=>$ext_id, offset=>0, positions=>[], size=>$length};
			push(@{$frag->{positions}}, @{$ctg->{Position}}[0..($#{$ctg->{Position}}-1)]);
			push(@fragments, $frag);
		}
	}
	return \@fragments;
=begin
distance	labels		length
262,980	87	(11, 41)	690,669
229,483	127	(16, 75)	940,814
336,367	123	(21, 69)	1,334,124
117,817	51	(12, 22)	400,018
191,827	50	(12, 19)	441,939
217,457	24	(4,6)		295,016
160,354	21	(1,5)		231,155
249,518	29	(8,39)		918,017
=cut
}

sub matchFragments
{
	my ($gappedFragments, $contigs, $nthreads) = @_;

	my %matches = ();
	my @sorted = sort{$contigs->[$a]->{size} <=> $contigs->[$b]->{size}} 0..$#{$contigs};
	my ($cnt, $total) = (0, scalar(@{$gappedFragments}));
	my %gapped = ();
	foreach (@{$gappedFragments}){
		$gapped{$_->{cmapid}} = 1;
	}
	my $nFrag = scalar(keys %gapped);
	print "found $total gaps within $nFrag cmaps.\n";
	print "found " . scalar(@sorted) . " small cmaps for gap filling.\n";
	my $pm = Parallel::ForkManager->new($nthreads);
	$pm->run_on_finish(
		sub{
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
			my $id = $data_structure_reference->{input};
			my $match = $data_structure_reference->{result};
			if(defined $match){
				if(!defined $matches{$id}){
					$matches{$id} = [];
				}
				push(@{$matches{$id}}, $match);
				print "+";
			}
			else{
				print ".";
			}
			if( (++$cnt) % 10 == 0){
				print " ";
			}
		}
	);
# for DEBUG
=begin
	my %selected = ( # false negtives
#		18994=>[1004127],
#		19032=>[1000747], # N
#		245853=>[1004025], #
		359170=>[132409],
		0=>0
	);
	my %selected = (
		11381=>[1001525],
		119=>[1001145],
		132354=>[1001677],
		132372=>[1001508],
		132465=>[1000920,1002188],
		135=>[1000770], #
		139=>[189115], #
		150=>[1000924],
		156=>[1000730],
		160=>[1002619,1003521],
		163=>[1000821],
		18963=>[1000973,1003254],
		18965=>[1001052],
		18994=>[1004127], #o
		19005=>[1002447],
		19018=>[1000922],
		19024=>[13539], #
		19032=>[1000747], # N
		19069=>[1002667],
		245774=>[1000182], #
		245853=>[1003066,1004025], #o
		28=>[1001070], #
		47=>[1003060],
		6=>[1002877], #
		7=>[1000458],
		75651=>[1001572], #
		75703=>[1001336,1001735,12095],
		75710=>[1001581],
		75714=>[13133], #
		76=>[1000974], #
		0=>0
	);
=cut
=begin
	my %selected = ( # false positives
		124=>[1001580],
		150=>[1004025], #
		191=>[1004156],
		12469=>[12095],
		135=>[1003380],
		18908=>[1003380],
		19041=>[1004156],
		359170=>[132409], #
		41=>[1004156],
		78=>[1003547],
		0=>0
	);
=cut
	foreach my $gappedFrag (@{$gappedFragments}){
		my $id = $gappedFrag->{cmapid};
		my $gapSize = $gappedFrag->{gapsize};
		my $pid = $pm->start and next;
		# forked thread
		my ($match, $optimalScore);
#		next if(!defined $selected{$id}); # for DEBUG
		foreach my $idx (@sorted){
			my $smallFrag = $contigs->[$idx];
			next if("$smallFrag->{cmapid}" eq "$id");
#			next if(!grep {/^$smallFrag->{cmapid}$/} @{$selected{$id}}); # for DEBUG
			my ($alignment, $score, $gapidx) = alignGapped($gappedFrag, $smallFrag);
			if(scalar(@{$alignment}) > 0){
				if(!defined $optimalScore or compareScore($score, $optimalScore) > 0){
					$optimalScore = $score;
					$match = {id=>$smallFrag->{cmapid}, align=>$alignment, gapidx=>$gapidx};
				}
			}
		}
=begin
			if(defined $match){
				if(!defined $matches{$id}){
					$matches{$id} = [];
				}
				push(@{$matches{$id}}, $match);
				print "+";
			}
			else{
				print ".";
			}
			if( (++$cnt) % 10 == 0){
				print " ";
			}
=cut
		$pm->finish(0, {result => $match, input => $id});
	}
	$pm->wait_all_children;
	print "\n";
	return \%matches;
}

sub alignGapped
{
	my ($gappedFrag, $smallFrag) = @_;
	my ($alignment1, $score1) = alignGappedDP($smallFrag->{positions}, $gappedFrag->{positions}, $gappedFrag->{offset2} - $gappedFrag->{offset});
	my @positions2 = reverse(@{$gappedFrag->{positions}});
	my $firstPos = $positions2[0];
	foreach (@positions2){
		$_ = $firstPos - $_;
	}
	my ($alignment2, $score2) = alignGappedDP($smallFrag->{positions}, \@positions2, scalar(@positions2) - ($gappedFrag->{offset2} - $gappedFrag->{offset}) );
	my ($alignment, $score, $offset_small, $offset_gapped, $gapidx);
	if(compareScore($score1, $score2) >= 0){
		($alignment, $score) = ($alignment1, $score1);
		($offset_small, $offset_gapped) = ($smallFrag->{offset}, $gappedFrag->{offset});
		if(scalar(@{$alignment}) > 0){
			foreach (@{$alignment}){
				my $tmp = $_->[0];
				$_->[0] = $_->[1] + $offset_gapped;
				$_->[1] = $tmp + $offset_small;
				if(!defined $gapidx and $_->[0] >= $gappedFrag->{offset2}){
					$gapidx = $_->[0];
				}
			}
		}
	}
	else{
		($alignment, $score) = ($alignment2, $score2);
		($offset_small, $offset_gapped) = ($smallFrag->{offset}, $gappedFrag->{offset});
		my $n2 = $#positions2;
		if(scalar(@{$alignment}) > 0){
			$alignment = [reverse(@{$alignment})];
			foreach (@{$alignment}){
				my $tmp = $_->[0];
				$_->[0] = $n2 - $_->[1] + $offset_gapped;
				$_->[1] = $tmp + $offset_small;
				if(!defined $gapidx and $_->[0] >= $gappedFrag->{offset2}){
					$gapidx = $_->[0];
				}
			}
		}
	}
	return ($alignment, $score, $gapidx);
}

sub fillCMap
{
	my ($cmaps, $matches) = @_;
	my %newCmap = ();
	my $baseCmap = $cmaps->[0];
	for my $field ('headers', 'channels', 'version', 'dataType', 'dataName'){
		if(ref($baseCmap->{$field}) eq ""){
			$newCmap{$field} = $baseCmap->{$field};
		}
		elsif(ref($baseCmap->{$field}) eq "ARRAY"){
			@{$newCmap{$field}} = @{dclone($baseCmap->{$field})}
		}
		else{
			$newCmap{$field} = dclone($baseCmap->{$field});
		}
	}
	my $baseContigs = $baseCmap->{contigs};
	my $newContigs = $newCmap{contigs} = {};
	my $data_name = $baseCmap->{dataName};
	my %used = ();
	my $cnt = 0;
	my $lengthFilled = 0;
	my ($len1, $len2);
	for my $hostid (keys %{$matches}){
		my $ctg = $newContigs->{$hostid} = {};
		my @final_positions = @{$baseContigs->{$hostid}->{Position}};
		my $offset;
		my @ins_positions = ();
		my $host_num_sites = $baseContigs->{$hostid}->{NumSites};
		$len1 = $baseContigs->{$hostid}->{ContigLength};
		foreach my $match ( sort{ $a->{gapidx} <=> $b->{gapidx} } @{$matches->{$hostid}} ){
			my $guestid = $match->{id};
			my $ctg2;
			if($guestid =~ /(\d+)_(\d+)/){
				$ctg2 = $cmaps->[$1]->{contigs}->{$2};
			}
			else{
				$ctg2 = $baseContigs->{$guestid};
			}
			$len2 = $ctg2->{ContigLength};
			$lengthFilled += $len2;
			my $positions2 = $ctg2->{Position};
			my $alignment = $match->{align};
			my $gapidx = $match->{gapidx};
			my ($align1, $align2) = splitAlign($alignment, $gapidx);
			next if(!defined $align1 or !defined $align2); # usually this can't happen
			if($guestid !~ /_/){
				my $guest_num_sites = $ctg2->{NumSites};
				my ($host_left, $host_right) = ($alignment->[0]->[0], $host_num_sites - $alignment->[-1]->[0] - 1);
				my ($guest_left, $guest_right);
				if($alignment->[0]->[1] <= $alignment->[-1]->[1]){
					($guest_left, $guest_right) = ($alignment->[0]->[1], $guest_num_sites - $alignment->[-1]->[1] - 1);
				}
				else{
					($guest_left, $guest_right) = ($guest_num_sites - $alignment->[0]->[1] - 1, $alignment->[-1]->[1]);
				}
				my $none_covered = max($guest_left - $host_left, 0) + max($guest_right - $host_right, 0);
				$used{$guestid} = 1 if($none_covered < max(10, $guest_num_sites / 4));
			}
			$cnt++;
			($match->{insert}, $offset) = pickOutInsert(\@final_positions, $positions2, $align1, $align2);
			push(@ins_positions, @{$match->{insert}});
			shiftPositions(\@final_positions, $gapidx, $offset);
			print "$cnt\t$hostid\t$len1\t$guestid\t$len2\n";
			foreach (@{$align1}){
				print "($_->[0],$_->[1])";
			}
			print "\n";
			foreach (@{$align2}){
				print "($_->[0],$_->[1])";
			}
			print "\n";
		}
		push(@final_positions, @ins_positions);
		my @indices = sort {$final_positions[$a] <=> $final_positions[$b]} 0..$#final_positions;
		my $nsites = scalar(@indices);
		my $end_position = sprintf("%.1f", $final_positions[$indices[-1]] + 20);
		$ctg->{$data_name->[1]} = $end_position; # ContigLength
		$ctg->{$data_name->[2]} = $nsites; # NumSites
		$ctg->{$data_name->[3]} = [1..$nsites]; # SiteID
		$ctg->{$data_name->[4]} = [(1) x $nsites]; # LabelChannel
		$ctg->{$data_name->[5]} = [@final_positions[@indices]]; # Position
		$ctg->{$data_name->[6]} = [("0.0") x $nsites]; # StdDev
		for(my $j=7; $j<scalar(@{$data_name}); $j++){
			$ctg->{$data_name->[$j]} = [("1.0") x $nsites]; # Coverage, Occurrence, etc.
		}
# the last one
		push(@{$ctg->{$data_name->[3]}}, $nsites + 1); # SiteID
		push(@{$ctg->{$data_name->[4]}}, 0); # LabelChannel
		push(@{$ctg->{$data_name->[5]}}, $end_position); # Position
		push(@{$ctg->{$data_name->[6]}}, "0.0"); # StdDev
		push(@{$ctg->{$data_name->[7]}}, "1.0"); # Coverage
		push(@{$ctg->{$data_name->[8]}}, "0.0"); # Occurrence
		for(my $k=9; $k<scalar(@{$data_name}); $k++){
			push(@{$ctg->{$data_name->[$k]}}, "1.0");
		}
	}
	my $nCmap = scalar(keys %{$matches});
	$lengthFilled = sprintf("%.1f", $lengthFilled);
	print "\n $cnt gaps in $nCmap cmaps have been filled, summing up to $lengthFilled bp.\n";
	foreach (keys %{$baseContigs}){
		next if(defined $newContigs->{$_} or defined $used{$_});
		$newContigs->{$_} = dclone($baseContigs->{$_});
	}
	$newCmap{nContigs} = scalar(keys %{$newContigs});

	return \%newCmap;
}

sub splitAlign
{
	my ($alignment, $gapidx) = @_;
	my $i;
	for($i=0; $i<scalar(@{$alignment}); $i++){
		last if($alignment->[$i]->[0] >= $gapidx);
	}
	my ($align1, $align2);
	$align1 = [@{$alignment}[0..($i-1)]] if($i >= 1);
	$align2 = [@{$alignment}[$i..$#{$alignment}]] if($i<=$#{$alignment});

	return ($align1, $align2);
}

sub pickOutInsert
{
	my ($positions1, $positions2, $align1, $align2) = @_;

	my ($ys_idx, $ye_idx) = ($align1->[0]->[0], $align1->[-1]->[0]);
	my ($xs_idx, $xe_idx) = ($align1->[0]->[1], $align1->[-1]->[1]);
	my ($ystart, $yend) = ($positions1->[$ys_idx], $positions1->[$ye_idx]);
	my ($xstart, $xend) = ($positions2->[$xs_idx], $positions2->[$xe_idx]);

# calculate the linear fitting parameter
	my ($slope, $intercept) = getPerfectLinearParams($xstart, $xend, $ystart, $yend);
# adjusted the intercept
	my ($i, $pair, $diff);
	for($diff=0,$i=0; $i<scalar(@{$align1}); $i++){
		$pair = $align1->[$i];
		$diff += ($slope * $positions2->[$pair->[1]] + $intercept) - $positions1->[$pair->[0]];
	}
	$diff /= scalar(@{$align1});
	$intercept -= $diff;

	my @positions_ins = ();
	my ($start, $end) = (min($align1->[-1]->[1], $align2->[0]->[1]) + 1, max($align1->[-1]->[1], $align2->[0]->[1]) - 1);
	for(my $k=$start; $k<=$end; $k++){
		push(@positions_ins, sprintf("%.1f", $slope * $positions2->[$k] + $intercept) );
	}

	for($diff=0,$i=0; $i<scalar(@{$align2}); $i++){
		$pair = $align2->[$i];
		$diff += ($slope * $positions2->[$pair->[1]] + $intercept) - $positions1->[$pair->[0]];
	}
	$diff /= scalar(@{$align2});

	return (\@positions_ins, $diff);
}

sub shiftPositions
{
	my ($positions, $idx, $offset) = @_;

	for(my $i=$idx; $i<scalar(@{$positions}); $i++){
		$positions->[$i] = sprintf("%.1f", $positions->[$i] + $offset);
	}
}
