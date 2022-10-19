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
	unshift(@INC, $progpath);
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

use BNG::Utility;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for merging a pair of CMAP files with different signal channels
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-r, --reference <str>  The reference CMAP file(s) for merging (REQUIRED)
	-l, --layout <str>   The layout file in TSV format (REQUIRED)
	-o, --output <str>   The output path (REQUIRED)
	-h, --help           Help

Example:
	$program -r BSPQI_r.cmap -r BSSSI_r.cmap -l paths.tsv -o multicolor
USAGE

our $cml = "$program " . join(" ", @ARGV);
use vars qw(@opt_r $opt_l $opt_o $opt_h);
GetOptions( "r|reference=s" => \@opt_r,
			"l|layout=s" => \$opt_l,
			"o|output=s" => \$opt_o,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

die("**ERROR: at least two -r options must be specified\n") unless(scalar(@opt_r)>=2);
die("**ERROR: -l option must be specified\n") unless(defined $opt_l);
die("**ERROR: -o option must be specified\n") unless(defined $opt_o);

my ($outdir, $outpre);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = ($opt_o eq "") ? "/" : $opt_o;
	$outpre = basename($opt_l); $outpre =~ s/\.tsv$//;
}
else{
	$outdir = dirname($opt_o);
	$outpre = basename($opt_o);
}

my ($cmd, $retCode);
$cmd = "mkdir -p $outdir";
$retCode = system($cmd);
die("**ERROR: Can not create directory \"$outdir\"\n") if($retCode != 0);

my @cmaps = ();
for(my $i=0; $i<scalar(@opt_r); $i++){
	($cmaps[$i+1]) = readCMap($opt_r[$i]);
}
my $paths = &readPaths($opt_l);

my $combined_cmap = &combineCMaps($paths, \@cmaps);
my $outfile = "$outdir/$outpre.cmap";
writeCMapFile($combined_cmap, $outfile, 1);

exit 0;

sub readPaths
{
	my ($infile) = @_;
	my @paths = ();
	my @mustFields = ("path_id", "node_id", "ori_id", "type", "start", "end", "dir", "xstart", "xend", "xlength");
	my $in;
	if($infile eq "-"){
		$in = \*STDIN;
		$infile = "STDIN";
	}
	else{
		open(TSV, "$infile") or die("Can not open \"$infile\" for reading\n");
		$in = \*TSV;
	}
	my $line;
	my @fields = ();
	my @indices = ();
	my ($nfields, $i, $idx);
	my @columns;
	my $row;
	my $old_path_id = 0;
	my $path_id;
	my $no = 0;
	while($line = <$in>){
		chomp($line);
		if(scalar(@fields) == 0){
			@fields = split(/\t/, $line);
			$nfields = scalar(@fields);
			@indices = &getIndices(\@fields, \@mustFields);
			if(scalar(@indices) != scalar(@mustFields)){
				last;
			}
			next;
		}
		@columns = split(/\t/, $line);
		next if(scalar(@columns != $nfields));
		$path_id = $columns[$indices[0]];
		if($path_id != $old_path_id){
			$no++;
			die("invalid content found in \"$infile\"\n") unless($no == $path_id);
			$old_path_id = $path_id;
			push(@paths, []);
		}
		$row = undef;
		for($i=1; $i<scalar(@indices); $i++){
			$idx = $indices[$i];
			$row->{$fields[$idx]} = $columns[$idx];
		}
		push($paths[-1], $row);
	}
	close $in;

	die("**ERROR: some of the required fields is missing in \"$infile\"\n") unless(scalar(@indices) == scalar(@mustFields));

	return \@paths;
}

sub combineCMaps
{
	my ($paths, $cmaps) = @_;
	our $cml;
	my $cmap;
	my @header = ();
	my ($data_name, $data_type) = &commonDataNameType($cmaps);
	my ($channels, $codes) = &unionOfChannels($cmaps);
	$cmap->{"version"} = &getMinVersion($cmaps);
	$cmap->{"dataName"} = $data_name;
	$cmap->{"dataType"} = $data_type;
	$cmap->{"channels"} = $channels;
	my $nmaps = scalar(@{$paths});
	my $hostname  = `hostname`; chomp($hostname);
	push(@header, "# hostname=$hostname");
	push(@header, "# \$ $cml");
	push(@header, "# CMAP File Version:\t$cmap->{version}");
	push(@header, "# Label Channels:\t" . scalar(keys $channels));
	foreach my $id  (sort { $a <=> $b } keys %{$channels}){
		my $seq = $channels->{$id};
		push(@header, "# Nickase Recognition Site $id:\t$seq");
	}
	push(@header, "# Number of Consensus Maps:\t$nmaps");
	push(@header, "#h " . join("\t", @{$data_name}[0..$#{$data_name}]));
	push(@header, "#f " . join("\t", @{$data_type}[0..$#{$data_type}]));
	$cmap->{"headers"} = \@header;

	for(my $i=0; $i<scalar(@{$paths}); $i++){
		$cmap->{contigs}->{$i+1} = &makeCmap($paths->[$i], $cmaps, $data_name, $codes, scalar(keys %$channels));
	}

	return $cmap;
}

sub makeCmap
{
	my ($path, $cmaps, $data_name, $codes, $num_channels) = @_;
	my %contigs = ();
	
	my @final_positions  = ();
	my @final_channels = ();
	$path = [sort {$a->{start} <=> $b->{start}} @{$path}];
	my @settled_positions = (undef) x ($num_channels + 1);
	my @pending_positions = (undef) x ($num_channels + 1);
	my @rear_graces = (undef) x ($num_channels + 1);
	my ($position, $code);
	foreach my $region (@{$path}){
		my $ctg = $cmaps->[$region->{type}]->{contigs}->{$region->{ori_id}};
		my $code_tbl = $codes->[$region->{type}];
		my ($start, $end) = ($region->{start}, $region->{end});
		my ($slope, $intercept) = getPerfectLinearParams($region->{xstart}, $region->{xend}, $start, $end);
		my @prev_positions = (undef) x ($num_channels + 1);
		my @core_positions = (undef) x ($num_channels + 1);
		my @post_positions = (undef) x ($num_channels + 1);
		my $num_sites = $ctg->{$data_name->[2]}; # NumSites
		my %update = ();
		if($region->{dir} eq "+"){
			for(my $j=0; $j<$num_sites; $j++){
				$position = $ctg->{$data_name->[5]}->[$j];
				$code = $code_tbl->{$ctg->{$data_name->[4]}->[$j]}; # LabelChannel
				$position = sprintf("%.1f", $position * $slope + $intercept); # new position in the combined cmap
				if($position < $start){
					push(@{$prev_positions[$code]}, $position);
				}
				elsif($position > $end){
					push(@{$post_positions[$code]}, $position);
				}
				else{
					push(@{$core_positions[$code]}, $position);
				}
				$update{$code} = 1;
			}
		}
		else{
			for(my $j=$num_sites-1; $j>=0; $j--){
				$position = $ctg->{$data_name->[5]}->[$j];
				$code = $code_tbl->{$ctg->{$data_name->[4]}->[$j]}; # LabelChannel
				$position = sprintf("%.1f", $position * $slope + $intercept); # new position in the combined cmap
				if($position < $start){
					push(@{$prev_positions[$code]}, $position);
				}
				elsif($position > $end){
					push(@{$post_positions[$code]}, $position);
				}
				else{
					push(@{$core_positions[$code]}, $position);
				}
				$update{$code} = 1;
			}
		}
		my ($front_grace, $rear_grace);
		if($slope > 0){
			($front_grace, $rear_grace) = ( $region->{xstart} * (1 - $slope), ($region->{xlength} - $region->{xend}) * (1 - $slope) );
		}
		else{
			($front_grace, $rear_grace) = ( ($region->{xlength} - $region->{xstart}) * (1 + $slope), $region->{xend} * (1 + $slope) );
		}
		my $offset;
		foreach my $k (keys %update){
			($settled_positions[$k], $offset) = mergePositions($settled_positions[$k], $pending_positions[$k], $prev_positions[$k], $core_positions[$k], $rear_graces[$k], $front_grace);
			$rear_graces[$k] = $rear_grace;
			if(defined $offset){
				foreach (@{$post_positions[$k]}){
					$_ -= $offset;
				}
			}
			$pending_positions[$k] = $post_positions[$k];
		}
	}
	my $nsites = 0;
	for(my $k=1; $k<=$num_channels; $k++){
		if(defined $pending_positions[$k]){
			$pending_positions[$k] = cutFromFront($pending_positions[$k]);
			push(@{$settled_positions[$k]}, @{$pending_positions[$k]});
		}
		next if(!defined $settled_positions[$k]);
		push(@final_positions, @{$settled_positions[$k]}); # Position
		push(@final_channels, ($k) x scalar(@{$settled_positions[$k]})); # LabelChannel
		$nsites += scalar(@{$settled_positions[$k]});
	}
	my @indices = sort {$final_positions[$a] <=> $final_positions[$b]} 0..$#final_positions;
	my $shift = $final_positions[$indices[0]] - 20;
	foreach (@final_positions){
		$_ = sprintf("%.1f", $_ - $shift);
	}
	my $end_position = sprintf("%.1f", $final_positions[$indices[-1]] + 20);

	$contigs{$data_name->[1]} = $end_position; # ContigLength
	$contigs{$data_name->[2]} = $nsites; # NumSites
	$contigs{$data_name->[3]} = [1..scalar(@indices)]; # SiteID
	$contigs{$data_name->[4]} = [@final_channels[@indices]]; # LabelChannel
	$contigs{$data_name->[5]} = [@final_positions[@indices]]; # Position
	$contigs{$data_name->[6]} = [("0.0") x $nsites]; # StdDev
	for(my $j=7; $j<scalar(@{$data_name}); $j++){
		$contigs{$data_name->[$j]} = [("1.0") x $nsites]; # Coverage, Occurrence, etc.
	}
	# the last one
	push(@{$contigs{$data_name->[3]}}, scalar(@indices)+1); # SiteID
	push(@{$contigs{$data_name->[4]}}, 0); # LabelChannel
	push(@{$contigs{$data_name->[5]}}, $end_position); # Position
	push(@{$contigs{$data_name->[6]}}, "0.0"); # StdDev
	push(@{$contigs{$data_name->[7]}}, "1.0"); # Coverage
	push(@{$contigs{$data_name->[8]}}, "0.0"); # Occurrence
	for(my $k=9; $k<scalar(@{$data_name}); $k++){
		push(@{$contigs{$data_name->[$k]}}, "1.0");
	}

	return \%contigs;
}

sub commonDataNameType
{
	my ($cmaps) = @_;
	my $data_name1 = $cmaps->[1]->{"dataName"};
	my $data_type1 = $cmaps->[1]->{"dataType"};
	my $count = min(scalar(@{$data_name1}), scalar(@{$data_type1}));
	my ($i, $k);
	my (@data_name, @data_type) = ((), ());
	for($i=0; $i<$count; $i++){
		for($k=$#{ $cmaps }; $k>=2; $k--){
			last if($cmaps->[$k]->{"dataName"}->[$i] ne $data_name1->[$i]);
			last if($cmaps->[$k]->{"dataType"}->[$i] ne $data_type1->[$i]);
		}
		last if($k>=2);
		push(@data_name, $data_name1->[$i]);
		push(@data_type, $data_type1->[$i]);
	}

	return (\@data_name, \@data_type);
}

sub unionOfChannels
{
	my ($cmaps) = @_;
	my %channels = ();
	my @codes = ();
	my $num = 0;
	my %ids = ();
	for(my $i=1; $i<scalar(@{$cmaps}); $i++){
		while( my ($id, $seq) = each %{$cmaps->[$i]->{"channels"}} ){
			if(!defined $ids{$seq}){
				$channels{++$num} = $seq;
				$ids{$seq} = $num;
			}
			$codes[$i]{$id} = $ids{$seq};
		}
	}

	return (\%channels, \@codes);
}

sub getMinVersion
{
	my ($cmaps) = @_;
	my $min_ver = $cmaps->[1]->{"version"};
	for(my $i=2; $i<scalar(@{$cmaps}); $i++){
		my $ver = $cmaps->[$i]->{"version"};
		if($ver < $min_ver){
			$min_ver = $ver;
		}
	}
	return $min_ver;
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

sub cutFromRear
{
	my ($positions) = @_;
	my @cutted = ();
	my ($firstPos, $curPos);
	for(my $i=$#{$positions}; $i>=0; $i--){
		$curPos = $positions->[$i];
		if(defined $firstPos){
			# if the distance between 2 adjacent positions is greater than 100k,
			#  then break
			last if($curPos - $firstPos < -100000);
		}
		unshift(@cutted, $curPos);
		$firstPos = $curPos;
	}
	return \@cutted;
}

sub cutFromFront
{
	my ($positions) = @_;
	my @cutted = ();
	my ($lastPos, $curPos);
	for(my $i=0; $i<scalar(@{$positions}); $i++){
		$curPos = $positions->[$i];
		if(defined $lastPos){
			# if the distance between 2 adjacent positions is greater than 100k,
			#  then break
			last if($curPos - $lastPos > 100000);
		}
		push(@cutted, $curPos);
		$lastPos = $curPos;
	}
	return \@cutted;
}

sub mergePositions
{
	my ($mergedA, $partA, $partB, $mergedB, $graceA, $graceB) = @_;

	my ($merged, $offset);
	push(@{$mergedA}, @{$partA}) if(defined $partA);
	if(!defined $mergedA){
		$merged = cutFromRear($partB); # the boundary part
		push(@{$merged}, @{$mergedB});
		return ($merged, $offset);
	}
	unshift(@{$mergedB}, @{$partB}) if(defined $partB);
	($merged, $offset) = mergePositionsGracefully($mergedA, $mergedB, $graceA, $graceB);
	return ($merged, $offset);
}

sub mergePositionsGracefully
{
	my ($partA, $partB, $graceA, $graceB) = @_;

	my @merged = ();

	$graceA = 50000 if($graceA < 50000);
	$graceB = 50000 if($graceB < 50000);

# align the overlapped regions
	my @positions1 = ();
	my @positions2 = ();

	my ($idx1, $idx2);
	my $lower_pos = $partB->[0] - $graceB;
	$idx1 = locate($partA, $lower_pos);
	push(@positions1, @{$partA}[$idx1..$#{$partA}]);

	my $upper_pos = $partA->[-1] + $graceA;
	$idx2 = locate($partB, $upper_pos);
	push(@positions2, @{$partB}[0..($idx2-1)]);

	my $alignment = localAlign(\@positions1, \@positions2);
	my ($positions, $offset) = combinePositions($alignment, \@positions1, \@positions2);
	push(@merged, @{$partA}[0..($idx1-1)]);
	push(@merged, @{$positions});
	if(defined $offset){
		for(my $k=0; $k<scalar(@{$partB}); $k++){
			$partB->[$k] -= $offset;
		}
	}
	push(@merged, @{$partB}[$idx2..$#{$partB}]);

	return (\@merged, $offset);
}

sub combinePositions
{
	my ($alignment, $positions1, $positions2) = @_;

	my @combined = ();
	my $offset;
	if(scalar(@{$alignment}) == 0){ # no significant overlap found
		if(scalar(@{$positions1}) == 0){
			return $positions2;
		}
		if(scalar(@{$positions2}) == 0){
			return $positions1;
		}
		my $idx = locate($positions2, $positions1->[-1]);
		push(@combined, @{$positions1});
		push(@combined, @{$positions2}[$idx..$#{$positions2}]);
	}
	else{
		my ($idx1, $idx2) = ($alignment->[0]->[0], $alignment->[-1]->[1]);
		push(@combined, @{$positions1}[0..($idx1-1)]);

		my ($i1, $i2);
		$offset = 0;
		foreach my $aln (@{$alignment}){
			($i1, $i2) = @{$aln};
			$offset += $positions2->[$i2] - $positions1->[$i1];
		}
		$offset = sprintf("%.1f", $offset / scalar(@{$alignment}));
		foreach (@$positions2){
			$_ -= $offset;
		}
		my ($l1, $l2);
		my $pos;
		foreach my $aln (@{$alignment}){
			($i1, $i2) = @{$aln};
			if(defined $l1){
				for(my $k=$l1+1; $k<$i1; $k++){
					$pos = sprintf("%.1f", $positions1->[$k]);
					push(@combined, $pos);
				}
				for(my $k=$l2+1; $k<$i2; $k++){
					$pos = sprintf("%.1f", $positions2->[$k]);
					push(@combined, $pos);
				}
			}
			$pos = sprintf("%.1f", ($positions1->[$i1] + $positions2->[$i2]) / 2); 
			push(@combined, $pos);
			($l1, $l2) = ($i1, $i2);
		}

		push(@combined, @{$positions2}[($idx2+1)..$#{$positions2}]);
	}

	return (\@combined, $offset);
}

sub localAlign
{
	my ($positions1, $positions2) = @_;
	my ($m, $n) = (scalar(@{$positions1}), scalar(@{$positions2}));
	my $maxSeedSize = 50;
	if($m <=1 or $n <= 1){
		return [];
	}
	if($m <= $maxSeedSize and $n <= $maxSeedSize){
		my ($alignment) = alignDP($positions1, $positions2);
		return $alignment;
	}
	my $selSize;
	my @seeds;
	$selSize = ($m > $maxSeedSize) ? $maxSeedSize : $m;
	@seeds = @{$positions1}[($m-$selSize)..($m-1)];
	my ($alignment1, $score1) = alignDP(\@seeds, $positions2);
	if($selSize != $m){
		my $offset = $m - $selSize;
		foreach (@{$alignment1}){
			$_->[0] += $offset;
		}
	}

	$selSize = ($n > $maxSeedSize) ? $maxSeedSize : $n;
	@seeds = @{$positions2}[0..($selSize-1)];
	my ($alignment2, $score2) = alignDP($positions1, \@seeds);

	return (scalar(@{$alignment1}) > scalar(@{$alignment2}) or
			scalar(@{$alignment1}) == scalar(@{$alignment2}) and $score1 < $score2) ?
			$alignment1 : $alignment2;
}
