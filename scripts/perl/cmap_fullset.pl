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
my $usage =<< "USAGE";
$program: A perl script for merging cmaps in different channels into a multi-channel CMAP file
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-i, --input  The input single channel CMAP(s) (REQUIRED)
	-o, --output The output path (REQUIRED)
	-h, --help   Help

Example:
	$program -i BSPQI.cmap -i BSSSI.cmap -i DLE1.cmap -o ./multicolors
USAGE

our $cml = "$program " . join(" ", @ARGV);
use vars qw(@opt_i $opt_o $opt_h);
GetOptions( "i|input=s" => \@opt_i,
			"o|output=s" => \$opt_o,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

die("**ERROR: at leat two -i options must be specified\n") unless(scalar(@opt_i) >= 2);
die("**ERROR: -o option must be specified\n") unless(defined $opt_o);

my ($outdir, $outpre);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = ($opt_o eq "") ? "/" : $opt_o;
	$outpre = "multicolor";
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
for(my $i=0; $i<scalar(@opt_i); $i++){
	($cmaps[$i]) = readCMap($opt_i[$i]);
}

my $combined_cmap = &combineCMaps(\@cmaps);
my $outfile = "$outdir/$outpre.cmap";
writeCMapFile($combined_cmap, $outfile, 1);

exit 0;

sub combineCMaps
{
	my ($cmaps) = @_;
	our $cml;
	my $cmap;
	my @header = ();
	my ($data_name, $data_type) = &commonDataNameType($cmaps);
	my ($channels, $codes) = &unionOfChannels($cmaps);
	$cmap->{"version"} = &getMinVersion($cmaps);
	$cmap->{"dataName"} = $data_name;
	$cmap->{"dataType"} = $data_type;
	$cmap->{"channels"} = $channels;
	my $hostname  = `hostname`; chomp($hostname);
	push(@header, "# hostname=$hostname");
	push(@header, "# \$ $cml");
	push(@header, "# CMAP File Version:\t$cmap->{version}");
	push(@header, "# Label Channels:\t" . scalar(keys $channels));
	foreach my $id  (sort { $a <=> $b } keys %{$channels}){
		my $seq = $channels->{$id};
		push(@header, "# Nickase Recognition Site $id:\t$seq");
	}
	my $lengths = &getMaxLengths($cmaps);
	my $nmaps = scalar(keys %$lengths);
	push(@header, "# Number of Consensus Maps:\t$nmaps");
	push(@header, "#h " . join("\t", @{$data_name}[0..$#{$data_name}]));
	push(@header, "#f " . join("\t", @{$data_type}[0..$#{$data_type}]));
	$cmap->{"headers"} = \@header;

	while( my ($id, $length) = each %$lengths ){
		$cmap->{contigs}->{$id} = &makeCmap($cmaps, $id, $length, $data_name, $codes, scalar(keys %$channels));
	}

	return $cmap;
}

sub makeCmap
{
	my ($cmaps, $id, $length, $data_name, $codes, $num_channels) = @_;
	my %contigs = ();
	
	my @positions  = ();
	my @channels = ();
	for(my $i=0; $i<scalar(@$cmaps); $i++){
		my $cmap = $cmaps->[$i];
		my $ctg = $cmap->{contigs}->{$id};
		next if(!defined $ctg);
		my $num_sites = $ctg->{$data_name->[2]}; # NumSites
		next if($num_sites + 1 != scalar(@{$ctg->{Position}}));
		push(@positions, @{$ctg->{Position}}[0..($num_sites-1)]);
		for(my $j=0; $j<$num_sites; $j++){
			push(@channels, $codes->[$i]{$ctg->{LabelChannel}->[$j]});
		}
	}
	my @indices = sort {$positions[$a] <=> $positions[$b]} 0..$#positions;
	my $nsites = scalar(@indices);
	my $end_position = sprintf("%.1f", $length);

	$contigs{$data_name->[1]} = $end_position; # ContigLength
	$contigs{$data_name->[2]} = $nsites; # NumSites
	$contigs{$data_name->[3]} = [1..scalar(@indices)]; # SiteID
	$contigs{$data_name->[4]} = [@channels[@indices]]; # LabelChannel
	$contigs{$data_name->[5]} = [@positions[@indices]]; # Position
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
	my $data_name1 = $cmaps->[0]->{"dataName"};
	my $data_type1 = $cmaps->[0]->{"dataType"};
	my $count = min(scalar(@{$data_name1}), scalar(@{$data_type1}));
	my ($i, $k);
	my (@data_name, @data_type) = ((), ());
	for($i=0; $i<$count; $i++){
		for($k=$#{ $cmaps }; $k>=1; $k--){
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
	for(my $i=0; $i<scalar(@{$cmaps}); $i++){
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

sub getMaxLengths
{
	my ($cmaps) = @_;
	my %lengths = ();
	foreach my $cmap (@$cmaps){
		my $contigs = $cmap->{contigs};
		foreach my $id (keys %$contigs){
			if(!defined $lengths{$id} or $contigs->{$id}->{ContigLength} > $lengths{$id}){
				$lengths{$id} = $contigs->{$id}->{ContigLength};
			}
		}
	}

	return \%lengths;
}
