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
$program: A perl script for joining a pair of CMAP files
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-f, --first <str>    The first CMAP file for joining (REQUIRED)
	-s, --second <str>   The second CMAP file for joining (REQUIRED)
	-o, --output <str>   The output path (REQUIRED)
	-h, --help           Help
USAGE

our $cml = "$program " . join(" ", @ARGV);
use vars qw($opt_f $opt_s $opt_o $opt_h);
GetOptions( "f|first=s" => \$opt_f,
			"s|second=s" => \$opt_s,
			"o|output=s" => \$opt_o,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

# check the required arguments
die("**ERROR: -f option must be specified\n") if(!defined $opt_f);
die("**ERROR: -s option must be specified\n") if(!defined $opt_s);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);

# get the arguments
my ($first_cmap_in, $second_cmap_in) = ($opt_f, $opt_s);

die("**ERROR: the 1st cmap \"$first_cmap_in\" does not exist\n") unless(-f $first_cmap_in);
die("**ERROR: the 2nd cmap \"$second_cmap_in\" does not exist\n") unless(-f $second_cmap_in);


my ($outdir, $outpre);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = ($opt_o eq "") ? "/" : $opt_o;
	$outpre = basename($opt_f);
	$outdir = "$outdir/$outpre";
}
else{
	$outdir = dirname($opt_o);
	$outpre = basename($opt_o);
}

my ($cmd, $retCode);

$cmd = "mkdir -p $outdir";
$retCode = system($cmd);
die("**ERROR: can not create directory \"$outdir\"\n") if($retCode != 0);
$outdir = abs_path("$outdir");

my ($main_cmap) = readCMap($first_cmap_in);
my ($add_cmap) = readCMap($second_cmap_in);
if( !&isCMapConsistent($main_cmap, $add_cmap) ){
	die("ERROR: The field definitions of the input files are not consistent");
}
my $final_cmap = &concatenateCMaps($main_cmap, $add_cmap);

writeCMapFile($final_cmap, "$outdir/$outpre.cmap");

exit 0;

#### subroutines
sub concatenateCMaps
{
	my ($cmap1, $cmap2) = @_;
	our $cml;
	my ($data_name, $data_type) = &commonDataNameType(\@_);
	my $combinedCmap = {};
	$combinedCmap->{dataName} = $data_name;
	$combinedCmap->{dataType} = $data_type;
	$combinedCmap->{channels} = $cmap1->{channels};
	$combinedCmap->{"version"} = &getMinVersion(\@_);
	my $contigs = {};
	for my $cmap ($cmap1, $cmap2){
		while( my ($id, $value) = each %{ $cmap->{contigs} } ){
			$contigs->{$id} = $value;
		}
	}
	$combinedCmap->{contigs} = $contigs;
	my @header = ();
	my $hostname  = `hostname`; chomp($hostname);
	push(@header, "# hostname=$hostname");
	push(@header, "# \$ $cml");
	push(@header, "# CMAP File Version:\t$combinedCmap->{version}");
	my $channels = $combinedCmap->{channels};
	push(@header, "# Label Channels:\t" . scalar(keys %$channels));
	foreach my $id  (sort { $a <=> $b } keys %{$channels}){
		my $seq = $channels->{$id};
		push(@header, "# Nickase Recognition Site $id:\t$seq");
	}
	my $nmaps = scalar(keys %$contigs);
	push(@header, "# Number of Consensus Maps:\t$nmaps");
	push(@header, "#h " . join("\t", @{$data_name}[0..$#{$data_name}]));
	push(@header, "#f " . join("\t", @{$data_type}[0..$#{$data_type}]));
	$combinedCmap->{"headers"} = \@header;

	return $combinedCmap;
}

sub isCMapConsistent
{
	my ($cmap1, $cmap2) = @_;
	if(!&isChannelsConsistent($cmap1->{channels}, $cmap2->{channels})){
		return 0;
	}
	return 1;
}

sub isChannelsConsistent
{
	my ($channels, $channels2) = @_;
	if(!defined $channels or (scalar(keys %$channels) == 0) ){
		return 0;
	}
	if(!defined $channels2 or (scalar(keys %$channels2) == 0) ){
		return 1;
	}
	if( scalar(keys %$channels) != scalar(keys %$channels2) ){
		return 0;
	}
	foreach my $id (keys %$channels){
		if($channels->{$id} ne $channels2->{$id}){
			return 0;
		}
	}
	return 1;
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
			last if(!defined $cmaps->[$k]->{"dataName"}->[$i]);
			last if(!defined $cmaps->[$k]->{"dataType"}->[$i]);
			last if($cmaps->[$k]->{"dataName"}->[$i] ne $data_name1->[$i]);
			last if($cmaps->[$k]->{"dataType"}->[$i] ne $data_type1->[$i]);
		}
		last if($k>=1);
		push(@data_name, $data_name1->[$i]);
		push(@data_type, $data_type1->[$i]);
	}

	return (\@data_name, \@data_type);
}

sub getMinVersion
{
	my ($cmaps) = @_;
	my $min_ver = $cmaps->[0]->{"version"};
	for(my $i=1; $i<scalar(@{$cmaps}); $i++){
		my $ver = $cmaps->[$i]->{"version"};
		next if(!defined $ver or $ver eq "");
		if($ver < $min_ver){
			$min_ver = $ver;
		}
	}
	return $min_ver;
}
