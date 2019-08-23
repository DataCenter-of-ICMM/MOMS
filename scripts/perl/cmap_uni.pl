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

BEGIN{
	my $progpath = abs_path(dirname(abs_path($0)));
	unshift(@INC, $progpath);
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

use BNG::Utility;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for transform all the channels in a multi-channel CMAP file to a unified channel to form a new CMAP file
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-i, --input   The input multi-channel cmap file (REQUIRED)
	-o, --output  The output path (REQUIRED)
	-h, --help    Help

Example:
	$program -i multicolor.cmap -o ./BSPQI-BSSSI
USAGE

our $cml = "$program " . join(" ", @ARGV);
use vars qw($opt_i $opt_o $opt_h);
GetOptions( "i|input=s" => \$opt_i,
			"o|output=s" => \$opt_o,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

die("**ERROR: -i option must be specified\n") unless(defined $opt_i);
die("**ERROR: -o option must be specified\n") unless(defined $opt_o);

my ($outdir, $outpre);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = ($opt_o eq "") ? "/" : $opt_o;
	$outpre = basename($opt_i); $outpre =~ s/\.cmap$//;
	$outpre .= "_unified";
}
else{
	$outdir = dirname($opt_o);
	$outpre = basename($opt_o);
}
my $cmd = "mkdir -p $outdir";
my $retCode = system($cmd);
die("**ERROR: can not create directory $outdir") if($retCode != 0);

my ($cmap) = readCMap($opt_i);

my $uni_cmap = &unifyCMap($cmap);
my $outfile = "$outdir/$outpre.cmap";

writeCMapFile($uni_cmap, $outfile, 1);

exit 0;

sub unifyCMap
{
	my ($cmap) = @_;
	our $cml;
	my ($data_name, $data_type) = ($cmap->{"dataName"}, $cmap->{"dataType"});
	my $channels = &condenseChannels($cmap);
	$cmap->{"channels"} = &condenseChannels($cmap);
	my $hostname = `hostname`; chomp($hostname);
	my @header = ();
	push(@header, "# hostname=$hostname");
	push(@header, "# \$ $cml");
	push(@header, "# CMAP File Version:\t$cmap->{version}");
	push(@header, "# Label Channels:\t" . scalar(keys $channels));
	foreach my $id  (sort { $a <=> $b } keys %{$channels}){
		my $seq = $channels->{$id};
		push(@header, "# Nickase Recognition Site $id:\t$seq");
	}
	push(@header, "# Number of Consensus Maps:\t$cmap->{nContigs}");
	push(@header, "#h " . join("\t", @{$data_name}[0..$#{$data_name}]));
	push(@header, "#f " . join("\t", @{$data_type}[0..$#{$data_type}]));
	$cmap->{"headers"} = \@header;
	my @cmapIds = sort {$a <=> $b} keys %{ $cmap->{contigs} };
	my $no = scalar(@cmapIds);
	my $ichannel;
	my ($i, $j, $k);
	for ($i=0; $i < $no; $i++) {
		my $ctg = $cmap->{contigs}->{$cmapIds[$i]};
		my $contig_len = $ctg->{$data_name->[1]};
		my $numsites = $ctg->{$data_name->[2]};
		$k = 0;
		for ($j=0; $j <= $numsites; $j++) {
			$ichannel = $ctg->{$data_name->[4]}->[$j];
			for (my $m=5; $m<=$#{$data_name}; $m++){
				$ctg->{$data_name->[$m]}->[$k] = $ctg->{$data_name->[$m]}->[$j];
			}
			$ctg->{$data_name->[3]}->[$k] = $k+1; # siteID
			$ctg->{$data_name->[4]}->[$k] = ($ichannel > 0) ? 1 : 0; # unified LabelChannel
			$k++;
		}
		$k--; # minus the last label in channel 0
		if($k > 0){
			$ctg->{$data_name->[5]}->[$k] = $contig_len;
		}
		for(my $m=3; $m<=$#{$data_name}; $m++){
			$j=$numsites;
			while($j>$k){
				pop(@{ $ctg->{$data_name->[$m]}});
				$j--;
			}
		}
		$ctg->{$data_name->[2]} = $k; # numSites
	}

	return $cmap;
}

sub condenseChannels
{
	my ($cmap) = @_;
	my %channels = ();
	my $len = 0;
	while( my ($id, $seq) = each %{$cmap->{"channels"}} ){
		if(length($seq) > $len){
			$len = length($seq);
		}
	}
	$channels{1} = join('', ('N') x $len);

	return \%channels;
}
