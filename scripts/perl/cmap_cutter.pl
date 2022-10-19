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
	unshift(@INC, $progpath);
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

use BNG::Utility;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for cutting a CMAP file at breakpoints according to a coordinate mapping file
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-i, --input <str>  The input CMAP file (REQUIRED)
	-t, --trans <str>  The file containing the coordinate transformation information during cutting (REQRUIED)
	-o, --output <str> The output CMAP file (REQUIRED)
	-h, --help         Help

Example:
	$program -i BSPQI.cmap -t BSPQI_coord_translation.txt -o BSPQI_cut.cmap
USAGE

use vars qw($opt_i $opt_t $opt_o $opt_h);
GetOptions( "i|input=s" => \$opt_i,
			"t|trans=s" => \$opt_t,
			"o|output=s" => \$opt_o,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

# check the required arguments
die("**ERROR: -i option must be specified\n") if(!defined $opt_i);
die("**ERROR: -t option must be specified\n") if(!defined $opt_t);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);

# get the arguments
my ($cmap_in, $cut_coord_file) = ($opt_i, $opt_t);
die("**ERROR: the input cmap file \"$cmap_in\" does not exist\n") unless(-f $cmap_in);
die("**ERROR: the coordinate transformation file \"$cut_coord_file\" does not exist\n") unless(-f $cut_coord_file);
my $outdir = dirname($opt_o);
my $outfile = basename($opt_o);
my ($cmd, $retCode);
$cmd = "mkdir -p $outdir";
$retCode = system($cmd);
die("**ERROR: can not create directory \"$outdir\"\n") if($retCode != 0);
$outdir = abs_path("$outdir");
my $cmap_out = "$outdir/$outfile";

my ($cmap) = readCMap($cmap_in);
my $segments = readCoordTransFile($cut_coord_file);
my $new_cmap = cutCMap($cmap, $segments);
writeCMapFile($new_cmap, "$outdir/$outfile", 1);

exit(0);

sub readCoordTransFile
{
	my ($in_file) = @_;
	my %segments = ();
	my $line;
	my @columns = ();
	my ($oldId, $oldStart, $oldEnd, $newId);
	my $segment;
	open(IN, "<$in_file") || die("ERROR: can not open cut coordinates file $in_file for reading\n");
	while( $line = <IN> ){
		chomp($line);
		@columns = split(/\t/, $line);
		die("ERROR: invalid data format detected for cut coordinates file $in_file\n") if(scalar(@columns) != 6);
		if($line =~ /^\D+/){
			next; # skipping header
		}
		($oldId, $oldStart, $oldEnd, $newId) = @columns[0..3];
		die("ERROR: invalid data format detected for cut coordinates file $in_file\n") if($oldStart >= $oldEnd);
		if(!defined $segments{$oldId}){
			$segments{$oldId} = [{start=>$oldStart, end=>$oldEnd, id=>$newId}];
			next;
		}
		$segment = $segments{$oldId}->[-1];
		die("ERROR: invalid data format detected for cut coordinates file $in_file\n") if($segment->{end} + 1 != $oldStart);
		push(@{$segments{$oldId}}, {start=>$oldStart, end=>$oldEnd, id=>$newId});
	}
	close IN;

	return \%segments;
}

sub cutCMap
{
	my ($cmap, $segments) = @_;
	
	my $new_cmap = {};
	my @headers = @{ $cmap->{headers} };
	my @new_headers = ();
	my @data_name = ();
	my @data_type = ();
	my $line;
	my $NUM_C_LIMITED = 9;
	for (my $i=0; $i<=$#headers; $i++) {
		$line = $headers[$i];
		if($line =~ /^#h\s+/){ # data name
			my $headTag = ""; # the #h
			($headTag, @data_name) = split(/\s+/, $line);
			splice @data_name, $NUM_C_LIMITED;
			push(@new_headers, join("\t", "$headTag $data_name[0]", @data_name[1..$#data_name]));
			next;
		}
		if ($line =~ /^#f\s+/){ # data type
			my $headTag = ""; # the #f
			($headTag, @data_type) = split(/\s+/, $line);
			splice @data_type, $NUM_C_LIMITED;
			push(@new_headers, join("\t", "$headTag $data_type[0]", @data_type[1..$#data_type]));
			next;
		}
		push(@new_headers, $line);
	}
	my $segs;
	foreach my $id (keys %{$cmap->{contigs}}){
		my $ctg = $cmap->{contigs}->{$id};
		$segs = $segments->{$id};
		if( scalar(@{$segs}) <= 1 ){
			if( scalar(@{$segs}) == 0 ){
				next;
			}
			if( (scalar(@{$segs}) == 1) and
				(($segs->[0]->{end} - $segs->[0]->{start} + 1) == int($ctg->{$data_name[1]})) ){
				$new_cmap->{contigs}->{$segs->[0]->{id}} = $ctg;
				next;
			}
		}
		my $numsites = $ctg->{$data_name[2]}; # NumSites
		my ($i, $j) = (0, 0);
		my ($lower, $upper);
		$lower = $segs->[0]->{start};
		while($i<scalar(@{$segs})){
			$upper = $segs->[$i]->{end};
			my %contig = ();
			my $k;
			my $position;
			for($k=0; $j<=$numsites; $j++,$k++){
				$position = $ctg->{$data_name[5]}->[$j]; # Position
				next if($position <= $lower);
				last if($position > $upper);
				$contig{$data_name[3]}->[$k] = ($k+1); # SiteID
				$contig{$data_name[4]}->[$k] = $ctg->{$data_name[4]}->[$j]; # LabelChannel
				$contig{$data_name[5]}->[$k] = sprintf("%.1f", $position - $lower); # Position
				$contig{$data_name[6]}->[$k] = "1.0"; # StdDev
				$contig{$data_name[7]}->[$k] = 1; # Coverage
				$contig{$data_name[8]}->[$k] = 1; # Occurrence
			}
			$contig{$data_name[1]} = sprintf("%.1f", $upper - $lower); # ContigLength
			$contig{$data_name[2]} = $k; # NumSites
			$contig{$data_name[3]}->[$k] = ($k+1); # SiteID
			$contig{$data_name[4]}->[$k] = 0; # LabelChannel
			$contig{$data_name[5]}->[$k] = sprintf("%.1f", $upper - $lower); # Position
			$contig{$data_name[6]}->[$k] = "0.0"; # StdDev
			$contig{$data_name[7]}->[$k] = 1; # Coverage
			$contig{$data_name[8]}->[$k] = 0; # Occurrence
			$lower = $upper;
			$new_cmap->{contigs}->{$segs->[$i]->{id}} = \%contig;
			$i++;
		}
	}
	$new_cmap->{headers} = \@new_headers;
	$new_cmap->{dataName} = \@data_name;
	$new_cmap->{dataType} = \@data_type;
	$new_cmap->{channels} = $cmap->{channels};
	$new_cmap->{version} = $cmap->{version};
	$new_cmap->{nContigs} = scalar(keys %{$new_cmap->{contigs}});

	return $new_cmap;
}
