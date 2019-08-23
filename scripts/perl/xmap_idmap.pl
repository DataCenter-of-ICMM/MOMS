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
	unshift @INC, $progpath;
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

use BNG::Utility;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for map query IDs in a XMAP file and output a new XMAP file
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-i, --input <str>   The input XMAP file (Required)
	-o, --output <str>  The output XMAP file (Required)
	-k, --keys <str>    The key mapping file output by fa2cmap (Required)
	-h, --help          Help

Exmaple:
	$program -i input.xmap -k contigs_BSPQI_BSSSI_key.txt -o final.xmap
USAGE

use vars qw($opt_i $opt_o $opt_k $opt_h);
GetOptions( "i|input=s" => \$opt_i,
			"o|output=s" => \$opt_o,
			"k|keys=s" => \$opt_k,
			"h|help" => \$opt_h);

die($usage) if($opt_h);
die("**ERROR: -i option must be specified\n") unless(defined $opt_i);
die("**ERROR: -o option must be specified\n") unless(defined $opt_o);
die("**ERROR: -k option must be specified\n") unless(defined $opt_k);

my ($xmap_file, $key_file, $new_xmap_file) = ($opt_i, $opt_k, $opt_o);
my $xmap = readXMap($xmap_file);
die("**ERROR: some of the required field is missing in \"$xmap_file\"\n") unless(&checkXmapFields($xmap->{dataName}));
my $map = &readKeyFile($key_file);
$xmap = &mapKeys($xmap, $map);
writeXMapFile($xmap, $new_xmap_file);

exit 0;

sub checkXmapFields
{
	my ($data_name) = @_;
	my %hasField = ();
	for(my $i=0; $i<scalar(@{$data_name}); $i++){
		$hasField{$data_name->[$i]} = 1;
	}
	my @required_common = ("QryContigID", "QryLen");
	foreach(@required_common){
		if(!defined $hasField{$_}){
			return 0;
		}
	}
	my @required_features = ("RefContigID", "RefLen", "Confidence", "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", "Orientation");
	foreach(@required_features){
		if(!defined $hasField{$_}){
			return 0;
		}
	}
	return 1;
}

sub readKeyFile
{
	my ($key_file) = @_;
	my %map = ();
	open(IN, "$key_file") or die("Can not open \"$key_file\" for reading.\n");
	my $line;
	my @fields;
	my ($id, $newId);
	while($line = <IN>){
		chomp($line);
		if($line !~ /^\d/){
			next if($line =~ /^#/);
			@fields = split(/\t/, $line);
			if(scalar(@fields) != 3 or $fields[0] ne "CompntId" or $fields[1] ne "CompntName" or $fields[2] ne "CompntLength"){
				die("invalid key file: $key_file\n");
			}
			next;
		}
		@fields = split(/\t/, $line);
		next if(scalar(@fields) < 3);
		next if($fields[0] !~ /^\d+$/);

		my @features = split(/\|/, $fields[1]);
		next if(scalar(@features) < 3);
		next if($features[-2] !~ /^\d+$/);

		($id, $newId) = ($fields[0], $features[-2]);
		$map{$id} = $newId;
	}
	close IN;

	return \%map;
}

sub mapKeys
{
	my ($xmap, $map) = @_;

	my $QryContigIDs = $xmap->{hits}->{QryContigID};
	my $nHits = $xmap->{totalHits};
	my ($qryId, $newQryId);
	for(my $i=0; $i<$nHits; $i++){
		$qryId = $QryContigIDs->[$i];
		$newQryId = $map->{$qryId};
		if(!defined $newQryId){
			warn "mapKeys: new ID for $qryId is undefined\n";
			next;
		}
		$QryContigIDs->[$i] = $newQryId;
	}
	return $xmap;
}
