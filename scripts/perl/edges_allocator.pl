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
use File::Basename;
use Cwd 'abs_path';

BEGIN{
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for allocating edges 
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program edges_ENZYME1_ENZYME2.tsv [edges_ENZYME3_ENZYME4.tsv ...] glues_ENZYME1_ENZYME2.tsv [glues_ENZYMES_ENZYME4.tsv ...] outdir
USAGE

if(scalar(@ARGV) < 3){
	print $usage;
	exit(1);
}

my %enzymeINDEX = (
	"BSPQI" => 1,
	"BSSSI" => 2,
	"BBVCI" => 3,
	"BSMI" => 4,
	"BSRDI" => 5,
	"BSECI" => 6,
	"BAMHI" => 7,
	"DLE1" => 8
);

my $file_type;
my @edges = ();
my @glues = ();
my @edgeMustFields = ("id_str", "nglues", "overlap", "intercept", "slope", "id1", "length1", "id2", "length2", "glues");
my @glueMustFields = ("ID", "QryContigID", "QryLen", "RefContigID.x", "RefLen.x", "RefContigID.y", "RefLen.y", "Orientation.x", "Orientation.y", "Overlap", "Confidence.x", "Confidence.y", "QryStartPos.x", "QryEndPos.x", "RefStartPos.x", "RefEndPos.x", "QryStartPos.y", "QryEndPos.y", "RefStartPos.y", "RefEndPos.y", "OriContigID.x", "offset.x", "len.x", "OriContigID.y", "offset.y", "len.y");
my ($name1, $name2);
my ($type1, $type2);
my %usedType = ();
my $outdir;
foreach my $file (@ARGV){
	if($file =~ /edges_([^_]*)_([^_]*)\.tsv$/){
		($name1, $name2) = ($1, $2);
		$file_type = "edge";
	}
	elsif($file =~ /glues_([^_]*)_([^_]*)\.tsv$/){
		($name1, $name2) = ($1, $2);
		$file_type = "glue";
	}
	else{
		if(!defined $outdir){
			$outdir = $file;
			my $cmd = "mkdir -p $outdir";
			my $retCode = system($cmd);
			die("**ERROR: can not create directory \"$outdir\"\n") if($retCode != 0);
			next;
		}
		die("$file does not match the pattern \"edges_PAT1_PAT2.tsv\" or \"glues_PAT1_PAT2.tsv\"")
	}
	$type1 = $enzymeINDEX{$name1};
	die("\"$name1\" is not recognizable") if(!defined $type1);
	$usedType{$type1} = 1;
	$type2 = $enzymeINDEX{$name2};
	die("\"$name2\" is not recognizable") if(!defined $type2);
	$usedType{$type2} = 1;
	if($file_type eq "edge"){
		&readTable(\@edges, \@edgeMustFields, $type1, $type2, $file);
	}
	else{ # "glue"
		&readTable(\@glues, \@glueMustFields, $type1, $type2, $file);
	}
}
if(!defined $outdir){
	$outdir = ".";
}
my $glueMap = &numberGlues(\@glues);
&renumberNodes(\@edges, [keys %usedType], $glueMap, \@glues);

&writeTable(\@edges, ["id_str", "nglues", "overlap", "intercept", "slope", "type1", "id1", "length1", "type2", "id2", "length2", "glues"], "$outdir/edges.tsv");
&writeTable(\@glues, ["ID", "QryContigID", "QryLen", "type1", "RefContigID.x", "RefLen.x", "type2", "RefContigID.y", "RefLen.y", "Orientation.x", "Orientation.y", "Overlap", "Confidence.x", "Confidence.y", "QryStartPos.x", "QryEndPos.x", "RefStartPos.x", "RefEndPos.x", "QryStartPos.y", "QryEndPos.y", "RefStartPos.y", "RefEndPos.y", "OriContigID.x", "offset.x", "len.x", "OriContigID.y", "offset.y", "len.y"], "$outdir/glues.tsv");

exit(0);

sub readTable
{
	my ($table, $mustFields, $type1, $type2, $infile) = @_;
	open(TSV, "$infile") or die("Can not open \"$infile\" for reading\n");
	my $line;
	my @fields = ();
	my @indices = ();
	my ($nfields, $i, $idx);
	my @columns;
	my $row;
	while($line = <TSV>){
		chomp($line);
		if(scalar(@fields) == 0){
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
		($row->{type1}, $row->{type2}) = ($type1, $type2);
		push(@{$table}, $row);
	}
	close TSV;

	die("**ERROR: some of the required fields is missing in \"$infile\"\n") unless(scalar(@indices) == scalar(@$mustFields));

	return $table;
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

sub numberGlues
{
	my ($table) = @_;
	my %glueMap = ();
	my $type;
	for(my $i=0; $i<scalar(@{$table}); $i++){
		$type = $table->[$i]->{type1} . "_" . $table->[$i]->{type2};
		$glueMap{$type}{$table->[$i]->{ID}} = ($i + 1);
		$table->[$i]->{ID} = ($i + 1);
	}

	return \%glueMap;
}

sub renumberNodes
{
	my ($table, $types, $glueMap, $glueTable) = @_;
	my @nodes = (undef) x (scalar(keys %enzymeINDEX) + 1);
	my ($type1, $type2, $id1, $id2);
	my $row;
	foreach my $type (@$types){
		$nodes[$type] = {};
	}
# renumber IDs of nodes of each type
	for(my $i=0; $i<scalar(@{$table}); $i++){
		$row = $table->[$i];
		($type1, $type2) = ($row->{type1}, $row->{type2});
		($id1, $id2) = split(/_/, $row->{id_str});
		if(!defined $nodes[$type1]{$id1}){
			$nodes[$type1]{$id1} = scalar(keys %{$nodes[$type1]}) + 1;
		}
		if(!defined $nodes[$type2]{$id2}){
			$nodes[$type2]{$id2} = scalar(keys %{$nodes[$type2]}) + 1;
		}
	}
# solve ID conflicts between types
	my $offset;
	my $type = 0;
	my %finalType = ();
	for(my $k=1; $k<scalar(@nodes); $k++){
		if(defined $nodes[$k]){
			$finalType{$k} = ++$type;
			if(defined $offset){
				for my $id (keys %{$nodes[$k]}){
					$nodes[$k]{$id} += $offset;
				}
			}
			else{
				$offset = 0;
			}
			$offset += scalar(keys %{$nodes[$k]});
		}
	}
	my $id_str;
	my ($ori_id1, $ori_id2);
	for(my $i=0; $i<scalar(@{$table}); $i++){
		$row = $table->[$i];
		($type1, $type2) = ($row->{type1}, $row->{type2});
		my $type = "${type1}_${type2}";
		my @glues = split(/,/, $row->{glues});
		for(my $j=0; $j<@glues; $j++){
			$glues[$j] = $glueMap->{$type}{$glues[$j]};
		}
		$row->{glues} = join(",", @glues);
		($ori_id1, $ori_id2) = split(/_/, $row->{id_str});
		($id1, $id2) = ($nodes[$type1]{$ori_id1}, $nodes[$type2]{$ori_id2});
		$id_str = join("_", ($id1, $id2));
		$row->{id_str} = $id_str;
		$row->{type1} = $finalType{$type1};
		$row->{id1} = $ori_id1;
		$row->{type2} = $finalType{$type2};
		$row->{id2} = $ori_id2;
	}
	for(my $i=0; $i<scalar(@{$glueTable}); $i++){
		$row = $glueTable->[$i];
		$row->{type1} = $finalType{$row->{type1}};
		$row->{type2} = $finalType{$row->{type2}};
	}
}

sub writeTable
{
	my ($table, $fields, $outfile) = @_;

	open(OUT, ">$outfile") or die("Can not open \"$outfile\" for writing.\n");
	print OUT join("\t", @{$fields}) . "\n";
	my ($i, $j);
	my $no = scalar(@{$table});
	for($i=0; $i<$no; $i++){
		for($j=0; $j<scalar(@{$fields}); $j++){
			print OUT $table->[$i]->{$fields->[$j]};
			if($j != $#{ $fields }){
				print OUT "\t";
			} else {
				print OUT "\n";
			}
		}
	}
	close OUT;
}
