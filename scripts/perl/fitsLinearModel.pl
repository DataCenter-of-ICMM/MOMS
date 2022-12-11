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
use Statistics::LineFit;
use List::Util qw[sum min max];

BEGIN{
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script which fits linear model using pairing information of contigs from two single-enzyme assemblies
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-i, --input <str>     The TSV file containing the pairing information of contigs from different assemblies (REQUIRED)
	-n, --nglue <int>     The minimum number of glues for a contig pair to be adopted (default: 1)
	-c, --coef <float>    The minimum absolute value of correlation coefficient for a contig pair to be adopted (default: 0.9)
	-o, --output <str>    The path of the output table (REQUIRED)
	-h, --help            Help

Example:
	$program -i glue_of_pairs.tsv -o edges
USAGE

use vars qw($opt_i $opt_n $opt_c $opt_o $opt_h);
GetOptions( "i|input=s" => \$opt_i,
			"n|nglue=i" => \$opt_n,
			"c|coef=f" => \$opt_c,
			"o|output=s" => \$opt_o,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

die("**ERROR: -i option must be specified\n") unless(defined $opt_i);
die("**ERROR: -o option must be specified\n") unless(defined $opt_o);

my $min_glues = (defined $opt_n) ? $opt_n : 1;
my $min_coef = (defined $opt_c) ? $opt_c : 0.9;

my ($outdir, $outpre);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = ($opt_o eq "") ? "/" : $opt_o;
	$outpre = basename($opt_i); $outpre =~ s/\.tsv$//;
}
else{
	$outdir = dirname($opt_o);
	$outpre = basename($opt_o);
}

my ($cmd, $retCode);
$cmd = "mkdir -p $outdir";
$retCode = system($cmd);
die("**ERROR: Can not create directory \"$outdir\"\n") if($retCode != 0);

my $glueTable = &readAssociation($opt_i);
&appendIDstr($glueTable);
&filterRowByGlues($glueTable, $min_glues);

my ($fitTable, $fitFields) = &makeFitTable($glueTable, $min_coef);

# output the processed table
my $outfile  = "$outdir/$outpre.tsv";
&outputTable($fitTable, $fitFields, $outfile);

exit 0;

sub readAssociation
{
	my ($infile) = @_;
	my @table = ();
	my @mustFields = ("ID", "QryContigID", "QryLen", "RefContigID.x", "RefLen.x", "RefContigID.y", "RefLen.y", "Orientation.x", "Orientation.y", "Confidence.x", "Confidence.y", "QryStartPos.x", "QryEndPos.x", "RefStartPos.x", "RefEndPos.x", "QryStartPos.y", "QryEndPos.y", "RefStartPos.y", "RefEndPos.y", "Overlap");
	my $in;
	if($infile eq "-"){
		$in = \*STDIN;
	} else {
		open(TSV, "$infile") or die("Can not open \"$infile\" for reading\n");
		$in = \*TSV;
	}
	my $line;
	my @fields = ();
	my @indices = ();
	my ($nfields, $i, $idx);
	my @columns;
	my $no = 0;
	my $row;
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
		next if(scalar(@columns) != $nfields);
		$row = undef;
		for($i=0; $i<scalar(@indices); $i++){
			$idx = $indices[$i];
			$row->{$fields[$idx]} = $columns[$idx];
		}
		$table[$no++] = $row;
	}

	close $in;

	die("**ERROR: some of the required field is missing in \"$infile\"\n") unless(scalar(@indices) == scalar(@mustFields));
	return \@table;
}

sub appendIDstr
{
	my ($table) = @_;

	for(my $i=0; $i<scalar(@{$table}); $i++){
		$table->[$i]->{"id_str"} = $table->[$i]->{"RefContigID.x"} . "_" . $table->[$i]->{"RefContigID.y"};
	}
}

sub filterRowByGlues
{
	my ($table, $min_glues) = @_;

	unless($min_glues > 1){
		return;
	}

	my $no = scalar(@{$table});
	my %nGlues = ();
	my $id_str;
	for(my $i=0; $i<$no; $i++){
		$id_str = $table->[$i]->{"id_str"};
		if(defined $nGlues{$id_str}){
			$nGlues{$id_str}++;
		}
		else{
			$nGlues{$id_str} = 1;
		}
	}
	my $k = 0;
	for(my $i=0; $i<$no; $i++){
		$id_str = $table->[$i]->{"id_str"};
		next if($nGlues{$id_str} < $min_glues);
		$table->[$k++] = $table->[$i];
	}
	while($k < $no){
		pop @{$table};
		$k++;
	}
}

sub makeFitTable
{
	my ($table, $min_coef) = @_;
	my @fitTable = ();
	my %offsets = ();

	my $refLengths = &getRefLengths($table);

# calculate relative direction and relative offset of contig pairs
	my ($rel_dir, $rel_offset, $offset_x, $offset_y);

	my $idx;
	my $id_str;
	my $no = scalar(@{$table});
	for(my $i=0; $i<$no; $i++){
		$id_str = $table->[$i]->{"id_str"};
		($rel_dir, $rel_offset, $offset_x, $offset_y) = &calculateOffset($table->[$i]);
		if(!defined $offsets{$id_str}){
			$offsets{$id_str}{x} = [];
			$offsets{$id_str}{y} = [];
			$offsets{$id_str}{w} = [];
			$offsets{$id_str}{i} = [];
			$offsets{$id_str}{ex} = sprintf("%.1f", $table->[$i]->{"RefLen.x"} / 2);
			$offsets{$id_str}{ey} = sprintf("%.1f", $table->[$i]->{"RefLen.y"} / 2);
			$offsets{$id_str}{slope} = ($rel_dir eq '+') ? 1 : -1;
			$offsets{$id_str}{intercept} = $rel_offset;
		}
		else{
			$idx = findIndex($offsets{$id_str}{x}, $offset_x);
			if($idx >= 0){
				if($offsets{$id_str}{w}->[$idx] < $table->[$i]->{"Overlap"}){
					$offsets{$id_str}{x}->[$idx] = $offset_x;
					$offsets{$id_str}{y}->[$idx] = $offset_y;
					$offsets{$id_str}{w}->[$idx] = $table->[$i]->{"Overlap"};
					$offsets{$id_str}{i}->[$idx] = $table->[$i]->{"ID"};
				}
				next;
			}
			$idx = findIndex($offsets{$id_str}{y}, $offset_y);
			if($idx >= 0){
				if($offsets{$id_str}{w}->[$idx] < $table->[$i]->{"Overlap"}){
					$offsets{$id_str}{x}->[$idx] = $offset_x;
					$offsets{$id_str}{y}->[$idx] = $offset_y;
					$offsets{$id_str}{w}->[$idx] = $table->[$i]->{"Overlap"};
					$offsets{$id_str}{i}->[$idx] = $table->[$i]->{"ID"};
				}
				next;
			}
		}
		push(@{$offsets{$id_str}{x}}, $offset_x);
		push(@{$offsets{$id_str}{y}}, $offset_y);
		push(@{$offsets{$id_str}{w}}, $table->[$i]->{"Overlap"}); # the weights
		push(@{$offsets{$id_str}{i}}, $table->[$i]->{"ID"});
	}

# renumber the contig IDs
	my %inX = ();
	my %inY = ();
	my ($id_x, $id_y);
	my ($j, $k) = (0, 0);
	for $id_str (keys %offsets){
		($id_x, $id_y) = split(/_/, $id_str);
		if(!defined $inX{$id_x}){
			$inX{$id_x} = ++$j;
		}
		if(!defined $inY{$id_y}){
			$inY{$id_y} = ++$k;
		}
	}
	for $id_y (keys %inY){
		$inY{$id_y} += $j;
	}

# calculate fitting parameters of contig pairs
	$no = 0;
	my $lineFit = Statistics::LineFit->new();
	my ($coef, $intercept, $slope, $weight);
	while (my ($key, $value) = (each %offsets) ){
		$coef = ((scalar(@{$value->{x}}) > 1) ? &calculateCorrelation($value->{x}, $value->{y}) : "1.000000");
		next if(abs($coef) < $min_coef); # filter by correlation coefficient
		$value->{id_str} = $key;
		($id_x, $id_y) = split(/_/, $key);
		$value->{nglues} = scalar(@{$value->{x}});
		$value->{glues} = join(",", sort{$a <=> $b}(@{$value->{i}}));
		$weight = sum(@{$value->{w}});
		$value->{overlap} = $weight;
		$value->{coef} = $coef;
		if(scalar(@{$value->{x}}) > 1){
			$lineFit->setData($value->{x}, $value->{y}, $value->{w});
			($intercept, $slope) = $lineFit->coefficients();
		}
		else{
			($intercept, $slope) = ($value->{intercept}, $value->{slope});
		}
		$intercept += $value->{ey} - $slope * $value->{ex};
		$value->{intercept} = sprintf("%.1f", $intercept);
		$value->{slope} = sprintf("%.6f", $slope);
		$value->{id1} = $inX{$id_x};
		$value->{id2} = $inY{$id_y};
		$value->{length1} = $refLengths->[0]->{$id_x};
		$value->{length2} = $refLengths->[1]->{$id_y};
		$fitTable[$no++] = $value;
	}

	return (\@fitTable, ["id_str", "nglues", "overlap", "intercept", "slope", "coef", "id1", "length1", "id2", "length2", "glues"]);
}

sub findIndex
{
	my ($arr, $element) = @_;
	my $idx;
	for($idx=$#{$arr}; $idx>=0; $idx--){
		if($arr->[$idx] == $element){
			last;
		}
	}
	return $idx;
}

sub outputTable
{
	my ($table, $fields, $outfile) = @_;
	open(OUT, ">$outfile") or die("Can not open \"$outfile\" for writing\n");

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

sub getRefLengths
{
	my ($table) = @_;
	my @refLengths = ();
	my $row;
	my @ids;
	my @lengths;
	my $id;
	for(my $i=0; $i<scalar(@{$table}); $i++){
		$row = $table->[$i];
		@ids = ($row->{"RefContigID.x"}, $row->{"RefContigID.y"});
		@lengths = ($row->{"RefLen.x"}, $row->{"RefLen.y"});
		for(my $k=0; $k<2; $k++){
			$id = $ids[$k];
			if(!defined $refLengths[$k]->{$id}){
				$refLengths[$k]->{$id} = $lengths[$k];
			}
			else{
				if($refLengths[$k]->{$id} != $lengths[$k]){
					die("RefLength." . (($k==0) ? "x" : "y") . " of contig $id has conflicting values\n");
				}
			}
		}
	}
	return \@refLengths;
}

sub calculateOffset
{
	my ($row) = @_;

	my @orientation = ($row->{"Orientation.x"}, $row->{"Orientation.y"});
	my @qry_start = ($row->{"QryStartPos.x"}, $row->{"QryStartPos.y"});
	my @qry_end = ($row->{"QryEndPos.x"}, $row->{"QryEndPos.y"});
	my $qry_len = $row->{"QryLen"};
	my @ref_start = ($row->{"RefStartPos.x"}, $row->{"RefStartPos.y"});
	my @ref_end = ($row->{"RefEndPos.x"}, $row->{"RefEndPos.y"});
	my @ref_len = ($row->{"RefLen.x"}, $row->{"RefLen.y"});
	my $overlap_start = max(($orientation[0] eq '+') ? $qry_start[0] : $qry_end[0], ($orientation[1] eq '+') ? $qry_start[1] : $qry_end[1]);
	my $overlap_center = $overlap_start + ($row->{"Overlap"} - 1) / 2;
	my $overlap_offset = ($overlap_center - $qry_len/2);

	my @signs = ();
	my @offset = ();
	for(my $i=0; $i<2; $i++){
		$signs[$i] = ($orientation[$i] eq '+') ? 1 : -1;
		$offset[$i] = sprintf("%.1f", (($ref_start[$i] + $ref_end[$i] - $ref_len[$i]) - $signs[$i] * ($qry_start[$i] + $qry_end[$i] - $qry_len)) / 2 + $signs[$i] * $overlap_offset);
	}
	my $sign = $signs[0] * $signs[1];
	my $rel_offset = $offset[1] - $sign * $offset[0];
	my $rel_dir = ($sign == 1) ? '+' : '-';

	return ($rel_dir, $rel_offset, $offset[0], $offset[1]);
}

sub calculateCorrelation
{
	my ($offset_x, $offset_y) = @_;
	my ($mean_x, $mean_y) = (&mean($offset_x), &mean($offset_y));
	my $ssxx = &variance($offset_x, $offset_x, $mean_x, $mean_x);
	my $ssyy = &variance($offset_y, $offset_y, $mean_y, $mean_y);
	my $ssxy = &variance($offset_x, $offset_y, $mean_x, $mean_y);

	my $sign = $ssxy/abs($ssxy);
	my $correl = sprintf("%.6f", $sign * sqrt($ssxy * $ssxy / ($ssxx * $ssyy)));

	return $correl;
}

sub mean
{
	my ($elements) = @_;

	die unless (scalar(@{$elements}) > 0);

	my $sum = 0;
	for(my $i=$#{$elements}; $i>=0; $i--){
		$sum += $elements->[$i];
	}
	return ($sum / scalar(@{$elements}));
}

sub variance
{
	my ($elements1, $elements2, $mean1, $mean2) = @_;

	die unless (scalar(@{$elements1}) > 0 and  scalar(@{$elements2}) > 0);
	die unless (scalar(@{$elements1}) == scalar(@{$elements2}));
	$mean1 = &mean($elements1) if(!defined $mean1);
	$mean2 = &mean($elements2) if(!defined $mean2);
	my $sum = 0;
	for(my $i=$#{$elements1}; $i>=0; $i--){
		$sum += ($elements1->[$i] - $mean1) * ($elements2->[$i] - $mean2);
	}
	return $sum;
}
