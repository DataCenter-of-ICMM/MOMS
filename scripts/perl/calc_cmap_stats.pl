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
use Getopt::Long qw(:config no_ignore_case bundling);
use Scalar::Util qw(looks_like_number);
use Cwd qw(abs_path);

my $dir=dirname($0);
require "$dir/msg.pm";

my $program=basename($0);
my $usage = << "USAGE";
$program: A perl script for calculate CMAP file statistics
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences
Usage: $program file.cmap [file2.cmap [...]]

Required options:
  file           One or more CMAP files

Optional options:
  -f             Output format for a table or list (default: table)
  -o <str>       The output file (default: STDOUT)
  -h             Help

Example:
	$program sample1.cmap -o cmap_stats.tsv

USAGE

use vars qw($opt_o $opt_f $opt_h);
GetOptions( "o=s" => \$opt_o,
			"f"   => \$opt_f,
			"h"   => \$opt_h) or exit;

die($usage) if($opt_h or scalar(@ARGV)==0);
select STDERR;

my $fmt = defined($opt_f);

foreach my $file ( @ARGV ){
	&dieMsg("File '$file' does not exist.") unless( -s $file );
}

my $fh;
if( defined($opt_o) ){
	open(OUT, '>', $opt_o) or &dieMsg("Can not open '$opt_o' for writing.");
	$fh = \*OUT;
}
else{
	$fh = \*STDOUT;
}

# header
if( !$fmt ){
	print $fh "file\tCount \tMin length (Mbp)\tMedian length (Mbp)\tMean length (Mbp)\tN50 length (Mbp)\tMax length (Mbp)\tTotal length (Mbp)\n";
}

my $null = 0;
foreach my $file ( @ARGV ){
	$file = abs_path($file);
	my $name = $file;
	$name =~ s/^.*\///;

	my ($idCol, $lenCol);
	my $line = `grep -m 1 '^#h' $file`;
	$line =~ s/^#h\s+//;
	chomp($line);
	my @columns = split(/\t/, $line);
	for( my $i=0; $i<scalar(@columns); $i++ ){
		if( $columns[$i] eq 'CMapId' ){
			$idCol = $i+1;
		}
		elsif( $columns[$i] eq 'ContigLength' ){
			$lenCol = $i+1;
		}
	}
	&dieMsg("File '$file' format is incorrect.") unless( $idCol and $lenCol );

	&infoMsg($name, 'Loading file', 32);
	open(IN, "grep -v '^#' $file | cut -f$idCol,$lenCol | uniq | cut -f2 | ") or die $!;
	my @len = <IN>;
	chomp(@len);
	close IN;

	my @ret = &stats(@len);

	if( $fmt ){
		print $fh "\n" if($null);
		$null = 1;
		print $fh "Maps Stats ($file):\n";
		print $fh "Count = $ret[0]\n";
		print $fh "Min length (Mbp) = $ret[1]\n";
		print $fh "Median length (Mbp) = $ret[2]\n";
		print $fh "Mean length (Mbp) = $ret[3]\n";
		print $fh "N50 length (Mbp) = $ret[4]\n";
		print $fh "Max length (Mbp) = $ret[5]\n";
		print $fh "Total length (Mbp) = $ret[6]\n";
	}
	else{
		print $fh "$name\t";
		print $fh join("\t", @ret) . "\n";
	}
}

close $fh;
print "\n";
exit(0);

# =============================================================================================================
sub stats
{
	my @len = sort {$b<=>$a} @_;

	my $qty = scalar(@len);
	if( $qty == 0 ){
		return ('NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA');
	}

	my $ratio = 10 ** 6;

	my $median = 0;
	my $m = int($qty/2);
	if( $qty%2 ){
		$median = $len[$m];
	}
	else{
		$median = ($len[$m-1] + $len[$m]) / 2;
	}
	$median = sprintf("%.3f", $median/$ratio);

	my $min = sprintf("%.3f", $len[-1]/$ratio);
	my $max = sprintf("%.3f", $len[0]/$ratio);

	my $total = eval join('+', @len);
	my $mean = sprintf("%.3f", $total/($qty*$ratio));

	my $half = $total / 2;
	my ($n50, $n50val) = (0, 0);
	for(my $i=0; $i<$qty; $i++){
		$n50 += $len[$i];
		if( $n50 >= $half ){
			$n50val = $len[$i];
			last;
		}
	}

	$total = sprintf("%.3f", $total/$ratio);
	$n50val = sprintf("%.3f", $n50val/$ratio);

	return ($qty, $min, $median, $mean, $n50val, $max, $total);
}
