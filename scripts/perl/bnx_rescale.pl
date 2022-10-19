#!/usr/bin/env  perl
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
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for rescaling molecules in a BNX file given bpp
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-i, --input <str>   The input BNX file for rescaling (REQUIRED)
	-b, --bpp <float>   The detected bpp (bases per pixel) (REQUIRED)
	-o, --output <str>  The output path (default: STDOUT)
	-h, --help          Help

Example:
	$program -i molecules-bspqi.bnx -b 472.2 -o output/bspqi-rescaled
USAGE

use vars qw($opt_i $opt_b $opt_o $opt_h);
GetOptions( "i|input=s" => \$opt_i,
			"b|bpp=f" => \$opt_b,
			"o|output=s" => \$opt_o,
			"h|help" => \$opt_h);
die($usage) if($opt_h);

# check the required arguments
die("**ERROR: -i option must be specified\n") if(!defined $opt_i);
die("**ERROR: -b option must be specified\n") if(!defined $opt_b);

my $infile = $opt_i;
my $bpp = $opt_b;
my $outfile;
if(defined $opt_o){
	my ($outdir, $outpre);
	if($opt_o =~ /\/$/ or  $opt_o !~ /\//){
		$opt_o =~ s/\/$//;
		$outdir = $opt_o;
		$outpre = basename($opt_i);
		$outpre =~ s/\.bnx$//;
		$outpre .= "-rescaled.bnx";
	}
	else{
		$outdir = abs_path(dirname($opt_o));
		$outpre = basename($opt_o);
	}
	my $cmd = "mkdir -p $outdir";
	my $retCode = system($cmd);
	die("**ERROR: Can not create directory \"$outdir\"\n") if($retCode != 0);
	$outfile = "$outdir/$outpre.bnx";
}


my $scale_factor = ($bpp / 500.0);
processBNXfile($infile, $scale_factor, $outfile);
exit 0;

sub processBNXfile
{
	my ($infile, $scale_factor, $outfile) = @_;
	my $out;
	if(defined $outfile){
		open(OUT, ">$outfile") or die("Can not create $outfile for writing.\n");
		$out = \*OUT;
	}
	else{
		$out = \*STDOUT;
	}
	open(IN, "<$infile") or die("Can not open $infile for reading.\n");
	my @columns;
	my $line;
	while($line = <IN>){
		chomp($line);
		if($line =~ /^[01]/){
			@columns = split(/\t/, $line);
			if($line =~ /^0/){
				$columns[2] = sprintf("%.2f", $columns[2] * $scale_factor);
			}
			else{
				for(my $i=1; $i<scalar(@columns); $i++){
					$columns[$i] = sprintf("%.2f", $columns[$i] * $scale_factor);
				}
			}
			$line = join("\t", @columns);
		}
		print $out "$line\n";
	}
	close IN;
	close $out;
}
