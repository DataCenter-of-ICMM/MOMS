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
use POSIX qw/ceil/;

BEGIN{
	my $progpath = abs_path(dirname($0));
	unshift(@INC, $progpath);
	select(STDERR); $|= 1;
	select(STDOUT); $|= 1;
}

use BNG::Utility;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for filtering a CMAP file
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-c, --cmap <str>    The (CMAP) file for filtering (REQUIRED)
	-i, --ids <str>     A text file lists the IDs of contigs to be output (default: STDIN)
	-e, --exclude       Exclude the contigs with IDs listed in 'ids' (default: include)
	-o, --output <str>  The output path (REQUIRED)
	-h, --help          Help

Example:
	$program -c NGS.cmap -i ids.txt -o NGS_filtered
USAGE

use vars qw($opt_c $opt_i $opt_e $opt_o $opt_h);
GetOptions( "c|cmap=s" => \$opt_c,
			"i|ids=s" => \$opt_i,
			"e|exclude" => \$opt_e,
			"o|output=s" => \$opt_o,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

# check the required arguments
die("**ERROR: -c option must be specified\n") if(!defined $opt_c);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);

# get the arguments
my ($outdir, $outfile);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = ($opt_o eq "") ? "/" : $opt_o;
	$outfile = basename($opt_c) . "_filtered";
}
else{
	$outdir = dirname($opt_o);
	$outfile = basename($opt_o);
}
$outfile .= ".cmap";

# read in the cmap file as well as the ID list file
my ($ngs_cmap) = readCMap($opt_c);
print "Read in $opt_c completed with $ngs_cmap->{nContigs} cmaps.\n";

my $in;
if(!defined($opt_i)){
	$in = \*STDIN;
} else {
	open(IN, "$opt_i") or die("Can not open \"$opt_i\" for reading\n");
	$in = \*IN;
}
my %listed = ();
my $line;
while($line = <$in>){
	chomp($line);
	$listed{$line} = 1;
}
close $in if (defined $opt_i);

# output

my ($cmd, $retCode);
$cmd = "mkdir -p $outdir";
$retCode = system($cmd);
die("**ERROR: can not create directory \"$outdir\"\n") if($retCode != 0);

my $num = &writeCMap($ngs_cmap, "$outdir/$outfile", \%listed, !defined($opt_e));
print "Output to $outdir/$outfile with $num cmaps.\n";

exit 0;

sub writeCMap
{
	my ($cmap, $cmap_file, $list, $include) = @_;
	my @data_name = @{ $cmap->{dataName} };
	my @cmapIds = sort {$a <=> $b} keys %{ $cmap->{contigs} };
	my @headers = @{ $cmap->{headers} };
	open(CMF, ">$cmap_file") or die ("ERROR: Unable to write to file \"$cmap_file\"\n");
	# print headers
	print CMF "$headers[-2]\n";
	print CMF "$headers[-1]\n";

	# print cmaps
	my $num = 0;
	for (my $i=0; $i < scalar(@cmapIds); $i++) {
		my $c = $cmapIds[$i];
		next if(defined($list->{$c}) ^ $include);
		my $numsites = $cmap->{contigs}->{$c}->{$data_name[2]};
		for (my $j=0; $j <= $numsites; $j++) {
			print CMF $c . "\t" . $cmap->{contigs}->{$c}->{$data_name[1]} . "\t" . $numsites;
			for (my $m=3; $m <= $#data_name; $m++) {
				print CMF  "\t" . $cmap->{contigs}->{$c}->{$data_name[$m]}->[$j];
			}
			print CMF "\n";
		}
		$num++;
	}
	close CMF;

	return $num;
}
