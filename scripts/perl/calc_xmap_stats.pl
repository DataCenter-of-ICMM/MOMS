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
use File::Path;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case bundling);
use Scalar::Util qw(looks_like_number);
use Cwd qw(abs_path getcwd);
use POSIX qw(ceil);

my $dir=dirname($0);
require "$dir/msg.pm";

my $program=basename($0);
my $usage = << "USAGE";
$program: A perl script for calculate XMAP file statistics
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences
Usage: $program -q file -r file [options]

Required options:
  -q             The query files of XMAP format
  -r             The reference file of CMAP format

Optional options:
  -o <str>       The output file of separated by TAB (default: STDOUT)
  -h             Help

Example:
	$program -q file.xmap -r file.cmap -o xmap_stats.tsv

USAGE

use vars qw($opt_q $opt_r $opt_o $opt_f $opt_h);
GetOptions( "q=s" => \$opt_q,
			"r=s" => \$opt_r,
			"o=s" => \$opt_o,
			"h"   => \$opt_h) or exit;

die($usage) if($opt_h);
die($usage) unless( $opt_q and $opt_r );
select STDERR;

&dieMsg("File '$opt_q' does not exist.") unless( -s $opt_q );
&dieMsg("File '$opt_r' does not exist.") unless( -s $opt_r );

my $line;
my @columns;

my %cmapLengthRef = ();
&infoMsg($opt_r, 'Loading file', 32);
open(IN, '<', $opt_r) or die("Can not open '$opt_r' for reading.\n");
while( defined($line=<IN>) ){
	next if( $line =~ /^#/ );
	chomp($line);
	@columns = split(/\t/, $line);

	if( scalar(@columns)>4 and $columns[4]==0 ){
		$cmapLengthRef{$columns[0]} = {
			size				=> $columns[1],
			overlapRef 			=> [],
			sumCoverage			=> 0,
			nrCoverage			=> 0,
			nrPercentCoverage	=> 0
		};
	}
}
close IN;

&infoMsg($opt_q, 'Loading file', 32);
open(IN, '<', $opt_q) or die("Can not open '$opt_q' for reading.\n");
while( defined($line=<IN>) ){
	next if( $line =~ /^#/ );
	chomp($line);
	$line =~ s/^\s+//;
	$line =~ s/\s+/\t/g;
	@columns = split(/\t/, $line);
	my ($qId, $refId, $qStart, $qEnd, $start, $end) = @columns[1..6];
	&dieMsg("Cannot get reference contig '$refId' information") unless( $cmapLengthRef{$refId} );
	push(@{$cmapLengthRef{$refId}{overlapRef}}, {start => $start, end => $end});
}
close IN;

&infoMsg("Calculating coverage ... ");
foreach my $id (keys %cmapLengthRef){
	next unless(scalar(@{$cmapLengthRef{$id}{overlapRef}}));

	$cmapLengthRef{$id}{overlapRef} = &sortByCoord($cmapLengthRef{$id}{overlapRef}, 'start');
	my $oStart = $cmapLengthRef{$id}{overlapRef}[0]{start};
	my $oEnd = $cmapLengthRef{$id}{overlapRef}[0]{end};
	my $oSizeSum = 0;
	my $redundantSizeSum = $oEnd - $oStart;
	for( my $i=1; $i<scalar(@{$cmapLengthRef{$id}{overlapRef}}); $i++ ){
		my $iRef = $cmapLengthRef{$id}{overlapRef}[$i];
		$redundantSizeSum += ($iRef->{end} - $iRef->{start});
		if( $oStart <= $iRef->{start} and $iRef->{start} <= $oEnd ){
			$oEnd = ($oEnd < $iRef->{end}) ? $iRef->{end} : $oEnd;
		}
		else{
			$oSizeSum += ($oEnd - $oStart);
			$oStart = $iRef->{start};
			$oEnd = $iRef->{end};
		}
	}
	$oSizeSum += ($oEnd - $oStart);
	$cmapLengthRef{$id}{sumCoverage} = $redundantSizeSum;
	$cmapLengthRef{$id}{nrCoverage} = $oSizeSum;
	$cmapLengthRef{$id}{nrPercentCoverage} = ($oSizeSum / $cmapLengthRef{$id}{size} * 100);
}
print "done\n";

my $fh;
if( defined $opt_o ){
	open(OUT, '>', $opt_o) or &dieMsg("Can not open '$opt_o' for writing.");
	$fh = \*OUT;
}
else{
	$fh = \*STDOUT;
}

print $fh "refContig\trefContig_size\tquery_cov_bp\tquery_unique_cov_bp\tquery_unique_cov_percent\n";
foreach my $id ( keys %cmapLengthRef ){
	my $cRef = $cmapLengthRef{$id};
	my $sNrPercentCoverage = sprintf("%.2f", $cRef->{nrPercentCoverage});
	print $fh join("\t", $id, $cRef->{size}, $cRef->{sumCoverage}, $cRef->{nrCoverage}, $sNrPercentCoverage)."\n";
}

close $fh;
print "\n";
exit(0);

# =====================================================================================================================
sub sortByCoord
{
	my ($hRef, $coord) = @_;

	@$hRef = sort {$a->{$coord} <=> $b->{$coord}} @$hRef;
	return $hRef;
}
