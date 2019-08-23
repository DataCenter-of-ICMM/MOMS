#!/usr/bin/env  perl
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
use List::Util qw[min max];

BEGIN{
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for collecting evidences given alignment information
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-i --input <str>    The file containing alingment information (REQUIRED)
		it's a <TAB> delimited file containing the following fields:
		1.no; 2.qryid; 3.refid; 4.qryStart; 5.qryEnd; 6.refStart; 7.refEnd; 8.orientation;
		9.qryLength; 10.refLength; 11.qrycnt; 12.refcnt; 13.channel; 44.alignment
	--min <int>     Minimum number of labels required (default: 5)
	--max <int>     Maximum number of labels allowed (default: none)
	-o, --output <str>  The output path (default: STDOUT)
	-h, --help          Help

Example:
	$program -i alignment/bspqi-alignment.tsv
USAGE

use vars qw($opt_i $opt_o $opt_m $opt_M $opt_h);
GetOptions( "i|input=s" => \$opt_i,
			"o|output=s" => \$opt_o,
			"min=i" => \$opt_m,
			"max=i" => \$opt_M,
			"h|help" => \$opt_h);
die($usage) if($opt_h);

# check the required arguments
die("**ERROR: -i option must be specified\n") if(!defined $opt_i);

my $minLabel = (defined $opt_m) ? $opt_m : 5;
my $maxLabel = (defined $opt_M) ? $opt_M : 9 ** 99;
my $infile = $opt_i;
my $outfile;
if(defined $opt_o){
	my ($outdir, $outpre);
	if($opt_o =~ /\/$/ or  $opt_o !~ /\//){
		$opt_o =~ s/\/$//;
		$outdir = $opt_o;
		$outpre = "evidences";
	}
	else{
		$outdir = abs_path(dirname($opt_o));
		$outpre = basename($opt_o);
	}
	my $cmd = "mkdir -p $outdir";
	my $retCode = system($cmd);
	die("**ERROR: Can not create directory \"$outdir\"\n") if($retCode != 0);
	$outfile = "$outdir/$outpre.tsv";
}

my $result = readAlignments($infile);
my $evidences = collectEvidences($result, $minLabel, $maxLabel);
outputEvidences($evidences, $outfile);

exit 0;

sub readAlignments
{
	my ($infile) = @_;
	my %qryInfo = ();
	my %refInfo = ();
	my @hits = ();
	my %result = (qry=>\%qryInfo, ref=>\%refInfo, hits=>\@hits);
	my $line;
	my @columns;
	open(IN, "<$infile") or die("Can not open \"$infile\" for reading.\n");
	while($line = <IN>){
		chomp($line);
		if($line =~ /^#/){
			next;
		}
		@columns = split(/\t/, $line);
		my ($qryid, $refid, $qryStart, $qryEnd, $refStart, $refEnd) = @columns[1..6];
		my ($dir, $qrylen, $reflen, $qrycnt, $refcnt, $iChannel, $alignstr) = @columns[7..13];
		if(!defined $qryInfo{$qryid}){
			$qryInfo{$qryid} = {len=>$qrylen, num=>$qrycnt};
		}
		if(!defined $refInfo{$refid}){
			$refInfo{$refid} = {len=>$reflen, num=>$refcnt};
		}
		my @alignment;
		$alignstr =~ s/^\(|\)$//g;
		my @pairs = split(/\)\(/, $alignstr);
		foreach (@pairs){
			my ($x, $y) = split(/,/, $_);
			push(@alignment, [$x-1, $y-1]);
		}
		my $end = (($alignment[0]->[0] + $alignment[-1]->[0]) < $qrycnt) ? '5': '3';
		my $ori = (($dir eq "+") ^ ($end eq "3")) ? '-' : '+';
		my $linkage = {qryStart=>$qryStart, qryEnd=>$qryEnd, refStart=>$refStart, refEnd=>$refEnd};

		my $hit = {ref=>$refid, qry=>$qryid, dir=>$dir, end=>$end, ori=>$ori, channel=>$iChannel, linkage=>$linkage, alignment=>\@alignment};
		push(@hits, $hit);
	}
	close IN;
	return \%result;
}

sub collectEvidences
{
	my ($result, $minLabel, $maxLabel) = @_;

	my @evidences = ();

	my $hits = $result->{hits};
	my $refInfo = $result->{ref};
	my $qryInfo = $result->{qry};

	my ($ctg1, $ctg2, $end1, $end2, $molid, $distance);
	my ($refid, $qryid1, $qryid2, $refItem, $qryItem1, $qryItem2, $alignment);
	my ($hit1, $hit2);
	my $lastHit;
	my %clustered = ();
	foreach my $hit (@$hits) {
		$refid = $hit->{ref};
		if(!defined $clustered{$refid}){
			$clustered{$refid} = [];
		}
		push(@{$clustered{$refid}}, $hit);
	}
	my $evi;
	foreach $refid (keys %clustered){
		$refItem = $refInfo->{$refid};
		my $clusteredHits = $clustered{$refid};
		for(my $i=0; $i<$#$clusteredHits; $i++){
			$hit1 = $clusteredHits->[$i];
			$qryid1 = $hit1->{qry};
			$qryItem1 = $qryInfo->{$qryid1};
			for(my $j=$i+1; $j<scalar(@$clusteredHits); $j++){
				$hit2 = $clusteredHits->[$j];
				$qryid2 = $hit2->{qry};
				$qryItem2 = $qryInfo->{$qryid2};
				$distance = calculateDistance($hit1->{linkage}, $hit1->{ori}, $hit2->{linkage}, $hit2->{ori});
				next unless(defined $distance);
				if($hit1->{ori} eq '+'){
					$distance -= ($hit1->{dir} eq '+') ? ($qryItem1->{len} - $hit1->{linkage}->{qryEnd}) : $hit1->{linkage}->{qryEnd};
					$distance -= ($hit2->{dir} eq '+') ? $hit2->{linkage}->{qryStart} : ($qryItem2->{len} - $hit2->{linkage}->{qryStart});
				}
				else{
					$distance -= ($hit2->{dir} eq '+') ? ($qryItem2->{len} - $hit2->{linkage}->{qryEnd}) : $hit2->{linkage}->{qryEnd};
					$distance -= ($hit1->{dir} eq '+') ? $hit1->{linkage}->{qryStart} : ($qryItem1->{len} - $hit1->{linkage}->{qryStart});
				}
				next if ($distance > 1000000 or $distance < -1000000);
				my $nlabel = min(scalar(@{$hit1->{alignment}}), scalar(@{$hit2->{alignment}}));
				next if($nlabel < $minLabel or $nlabel > $maxLabel);
				$evi = {id1=>$qryid1, id2=>$qryid2, end1=>$hit1->{end}, end2=>$hit2->{end}, mid=>$refid, channel=>$hit1->{channel}, nlabels=>$nlabel, distance=>sprintf("%.1f", $distance)};
				push(@evidences, $evi);
			}
		}
	}
	return \@evidences;
}

sub calculateDistance
{
	my ($linkage1, $ori1, $linkage2, $ori2) = @_;
	if($ori1 eq $ori2){
		return undef;
	}
	my $dist = ($ori1 eq '+') ? ($linkage2->{refStart} - $linkage1->{refEnd}) : ($linkage1->{refStart} - $linkage2->{refEnd});
	my ($x1, $x2) = (abs($linkage1->{refEnd} - $linkage1->{refStart}), abs($linkage2->{refEnd} - $linkage2->{refStart}));
	my ($y1, $y2) = (abs($linkage1->{qryEnd} - $linkage1->{qryStart}), abs($linkage2->{qryEnd} - $linkage2->{qryStart}));
	if(abs($x1 - $x2) < 1e-6){
		return undef;
	}
	my $k = ($y2 - $y1) / ($x2 - $x1);
	my $b = $y1 - $k * $x1;
	$dist = sprintf("%.1f", $k * $dist + $b);
	return $dist;
}

sub outputEvidences
{
	my ($evidences, $outfile) = @_;
	my $out;
	if(defined $outfile){
		open(OUT, ">$outfile") or die("Can not open \"$outfile\" for writing.\n");
		$out = \*OUT;
	}
	else{
		$out = \*STDOUT;
	}
	print $out "#no\tctg1\tctg2\tend1\tend2\tmol\tminlabels\tchannel\tdistance\n";
	my $no = 0;
#	foreach my $evi (sort{ $a->{id1} <=> $b->{id1} || $a->{id2} <=> $b->{id2} || $a->{mid} <=> $b->{mid} } @$evidences){
	foreach my $evi (@$evidences){
		$no++;
		print $out "$no\t$evi->{id1}\t$evi->{id2}\t$evi->{end1}\'\t$evi->{end2}\'\t$evi->{mid}\t$evi->{nlabels}\t$evi->{channel}\t$evi->{distance}\n";
	}
	if(defined $outfile){
		close $out;
	}
}
