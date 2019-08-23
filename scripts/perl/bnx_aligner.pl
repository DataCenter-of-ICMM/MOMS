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
use Parallel::ForkManager;
use Cwd 'abs_path';
use List::Util qw[min max];

BEGIN{
	my $progpath = dirname(abs_path($0));
	unshift(@INC, $progpath);
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

use BNG::Utility;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for aligning molecules in a BNX file with cmaps
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-r, --reference <str>   The reference BNX file for alignment (REQUIRED)
	-q, --query <str>       The query CMAP file for alignment (REQUIRED)
	-o, --output <str>      The output path (REQUIRED)
	-c, --channel <int>     The channel # (default: 1)
	-m, --min <int>,<int>   Minimum numbers of labels for REF and QRY respectively (default: 30,15)
	-t, --threads <int>     The number of threads for parallel processing (default: 8)
	-h, --help              Help

Example:
	$program -r molecules-bspqi.bnx -q bspqi.cmap -o alignment/bspqi -t 4
USAGE

use vars qw($opt_r $opt_q $opt_o $opt_c $opt_m $opt_t $opt_h);
GetOptions( "r|reference=s" => \$opt_r,
			"q|query=s" => \$opt_q,
			"o|output=s" => \$opt_o,
			"c|channel=i" => \$opt_c,
			"m|min=s" => \$opt_m,
			"t|threads=i" => \$opt_t,
			"h|help" => \$opt_h);
die($usage) if($opt_h);

# check the required arguments
die("**ERROR: -r option must be specified\n") if(!defined $opt_r);
die("**ERROR: -q option must be specified\n") if(!defined $opt_q);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);

my $refbnx = $opt_r;
my $qrycmap = $opt_q;
my ($outdir, $outpre);
if($opt_o =~ /\/$/ or  $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = $opt_o;
	$outpre = basename($opt_q);
	$outpre =~ s/\.cmap$//;
	$outpre .= "_v_" . basename($opt_r);
	$outpre =~ s/\.bnx$//;
}
else{
	$outdir = abs_path(dirname($opt_o));
	$outpre = basename($opt_o);
}
my $ichannel = (defined $opt_c && ($opt_c =~ /^\d+$/)) ? $opt_c : 1; 
my ($minref, $minqry) = (defined $opt_m) ? split(/,/, $opt_m) : (30, 15);
my $nthreads = (defined $opt_t && ($opt_t =~ /^\d+$/)) ? (($opt_t < 1) ? 1 : $opt_t) : 8;

my $cmd = "mkdir -p $outdir";
my $retCode = system($cmd);
die("**ERROR: Can not create directory \"$outdir\"\n") if($retCode != 0);

my ($bnxInfo, $bnxPositions) = readBNXPositions($opt_r, $minref);
my ($cmapInfo, $cmapPositions) = readCMAPPositions($opt_q, $minqry);

my $hits = alignBNXwithCMAP($bnxPositions, $cmapPositions, $nthreads);
processAlignments($bnxInfo, $bnxPositions, $cmapInfo, $cmapPositions, $hits);
outputAlignments($bnxInfo, $cmapInfo, $hits, $ichannel, "$outdir/$outpre-alignments.tsv");
collectEvidences($bnxInfo, $cmapInfo, $hits, $ichannel, "$outdir/$outpre-evidences.tsv");
exit 0;

sub readBNXPositions
{
	my ($infile, $minLabels) = @_;
	my %info = ();
	my %positions = ();
	open(IN, "<$infile") or die("Can not open \"$infile\" for reading\n");
	my $line;
	my ($idx_id, $idx_len, $idx_num, $idx_max);
	my ($id, $len, $num);
	my @columns;
	my @names = ();
	while($line = <IN>){
		chomp($line);
		if($line !~ /^\d\t/){
			if($line =~ /^#0h /){
				@names = split(/\t/, $line);
				$names[0] =~ s/^#0h //;
				$idx_id = getIndex(\@names, "MoleculeID");
				if(!defined $idx_id){
					print STDERR "wrong format of \"$infile\": MoleculeID is undefined";
					last;
				}
				$idx_len = getIndex(\@names, "Length");
				if(!defined $idx_len){
					print STDERR "wrong format of \"$infile\": Length is undefined";
					last;
				}
				$idx_num = getIndex(\@names, "NumberofLabels");
				if(!defined $idx_num){
					print STDERR "wrong format of \"$infile\": NumberofLabels is undefined";
					last;
				}
				$idx_max = max(max($idx_id, $idx_len), $idx_num);
			}
			next;
		}
		@columns = split(/\t/, $line);
		next if(scalar(@columns) < 2);
		if($columns[0] eq 0){
			if(scalar(@columns) <= $idx_max){
				print STDERR "irregular data found in \"$infile\": $line";
				last;
			}
			($id, $len, $num) = ($columns[$idx_id], $columns[$idx_len], $columns[$idx_num]);
			$info{$id} = {id=>$id, len=>$len, num=>$num};
		}
		elsif($columns[0] eq 1){
			if(scalar(@columns) != $num + 2){
				print STDERR "irregular data found in \"$infile\": $line";
				last;
			}
			if(defined $id){
				if($num >= $minLabels){
					$positions{$id} = [@columns[1..$num]];
				}
				$id = undef;
			}
		}
	}
	close IN;

	return (\%info, \%positions);
}

sub readCMAPPositions
{
	my ($infile, $minLabels) = @_;
	my %positions = ();
	open(IN, "<$infile") or die("Can not open \"$infile\" for reading\n");
	my $line;
	my ($idx_id, $idx_len, $idx_num, $idx_no, $idx_pos, $idx_max);
	my ($id, $len, $numSites, $siteNo, $pos);
	my @columns;
	my %info = ();
	my @names = ();
	while($line = <IN>){
		chomp($line);
		if($line !~ /^\d/){
			if($line =~ /^#h /){
				@names = split(/\t/, $line);
				$names[0] =~ s/^#h //;
				$idx_id = getIndex(\@names, "CMapId");
				if(!defined $idx_id){
					print STDERR "wrong format of \"$infile\": CMapId is undefined";
					last;
				}
				$idx_len = getIndex(\@names, "ContigLength");
				if(!defined $idx_len){
					print STDERR "wrong format of \"$infile\": ContigLength isundefined";
					last;
				}
				$idx_num = getIndex(\@names, "NumSites");
				if(!defined $idx_num){
					print STDERR "wrong format of \"$infile\": NumSites is undefined";
					last;
				}
				$idx_no = getIndex(\@names, "SiteID");
				if(!defined $idx_no){
					print STDERR "wrong format of \"$infile\": SiteID is undefined";
					last;
				}
				$idx_pos = getIndex(\@names, "Position");
				if(!defined $idx_pos){
					print STDERR "wrong format of \"$infile\": Position is undefined";
					last;
				}
				$idx_max = max($idx_len, max(max($idx_id, $idx_num), max($idx_no, $idx_pos)));
			}
			next;
		}
		@columns = split(/\t/, $line);
		if(scalar(@columns) <= $idx_max){
			print STDERR "irregular data found in \"$infile\": $line";
			last;
		}
		($id, $len, $numSites, $siteNo, $pos) = ($columns[$idx_id], $columns[$idx_len],  $columns[$idx_num], $columns[$idx_no], $columns[$idx_pos]);
		if(!defined $info{$id}){
			$info{$id} = {id=>$id, len=>$len, num=>$numSites};
		}
		next if($info{$id}->{num} < $minLabels);
		if(!defined $positions{$id}){
			$positions{$id} = [(0) x $numSites];
		}
		if($siteNo <= $numSites){
			$positions{$id}->[$siteNo-1] = $pos;
		}
	}
	close IN;

	return (\%info, \%positions);
}

sub getIndex
{
	my ($names, $name) = @_;
	for(my $i=0; $i<scalar(@$names); $i++){
		if($names->[$i] eq $name){
			return $i;
		}
	}
	return undef;
}

sub alignBNXwithCMAP
{
	my ($bnxPositions, $cmapPositions, $nthreads) = @_;

	my @hits = ();
	my $alignment;
	my $cnt = 0;
	my $ten = 0;
	my $hun = 0;
	my $bFoundInTen = 0;
	my $pm = Parallel::ForkManager->new(min($nthreads, scalar(keys %$bnxPositions)));
	$pm->run_on_finish(
		sub{
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
			my $localHits = $data_structure_reference->{result};
			my $refid = $data_structure_reference->{id};
			if(scalar(@$localHits) > 0){
				push(@hits, @{$localHits});
				$bFoundInTen = 1;
			}
			else{
				print STDERR "$refid\n";
			}
			$cnt++;
			if($cnt % 10 == 0){
				if($bFoundInTen){
					print "+";
					$bFoundInTen = 0;
				}
				else{
					print ".";
				}
				$ten++;
				if($ten % 10 == 0){
					$hun++;
					print (($hun % 10 == 0) ? "\n" : " ");
				}
			}
		}
	);
	foreach my $refid (keys %$bnxPositions){
		my $reflocs = $bnxPositions->{$refid};
		my $pid = $pm->start and next;
		# forked thread
		my @localHits = ();
		foreach my $qryid (keys %$cmapPositions){
			my $qrylocs = $cmapPositions->{$qryid};
			$alignment = localAlign($reflocs, $qrylocs);
			if(scalar(@{$alignment}) > 0){ # found significant overlap
				push(@localHits, {ref=>$refid, qry=>$qryid, alignment=>$alignment});
			}
		}
		$pm->finish(0, {result => \@localHits, id=>$refid});
	}
	$pm->wait_all_children;
	print "\n";
	return \@hits;
}

sub packedAlignDP
{
	my ($qrylocs, $reflocs) = @_;

	my $revqry = getReversePositions($qrylocs);
	my ($alignment1, $score1) = alignDP($qrylocs, $reflocs, 3);
	my ($alignment2, $score2) = alignDP($revqry, $reflocs, 3);
	if( scalar(@$alignment1) > scalar(@$alignment2) or 
		scalar(@$alignment1) == scalar(@$alignment2) and ($score1 < $score2) ){
		return $alignment1;
	}
	my $lastIdx = $#$revqry;
	foreach (@{$alignment2}){
		$_->[0] = $lastIdx - $_->[0];
	}
	return $alignment2;
}

sub combineAlign
{
	my ($qrylocs, $reflocs, $nSeed) = @_;

	my $m = scalar(@{$qrylocs});
	my @seeds;
	@seeds = @{$qrylocs}[0..($nSeed-1)];
	my ($alignment1, $score1) = alignDP(\@seeds, $reflocs, 1);

	@seeds = @{$qrylocs}[($m-$nSeed)..($m-1)];
	my ($alignment2, $score2) = alignDP(\@seeds, $reflocs, 2);
	my $offset = $m - $nSeed;
	foreach (@{$alignment2}){
		$_->[0] += $offset;
	}

	if(scalar(@{$alignment1}) > 0){
		if(scalar(@{$alignment2}) > 0){ # both aligned
			foreach (@{$alignment2}){
				if( $_->[0] > $alignment1->[-1]->[0] ||
					$_->[0] == $alignment1->[-1]->[0] &&
					$_->[1] > $alignment1->[-1]->[1] ){
					push(@$alignment1, $_);
				}
			}
			$score1 += $score2;
		}
		return ($alignment1, $score1);
	}
	if(scalar(@{$alignment2}) > 0){
		return ($alignment2, $score2);
	}

	return ([], 0);
}

sub packedCombineAlign
{
	my ($qrylocs, $reflocs, $nSeed) = @_;

	my $revqry = getReversePositions($qrylocs);
	my ($alignment1, $score1) = combineAlign($qrylocs, $reflocs, $nSeed);
	my ($alignment2, $score2) = combineAlign($revqry, $reflocs, $nSeed);
	if( scalar(@$alignment1) > scalar(@$alignment2) or 
		scalar(@$alignment1) == scalar(@$alignment2) and ($score1 < $score2) ){
		return $alignment1;
	}
	my $lastIdx = $#$revqry;
	foreach (@{$alignment2}){
		$_->[0] = $lastIdx - $_->[0];
	}
	return $alignment2;
}

sub localAlign
{
	my ($reflocs, $qrylocs) = @_;

	my $alignment;
	my ($m, $n) = (scalar(@{$qrylocs}), scalar(@{$reflocs}));
	my $nSeed = 14;
	if($m <=1 or $n <= 1){
		return [];
	}
	if($m <= $nSeed * 1.5){
		$alignment = packedAlignDP($qrylocs, $reflocs);
	}
	else{
		$alignment = packedCombineAlign($qrylocs, $reflocs, $nSeed);
	}
	return $alignment;
}

sub getReversePositions
{
	my ($positions) = @_;
	my @newpositions = ();
	my $len = $positions->[-1] + 20;
	for(my $i=$#$positions; $i>=0; $i--){
		push(@newpositions, sprintf("%.1d", $len - $positions->[$i]));
	}
	return \@newpositions;
}

sub processAlignments
{
	my ($refInfo, $refPositions, $qryInfo, $qryPositions, $hits, $nthreads) = @_;
	my ($qryid, $refid, $qrylen, $reflen, $qrycnt, $refcnt);
	my ($alignment, $Qpositions, $Rpositions, $adjust);
	my ($start, $end, $mid, $ori);
	my $pm = Parallel::ForkManager->new($nthreads);
	foreach (sort {$a->{qry} <=> $b->{qry}} @$hits){
		my $pid = $pm->start and next;
		($qryid, $refid) = ($_->{qry}, $_->{ref});
		$Qpositions = $qryPositions->{$qryid};
		$Rpositions = $refPositions->{$refid};
		my ($qryItem, $refItem) = ($qryInfo->{$qryid}, $refInfo->{$refid});
		($qrylen, $reflen, $qrycnt, $refcnt) = ($qryItem->{len}, $refItem->{len}, $qryItem->{num}, $refItem->{num});
		$alignment = $_->{alignment};
		$_->{end} = ($alignment->[0]->[0] + $alignment->[-1]->[0] < $qryItem->{num}) ? "5'" : "3'"; # aligned in which end
		$_->{dir} = (($alignment->[0]->[0] <=$alignment->[-1]->[0]) ? "+" : "-");
# to calculate the region that can be used for retrieving linkage information
		if(($_->{dir} eq "+") ^ ($_->{end} eq "3'")){
			$adjust = ($_->{dir} eq "+") ? $Qpositions->[$alignment->[0]->[0]] : $qrylen - $Qpositions->[$alignment->[0]->[0]];
			($start, $mid, $end, $ori) = (0, max($Rpositions->[$alignment->[0]->[1]] - $adjust, 0), $Rpositions->[$alignment->[-1]->[1]], "-");
		}
		else{
			$adjust = ($_->{dir} eq "+") ? ($qrylen - $Qpositions->[$alignment->[-1]->[0]]) : $Qpositions->[$alignment->[-1]->[0]];
			($start, $mid, $end, $ori) = ($Rpositions->[$alignment->[0]->[1]], min($Rpositions->[$alignment->[-1]->[1]] + $adjust, $refItem->{len}), $refItem->{len}, "+");
		}
		$_->{linkage} = {start=>$start, mid=>$mid, end=>$end, ori=>$ori};
		$pm->finish;
	}
	$pm->wait_all_children;
}

sub outputAlignments
{
	my ($refInfo, $qryInfo, $hits, $ichannel, $outfile) = @_;
	open(OUT, ">$outfile") or die("Can not open \"$outfile\" for writing.\n");
	print OUT "#no\tqryid\trefid\tqryLength\trefLength\tqrycnt\trefcnt\tchannel\torientation\t[S,M,T]\tend\talignment\n";
	my ($qryid, $refid, $qrylen, $reflen, $qrycnt, $refcnt);
	my $no = 0;
	foreach (sort {$a->{qry} <=> $b->{qry} || $a->{ref} <=> $b->{ref}} @$hits){
		$no++;
		($qryid, $refid) = ($_->{qry}, $_->{ref});
		my ($qryItem, $refItem) = ($qryInfo->{$qryid}, $refInfo->{$refid});
		($qrylen, $reflen, $qrycnt, $refcnt) = ($qryItem->{len}, $refItem->{len}, $qryItem->{num}, $refItem->{num});
		print OUT "$no\t$qryid\t$refid\t$qrylen\t$reflen\t$qrycnt\t$refcnt\t$ichannel\t";
		print OUT $_->{dir} . "\t";
		print OUT "[" . $_->{linkage}->{start} . ", ", $_->{linkage}->{mid}, ", ", $_->{linkage}->{end} . "]\t";
		print OUT $_->{end} ."\t";
		foreach (@{$_->{alignment}}){
			print OUT "(" . ($_->[0]+1) . "," . ($_->[1]+1) . ")";
		}
		print OUT "\n";
	}
	close OUT;
}

sub collectEvidences
{
	my($refInfo, $qryInfo, $hits, $ichannel, $outfile) = @_;
	open(OUT, ">$outfile") or die("Can not open \"$outfile\" for writing.\n");
	print OUT "#no\tctg1\tctg2\tend1\tend2\tmol\tdistance\n";
	my ($ctg1, $ctg2, $end1, $end2, $molid, $distance);
	my ($refid, $qryid1, $qryid2, $refItem, $qryItem1, $qryItem2, $alignment);
	my ($hit1, $hit2, $intersect);
	my $lastHit;
	my %clustered = ();
	foreach my $hit (@$hits) {
		$refid = $hit->{ref};
		if(!defined $clustered{$refid}){
			$clustered{$refid} = [];
		}
		push(@{$clustered{$refid}}, $hit);
	}
	my $no = 0;
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
				$intersect = intersect($hit1->{linkage}, $hit2->{linkage});
				next unless(%$intersect);
				$no++;
				$distance = distance($hit1->{linkage}, $hit2->{linkage});
				print OUT "$no\t$qryid1\t$qryid2\t$hit1->{end}\t$hit2->{end}\t$refid\t$distance\n";
			}
		}
	}
	close OUT;
}

sub intersect
{
	my ($rgn1, $rgn2) = @_;
	if($rgn1->{ori} eq $rgn2->{ori}){
		return {};
	}
	my $start = max($rgn1->{start}, $rgn2->{start});
	my $end = min($rgn1->{end}, $rgn2->{end});
	return ($start <= $end) ? {start=>$start, end=>$end} : {};
}

sub distance
{
	my ($rgn1, $rgn2) = @_;
	my ($ori1, $ori2) = ($rgn1->{ori}, $rgn2->{ori});
	if($ori1 eq $ori2){
		return undef;
	}
	return sprintf("%.1f", ($ori1 eq "+") ? ($rgn2->{mid} - $rgn1->{mid}) : ($rgn1->{mid} - $rgn2->{mid}));
}
