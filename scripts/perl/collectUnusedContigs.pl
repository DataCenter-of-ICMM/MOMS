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
use Cwd qw(abs_path);
use List::Util qw(min);
use POSIX;

BEGIN{
	my $progpath = abs_path(dirname(abs_path($0)));
	unshift(@INC, $progpath);
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

my $program = basename($0);
my $usage = << "USAGE";
$program A perl script for outputing sequences of unused contigs given IDs of used contigs
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]:
	-i, input <str>   The input FASTA file from NGS assembly (REQUIRED)
	-m, --map <str>   The file containing the mapping information from ID to sequence name (REQUIRED)
	-t, --trans <str> The file containing the coordinate transformation information when resolving conflicts (REQUIRED)
	-o, output <str>  The output FASTA file (REQUIRED)
	-u, used <str>    The IDs for used contigs (default: NONE)
	-minlen <int>     The minimum allowed length for outputing (default: 0)
	-maxlen <int>     The maximum allowed length for outputing (default: ~0)
	-h, help          Help
USAGE

use vars qw($opt_i $opt_m $opt_t $opt_o $opt_u $opt_min $opt_max $opt_h);
if(!GetOptions( "i|input=s" => \$opt_i,
	"m|map=s" => \$opt_m,
	"t|trans=s" => \$opt_t,
	"o|output=s" => \$opt_o,
	"u|used=s" => \$opt_u,
	"minlen=i" => \$opt_min,
	"maxlen=i" => \$opt_max,
	"h|help" => \$opt_h)){
	print ("Please try -h for more details\n");
	exit 1;
}
die($usage) if($opt_h);

# check the required arguments
die("**ERROR: -i option must be specified\n") if(!defined $opt_i);
die("**ERROR: -m option must be specified\n") if(!defined $opt_m);
die("**ERROR: -t option must be specified\n") if(!defined $opt_t);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);
my $minlen = (defined $opt_min) ? $opt_min : 0;
my $maxlen = (defined $opt_max) ? $opt_max : ~0;

# set the input files
my $fasta_file = $opt_i;
my $ngs_namemap_file = $opt_m;
my $cut_coord_file = $opt_t;

# set the output directory and prefix
my ($outdir, $outpre);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = $opt_o;
	$outpre = basename($opt_i);
	$outpre =~ s/\.[^.]+$//;
}
else{
	$outdir = dirname($opt_o);
	$outpre = basename($opt_o);
}
my ($cmd, $retCode);
$cmd = "mkdir -p $outdir";
$retCode = system($cmd);
die("**ERROR: Can not create directory \"$outdir\"\n") if($retCode != 0);

my $seqTable = readFasta($ngs_namemap_file, $fasta_file);
my $contigInfo = readContigInfo($seqTable, $cut_coord_file);
my $usedIds = readIds($opt_u);
$usedIds = expandIds($usedIds, $contigInfo);
$contigInfo = filterContig($contigInfo, $usedIds, $minlen, $maxlen);

printFasta($contigInfo, $seqTable, "$outdir/$outpre.fa");
printFasta($contigInfo, $seqTable, "$outdir/$outpre-small.fa", -1);
printFasta($contigInfo, $seqTable, "$outdir/$outpre-large.fa", 1);

exit 0;

sub readIds
{
	my ($infile) = @_;
	my $in;
	if(!defined $infile){
		return {};
	}
	if($infile ne "-"){
		open(IN, "$infile") or die("Can not open \"$infile\" for reading\n");
		$in = \*IN;
	}
	else{
		$in = \*STDIN;
	}
	my %usedIDs = ();
	my $line;
	while($line = <$in>){
		chomp($line);
		if($line =~ /^\d+$/){
			$usedIDs{$line} = 1;
		}
	}
	close $in if($infile ne "-");

	return \%usedIDs;
}

sub expandIds
{
	my ($usedIDs, $contigInfo) = @_;
	my %added = ();
	for my $id (keys %$usedIDs){
		while ( my ($sibling_id, $sibling) = each (%{$contigInfo->{$id}->{sibling}}) ){
			next if(defined $usedIDs->{$sibling_id});
			next unless ($sibling->{overlap} / $contigInfo->{$sibling_id}->{len} > 0.75);
			$added{$sibling_id} = 1;
		}
	}
	for my $id (keys %added){
		$usedIDs->{$id} = 1;
	}

	return $usedIDs;
}

sub readFasta
{
	my ($namemap_file, $fasta_file) = @_;
	my %seqTable = ();
	my ($seqid, $name);
	my %rev_map = (); # temporary map for retrieving ID given NAME
	my $line;
	open(NM, "<$namemap_file") or die("Can not open \"$namemap_file\" for reading\n");
	while($line = <NM>){
		chomp($line);
		next if($line !~ /^\d/);
		my @tokens = split(/\t/, $line);
		next if(scalar(@tokens) < 3);
		$seqid = $tokens[0];
		$name = (split(/\s+/, $tokens[1]))[0];
		$seqTable{$seqid}{name} = $name;
		$seqTable{$seqid}{len} = $tokens[2];
		$rev_map{$name} = $seqid;
	}
	close NM;
	my $seq = "";
	open(FA, "<$fasta_file") or die("Can not open \"$fasta_file\" for reading\n");
	while($line = <FA>){
		chomp($line);
		if($line =~ /^>(\S+)/){
			if($seq ne ""){
				$seqid = $rev_map{$name};
				if(defined $seqid){
					$seqTable{$seqid}{seq} = $seq;
				}
				else{
					warn "Warning: $name is not recognized to be a valid NAME, check $namemap_file or $fasta_file\n";
				}
				$seq = "";
			}
			$name = $1;
			next;
		}
		$seq .= $line;
	}
	close FA;
	if($seq ne ""){
		$seqid = $rev_map{$name};
		if(defined $seqid){
			$seqTable{$seqid}{seq} = $seq;
		}
		else{
			warn "Warning: $name is not recognized to be a valid NAME, check $namemap_file or $fasta_file\n";
		}
	}
	return \%seqTable;
}

sub readContigInfo
{
	my ($seqTable, $coord_file) = @_;
	my $line;
	my %contigInfo = ();
	my %subseqs = ();
	my ($id, $name, $newId, $start, $end, $len);
	open(CR, "<$coord_file") or die ("Can not open \"$coord_file\" for reading\n");
	while($line = <CR>){
		chomp($line);
		next if($line =~ /^oldId/); # skipping header
		my @tokens = split(/\t/, $line);
		next if(scalar(@tokens) < 6);
		($id, $newId) = ($tokens[0], $tokens[3]);
		($start, $end) = ($tokens[1], $tokens[2]);
		$len = $end - $start + 1;
		$contigInfo{$newId} = {seqid=>$id, start=>$start, end=>$end, len=>$len, length=>$seqTable->{$id}->{len}, sibling=>{}};
		if(!exists $subseqs{$id}){
			$subseqs{$id} = [];
		}
		push(@{$subseqs{$id}}, $newId);
	}
	close CR;

	foreach my $group (values %subseqs){
		next unless (scalar(@$group) > 1);
		my @sorted_group = sort { $contigInfo{$a}->{start} <=> $contigInfo{$b}->{start} || $contigInfo{$b}->{end} <=> $contigInfo{$a}->{end} } @{$group};
		my ($i, $j, $id1, $id2);
		my ($start1, $end1, $start2, $end2);
		my ($offset, $overlap);
		my ($sibling1, $sibling2);
		for($i=0; $i<$#sorted_group; $i++){
			$id1 = $sorted_group[$i];
			$sibling1 = $contigInfo{$id1}->{sibling};
			($start1, $end1) = ($contigInfo{$id1}->{start}, $contigInfo{$id1}->{end});
			for($j=$i+1; $j<scalar(@sorted_group); $j++){
				$id2 = $sorted_group[$j];
				$start2 = $contigInfo{$id2}->{start};
				last unless ($start2 <= $end1);
				# found an overlap
				$end2 = $contigInfo{$id2}->{end};
				($offset, $overlap) = ($start2 - $start1, min($end1, $end2) - $start2 + 1);
				$sibling2 = $contigInfo{$id2}->{sibling};
				$sibling1->{$id2} = {offset => $offset, overlap => $overlap};
				$sibling2->{$id1} = {offset => -$offset, overlap => $overlap};
			}
		}
	}
	return \%contigInfo;
}

sub filterContig
{
	my ($contigInfo, $usedIds, $minlen, $maxlen) = @_;
	my ($ndeleted, $nlarge, $nsmall, $nmiddle, $nused) = (0) x 5;
	printf("Statistics of contig filtering:\n");
	my $ntotal = scalar(keys %$contigInfo);
	my ($nsum_used_large, $nsum_used_small, $nsum_used_middle, $nsum_unused_large, $nsum_unused_small, $nsum_unused_middle) = (0) x 6;
	my ($n_used_large, $n_used_small, $n_used_middle, $n_unused_large, $n_unused_small, $n_unused_middle) = (0) x 6;
	foreach my $id (keys %$contigInfo){
		my $len = $contigInfo->{$id}->{len};
		if(defined $usedIds->{$id} or $len < $minlen or $len > $maxlen){
			$nused++ if(defined $usedIds->{$id});
			if($len > $maxlen){
				$nlarge++;
			}
			elsif($len >= $minlen){
				$nmiddle++;
			}
			else{
				$nsmall++;
			}
			$ndeleted++;
			if(defined $usedIds->{$id}){
				if($len < $minlen){
					$nsum_used_small += $len;
					$n_used_small++;
				}
				elsif($len > $maxlen){
					$nsum_used_large += $len;
					$n_used_large++;
				}
				else{
					$nsum_used_middle += $len;
					$n_used_middle++;
				}
				delete $contigInfo->{$id};
			}
			else{
				if($len < $minlen){
					$nsum_unused_small += $len;
					$n_unused_small++;
					$contigInfo->{$id}->{flag} = -1;
				}
				else{ # $len > $maxlen
					$nsum_unused_large += $len;
					$n_unused_large++;
					$contigInfo->{$id}->{flag} = 1;
				}
			}
		}
		else{
			$nsum_unused_middle += $len;
			$n_unused_middle++;
			$contigInfo->{$id}->{flag} = 0;
		}
	}
	printf("%8s contigs were filtered from %8s contigs (%8s contigs left), where\n", commify($ndeleted), commify($ntotal), commify($ntotal - $ndeleted));
	printf("%8s contigs were large contigs (> %sk bp)\n", commify($nlarge), commify(int($maxlen/1000))) if($maxlen < ~0);
	printf("%8s contigs were middle contigs\n", commify($nmiddle));
	printf("%8s contigs were small contigs (< %sk bp)\n", commify($nsmall), commify(ceil($minlen/1000))) if($minlen > 0);
	printf("%8s contigs were used contigs\n", commify($nused));
	printf("Statistics of total lengths:\n");
	printf("%15s bp for %8s used large contigs (> %sk bp).\n", commify($nsum_used_large), commify($n_used_large), commify(int($maxlen/1000))) if($maxlen < ~0);
	printf("%15s bp for %8s used middle contigs.\n", commify($nsum_used_middle), commify($n_used_middle));
	printf("%15s bp for %8s used small contigs (< %sk bp).\n", commify($nsum_used_small), commify($n_used_small), commify(ceil($minlen/1000))) if($minlen > 0);
	printf("%15s bp for %8s unused large contigs (> %sk bp).\n", commify($nsum_unused_large), commify($n_unused_large), commify(int($maxlen/1000))) if($maxlen < ~0);
	printf("%15s bp for %8s unused middle contigs.\n", commify($nsum_unused_middle), commify($n_unused_middle));
	printf("%15s bp for %8s unused small contigs (< %sk bp).\n", commify($nsum_unused_small), commify($n_unused_small), commify(ceil($minlen/1000))) if($minlen > 0);

	return $contigInfo;
}

sub commify
{
	my $text = reverse $_[0];
	$text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text;
}

sub printFasta
{
	my ($contigInfo, $seqTable, $outfile, $flag) = @_;
	$flag = 0 if(!defined $flag);
	open (OUT, ">$outfile") or die("Can not open \"$outfile\" for writing\n");
	my ($ctg, $seqid, $name, $start, $end, $seq, $len);
	for my $id (sort { $a <=> $b } keys %$contigInfo){
		$ctg = $contigInfo->{$id};
		next if($ctg->{flag} ne $flag);
		($seqid, $start, $end, $len) = ($ctg->{seqid}, $ctg->{start}, $ctg->{end}, $ctg->{len});
		$name = $seqTable->{$seqid}->{name};
		print OUT ">${name}_${start}_${end}|$id|$len\n";
		$seq = retrieveSequence($seqid, $seqTable, '+', $start, $len);
		print OUT "$seq\n";
	}
	close OUT;
}

sub retrieveSequence
{
	my ($seqid, $seqTable, $dir, $offset, $len) = @_;
	$len = $seqTable->{$seqid}{len} if(!defined $len);
	if($len <= 0){
		return "";
	}
	$dir = "+" if(!defined $dir);
	$offset = 0 if(!defined $offset or $offset < 0);
	my $seq = substr($seqTable->{$seqid}{seq}, $offset, $len);
	if($dir ne "+"){
		$seq = reverseComplement($seq);
	}
	return $seq;
}

sub reverseComplement
{
	my ($str) = @_;
	$str =~ tr/ACGTacgt/TGCAtgca/;
	return reverse($str);
}
