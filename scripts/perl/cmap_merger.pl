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
use XML::Simple;
use POSIX qw/ceil/;
use List::Util qw(min max);
use File::Path qw(make_path);
use File::Copy;

BEGIN{
	my $progpath = abs_path(dirname($0));
	unshift(@INC, $progpath);
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

use BNG::refAlignerRun;
use BNG::Utility;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for merging a pair of CMAP files
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-f, --first <str>    The first (NGS) CMAP file for merging (REQUIRED)
	-s, --second <str>   The second (BN) CMAP file for merging (REQUIRED)
	-e, --errbin <str>   The binary file containing error statistics of the alignment (REQUIRED)
	-o, --output <str>   The output path (REQUIRED)
	-t, --threads <int>  The number of threads for parallel processing (default: by aligner)
	-h, --help           Help

Options for invoking \$bionano/binary/RefAligner:
	-x, --xml <str>      Parameter file in XML format (default: \$bionano/xml/mergeArguments.xml)

USAGE

use vars qw($opt_f $opt_s $opt_e $opt_o $opt_t $opt_x $opt_h);
GetOptions( "f|first=s" => \$opt_f,
			"s|second=s" => \$opt_s,
			"e|errbin=s" => \$opt_e,
			"o|output=s" => \$opt_o,
			"t|threads=i" => \$opt_t,
			"x|xml=s" => \$opt_x,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

# get the bionano tool path
my $bionano = abs_path(dirname($0)) . "/../bionano";
my $toolpath = "$bionano/binary";
die("**ERROR: \"$toolpath\" directory is not found\n") unless(-d $toolpath);
$toolpath = abs_path($toolpath);

# get the BNG aligner
my $refAligner = "$toolpath/RefAligner";
die("**ERROR: Can not find RefAligner at $toolpath\n") unless(-f $refAligner);

# check the required arguments
die("**ERROR: -f option must be specified\n") if(!defined $opt_f);
die("**ERROR: -s option must be specified\n") if(!defined $opt_s);
die("**ERROR: -e option must be specified\n") if(!defined $opt_e);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);

# get the arguments
my ($ngs_cmap_in, $bng_cmap_in) = ($opt_f, $opt_s);
my $optxml = (defined $opt_x) ? $opt_x : "$bionano/xml/mergeArguments.xml";
my $errbin_in = $opt_e;

die("**ERROR: the 1st cmap \"$ngs_cmap_in\" does not exist\n") unless(-f $ngs_cmap_in);
die("**ERROR: the 2nd cmap \"$bng_cmap_in\" does not exist\n") unless(-f $bng_cmap_in);
die("**ERROR: the errbin file \"$errbin_in\" does not exist\n") unless(-f $errbin_in);
die("**ERROR: the parameter file \"$optxml\" does not exist\n") unless(-f $optxml);

my ($outdir, $outpre);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = ($opt_o eq "") ? "/" : $opt_o;
	$outpre = basename($opt_e);
	$outpre =~ s/\.errbin$//;
}
else{
	$outdir = dirname($opt_o);
	$outpre = basename($opt_o);
	$outdir = "$outdir/$outpre";
}

my ($cmd, $params, $retCode);

$cmd = "mkdir -p $outdir";
$retCode = system($cmd);
die("**ERROR: can not create directory \"$outdir\"\n") if($retCode != 0);
$outdir = abs_path("$outdir");

my $XML = new XML::Simple(KeyAttr=>[]);
my $configs = $XML->XMLin($optxml);
my $paraDict;
$paraDict = &parseConfig($configs, "global");
my $maxthreads = $paraDict->{'maxthreads'}{val};
my $nthreads =  (defined $opt_t && ($opt_t =~ /^\d+$/)) ? (($opt_t < 1) ? 1 : $opt_t) : $maxthreads;
$params = makeParams($paraDict);
$paraDict = parseConfig($configs, "mergeNGS_BN");
$params .= " " . makeParams($paraDict);
$params =~ s/ -maxthreads \S+//;
$params .= " -maxthreads $nthreads";

#####################
# Step 0. Initiation.
#####################
my ($ngs_cmap_mid, $bng_cmap_mid);

# 0.1 Read in the BioNano Assembled Cmap file
my ($bng_cmap) = readCMap($bng_cmap_in);
print "Read BioNano contig $bng_cmap_in completed with $bng_cmap->{nContigs} cmaps.\n";
# 0.2 Read in the NGS nick-encoded Cmap file
my ($ngs_cmap) = readCMap($ngs_cmap_in);
print "Read NGS contig $ngs_cmap_in completed with $ngs_cmap->{nContigs} cmaps.\n";
# 0.3 Shift BioNano ContigId by offset to prevent ContigID collision,
my $id_shift = 10 ** (ceil(log(max(keys $ngs_cmap->{contigs}) + 1000) / log(10)));
shiftCMapIds($bng_cmap, $id_shift);
$bng_cmap_mid = basename($bng_cmap_in);
$bng_cmap_mid =~ s/\.cmap$/_idshift.cmap/;
$bng_cmap_mid = "$outdir/$bng_cmap_mid";
writeCMapFile($bng_cmap, $bng_cmap_mid, 0);
print "Output BioNano contig with ID shift of $id_shift completed and output to $bng_cmap_mid.\n";

# 0.4 Separate the sequence entries into those with too few sites and those with normal number of sites
my $unMergeableMaxSites = 4;
my $mergeableMinSites = $unMergeableMaxSites + 1;
$ngs_cmap_mid = basename($ngs_cmap_in);
$ngs_cmap_mid =~ s/\.cmap$//;
$ngs_cmap_mid .= "_maxsites_$unMergeableMaxSites.cmap";
$ngs_cmap_mid = "$outdir/$ngs_cmap_mid";
writeCMapFile($ngs_cmap, $ngs_cmap_mid, 1, undef, $unMergeableMaxSites);

$ngs_cmap_mid = basename($ngs_cmap_in);
$ngs_cmap_mid =~ s/\.cmap$//;
$ngs_cmap_mid .= "_minsites_$mergeableMinSites.cmap";
$ngs_cmap_mid = "$outdir/$ngs_cmap_mid";
writeCMapFile($ngs_cmap, $ngs_cmap_mid, 1, $mergeableMinSites, undef);

# 0.5 end-mask to mark the entity of Bionano maps and sequences so that pairmerge does not merge between Bionano and Bionano, sequence with sequence
&maskCMap($bng_cmap, "0x8");
&maskCMap($ngs_cmap, "0x10");

my ($bng_cmap_mergeable, $ngs_cmap_mergeable) = ($bng_cmap_mid, $ngs_cmap_mid);
$bng_cmap_mergeable =~ s/\.cmap$/_masked.cmap/;
$ngs_cmap_mergeable =~ s/\.cmap$/_masked.cmap/;

writeCMapFile($bng_cmap, $bng_cmap_mergeable, 0);
writeCMapFile($ngs_cmap, $ngs_cmap_mergeable, 1, $mergeableMinSites, undef);

print "Step 0. Initiation Completed\n";

#####################
# Step 1. Merge
#####################
print "Step 1. Pair Merge Repeat between NGS\/BioNano started.\n";

my $mrg_output_prefix = "$outdir/Mrg";
$params = "-o $mrg_output_prefix -stdout -stderr -i $bng_cmap_mergeable -i $ngs_cmap_mergeable -readparameters $errbin_in $params";

$cmd = "$refAligner $params";
print "$cmd\n"; 
$retCode = system($cmd);
die("**ERROR: 1.0 In performing pairmerge $cmd.\n") if($retCode != 0);

my ($numMrg, $this_mrg_pairs) = parsingFastMrgStdOut("$mrg_output_prefix.stdout");
die("ERROR: 1.1. There are no merge possible between NGS and BioNano map. Please loosen parameters (e.g. -T or -pairmerge) and check input files.\n") if(!$numMrg);
# keep record of which pairs of fragments were merged
writeAllFastMrgPairs("$outdir/step1.merge.pairs.txt", $this_mrg_pairs, $numMrg);
print "1. We successfully merged $numMrg pairs.\n";

# now figure out which contigs are hybrid scaffolds, and sequences, and BioNano
my ($usedBioNanoRef, $usedSeqRef, $hybridRef) = determineParticipants($this_mrg_pairs, $numMrg, $id_shift);

# now iterates the output directory to figure out which file is hybrid, sequence leftover or bionano leftover
my ($hybridCmapsRef, $seqLeftOverCmapsRef, $bioNanoLeftOverCmapsRef) = categorizeCmapFiles($outdir, "Mrg_contig", $id_shift, $usedBioNanoRef, $usedSeqRef, $hybridRef);
my ($numHybridFiles, $numSeqLeftOverFiles, $numBioNanoLeftOverFiles) = (scalar(keys %$hybridCmapsRef), scalar(keys %$seqLeftOverCmapsRef), scalar(keys %$bioNanoLeftOverCmapsRef));
print "1.1 There are $numHybridFiles hybrid scaffold cmap files, $numSeqLeftOverFiles sequence left over cmap files, $numBioNanoLeftOverFiles BioNano left over cmap files.\n";

my $quickConcatCmap;
my ($outResults, $errResults, $job_status);
# now make a single hybrid cmap, a single sequence left over and a single BioNano left over cmap files
print "1.2 Concatenating into single hybrid, and single left over cmap files now.\n";
# print a test file
printIdFile("$outdir/step2.hybrid.id.txt", $hybridCmapsRef);
printIdFile("$outdir/step1.BN.id.txt", $bioNanoLeftOverCmapsRef);
printIdFile("$outdir/step1.NGS.id.txt", $seqLeftOverCmapsRef);

$quickConcatCmap = new BNG::refAlignerRun({
	binPath		=>	$refAligner,
	"if"		=>	"$outdir/step1.BN.id.txt",
	o			=>	"$outdir/step1.BN.naive",
	merge		=>	1,
	stdout		=>	1,
	f			=>	1
});
$cmd = $quickConcatCmap->getCMD();
print "$cmd\n";
($outResults, $errResults, $job_status) = $quickConcatCmap->runCMD();
if ($job_status != 0)	{
	# print out error
	$errResults = "ERROR: 1.2 In concatenating a single step1.BN.naive, error was encountered with exit code=$job_status; out info: $outResults; error info: $errResults";
	die($errResults);
}

$quickConcatCmap = new BNG::refAlignerRun({
	binPath		=>	$refAligner,
	"if"		=>	"$outdir/step1.NGS.id.txt",
	o			=>	"$outdir/step1.NGS.naive",
	merge		=>	1,
	stdout		=>	1,
	f			=>	1
});
$cmd = $quickConcatCmap->getCMD();
print "$cmd\n";
($outResults, $errResults, $job_status) = $quickConcatCmap->runCMD();
if ($job_status != 0)	{
	# print out error
	$errResults = "ERROR: 1.2 In concatenating a single step1.NGS.naive, error was encountered with exit code=$job_status; out info: $outResults; error info: $errResults";
	die($errResults);
}

$quickConcatCmap = new BNG::refAlignerRun({
	binPath		=>	$refAligner,
	"if"		=>	"$outdir/step2.hybrid.id.txt",
	o			=>	"$outdir/step2.hybrid",
	merge		=>	1,
	stdout		=>	1,
	f			=>	1
});
$cmd = $quickConcatCmap->getCMD();
print "$cmd\n";
($outResults, $errResults, $job_status) = $quickConcatCmap->runCMD();
if ($job_status != 0)	{
	# print out error
	$errResults = "ERROR: 1.2 In concatenating a single hybrid step2.hybrid, error was encountered with exit code=$job_status; out info: $outResults; error info: $errResults";
	die($errResults);
}

print "Step 1.2 file concatenation completed successfully\n";

print "$0 finished successfully.\n\n";

#####################################
# subroutines
sub printIdFile
{
	my ($file, $idsRef) = @_;
	open(OUT, ">$file") or die("ERROR: printIdFile: cannot write to $file: $!\n");
	foreach my $id (sort { $a <=> $b} keys %$idsRef)	{
		print OUT "$idsRef->{$id}\n";	# print the file name
	}
	close OUT;
} 

sub categorizeCmapFiles	{
	my ($dir, $filePrefix, $idShift, $usedBioNanoRef, $usedSeqRef, $hybridRef) = @_;

	my %hybridCmaps = ();
	my %seqLeftOverCmaps = ();
	my %bioNanoLeftOverCmaps = ();
	opendir(DIR, $dir) or die("ERROR: categorizeCmapFiles: cannot open dir $dir: $!\n");
	my @cmapFiles = grep {$_ =~ /^$filePrefix\d+\.cmap$/i} readdir DIR;
	closedir DIR;

	for (my $i = 0; $i < scalar(@cmapFiles); $i += 1)	{
		my $theId = $cmapFiles[$i];	$theId =~ s/^$filePrefix//;	$theId =~ s/\.cmap$//;
		if (exists $hybridRef->{$theId})	{
			# this file is a hybrid
			$hybridCmaps{$theId} = "$dir/$cmapFiles[$i]";
		} else	{
			# this file is a left over, check whether it is a sequence or a BioNano
			if ($theId < $idShift)	{
				$seqLeftOverCmaps{$theId} = "$dir/$cmapFiles[$i]";
			} else	{
				$bioNanoLeftOverCmaps{$theId} = "$dir/$cmapFiles[$i]";
			}
		}
	}
	return (\%hybridCmaps, \%seqLeftOverCmaps, \%bioNanoLeftOverCmaps);
}

sub determineParticipants
{
	my ($mrgPairsRef, $numMrg, $idShift) = @_;

	my %usedBioNano = ();
	my %usedSeq = ();
	my %hybrid = ();
	# stores hybrid ids, but would not store those that are not present in the output directory in the end, BECAUSE they have been merged with an entity either with a lower id or was totally encompassed

	for (my $i = 0; $i < $numMrg; $i += 1)	{
		my ($contig1, $contig2, $theHybrid) = ($mrgPairsRef->{ContigID1}[$i], $mrgPairsRef->{ContigID2}[$i], $mrgPairsRef->{ResultContigID}[$i]);
		if ($contig1 < $idShift){ # a sequence
			$usedSeq{$contig1} = 1;
		} else{ # a genome map
			$usedBioNano{$contig1} = 1;
		}

		if ($contig2 < $idShift){ # a sequence
			$usedSeq{$contig2} = 1;
		} else{ # a genome map
			$usedBioNano{$contig2} = 1;
		}

		$hybrid{$theHybrid} = 1;

		# now figure out if the particpant was a hybrid already AND that this merge resulted in an hybrid with a different id (either be a lower id, or that the participant was completely encompassed)
		delete $hybrid{$contig1} if (exists $hybrid{$contig1} && $theHybrid != $contig1);
		delete $hybrid{$contig2} if (exists $hybrid{$contig2} && $theHybrid != $contig2);
	}
	return (\%usedBioNano, \%usedSeq, \%hybrid);
}

# this subroutine marks the Bionano maps with a Mask bit, or marks the sequence with a different Mask bit 
sub maskCMap
{
	my ($cmap, $hexstr) = @_;
	my $tag = hex($hexstr);
	my $data_name = $cmap->{dataName};
	my @result = grep {/^Mask$/} @{$data_name};
	my @cmapIds = sort {$a <=> $b} keys %{ $cmap->{contigs} };
	my $nsitename = $data_name->[2];
	if(scalar(@result) == 0){
		my $data_type = $cmap->{dataType};
		my $headers = $cmap->{headers};
		push(@{$data_name}, "Mask");
		push(@{$data_type}, "Hex");
		for(my $i=0; $i<=$#{$headers}; $i++){
			if($headers->[$i] =~ /^#h /){
				$headers->[$i] =~ s/$/\tMask/;
			}
			elsif($headers->[$i] =~ /^#f /){
				$headers->[$i] =~ s/$/\tHex/;
			}
		}
		my $new_hexstr = sprintf("%x", $tag);
		for(my $i=0; $i<scalar(@cmapIds); $i++){
			my $ctg = $cmap->{contigs}->{$cmapIds[$i]};
			my $numsites = $ctg->{$nsitename};
			my @arr = (0) x ($numsites+1);
			$ctg->{"Mask"} = \@arr;
			$ctg->{"Mask"}->[0] = $ctg->{"Mask"}->[$numsites] = $new_hexstr;
		}
	}
	else{
		my $new_tag;
		for(my $i=0; $i<scalar(@cmapIds); $i++){
			my $ctg = $cmap->{contigs}->{$cmapIds[$i]};
			my $numsites = $ctg->{$nsitename};
			$new_tag = $tag | hex("0x" . $ctg->{"Mask"}->[0]);
			$ctg->{"Mask"}->[0] = sprintf("%x", $new_tag);
			$new_tag = $tag | hex("0x" . $ctg->{"Mask"}->[$numsites]);
			$ctg->{"Mask"}->[$numsites] = sprintf("%x", $new_tag);
		}
	}
}
