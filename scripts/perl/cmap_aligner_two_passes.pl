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
use Cwd qw(abs_path);
use XML::Simple;
use List::Util qw[min max];

BEGIN{
	my $progpath = abs_path(dirname(abs_path($0)));
	unshift(@INC, $progpath);
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

use BNG::Utility;
#use BNG::refAlignerRun;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for aligning a pair of CMAP files with two passes
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-r, reference <str> The scaffold CMAP file for alignment (REQUIRED)
	-q, query <str>     The query NGS CMAP file for alignment (REQUIRED)
	-o, --output <str>  The output path (REQUIRED)
	-s, --skip          Skip the second phase (default: no)
	-t, --threads <int> The number of threads for parallel processing (default: 8)
	-h, --help          Help

Options for invoking \$bionano/binary/RefAligner:
	-x, --xml <str>     Parameter file in XML format (default: \$bionano/xml/alignArguments.xml)

USAGE

use vars qw($opt_r $opt_q $opt_o $opt_s $opt_t $opt_h $opt_x);
GetOptions( "r|reference=s" => \$opt_r,
			"q|query=s" => \$opt_q,
			"o|output=s" => \$opt_o,
			"s|skip" => \$opt_s,
			"t|threads=i" => \$opt_t,
			"x|xml=s" => \$opt_x,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

# set the bionano tool path
my $progpath = dirname(abs_path($0));
my $bionano = abs_path($progpath . "/../bionano");
my $toolpath = abs_path("$bionano/binary");
die("**ERROR: \"$toolpath\" directory is not found\n") unless(-d $toolpath);

# set the BNG refaligner
my $refAligner = "$toolpath/RefAligner";
die("**ERROR: Can not find RefAligner at $toolpath\n") unless(-f $refAligner);

# set the filter
my $cmapfilter = abs_path($progpath . "/cmap_filter.pl");
die("**ERROR: Can not find cmap_filter.pl at " . abs_path($progpath) . "\n") unless(-f $cmapfilter);

# check the required arguments
die("**ERROR: -r option must be specified\n") if(!defined $opt_r);
die("**ERROR: -q option must be specified\n") if(!defined $opt_q);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);

# set the arguments
my $scaffoldCmapFile = $opt_r;
my $NGSCmapFile = $opt_q;
my $optxml = (defined $opt_x) ? $opt_x : "$bionano/xml/alignArguments.xml";

my ($outdir, $outpre);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = $opt_o;
	$outpre = basename($opt_q);
	$outpre =~ s/\.cmap$/_cmap/;
	$outpre .= "_v_" . basename($opt_r);
	$outpre =~ s/\.cmap$/_cmap/;
}
else{
	$outdir = dirname($opt_o);
	$outpre = basename($opt_o);
}

# check the input files
die("**ERROR: reference cmap \"$scaffoldCmapFile\" does not exist") unless(-f $scaffoldCmapFile);
die("**ERROR: query cmap \"$NGSCmapFile\" does not exist") unless(-f $NGSCmapFile);
die("**ERROR: parameter file \"$optxml\" does not exist") unless(-f $optxml);

# prepare output directory
my $cmd = "mkdir -p $outdir";
my $retCode = system($cmd);
die("**ERROR: can not create directory $outdir") if($retCode != 0);

# process input files
my $XML = new XML::Simple(KeyAttr=>[]);
my $configs = $XML->XMLin($optxml);
my $paraDict = parseConfig($configs, "global");
my $maxmem = $paraDict->{'maxmem'}{val};
my $maxthreads = $paraDict->{'maxthreads'}{val};
my $maxvirtmem = (defined $paraDict->{'maxvirtmem'}{val}) ? $paraDict->{'maxvirtmem'}{val} : 0;
my $nthreads =  (defined $opt_t && ($opt_t =~ /^\d+$/)) ? (($opt_t < 1) ? 1 : $opt_t) : $maxthreads;

my $paraStage1 = parseConfig($configs, "stage_NGS_1st");
my $paraStage2 = parseConfig($configs, "stage_NGS_2nd");
my $firstPassPrefix = "${outpre}_1st_pass";
my $secondPassPrefix = "${outpre}_2nd_pass";

### 1st pass ###
### align the input sequences (post-conflict resolved) to the hybrid scaffold 
### use the sequence as the reference, and the hybrid scaffold as the query in the 1st pass
my $params = makeParams($paraStage1);
$params .= " -RAmem " . $paraDict->{'RAmem'}{val} if(defined $paraDict->{'RAmem'}{val});
$cmd = "$refAligner -ref $NGSCmapFile -i $scaffoldCmapFile -o $outdir/$firstPassPrefix -stdout -stderr -maxmem $maxmem -maxthreads $nthreads -maxvirtmem $maxvirtmem $params;";
print "$cmd\n";
$retCode = system("$cmd");
die("ERROR: $retCode") if($retCode != 0);

my ($usedScaffold1stPassRef, $usedSeq1stPassRef) = &getUsedIDs("$outdir/$firstPassPrefix.xmap");
my $numTotalSeq = `grep -v '^#' $NGSCmapFile | cut -f1 | uniq | wc -l`;
my ($firstXmapFile, $secondXmapFile) = ("", "");
if (defined $opt_s or $numTotalSeq == scalar(keys %$usedSeq1stPassRef)){
	$firstXmapFile = "$outdir/$firstPassPrefix.xmap";
}
else{
	if(scalar(keys %$usedSeq1stPassRef) > 0){
		$firstXmapFile = "$outdir/$firstPassPrefix.xmap";
	}
	my $getNotUsedSeq = &getNotUsedSeqIDs($NGSCmapFile, $usedSeq1stPassRef);
	my $notUsedSeqPrefix = "${firstPassPrefix}_not_used_seq";
	my $notUsedSeqIdFile = "$outdir/${notUsedSeqPrefix}_id.txt";
	&writeIDs2File($notUsedSeqIdFile, $getNotUsedSeq);
	$cmd = "$cmapfilter -c $NGSCmapFile -i $notUsedSeqIdFile -o $outdir/$notUsedSeqPrefix";
	print "$cmd\n";
	$retCode = system("$cmd");
	die("ERROR: $retCode") if($retCode != 0);
	### 2nd pass ###
	$params = makeParams($paraStage2);
	$params .= " -RAmem " . $paraDict->{'RAmem'}{val} if(defined $paraDict->{'RAmem'}{val});
	$cmd = "$refAligner -ref $outdir/$notUsedSeqPrefix.cmap -i $scaffoldCmapFile -o $outdir/$secondPassPrefix -stdout -stderr -maxmem $maxmem -maxthreads $nthreads -maxvirtmem $maxvirtmem $params;";
	print "$cmd\n";
	$retCode = system("$cmd");
	die("ERROR: $retCode") if($retCode != 0);
	$secondXmapFile = "$outdir/$secondPassPrefix.xmap";
}

# switch the reference and query fields in the alignment file
my $outSeqFile = "$outdir/${outpre}_q.cmap";
my $outHybridFile = "$outdir/${outpre}_r.cmap";
my $noFilteredXmap = &combineSwitchReferenceQueryXmap($firstXmapFile, $secondXmapFile);
my $notFilteredXmapFile = "$outdir/${outpre}_not_filtered.xmap";
writeXMapFile($noFilteredXmap, $notFilteredXmapFile, 1, $outHybridFile, $outSeqFile);

# filter out non-first-place alignments for each query sequence
my $finalXmap = &getBestAlignment($noFilteredXmap);
my $finalXmapFile = "$outdir/${outpre}.xmap";
writeXMapFile($finalXmap, $finalXmapFile, 1, $outHybridFile, $outSeqFile);

my $firstCmapFile = $firstXmapFile;
my $secondCmapFile = $secondXmapFile; 
$firstCmapFile =~ s/\.xmap$/_r.cmap/;
$secondCmapFile =~ s/\.xmap$/_r.cmap/;
my $seqCmap = &combineReadCMaps($firstCmapFile, $secondCmapFile);
my ($hybridCmap) = readCMap($scaffoldCmapFile);
my ($usedSeqRef, $usedScaffoldRef) = &getUsedIDs($finalXmapFile);
$seqCmap = &filterCMap($seqCmap, $usedSeqRef);
$hybridCmap = &filterCMap($hybridCmap, $usedScaffoldRef);
writeCMapFile($seqCmap, $outSeqFile, 1);
writeCMapFile($hybridCmap, $outHybridFile, 1);

exit(0);

### subroutines ###
sub getUsedIDs
{
	my ($alignXmapFile) = @_;
	my $usedQryIDsRef = {};
	my $usedRefIDsRef = {};

	open(IN, $alignXmapFile) or die("ERROR: getUsedIDs: cannot open alignment XMAP file $alignXmapFile: $!\n");
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;

		next if ($line =~ /^#/);
		# data lines
		my @content = split(/\t/, $line);
		next if(scalar(@content) < 14);
		$usedQryIDsRef->{$content[1]} = 1; # QryContigID
		$usedRefIDsRef->{$content[2]} = 1; # RefContigID
	}
	close IN;
	return ($usedQryIDsRef, $usedRefIDsRef);
}

sub getNotUsedSeqIDs
{
	my ($file, $usedSeqRef) = @_;
	my %notUsedSeq = ();

	open(IN, "grep -v '^#' $file | cut -f1 |") or die ("ERROR: getNotUsedSeqIDs: cannot open original (post-conflict resolved) sequence file: $!\n");
	while (my $line = <IN>)	{
		chomp $line;
		next if($line !~ /^\d+$/);
		if(!defined $usedSeqRef->{$line}){
			$notUsedSeq{$line} = 1;
		}
	}
	close IN;
	return \%notUsedSeq;
}

sub writeIDs2File
{
	my ($idFile, $ids) = @_;
	open(OUT, ">$idFile") or die("ERROR: writeIDsFile cannot write to $idFile: $!\n");
	foreach my $id (sort {$a <=> $b} keys %$ids){
		print OUT "$id\n";
	}
	close OUT;
}

sub combineXmapHeaders
{
	my ($xmap1, $xmap2) = @_;
	my $xmap;
	my $headers = [];
	my $data_name = [];
	my $data_type = [];
	my $version;
	if(!defined $xmap1->{headers} || !defined $xmap2->{headers}){
		my $xmap = (defined $xmap1->{headers}) ? $xmap1 : $xmap2;
		$headers = \@{$xmap->{headers}};
		$data_name = \@{$xmap->{dataName}};
		$data_type = \@{$xmap->{dataType}};
		$version = $xmap->{version};
	}
	else{
		$headers = \@{$xmap1->{headers}};
		my @data_name1 = @{$xmap1->{dataName}};
		my @data_name2 = @{$xmap2->{dataName}};
		my $cnt1 = scalar(@data_name1);
		my $cnt2 = scalar(@data_name2);
		my @data_type1 =@{$xmap1->{dataType}};
		my @data_type2 =@{$xmap2->{dataType}};
		my $cnt = min($cnt1, $cnt2);
		for(my $i=0; $i<$cnt; $i++){
			if($data_name1[$i] ne $data_name2[$i]){
				last;
			}
			if(!defined $data_type1[$i] or !defined $data_type2[$i] or $data_type1[$i] ne $data_type2[$i]){
				last;
			}
			push(@$data_name, $data_name1[$i]);
			push(@$data_type, $data_type1[$i]);
		}
		$version = min($xmap1->{version}, $xmap2->{version});
	}
	$headers = &simplifyXMapHeaders($headers);

	return ($headers, $data_name, $data_type, $version);
}

sub simplifyXMapHeaders
{
	my ($headers) = @_;
	my @newHeaders = ();
	for my $line (@$headers){
		$line =~ s/\s+$//g;
		if( $line =~ /^# XMAP File Version:|^# Label Channels:|^# Reference Maps From:|^# Query Maps From:|^#h|^#f/ ){
			push(@newHeaders, $line);
		}
	}

	return \@newHeaders;
}

sub combineSwitchReferenceQueryXmap
{
	my ($xmapFile1, $xmapFile2) = @_;

	my $xmap1 = ($xmapFile1 ne "") ? readXMap($xmapFile1) : {};
	my $xmap2 = ($xmapFile2 ne "") ? readXMap($xmapFile2) : {};

	my $xmap = {};
	($xmap->{headers}, $xmap->{dataName}, $xmap->{dataType}, $xmap->{version}) = combineXmapHeaders($xmap1, $xmap2);
	combineSwitchXmapFields($xmap1, $xmap2, $xmap);

	return $xmap;
}

sub combineSwitchXmapFields
{
	my ($xmap1, $xmap2, $finalXmap) = @_;

	my @data_name = @{ $finalXmap->{dataName} };
	my $numc = scalar(@data_name);
	my %changed = ();
	for my $name ("QryContigID", "RefContigID", "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", "HitEnum", "QryLen", "RefLen", "Alignment"){
		$changed{$name} = 1;
	}
	my @names = ();
	for(my $j=1; $j<$numc; $j++){ # skip 0: XmapEntryID
		$finalXmap->{hits}->{$data_name[$j]} = [];
		if(!defined $changed{$data_name[$j]}){
			push(@names, $data_name[$j]);
		}
	}
	my $id = 0;
	my $finalHits = {};
	for my $xmap ($xmap1, $xmap2){
		next if(!defined $xmap->{headers});
		my $hits = $xmap->{hits};
		for(my $i=0; $i<$xmap->{totalHits}; $i++){
			# new XMAP id number
			my $orientation = $hits->{Orientation}->[$i];
			# swap query and reference id
			($finalHits->{QryContigID}->[$id], $finalHits->{RefContigID}->[$id]) = swap($hits->{QryContigID}->[$i], $hits->{RefContigID}->[$i]);
			# swap the co-ordinates
			($finalHits->{QryStartPos}->[$id], $finalHits->{QryEndPos}->[$id], $finalHits->{RefStartPos}->[$id], $finalHits->{RefEndPos}->[$id]) =
			 	&swapPositions($hits->{QryStartPos}->[$i], $hits->{QryEndPos}->[$i], $hits->{RefStartPos}->[$i], $hits->{RefEndPos}->[$i], $orientation);
			# cigar string
			$finalHits->{HitEnum}->[$id] = &alterCigarString($hits->{HitEnum}->[$i], $orientation);
			if(defined $hits->{QryLen} and defined $hits->{RefLen}){
				# swap query and reference length
				($finalHits->{QryLen}->[$id], $finalHits->{RefLen}->[$id]) = swap($hits->{QryLen}->[$i], $hits->{RefLen}->[$i]);
			}
			if(defined $hits->{Alignment}){
				# alignment string
				$finalHits->{Alignment}->[$id] = &alterAlignmentString($hits->{Alignment}->[$i], $orientation);
			}
			foreach my $name (@names){
				$finalHits->{$name}->[$id] = $hits->{$name}->[$i];
			}
			$finalHits->{XmapEntryID}->[$id] = $id+1;
			$id++;
		}
	}
	$finalXmap->{hits} = $finalHits;
	$finalXmap->{totalHits} = $id;
}

sub getBestAlignment
{
	my ($xmap) = @_;
	my $filteredXmap = {};
	my @data_name = @{$xmap->{dataName}};
	my @data_type = @{$xmap->{dataType}};
	$filteredXmap->{headers} = $xmap->{headers};
	$filteredXmap->{dataName} = \@data_name;
	$filteredXmap->{dataType} = \@data_type;
	$filteredXmap->{version} = $xmap->{version};
	my $filteredHits = {};
	my $numc = scalar(@data_name);
	for(my $j=0; $j<$numc; $j++){
		$filteredHits->{$data_name[$j]} = [];
	}
	my %data = ();
	my $hits = $xmap->{hits};
	for(my $i=0; $i<$xmap->{totalHits}; $i++){
		my $qryId = $hits->{QryContigID}->[$i];
		push(@{$data{$qryId}}, {"confidence" => $hits->{Confidence}->[$i], "id" => $i});
	}
	my $id = 0;
	foreach my $qryId (sort {$a <=> $b} keys %data)	{
		@{$data{$qryId}} = sort { $b->{confidence} <=> $a->{confidence} } @{$data{$qryId}};
		my $i = $data{$qryId}->[0]->{id};
		for(my $j=1; $j<$numc; $j++){
			$filteredHits->{$data_name[$j]}->[$id] = $hits->{$data_name[$j]}->[$i];
		}
		$filteredHits->{$data_name[0]}->[$id] = $id+1;
		$id++;
	}
	$filteredXmap->{hits} = $filteredHits;
	$filteredXmap->{totalHits} = $id;
	return $filteredXmap;
}

sub combineReadCMaps
{
	my ($cmapFile1, $cmapFile2) = @_;
	my ($cmap1, $cmap2);
	my $combinedCmap;
	if($cmapFile1 ne ""){
		($cmap1) = readCMap($cmapFile1);
	}
	if($cmapFile2 ne ""){
		($cmap2) = readCMap($cmapFile2);
	}
	if($cmapFile1 eq "" or $cmapFile2 eq ""){
		$combinedCmap = ($cmapFile1 ne "") ? $cmap1 : $cmap2;
	}
	else{
		$combinedCmap = {};
		$combinedCmap->{headers} = $cmap1->{headers};
		$combinedCmap->{dataName} = $cmap1->{dataName};
		$combinedCmap->{dataType} = $cmap1->{dataType};
		$combinedCmap->{channels} = $cmap1->{channels};
		$combinedCmap->{version} = $cmap1->{version};
		my $contigs = {};
		for my $cmap ($cmap1, $cmap2){
			while( my ($id, $value) = each %{ $cmap->{contigs} } ){
				$contigs->{$id} = $value;
			}
		}
		$combinedCmap->{contigs} = $contigs;
	}
	$combinedCmap->{headers} = &simplifyCMapHeaders($combinedCmap->{headers});
	return $combinedCmap;
}

sub simplifyCMapHeaders
{
	my ($headers) = @_;
	my @newHeaders = ();
	for my $line (@$headers){
		$line =~ s/\s+$//g;
		if( $line =~ /^# CMAP File Version:|^# Label Channels:|^# Nickase Recognition Site|^# Number of Consensus Maps:|^#h|^#f/ or
			$line =~ /^# Values corresponding to intervals/ ){
			push(@newHeaders, $line);
		}
	}

	return \@newHeaders;
}

sub filterCMap
{
	my ($cmap, $usedIds) = @_;
	for my $id (keys %{$cmap->{contigs}}){
		if(!defined $usedIds->{$id}){
			delete $cmap->{contigs}->{$id};
		}
	}
	return $cmap;
}

sub alterAlignmentString
{
	my ($alignmentString, $orientation) = @_;
	$alignmentString =~ s/^\(|\)$//g;
	my $bReverse = ($orientation =~ /^-$/) ? 1 : 0;
	my @pairs = split(/\)\(/, $alignmentString);
	my @indices = ();
	for(my $i=0; $i<scalar(@pairs); $i++){
		my ($refIdx, $qryIdx) = split(/,/, $pairs[$i]);
		if($bReverse){
			unshift(@indices, "$qryIdx,$refIdx");
		}
		else{
			push(@indices, "$qryIdx,$refIdx");
		}
	}
	return "(". join(")(", @indices) .")";
}

sub alterCigarString
{
	my ($cigarString, $orientation) = @_;
	my @tokens = ();
	my $bReverse = ($orientation =~ /^-$/) ? 1 : 0;
	while($cigarString =~ /(\d+M|\d+I|\d+D)/g){
		my $token = $1;
		$token =~ tr/ID/DI/;
		if($bReverse){
			unshift(@tokens, $token);
		}
		else{
			push(@tokens, $token);
		}
	}
	return join("", @tokens);
}

sub swapPositions
{
	my ($qryStart, $qryEnd, $refStart, $refEnd, $orientation) = @_;
	($qryStart, $refStart) = swap($qryStart, $refStart);
	($qryEnd, $refEnd) = swap($qryEnd, $refEnd);
	# reference start should always be less than reference end
	if ($orientation =~ /^-$/)	{
		# negative orientation alignment, needs to swap the start and end co-ordinates
		($qryStart, $qryEnd) = swap($qryStart, $qryEnd);
		($refStart, $refEnd) = swap($refStart, $refEnd);
	}
	return ($qryStart, $qryEnd, $refStart, $refEnd);
}

sub swap
{
	my ($data1, $data2) = @_;
	return ($data2, $data1);
}
