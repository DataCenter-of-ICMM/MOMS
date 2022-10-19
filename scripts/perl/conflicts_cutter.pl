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

BEGIN{
	my $progpath = abs_path(dirname($0));
	unshift @INC, $progpath;
	select(STDERR);	$| = 1;
	select(STDOUT);	$| = 1;
}

use BNG::Utility;

my $program=basename($0);
my $usage = << "USAGE";
$program: A perl script for cutting CMAPs at the detected conflict points
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-i, --xmap <str>        The XMAP file containing the alignment information (REQUIRED)
	-r, --reference <str>   The reference CMAP file for alignment (REQUIRED)
	-q, --query <str>       The query CMAP file for alignment (REQUIRED)
	-r0, --ref0  <str>      The original reference CMAP file for alignment (REQUIRED)
	-q0, --qry0  <str>      The original query CMAP file for alignment (REQUIRED)
	-c, --conflict <str     The conflict.txt file (REQUIRED)
	-o, --output <str>      The output path (REQUIRED)
	-s, --subset <str>      Comma separated string denote subset of a group (e.g. 1,3) (default: 1,1)
	-h, --help              Help

Options for invoking \$bionano/binary/RefAligner:
	-x, --xml <str>         Parameter file in XML format (default: \$bionano/xml/conflictsArguments.xml)

Example:
	$program -i align.xmap -r align_r.cmap -q align_q.cmap -r0 ngs_assembly.cmap -q0 bng_assembly.cmap -c conflicts.txt -o conflict_cut

USAGE

use vars qw($opt_i $opt_r $opt_q $opt_r0 $opt_q0 $opt_c $opt_o $opt_s $opt_x $opt_h);
GetOptions( "i|xmap=s" => \$opt_i,
			"r|reference=s" => \$opt_r,
			"q|query=s" => \$opt_q,
			"r0|ref0=s" => \$opt_r0,
			"q0|qry0=s" => \$opt_q0,
			"c|conflict=s" => \$opt_c,
			"o|output=s" => \$opt_o,
			"s|subset=s" => \$opt_s,
			"x|xml=s" => \$opt_x,
			"h|help"  => \$opt_h);

die($usage) if($opt_h);

my $progpath = abs_path(dirname($0));
my $bionano = abs_path(dirname($0) . "/../bionano");
my $cmap_cutter = "$progpath/cmap_cutter.pl";
my $refAligner  = "$bionano/binary/RefAligner";
my $optxml = (defined $opt_x) ? $opt_x : "$bionano/xml/conflictsArguments.xml";
die("**ERROR: Can not find RefAligner at $bionano/binary\n") unless(-f $refAligner);
die("**ERROR: Can not find cmap_cutter.pl at $progpath/\n") unless(-f $cmap_cutter);
die("**ERROR: parameter file \"$optxml\" does not exist") unless(-f $optxml);

# check the required arguments
die("**ERROR: -i option must be specified\n") if(!defined $opt_i);
die("**ERROR: -r option must be specified\n") if(!defined $opt_r);
die("**ERROR: -q option must be specified\n") if(!defined $opt_q);
die("**ERROR: -r0 option must be specified\n") if(!defined $opt_r0);
die("**ERROR: -q0 option must be specified\n") if(!defined $opt_q0);
die("**ERROR: -c option must be specified\n") if(!defined $opt_c);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);

my ($sub_i, $sub_n) = (defined $opt_s and $opt_s =~ /,/) ? split(/,/, $opt_s) : (1, 1);
die("**ERROR: invalid -s option specified\n") if($sub_i !~ /^\d+$/ or $sub_n !~ /^\d+$/ or ($sub_i < 0) or ($sub_i > $sub_n));

my ($xmap_in, $ref_cmap_in, $qry_cmap_in, $ref0_cmap_in, $qry0_cmap_in, $conflict_in) = ($opt_i, $opt_r, $opt_q, $opt_r0, $opt_q0, $opt_c);

# check the input files
die("**ERROR: XMAP file \"$xmap_in\" does not exist") unless(-f $xmap_in);
die("**ERROR: reference cmap \"$ref_cmap_in\" does not exist") unless(-f $ref_cmap_in);
die("**ERROR: query cmap \"$qry_cmap_in\" does not exist") unless(-f $qry_cmap_in);
die("**ERROR: original reference cmap \"$ref0_cmap_in\" does not exist") unless(-f $ref0_cmap_in);
die("**ERROR: original query cmap \"$qry0_cmap_in\" does not exist") unless(-f $qry0_cmap_in);
die("**ERROR: conflict file \"$conflict_in\" does not exist") unless(-f $conflict_in);

# set the output directory and prefix
my ($outdir, $outpre);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = $opt_o;
	$outpre = "classified";
}
else{
	$outdir = dirname($opt_o);
	$outpre = basename($opt_o);
}

# read in the parameter configuration
my $XML = new XML::Simple(KeyAttr=>[]);
my $configs = $XML->XMLin($optxml);
my $paraDict = &parseConfig($configs, "conflicts_cut");
my ($window_size, $min_quality, $min_coverage, $max_overhang) = ($paraDict->{'window_size'}{val}, $paraDict->{"min_quality_score_threshold"}{val}, $paraDict->{"min_coverage_threshold"}{val}, $paraDict->{"max_overhang"}{val});

###################################################################################
# define some of the outputs
my $ref0pre = basename($ref0_cmap_in);
$ref0pre =~ s/\..*$//;
my $qry0pre = basename($qry0_cmap_in);
$qry0pre =~ s/\..*$//;

my $mdldir = "$outdir/middle_files/$outpre";
system("mkdir -p $mdldir");

# 1) read in break point file
# extract the breakpoints (the query is the BioNano gm)
my ($headerLinesRef, $conflictsQueryRef, $conflictsAllRef) = getConflicts($conflict_in);
$conflictsQueryRef = extendByWindowSize($conflictsQueryRef, $window_size);
$conflictsQueryRef = sortByCoord($conflictsQueryRef, "start");

# 2) check the BioNano chimeric scores at the conflict loci
# now extract the chimeric quality score
my ($qScoresRef, $noQScoreFlag) = getQScores($qry0_cmap_in);
print "noQScore = $noQScoreFlag\n";
if ($noQScoreFlag == 0)	{
	# assumes that the cmap file is already sorted by the label position in the gm
	$conflictsQueryRef = findOverlap($conflictsQueryRef, $qScoresRef);
	$conflictsQueryRef = flagCut($conflictsQueryRef, $min_quality, $min_coverage);
	
	# update the conflictsAllRef hash
	$conflictsAllRef = updateConflictsAll($conflictsQueryRef, $conflictsAllRef);

# 3) for each conflict locus, determine whether it is the BioNano or the sequence assembly that is needed to be cut
	# check for alternate support beside the conflicting partners
	my $align1QueryCmapRef = readAlign1Cmap2Paint($qry_cmap_in);
	my $align1ReferenceCmapRef = readAlign1Cmap2Paint($ref_cmap_in);
	my ($align1Xmap4Ref2PaintRef, $align1Xmap4Qry2PaintRef) = readAlign1Xmap4Paint($xmap_in);
	$align1Xmap4Ref2PaintRef = sortXmap4PaintByConfScore($align1Xmap4Ref2PaintRef);	# sort the alignments by confidence score
	$align1Xmap4Qry2PaintRef = sortXmap4PaintByConfScore($align1Xmap4Qry2PaintRef);
	$align1QueryCmapRef = paintAlign1Cmap($align1QueryCmapRef, $align1Xmap4Qry2PaintRef);	# now paint the cmaps according to the score of their alignment partners
	$align1ReferenceCmapRef = paintAlign1Cmap($align1ReferenceCmapRef, $align1Xmap4Ref2PaintRef);		
	$conflictsAllRef = checkSupportingPaint($conflictsAllRef, $align1QueryCmapRef, $align1ReferenceCmapRef, $max_overhang);
} else	{
	### remember, if there is no quality column, print warning message
	print "WARNING: In the cut conflict stage, but there were no quality scores in the BioNano assembly\n";
} # if noQScoreFlag

# 4) write an updated break point status file
# print an updated breakpoint file
printUpdatedBreakPointFile("$mdldir/${outpre}_conflicts_cut_status.txt", $headerLinesRef, $conflictsAllRef, "auto");

# now read in the sequence file and the gm file to figure out the largest id in each file
my ($maxSeqId, $seqLengthsRef) = readCmapIdLength($ref0_cmap_in);
my ($maxGMId, $gMLengthsRef) = readCmapIdLength($qry0_cmap_in);

# print coordinate translations, unless the fragment has no conflict, is correct whenever there is conflict, and  not to be discarded
my ($ngs_coord_trans_file, $bn_coord_trans_file) = ("$outdir/${outpre}_auto_cut_NGS_coord_translation.txt", "$outdir/${outpre}_auto_cut_BN_coord_translation.txt");
my ($seqCutIdsRef, $seqDiscardIdsRef, $seqToCutIdFragments) = printCoordTranslation($ngs_coord_trans_file, $conflictsAllRef, "ref", $seqLengthsRef, $maxSeqId, [$sub_i, $sub_n]);
my ($gMCutIdsRef, $gMDiscardIdsRef, $gMToCutIdFragments) = printCoordTranslation($bn_coord_trans_file, $conflictsAllRef, "qry", $gMLengthsRef, $maxGMId, [1, 1]);

# bed file annotation window size, a small size to annotate where cuts were made
my $bedAnnotLength = 1000;	
# print bed file co-ordinates (ngs pre-cut ngs coordinates; BN pre-cut BN coordinates; BN pre-cut BN coordinate projected onto NGS coordinates)
printPreCutBedRespectiveCoords("$mdldir/${outpre}_ngs_pre_cut_annotations.bed", $seqCutIdsRef, $bedAnnotLength, "NGS", "255,102,0");
printPreCutBedRespectiveCoords("$mdldir/${outpre}_bn_pre_cut_annotations.bed", $gMCutIdsRef, $bedAnnotLength, "BN", "102,0,255");
printGMPreCutBedRefCoords("$mdldir/${outpre}_bn_pre_cut_projected_ngs_coord_annotations.bed", $conflictsAllRef, $seqDiscardIdsRef, $gMDiscardIdsRef, $bedAnnotLength, "102,0,255");

my ($outSeqCutPreDiscardCmap, $outGMCutPreDiscardCmap) = ("$mdldir/${ref0pre}_cut_pre_exclude.cmap", "$mdldir/${qry0pre}_cut_pre_exclude.cmap");

# now tell RefAligner to do the cut (remember that the co-ordinates fed into RefAligner MUST be in kilobase not bp)
callCmapCutter($cmap_cutter, $ref0_cmap_in, $ngs_coord_trans_file, $outSeqCutPreDiscardCmap);
callCmapCutter($cmap_cutter, $qry0_cmap_in, $bn_coord_trans_file, $outGMCutPreDiscardCmap);

### now this part deals with the fragments that need to be discarded
my $toKeepSeqIdsRef = distinguishKeepDiscards($seqLengthsRef, $seqToCutIdFragments, $seqDiscardIdsRef);
my $toKeepGMIdsRef = distinguishKeepDiscards($gMLengthsRef, $gMToCutIdFragments, $gMDiscardIdsRef);

# call RefAligner to make two set of files
my ($outSeqFile, $outKeepSeqIdList) = ("$outdir/${ref0pre}_cut.cmap", "$outdir/${ref0pre}_keep_ids.txt");
callRefAlignerKeepDiscard($refAligner, $outSeqCutPreDiscardCmap, $outSeqFile, $toKeepSeqIdsRef, $outKeepSeqIdList);
my ($outDiscardSeqFile, $outDiscardSeqIdList) = ("$mdldir/${ref0pre}_exclude.cmap", "$mdldir/${ref0pre}_exclude_ids.txt");
callRefAlignerKeepDiscard($refAligner, $outSeqCutPreDiscardCmap, $outDiscardSeqFile, $seqDiscardIdsRef, $outDiscardSeqIdList);

my ($outGMFile, $outKeepGMIdList) = ("$outdir/${qry0pre}_cut.cmap", "$outdir/${qry0pre}_keep_ids.txt");
callRefAlignerKeepDiscard($refAligner, $outGMCutPreDiscardCmap, $outGMFile, $toKeepGMIdsRef, $outKeepGMIdList);
my ($outDiscardGMFile, $outDiscardGMIdList) = ("$mdldir/${qry0pre}_exclude.cmap", "$mdldir/${qry0pre}_exclude_ids.txt");
callRefAlignerKeepDiscard($refAligner, $outGMCutPreDiscardCmap, $outDiscardGMFile, $gMDiscardIdsRef, $outDiscardGMIdList);

sub callRefAlignerKeepDiscard
{
	my ($refAligner, $inFile, $outFile, $idsRef, $idListFile) = @_;
	my $outFilePrefix = $outFile;	$outFilePrefix =~ s/\.cmap$//;

	# now make a list txt file, whose purpose is to feed into RefAligner with selectidf
	printIdListFile($idListFile, $idsRef);

	my $params = "-f -minsites 0 -i $inFile -o $outFilePrefix -merge -selectidf $idListFile -stdout 2>&1";
	my $cmd = "$refAligner $params";
	system($cmd);
	if($? != 0){
		die("ERROR: failed to execute: $!\n");
	}
}

sub printIdListFile
{
	my ($file, $idsRef) = @_;
	# this subroutine just output to a file with id on each line
	open(OUT, ">$file") or die "ERROR: printIdListFile: cannot write to $file: $!\n";
	foreach my $id (keys %$idsRef)	{
		print OUT "$id\n";
	} # foreach id
	close OUT;
}

sub distinguishKeepDiscards
{
	my ($allIdsRef, $cuttedIdsRef, $discardIdsRef) = @_;
	# remember allIdsRef has the lengths of all the pre-cut fragments
	# cuttedIdsRef has ids of those that went through cut, its old id and all new ids
	# discardIdsRef has ids of those need to be thrown away (MUST be mutually exclusive from the cuttedIdsRef)
	my %toKeepIds = ();
	
	# now deal with those fragments where no cut was involved
	foreach my $id (keys %$allIdsRef)	{
		next if (exists $discardIdsRef->{$id});	# skip those that have been flagged to be discard
		next if (exists $cuttedIdsRef->{$id});
		$toKeepIds{$id} = 1;			
	}

	# add in the new ids generated from performing cuts
	foreach my $id (keys %$cuttedIdsRef)	{
		for (my $i = 0; $i < scalar(@{$cuttedIdsRef->{$id}}); $i += 1)	{
			$toKeepIds{$cuttedIdsRef->{$id}[$i]{newId}} = 1;		
		}
	}

	return \%toKeepIds;
}

sub printPreCutBedRespectiveCoords
{
	my ($file, $toCutIdsRef, $bedAnnotLength, $dataType, $colourString) = @_;
	# this subroutine looks at an array of cut coordinates and determine the bed co-ordinates (in reference co-ordinates)
	my %toCutBedFragments = ();
	foreach my $id (keys %$toCutIdsRef)	{
		for (my $i = 0; $i < scalar(@{$toCutIdsRef->{$id}}); $i += 1)	{
			my $cRef = $toCutIdsRef->{$id}[$i];
			push(@{$toCutBedFragments{$id}}, {start => $cRef->{start} - $bedAnnotLength, end => $cRef->{start} + $bedAnnotLength, oriPosition => $cRef->{start},
				descString => "$dataType"."_proposed_cut_$dataType"."Id_$id"."_$dataType"."Position_$cRef->{start}"});
		}
	}

	# print to file
	open(OUT, ">$file") or die "ERROR: printNGSCutBedRefCoords: cannot write to $file: $!\n";
	my $count = 1;
	foreach my $id (sort { $a <=> $b } keys %toCutBedFragments)	{
		for (my $i = 0; $i < scalar(@{$toCutBedFragments{$id}}); $i += 1)	{
			my $aRef = $toCutBedFragments{$id}[$i];
			print OUT join("\t", $id, $aRef->{start}, $aRef->{end}, $aRef->{descString}, $count, "+", $aRef->{start}, $aRef->{end}, $colourString)."\n";
			$count += 1;
		}
	}
	close OUT;
}

sub printGMPreCutBedRefCoords
{
	my ($file, $conflictsAllRef, $toDiscardRefIdsRef, $toDiscardQryIdsRef, $bedAnnotLength, $colourString) = @_;
	my %toCutRefIds = ();		# it is okay for this array to store redundant reference coordinates, as multiple queries can conflict with the same reference at the same reference position
	foreach my $theKey (keys %$conflictsAllRef)	{
		my ($xId, $refId, $qryId) = split(/\t/, $theKey);	
		# if either refId or qryId is -1, then skip
		next if ($refId =~ /^-1$/ || $qryId =~ /^-1$/);
		# if either ref or qry needs to be discard, then skip
		next if (exists $toDiscardRefIdsRef->{$refId} || exists $toDiscardQryIdsRef->{$qryId});
		
		# now build an array of cut co-ordinates on the reference, but only if the qry needs to be cut
		for (my $i = 0; $i < scalar(@{$conflictsAllRef->{$theKey}}); $i += 1)	{
			my $aRef = $conflictsAllRef->{$theKey}[$i];
			
			if ($aRef->{qryLeftBkptStatus} =~ /^cut$/)	{
				push(@{$toCutRefIds{$refId}}, {refCoord => $aRef->{refLeftBkpt}, 
					refStart => $aRef->{refLeftBkpt} - $bedAnnotLength, refEnd => $aRef->{refLeftBkpt} + $bedAnnotLength,
					qryId => $qryId, qryCoord => $aRef->{qryLeftBkpt},
					descString => "BN_proposed_cut_BNId_$qryId"."_BNPosition_$aRef->{qryLeftBkpt}"."_projectedNGSId_$refId"."_projectedNGSPosition_$aRef->{refLeftBkpt}"});
			}
			
			if ($aRef->{qryRightBkptStatus} =~ /^cut$/)	{
				push(@{$toCutRefIds{$refId}}, {refCoord => $aRef->{refRightBkpt}, 
					refStart => $aRef->{refRightBkpt} - $bedAnnotLength, refEnd => $aRef->{refRightBkpt} + $bedAnnotLength, 
					qryId => $qryId, qryCoord => $aRef->{qryRightBkpt},
					descString => "BN_proposed_cut_BNId_$qryId"."_BNPosition_$aRef->{qryRightBkpt}"."_projectedNGSId_$refId"."_projectedNGSPosition_$aRef->{refRightBkpt}"});
			}

		}
	}
	
	# now print a bed file
	open(OUT, ">$file") or die "ERROR: printGMPreCutBedRefCoords: cannot write to $file: $!\n";
	my $count = 1;
	foreach my $refId (sort { $a <=> $b } keys %toCutRefIds)	{
		for (my $i = 0; $i < scalar(@{$toCutRefIds{$refId}}); $i += 1)	{
			my $aRef = $toCutRefIds{$refId}[$i];
			print OUT join("\t", $refId, $aRef->{refStart}, $aRef->{refEnd}, $aRef->{descString}, $count, "+", $aRef->{refStart}, $aRef->{refEnd}, $colourString)."\n";
			$count += 1;
		}
	}
	close OUT;
}

sub callCmapCutter
{
	my ($cmap_cutter, $cmap_in, $coord_trans_file, $cmap_out) = @_;
	my $cmd = "$cmap_cutter -i $cmap_in -t $coord_trans_file -o $cmap_out 2>&1";
	print "Running command: $cmd\n";
	system("$cmd");
	if($? != 0){
		die("ERROR: failed to execute: $!\n");
	}
}

sub callRefAlignerBreak
{
	my ($refAligner, $inFile, $outFile, $maxId, $cutIdsRef, $noQScoreFlag, $assemblyType) = @_;
	# actual call to RefAligner to do cut (remember that the co-ordinates must be in kilobase)
	# Aside: remember that AssignAlignType prints the breakpoints such that the conflict label is always kept with the aligned match group
	my $outFilePrefix = $outFile;	$outFilePrefix =~ s/\.cmap$//;

	my $params = "-f -minsites 0 -i $inFile -o $outFilePrefix -merge";
	if ($noQScoreFlag == 0 && scalar(keys %$cutIdsRef) > 0)	{
		# only break if there was chimeric score provided 
		$params .= " -break $maxId";
	}
	foreach my $id (keys %$cutIdsRef){
		$params .= " $id";
		for (my $i = 0; $i < scalar(@{$cutIdsRef->{$id}}); $i++){
			my $theCoord = $cutIdsRef->{$id}[$i]{start} / 1000;
			my $theCoordString = sprintf("%.4f", $theCoord);
			$params .= " $theCoordString";
		}
	}
	my $cmd = "$refAligner $params -stdout 2>&1";
	print "Running command: $cmd\n";
	system("$cmd");
	if($? != 0){
		die("ERROR: failed to execute: $!\n");
	}
}

sub printCoordTranslation
{
	my ($file, $conflictsAllRef, $assemblyType, $lengthsRef, $maxId, $subset) = @_;
	# this subroutine will print out the pre-cut and post-cut ids and co-ordinates in a file
	# and return a list of ids and the co-ordinates to be cut
	# but of course, do NOT mark a contig to be cut, if it is to be discarded 

	my $toDiscardIdsRef = determineDiscardFragments($conflictsAllRef, $assemblyType);

	my $toCutIdsRef = determineCutCoords($conflictsAllRef, $assemblyType, $toDiscardIdsRef);

	# sort by coord
	$toCutIdsRef = sortByCoord($toCutIdsRef, "start");	
	
	# determine the fragments
	my $toCutIdFragmentsRef = determineFragments($toCutIdsRef, $lengthsRef, $maxId, $subset);

	open(OUT, ">$file") or die ("ERROR printCoordTranslation: cannot write to $file");
	print OUT join("\t", "oldId", "oldStart", "oldEnd", "newId", "newStart", "newEnd")."\n";
	# first print out those that will be cut
	foreach my $id (sort { $a <=> $b } keys %$toCutIdFragmentsRef)	{
		for (my $i = 0; $i < scalar(@{$toCutIdFragmentsRef->{$id}}); $i += 1)	{
			my $aRef = $toCutIdFragmentsRef->{$id}[$i];
			print OUT join("\t", $id, $aRef->{oldStart}, $aRef->{oldEnd}, $aRef->{newId}, $aRef->{newStart}, $aRef->{newEnd})."\n";
		}
	}

	# then print out the ones that do not need to be cut
	foreach my $id (sort { $a <=> $b } keys %$lengthsRef)	{
		next if (exists $toCutIdFragmentsRef->{$id});	# skip those that needed to be cut
		my $endCoord = $lengthsRef->{$id} - 1;
		print OUT join("\t", $id, 0, $endCoord, $id, 0, $endCoord)."\n";
	}
	close OUT;	

	return ($toCutIdsRef, $toDiscardIdsRef, $toCutIdFragmentsRef);
}

sub determineFragments
{
	my ($toCutIdsRef, $lengthsRef, $maxId, $subset) = @_;
	# this subroutine looks at an array of cut coordinates and determine the fragment coordinates
	my ($sub_i, $sub_n) = @{$subset};
	my $groupId = $maxId * $sub_n;
	my $baseId = $maxId * $sub_i;
	my %toCutIdFragments = ();

	foreach my $id (keys %$toCutIdsRef)	{
		die "ERROR: cannot find the original length of the contig whose id=$id\n" if (! exists $lengthsRef->{$id});
		for (my $i = 0; $i < scalar(@{$toCutIdsRef->{$id}}); $i += 1)	{
			my $cRef = $toCutIdsRef->{$id}[$i];
			if ($i == 0)	{
				# first fragment
				push(@{$toCutIdFragments{$id}}, {oldStart => 0, oldEnd => $cRef->{start}, newStart => 0, newEnd => $cRef->{start}, newId => $baseId + $groupId * $i + $id});
			}
			if ($i == $#{$toCutIdsRef->{$id}})	{
				# last fragment
				push(@{$toCutIdFragments{$id}}, {oldStart => $cRef->{start} + 1, oldEnd => $lengthsRef->{$id} - 1, newStart => 0, newEnd => ($lengthsRef->{$id} - 1) - ($cRef->{start} + 1), newId => $baseId + $groupId * ($i + 1) + $id});
			} else	{
				# middle fragment
				push(@{$toCutIdFragments{$id}}, {oldStart => $cRef->{start} + 1, oldEnd => $toCutIdsRef->{$id}[$i + 1]{start}, newStart => 0, newEnd => $toCutIdsRef->{$id}[$i + 1]{start} - $cRef->{start} - 1, newId => $baseId + $groupId * ($i + 1) + $id});
			}
		}
	}
	return \%toCutIdFragments;
}

sub determineCutCoords
{
	my ($conflictsAllRef, $assemblyType, $toDiscardIdsRef) = @_;
	# this subroutines puts the cut coordinates into an array (after making a non-redundant list of co-ordinates first)
	# again, ignore those fragments that have been flagged as discard
	my %nrCutIds = ();
	foreach my $theKey (keys %$conflictsAllRef)	{
		my ($xId, $refId, $qryId) = split(/\t/, $theKey);
		for (my $i = 0; $i < scalar(@{$conflictsAllRef->{$theKey}}); $i += 1)	{
			my $aRef = $conflictsAllRef->{$theKey}[$i];
			if ($assemblyType eq "ref")	{
				# examine the sequence bkpt status
				$nrCutIds{$refId}{$aRef->{refLeftBkpt}} = 1 if ($aRef->{refLeftBkptStatus} =~ /^cut$/ && ! exists $toDiscardIdsRef->{$refId} && $refId !~ /^-1$/ && $aRef->{refLeftBkpt} != -1);
				$nrCutIds{$refId}{$aRef->{refRightBkpt}} = 1 if ($aRef->{refRightBkptStatus} =~ /^cut$/ && ! exists $toDiscardIdsRef->{$refId} && $refId !~ /^-1$/ && $aRef->{refRightBkpt} != -1);
			} else	{
				# examine the genome map bkpt status
				$nrCutIds{$qryId}{$aRef->{qryLeftBkpt}} = 1 if ($aRef->{qryLeftBkptStatus} =~ /^cut$/ && ! exists $toDiscardIdsRef->{$qryId} && $qryId !~ /^-1$/ && $aRef->{qryLeftBkpt} != -1);
				$nrCutIds{$qryId}{$aRef->{qryRightBkpt}} = 1 if ($aRef->{qryRightBkptStatus} =~ /^cut$/ && ! exists $toDiscardIdsRef->{$qryId} && $qryId !~ /^-1$/ && $aRef->{qryRightBkpt} != -1);
			}
		}
	}

	my %toCutIds = ();
	foreach my $theId (keys %nrCutIds)	{
		foreach my $theStart (keys %{$nrCutIds{$theId}})	{
			push(@{$toCutIds{$theId}}, {start => $theStart});
		}
	}

	return \%toCutIds;
}

sub determineDiscardFragments
{
	my ($conflictsAllRef, $assemblyType) = @_;
	my %toDiscardIds = ();
	foreach my $theKey (keys %$conflictsAllRef)	{
		my ($xId, $refId, $qryId) = split(/\t/, $theKey);
		for (my $i = 0; $i < scalar(@{$conflictsAllRef->{$theKey}}); $i += 1)	{
			my $cRef = $conflictsAllRef->{$theKey}[$i];
			$toDiscardIds{$refId} = 1 if ($cRef->{refStatus} =~ /^exclude$/ && $refId !~ /^-1$/ && $assemblyType =~ /^ref$/i);
			$toDiscardIds{$qryId} = 1 if ($cRef->{qryStatus} =~ /^exclude$/ && $qryId !~ /^-1$/ && $assemblyType =~ /^qry$/i);
		}
	}
	return \%toDiscardIds;
}

sub readCmapIdLength
{
	my ($file) = @_;
	# read in a cmap file, and determine the largest sequence id
	my $maxSeqId = -1;
	my %theLengths = ();
	open(IN, "$file") or die "ERROR: readCmapId: reading in $file: $!\n";
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;		
		next if ($line =~ /^#/);	# skip header lines
		my @content = split(/\t/, $line);
		my ($id, $aLength) = @content[0..1];
		$aLength = int($aLength);	# again change to integer
		$theLengths{$id} = $aLength;
		$maxSeqId = ($maxSeqId < $id) ? ($id) : ($maxSeqId);
	} # while line
	close IN;
	return ($maxSeqId, \%theLengths);	
}
=begin manual breakpoint
## read in subroutines
sub getManualBreakPointFile
{
	# this subroutine reads in manually modified break point status file, and return a hash containing that information
	my ($file, $headerLinesRef, $conflictsAllRef) = @_;
	open(IN, $file) or die "ERROR: getManualBreakPointFile: cannot read in $file\n: $!\n";
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^#/i)	{
			push(@$headerLinesRef, $line);
			next;
		}
		# store the conflict information into the conflictsAll hash
		$line =~ s/\s+/\t/g;	$line =~ s/\s+/\t/g;
		my @content = split(/\t/, $line);
		my $xId = $content[0];
		my ($refId, $refLeftBkpt, $refRightBkpt, $refAlignOrientation, $refLeftBkptStatus, $refRightBkptStatus, $refStatus) = @content[2..8];
		my ($qryId, $qryLeftBkpt, $qryRightBkpt, $qryAlignOrientation, $qryLeftBkptStatus, $qryRightBkptStatus, $qryStatus) = @content[10..16]; 
		
		breakPointFileCheck($line, $xId, $refId, $refLeftBkpt, $refRightBkpt, $refAlignOrientation, $refLeftBkptStatus, $refRightBkptStatus, $refStatus,
						$qryId, $qryLeftBkpt, $qryRightBkpt, $qryAlignOrientation, $qryLeftBkptStatus, $qryRightBkptStatus, $qryStatus);
		# change the bkpt to integer, and status to lower case
		($refLeftBkpt, $refRightBkpt, $qryLeftBkpt, $qryRightBkpt) = (int($refLeftBkpt), int($refRightBkpt), int($qryLeftBkpt), int($qryRightBkpt));
		($refLeftBkptStatus, $refRightBkptStatus, $qryLeftBkptStatus, $qryRightBkptStatus) = (lc($refLeftBkptStatus), lc($refRightBkptStatus), lc($qryLeftBkptStatus), lc($qryRightBkptStatus));
		($refStatus, $qryStatus) = (lc($refStatus), lc($qryStatus));
		# record
		push(@{$conflictsAllRef->{"$xId\t$refId\t$qryId"}}, {	refLeftBkpt => $refLeftBkpt, qryLeftBkpt => $qryLeftBkpt,
									refRightBkpt => $refRightBkpt, qryRightBkpt => $qryRightBkpt,
									refAlignOrientation => $refAlignOrientation, qryAlignOrientation => $qryAlignOrientation,
									refLeftBkptStatus => $refLeftBkptStatus, qryLeftBkptStatus => $qryLeftBkptStatus, 
									refRightBkptStatus => $refRightBkptStatus, qryRightBkptStatus => $qryRightBkptStatus,
									refStatus => $refStatus, qryStatus => $qryStatus
		});
	}
	close IN;
	return ($headerLinesRef, $conflictsAllRef);
}

sub breakPointFileCheck
{
	my ($line, $xId, $refId, $refLeftBkpt, $refRightBkpt, $refAlignOrientation, $refLeftBkptStatus, $refRightBkptStatus, $refStatus,
			$qryId, $qryLeftBkpt, $qryRightBkpt, $qryAlignOrientation, $qryLeftBkptStatus, $qryRightBkptStatus, $qryStatus) = @_;
	
	die ("ERROR: getManualBreakPointFile: Incorrect XMapId, must be a positive integer or -1. Line = $line\n") if ($xId !~ /^\d+$/ && $xId !~ /^-1$/);
	die ("ERROR: getManualBreakPointFile: Incorrect refId, must be a positive integer or -1. Line = $line\n") if ($refId !~ /^\d+$/ && $refId !~ /^-1$/);
	die ("ERROR: getManualBreakPointFile: Incorrect qryId, must be a positive integer or -1. Line = $line\n") if ($qryId !~ /^\d+$/ && $qryId !~ /^-1$/);
	die ("ERROR: getManualBreakPointFile: Incorrect conflict reference co-ordinate. Line = $line\n") if (!($refRightBkpt =~ /^\d+\.?\d*$/ || $refRightBkpt == -1) || !($refLeftBkpt =~ /^\d+\.?\d*$/ || $refLeftBkpt == -1));
	die ("ERROR: getManualBreakPointFile: Incorrect conflict query co-ordinate. Line = $line\n") if (!($qryRightBkpt =~ /^\d+\.?\d*$/ || $qryRightBkpt == -1) || !($qryLeftBkpt =~ /^\d+\.?\d*$/ || $qryLeftBkpt == -1));
	die ("ERROR: getManualBreakPointFile: Incorrect conflict reference conflict point status, must be okay, cut or -. Line = $line\n") if ($refRightBkptStatus !~ /^(okay|cut|-)$/i || $refLeftBkptStatus !~ /^(okay|cut|-)$/i);
	die ("ERROR: getManualBreakPointFile: Incorrect conflict query conflict point status, must be okay, cut or -. Line = $line\n") if ($qryRightBkptStatus !~ /^(okay|cut|-)$/i || $qryLeftBkptStatus !~ /^(okay|cut|-)$/i);
	die ("ERROR: getManualBreakPointFile: Incorrect reference status, must be okay or exclude. Line = $line\n") if ($refStatus !~ /^(okay|exclude|-)$/i);
	die ("ERROR: getManualBreakPointFile: Incorrect query status, must be okay or exclude. Line = $line\n") if ($qryStatus !~ /^(okay|exclude|-)$/i);
}
=end manual breakpoint
=cut

sub printUpdatedBreakPointFile
{
	my ($file, $headerLinesRef, $conflictsAllRef, $flag) = @_;
	open(OUT, ">$file") or die "ERROR printUpdatedBreakPointFile: cannot write to $file: $!\n";
	my @headerContent = split(/\t/, $headerLinesRef->[0]);

	my @headerIndecies = ();
	if ($flag eq "auto")	{
		# if flag is auto, that means that the input conflict file is from assignAlignType, thus not having the columns with cut and discard information
		push(@headerIndecies, (0, 5, 6, 10));
	} else{
		# if flag is manual, that means that the input conflict file is from previous cut_conflict_status file, thus already having the columns with cut and discard information
		push(@headerIndecies, (0, 5, 9, 13));
	}
	print OUT join("\t", @headerContent[$headerIndecies[0]..$headerIndecies[1]])."\t".join("\t", "ref_leftBkpt_toCut", "ref_rightBkpt_toCut", "ref_toDiscard")."\t";
	print OUT join("\t", @headerContent[$headerIndecies[$#headerIndecies - 1]..$headerIndecies[$#headerIndecies]])."\t".join("\t", "qry_leftBkpt_toCut", "qry_rightBkpt_toCut", "qry_toDiscard")."\n";

	print OUT join("\t", "# id/-1", 	"ref", "id/-1", "position/-1", "position/-1", "+/-", "okay/cut/-", "okay/cut/-", "okay/exclude/-")."\t";
	print OUT join("\t", 		"qry", "id/-1", "position/-1", "position/-1", "+/-", "okay/cut/-", "okay/cut/-", "okay/exclude/-")."\n";
	
	foreach my $theKey (keys %$conflictsAllRef)	{
		my ($xId, $refId, $qryId) = split(/\t/, $theKey);
		for (my $i = 0; $i < scalar(@{$conflictsAllRef->{$theKey}}); $i += 1)	{
			my $aRef = $conflictsAllRef->{$theKey}[$i];
			print OUT join("\t", $xId, 	"ref", $refId, $aRef->{refLeftBkpt}, $aRef->{refRightBkpt}, $aRef->{refAlignOrientation}, $aRef->{refLeftBkptStatus}, $aRef->{refRightBkptStatus}, $aRef->{refStatus})."\t";
			print OUT join("\t", 		"qry", $qryId, $aRef->{qryLeftBkpt}, $aRef->{qryRightBkpt}, $aRef->{qryAlignOrientation}, $aRef->{qryLeftBkptStatus}, $aRef->{qryRightBkptStatus}, $aRef->{qryStatus})."\n";
		}
	}
	close OUT;
}

sub checkSupportingPaint
{
	my ($conflictsAllRef, $paintedQryCmapRef, $paintedRefCmapRef, $max_overhang) = @_;

	my $shiftbp = 30;
	# d if there is no support @ seq junction and no support @ BN junction, then as already in conflictsAllRef
	# c if there is support @ seq junction and no support @ BN junction, then (if BN's current status is "okay" then turn both to "okay" <- pairmerge will ignore, but if current status is "cut" then leave BN as "cut")
	# b if there is no support @ seq junction and support @ BN junction, then (if BN's current status is "okay" then leave seq as "cut", but if current status "cut" then turn both to "okay" <- pairmerge will ignore) 
	# a if there is support @ seq junction and support @ BN junction, then ignore conflict, turn both to "okay" as pairmerge will ignore anyway

	foreach my $theKey (keys %$conflictsAllRef)	{
		my ($xId, $refId, $qryId) = split(/\t/, $theKey);
		for (my $i = 0; $i < scalar(@{$conflictsAllRef->{$theKey}}); $i += 1)	{
			die "ERROR: cut_conflicts: checkSupportingPaint: cannot find paint information for qry id = $qryId\n" if (! exists $paintedQryCmapRef->{$qryId});
			die "ERROR: cut_conflicts: checkSupportingPaint: cannot find paint information for ref id = $refId\n" if (! exists $paintedRefCmapRef->{$refId});
			my $cRef = $conflictsAllRef->{$theKey}[$i];
			my $sign = ($cRef->{qryAlignOrientation} eq '+') ? 1 : -1;
			my ($rRef, $qRef) = ($paintedRefCmapRef->{$refId}, $paintedQryCmapRef->{$qryId});
			my ($refLeftBkpt, $qryLeftBkpt, $refRightBkpt, $qryRightBkpt) = ($cRef->{refLeftBkpt}, $cRef->{qryLeftBkpt}, $cRef->{refRightBkpt}, $cRef->{qryRightBkpt});
			if (! ($cRef->{refLeftBkptStatus} eq "okay" && $cRef->{qryLeftBkptStatus} eq "okay"))	{
				# now check the paint colour at the left breakpoint (remember to shift up the left co-ordinates -- as it was shifted down in AssignAlignType)
				($cRef->{refLeftBkptStatus}, $cRef->{qryLeftBkptStatus}) = conflictResolveWithAltSupport($refId, $qryId, $refLeftBkpt + $shiftbp, $qryLeftBkpt + $sign * $shiftbp, $cRef->{refLeftBkptStatus}, $cRef->{qryLeftBkptStatus}, $rRef, $qRef, $max_overhang);
			} # if we need to deal with left bkpt
			if (! ($cRef->{refRightBkptStatus} eq "okay" && $cRef->{qryRightBkptStatus} eq "okay"))	{
				# remember to shift down the right co-ordinates -- as they were shifted up in AssignAlignType
				($cRef->{refRightBkptStatus}, $cRef->{qryRightBkptStatus}) = conflictResolveWithAltSupport($refId, $qryId, $refRightBkpt - $shiftbp, $qryRightBkpt - $sign * $shiftbp, $cRef->{refRightBkptStatus}, $cRef->{qryRightBkptStatus}, $rRef, $qRef, $max_overhang);
			} # if we need to deal with the right bkpt
		}
	}
	return $conflictsAllRef;
}

sub conflictResolveWithAltSupport
{
	my ($refId, $qryId, $refBkpt, $qryBkpt, $refBkptStatus, $qryBkptStatus, $paintedRefCmapIdRef, $paintedQryCmapIdRef, $max_overhang) = @_;
	my $refOIndeciesRef = searchForOverlap($refBkpt, $refBkpt, $paintedRefCmapIdRef);	# should have only 1 matching that co-ordinate

	die "ERROR: cut_conflicts conflictResolveWithAltSupport: cannot find the bkpt co-ordinate=$refBkpt for id=$refId in the paintedRefIdCmap array\n" if (scalar(@$refOIndeciesRef) == 0);
	my ($indexStart, $indexEnd) = (($refOIndeciesRef->[0] - $max_overhang < 0) ? (0) : ($refOIndeciesRef->[0] - $max_overhang), 
					($refOIndeciesRef->[0] + $max_overhang > $#$paintedRefCmapIdRef) ? ($#$paintedRefCmapIdRef) : ($refOIndeciesRef->[0] + $max_overhang));
	my ($altSupportForReferenceFlag, $altSupportForReferenceId) = findSupportingPaint($paintedRefCmapIdRef, $indexStart, $indexEnd, $qryId);	# find alt support spanning the conflict region +/- overhang label
	my $qryOIndeciesRef = searchForOverlap($qryBkpt, $qryBkpt, $paintedQryCmapIdRef);
	die "ERROR: cut_conflicts conflictResolveWithAltSupport: cannot find the bkpt co-ordinate=$qryBkpt for id=$qryId in the paintedQryIdCmap array\n" if (scalar(@$qryOIndeciesRef) == 0);
	($indexStart, $indexEnd) = (($qryOIndeciesRef->[0] - $max_overhang < 0) ? (0) : ($qryOIndeciesRef->[0] - $max_overhang),
					($qryOIndeciesRef->[0] + $max_overhang > $#$paintedQryCmapIdRef) ? ($#$paintedQryCmapIdRef) : ($qryOIndeciesRef->[0] + $max_overhang));
	my ($altSupportForQueryFlag, $altSupportForQueryId) = findSupportingPaint($paintedQryCmapIdRef, $indexStart, $indexEnd, $refId);	# find alt support spanning the conflict region +/- overhang label

	print "conflictResolveWithAltSupport: conflict between refId=$refId at bkpt=$refBkpt VS qryId=$qryId at bkpt=$qryBkpt: ";
	if ($altSupportForReferenceFlag == 1 && $altSupportForQueryFlag == 1)	{
		# case a
		($refBkptStatus, $qryBkptStatus) = ("okay", "okay");	# ignore conflict
		print "both have alt support: reference supported by altId=$altSupportForReferenceId, and query supported by altId=$altSupportForQueryId; so ignore conflict\n";
	} elsif ($altSupportForReferenceFlag == 0 && $altSupportForQueryFlag == 1)	{
		# case b
		print "reference has no alt support but query supported by altId=$altSupportForQueryId; ";
		if ($qryBkptStatus eq "okay")	{
			($refBkptStatus, $qryBkptStatus) = ("cut", "okay");
			print "query has molecule support, so reference is cut\n";
		} else	{
			($refBkptStatus, $qryBkptStatus) = ("okay", "okay");	# ignore conflict
			print "query has no molecule support, so ignore conflict\n";	
		}
	} elsif ($altSupportForReferenceFlag == 1 && $altSupportForQueryFlag == 0)	{
		# case c
		print "reference supported by altId=$altSupportForReferenceId; ";
		if ($qryBkptStatus eq "okay")	{
			($refBkptStatus, $qryBkptStatus) = ("okay", "okay");	# ignore conflict
			print "query has molecule support, so ignore conflict\n";
		} else	{
			($refBkptStatus, $qryBkptStatus) = ("okay", "cut");
			print "query has no molecule support, so query is cut\n";
		}
	} else	{
		# case d
		print "reference has no alt support and query has no alt support; leave reference as $refBkptStatus and query as $qryBkptStatus\n";
	}
	return ($refBkptStatus, $qryBkptStatus);
}

sub findSupportingPaint
{
	my ($paintedCmapRef, $indexStart, $indexEnd, $conflictPartnerId) = @_;
	my $foundAltSupport = 1;
	# look for a better alignment with another partner spanning the entire indexStart to indexEnd region (encompassing the conflict region, of course)
	my $firstLabelAlignPartnerId = $paintedCmapRef->[$indexStart]{topAlignPartnerId};
	$foundAltSupport = 0 if ($firstLabelAlignPartnerId == $conflictPartnerId || $firstLabelAlignPartnerId == -1);	# first label supported by an alternate partner?
	for (my $i = $indexStart + 1; $i <= $indexEnd && $foundAltSupport == 1; $i += 1)	{
		$foundAltSupport = 0 if ($paintedCmapRef->[$i]{topAlignPartnerId} != $firstLabelAlignPartnerId);	# still the same alternate partner?
	}
	return ($foundAltSupport, $firstLabelAlignPartnerId);
}

sub updateConflictsAll
{
	my ($conflictsQueryRef, $conflictsAllRef) = @_;
	foreach my $id (keys %$conflictsQueryRef)	{
		for (my $i = 0; $i < scalar(@{$conflictsQueryRef->{$id}}); $i += 1)	{
			my $cRef = $conflictsQueryRef->{$id}[$i];
			my ($xId, $refId, $qryId) = ($cRef->{xId}, $cRef->{refId}, $id);
			# now loop through conflictsAll 
			die "ERROR: updateConflictsAll: cannot find entry for xId=$xId; refId=$refId; qryId=$qryId in original breakpoint information\n" if (! exists $conflictsAllRef->{"$xId\t$refId\t$qryId"});
			for (my $a = 0; $a < scalar(@{$conflictsAllRef->{"$xId\t$refId\t$qryId"}}); $a += 1)	{
				my $aRef = $conflictsAllRef->{"$xId\t$refId\t$qryId"}[$a];
				# check if the query breakpoints are the same
				if (equal($aRef->{qryLeftBkpt}, $cRef->{oriStart}, 3))	{
					# left breakpoint
					if ($cRef->{toCut} == 1)	{
						# the left breakpoint was genome map at fault; label the leftBkpt Status with genome map being toCut
						$aRef->{qryLeftBkptStatus} = "cut";
					} else	{
						# the left breakpoint was sequence at fault; label the leftBkpt status with sequence being toCut
						$aRef->{refLeftBkptStatus} = "cut";
					}
				}
				if (equal($aRef->{qryRightBkpt}, $cRef->{oriStart}, 3))	{
					# right breakpoint
					if ($cRef->{toCut} == 1)	{
						$aRef->{qryRightBkptStatus} = "cut";	
					} else	{
						$aRef->{refRightBkptStatus} = "cut";
					}
				}
			}	
		}
	}
	return $conflictsAllRef;
}

sub equal
{
	# this subrountine checks whether two decimal numbers are the same (up to $decimalPlace precision)
	my ($num1, $num2, $decimalPlace) = @_;
	return sprintf("%.${decimalPlace}g", $num1) eq sprintf("%.${decimalPlace}g", $num2);
}

sub flagCut
{
	# this subroutine loops through all breakpoints and determine if there is any label whose support is below the cutThreshold percentage value
	my ($conflictsRef, $min_quality, $min_coverage) = @_;
	foreach my $id (sort { $a <=> $b } keys %$conflictsRef)	{
		for (my $i = 0; $i < scalar(@{$conflictsRef->{$id}}); $i += 1)	{
			my $toCut = 0;
			for (my $j = 0; $j < scalar(@{$conflictsRef->{$id}[$i]{neighbourQScoreRef}}) && $toCut == 0; $j += 1)	{
				print "flagCut: id=$id, position=$conflictsRef->{$id}[$i]{oriStart}, coverage=$conflictsRef->{$id}[$i]{neighbourCovRef}[$j], chim score=$conflictsRef->{$id}[$i]{neighbourQScoreRef}[$j]\n";
				if ($conflictsRef->{$id}[$i]{neighbourCovRef}[$j] < $min_coverage || $conflictsRef->{$id}[$i]{neighbourQScoreRef}[$j] < $min_quality)       {
					$toCut = 1;
				}
			}
			$conflictsRef->{$id}[$i]{toCut} = $toCut;
			print "\ttoCut=$toCut\n";
		}
	}
	return $conflictsRef;
}

sub findOverlap
{
	my ($data1Ref, $data2Ref) = @_;
	foreach my $id (keys %$data1Ref)        {
		die "ERROR: findOverlap: cannot find quality scores information for genome map = $id\n" if (! exists $data2Ref->{$id});
		for (my $i = 0; $i < scalar(@{$data1Ref->{$id}}); $i += 1)      {
			my $d1Ref = $data1Ref->{$id}[$i];
			my $oIndeciesRef = searchForOverlap($d1Ref->{start}, $d1Ref->{end}, $data2Ref->{$id});  
			for (my $o = 0; $o < scalar(@$oIndeciesRef); $o += 1)   {
			        # record the coverages and chimeric quality scores of the region
			        my $d2Ref = $data2Ref->{$id}[$oIndeciesRef->[$o]];
			        push(@{$d1Ref->{neighbourCovRef}}, $d2Ref->{coverage});
			        push(@{$d1Ref->{neighbourQScoreRef}}, $d2Ref->{chimQuality});
			}
		}
	}
	return $data1Ref;
}

sub searchForOverlap
{
	my ($qStart, $qEnd, $dataChromRef) = @_;
	my @oIndecies = ();
	my $bot = 0;
	my $top = $#{$dataChromRef};
	while($top > $bot + 1){
		my $mid = int(($top + $bot) / 2);
		if ($dataChromRef->[$mid]{start} < $qStart)       {
			$bot = $mid;
		} else  {
			$top = $mid;
		}
	}

	for (my $i = $bot; $i < scalar(@$dataChromRef); $i += 1)     {
		my $dRef = $dataChromRef->[$i];
		last if ($dRef->{start} > $qEnd);
		# any overlap
		if ( $qStart <= $dRef->{start} )   {
			push(@oIndecies, $i);
		}
	}
	return \@oIndecies;
}

sub getQScores
{
	my ($file) = @_;
	my %qScores = ();
	# read in the genome map cmap file, ASSUMING that the chimQuality is on column 10
	open(IN, "$file") or die "ERROR: getQScores: cannot open $file: $!\n";
	my $foundChimQuality = -1;
	while (my $line = <IN>)	{
		chomp $line;
		if ($line =~ /^#/)	{
			# assumes that data line comes after header lines
			if ($line =~ /^#h\s+/)	{
				# find that column storing the ChimQuality
				$line =~ s/^#h\s+//;
				my @headerContent = split(/\t/, $line);
				for (my $i = 0; $i < scalar(@headerContent); $i += 1)	{
					$foundChimQuality = $i if ($headerContent[$i] =~ /ChimQuality/i);
				}
			}
			next;
		}
		$line =~ s/^\s+/\t/;    $line =~ s/\s+/\t/g;
		my @content = split(/\t/, $line);
		if ($foundChimQuality == -1 || $content[$foundChimQuality] !~ /\d+/)	{
			warn "cut_conflicts.pl: getQScores: cmap file=$file does not have chimeric quality score\n";
			%qScores = (); # empty any qScores currently stored, and just flag this file as having no qScores
			return (\%qScores, 1);
		}
		my $labelChannel = $content[4];
		# skip if label channel is 0    (line indicating contig length)
		next if ($labelChannel == 0);
		my ($id, $start, $coverage, $chimQuality) = ($content[0], $content[5], $content[7], $content[$foundChimQuality]);
		push(@{$qScores{$id}}, {start => $start, coverage => $coverage, chimQuality => $chimQuality});
	}
	close IN;
	return (\%qScores, 0);
}

sub sortByCoord
{
	my ($dataRef, $coord) = @_;
	foreach my $id (keys %$dataRef) {
		@{$dataRef->{$id}} = sort       {
			$a->{$coord}    <=>     $b->{$coord}
		} @{$dataRef->{$id}}
	}
	return $dataRef;
}

sub extendByWindowSize
{
	my ($conflictsRef, $window_size) = @_;
	foreach my $id (keys %$conflictsRef)    {
		for (my $i = 0; $i < scalar(@{$conflictsRef->{$id}}); $i += 1)  {
			my $hRef = $conflictsRef->{$id}[$i];
			$hRef->{start} = ($hRef->{oriStart} - $window_size < 1) ? (1) : ($hRef->{oriStart} - $window_size);
			$hRef->{end} = $hRef->{oriStart} + $window_size;
		}
	}
	return $conflictsRef;
}

sub paintAlign1Cmap
{
	# this subroutine paints the labels on the cmap with the alignedPartner's id
	# the alignedPartner's with the lowest alignment confidence is painted first, then the better partner's id is painted next (potentially over-painting the initial partner's id)
	# and so on, the bext alignPartner's id would be painted last, and thus its id should be the "surface" paint
	my ($cmapRef, $xmapRef) = @_;
	foreach my $id (keys %$xmapRef)	{
		die "ERROR: cannot find any cmap information, yet alignment information is available for contig id = $id\n" if (! exists $cmapRef->{$id});
		for (my $i = 0; $i < scalar(@{$xmapRef->{$id}}); $i += 1)	{
			# now paint the corresponding labels
			my $xRef = $xmapRef->{$id}[$i];
			my ($alignmentStartLabel, $alignmentEndLabel, $partnerId) = ($xRef->{alignmentStartLabel}, $xRef->{alignmentEndLabel}, $xRef->{partnerId});
			my ($alignmentStartLabelIndex, $alignmentEndLabelIndex) = ($alignmentStartLabel - 1, $alignmentEndLabel - 1);	# since array indecies are 0-base
			die "ERROR: cut_conflicts, paintAlign1Cmap: alignment start and end label indecies are invalid, ie. out of bound of the cmap array\n\tcmap id = $id, alignmentStartLabel = $alignmentStartLabel, and alignmentEndLabel = $alignmentEndLabel" if (!(0 <= $alignmentStartLabelIndex && $alignmentStartLabelIndex <= $alignmentEndLabelIndex && $alignmentEndLabelIndex < scalar(@{$cmapRef->{$id}})));

			for (my $j = $alignmentStartLabelIndex; $j <= $alignmentEndLabelIndex; $j += 1)	{
				$cmapRef->{$id}[$j]{topAlignPartnerId} = $xRef->{partnerId};			# paint the label with its alignment partner id	
			}
		}
	}
	return $cmapRef;
}

sub sortXmap4PaintByConfScore
{
	# this subroutine sort the alignment match groups by the confScore (in ascending order)
	my ($xmapRef) = @_;
	foreach my $id (keys %$xmapRef)	{
		@{$xmapRef->{$id}}	=	sort	{
			$a->{confScore}	<=>	$b->{confScore}
		} @{$xmapRef->{$id}}
	}
	return $xmapRef;
}

sub readAlign1Xmap4Paint
{
	my ($file) = @_;
	my %xmap4Ref2Paint = ();	# remember ref is sequence
	my %xmap4Qry2Paint = ();	# remember qry is BioNano genome maps
	# xmap4Ref2Paint{$refId}[numAlignment] = {partnerId, confScore, alignmentStartLabel (wrt reference), alignmentEndLabel (wrt reference)}
	# xmap4Qry2Paint{$refId}[numAlignment] = {partnerId, confScore, alignemntStartLabel (wrt query), alignmentEndLabel (wrt query)}

	open(IN, $file) or die "ERROR: cannot open file $file for reading\n";
	while (my $line = <IN>)	{
		chomp $line;
		next if ($line =~ /^#/);	# skip header
		my @content = split(/\t/, $line);
		my ($qryId, $refId, $qStart, $qEnd, $rStart, $rEnd, $confScore, $alignmentString) = ($content[1], $content[2], $content[3], $content[4], $content[5], $content[6], $content[8], $content[13]);

		my ($refStartAlignLabel, $qryStartAlignLabel, $refEndAlignLabel, $qryEndAlignLabel) = (-1, -1, -1, -1);
		if ($alignmentString =~ /^\((\d+),(\d+)\)/)	{
			($refStartAlignLabel, $qryStartAlignLabel) = ($1, $2);
		}
		if ($alignmentString =~ /\((\d+),(\d+)\)$/)	{
			($refEndAlignLabel, $qryEndAlignLabel) = ($1, $2);
		}
		if ($refStartAlignLabel == -1 || $qryStartAlignLabel == -2 || $refEndAlignLabel == -1 || $qryEndAlignLabel == -1)	{
			die "ERROR: cannot parse alignment string: $alignmentString\n\tin line $line\n";
		}

		if ($qryEndAlignLabel < $qryStartAlignLabel)	{
			my $temp = $qryEndAlignLabel;
			$qryEndAlignLabel = $qryStartAlignLabel;
			$qryStartAlignLabel = $temp;
		}
		if ($refEndAlignLabel < $refStartAlignLabel)	{
			# this should not happen, that is, reference alignment start should ALWAYS be smaller than alignment end
			my $temp = $refEndAlignLabel;
			$refEndAlignLabel = $refStartAlignLabel;
			$refStartAlignLabel = $temp;
		}

		push(@{$xmap4Ref2Paint{$refId}}, {confScore => $confScore, alignmentStartLabel => $refStartAlignLabel, alignmentEndLabel => $refEndAlignLabel, partnerId => $qryId});
		push(@{$xmap4Qry2Paint{$qryId}}, {confScore => $confScore, alignmentStartLabel => $qryStartAlignLabel, alignmentEndLabel => $qryEndAlignLabel, partnerId => $refId});
	}
	close IN;
	return (\%xmap4Ref2Paint, \%xmap4Qry2Paint);
}

sub readAlign1Cmap2Paint
{
	my ($file) = @_;
	# assuming that cmap entries are sorted according to their id, and then their position
	my %cmap2Paint = ();
	open(IN, $file) or die "ERROR: cannot open file $file for reading.\n";
	while (my $line = <IN>)	{
		chomp $line;
		next if ($line =~ /^#/);
		my @content = split(/\t/, $line);
		my ($id, $siteId, $labelChannel, $position) = ($content[0], $content[3], $content[4], int($content[5]));
		# cmap length line
		next if ($labelChannel == 0);	# skip at labelChannel 0 lines
		# note that the array index should be $siteId - 1
		# cmap{id}[$labelNumIndex] = {start, end, labelNum, topAlignPartnerId} for both reference and query
		push(@{$cmap2Paint{$id}}, {start => $position, end => $position, siteId => $siteId, topAlignPartnerId => -1});
	}
	close IN;
	return \%cmap2Paint;
}

# this subroutine extracts the break point information from a given file
# it stores conflicts (wrt the qId) into conflictsQuery to then cross check with chimeric quality scores
# it stores conflicts (wrt to xMapId\trId\tqId triplet)
sub getConflicts
{
	my ($file) = @_;
	my @headerLines = ();
	my %conflictsQuery = ();
	my %conflictsAll = ();
	open(IN, "$file") or die "ERROR: cut_conflicts, getConflicts: cannot open file $file: $!\n";
	while (my $line = <IN>)	{
		chomp $line;
		if($line =~ /^#/){
			push(@headerLines, $line) unless($line !~ /^#\s+xMapId/);
			next;
		}
		# store the conflict information into the conflictsQuery hash and conflictsAll hash
		my @content = split(/\t/, $line);
		my $xId = $content[0];
		my ($refId, $refLeftBkpt, $refRightBkpt, $refAlignOrientation) = ($content[2], int($content[3]), int($content[4]), $content[5]);
		my ($qryId, $qryLeftBkpt, $qryRightBkpt, $qryAlignOrientation) = ($content[7], int($content[8]), int($content[9]), $content[10]);
		push(@{$conflictsAll{"$xId\t$refId\t$qryId"}}, {refLeftBkpt => $refLeftBkpt, qryLeftBkpt => $qryLeftBkpt, 
			refRightBkpt => $refRightBkpt, qryRightBkpt => $qryRightBkpt,
			refLeftBkptStatus => "okay", refRightBkptStatus => "okay", refStatus => "okay",
			qryLeftBkptStatus => "okay", qryRightBkptStatus => "okay", qryStatus => "okay",
			refAlignOrientation => $refAlignOrientation, qryAlignOrientation => $qryAlignOrientation});
		
		# only store conflicts in conflictsQuery
		if ($qryLeftBkpt != -1)	{
			push(@{$conflictsQuery{$qryId}}, {oriStart => $qryLeftBkpt, start => -1, end => -1, xId => $xId, refId => $refId, line => $line, neighbourCovRef => [], neighbourQScoreRef => [], toCut => 0});
		}
		if ($qryRightBkpt != -1)	{
			push(@{$conflictsQuery{$qryId}}, {oriStart => $qryRightBkpt, start => -1, end => -1, xId => $xId, refId => $refId, line => $line, neighbourCovRef => [], neighbourQScoreRef => [], toCut => 0});
		}
	}

	close IN;

	return (\@headerLines, \%conflictsQuery, \%conflictsAll);
}
