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
use Math::Round;
use Cwd qw(abs_path);

BEGIN{
	my $progpath = abs_path(dirname($0));
	unshift(@INC, $progpath);
	select(STDERR); $| = 1;
	select(STDOUT); $| = 1;
}

use BNG::Utility;
use Getopt::Std;
use Tie::File;

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for exporting AGP file and FASTA file, given relevant NGS and BNG information
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-i, --input <str>  The input XMAP file (REQUIRED)
	-o, --output <str> The output path containing the output directory and prefix (REQUIRED)
	-c, --cmap <str>   The input CMAP file from hybrid assembly (REQUIRED)
	-s, --seq <str>    The input FASTA file from NGS assembly (REQUIRED)
	-m, --map <str>    The file containing the mapping information from ID to sequence name (REQUIRED)
	-t, --trans <str>  The file containing the coordinate transformation information when resolving conflicts
	-g, --gap <int>    The length of gap to be inserted between overlapping NGS contigs in a scaffold (default: 13)
	-h, --help         Help
USAGE

#global varaibles

#enzyme cut site sequences
our %enzyme = (
	"BSPQI" => "GCTCTTC",
	"BSSSI" => "CACGAG",
	"BBVCI" => "CCTCAGC",
	"BSMI"  => "GAATGC",
	"BSRDI" => "GCAATG",
	"BSECI" => "ATCGAT",
	"BAMHI" => "GGATCC",
	"DLE1" => "CTTAAG"
	# You can add more enzymes here ...
	
);

use vars qw($opt_i $opt_o $opt_c $opt_s $opt_m $opt_t $opt_g $opt_h);
if(!GetOptions( "i|input=s" => \$opt_i,
			"o|output=s" => \$opt_o,
			"c|cmap=s" => \$opt_c,
			"s|seq=s" => \$opt_s,
			"m|map=s" => \$opt_m,
			"t|trans=s" => \$opt_t,
			"g|gap=i" => \$opt_g,
			"h|help" => \$opt_h)){
	print ("Please try -h for more details\n");
	exit(1);
}

die($usage) if($opt_h);

# check the required arguments
die("**ERROR: -i option must be specified\n") if(!defined $opt_i);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);
die("**ERROR: -c option must be specified\n") if(!defined $opt_c);
die("**ERROR: -s option must be specified\n") if(!defined $opt_s);
die("**ERROR: -m option must be specified\n") if(!defined $opt_m);

# get the XMAP file name
my $xmap_file = $opt_i;

# set the output directory and prefix
my ($outdir, $outpre);
if($opt_o =~ /\/$/ or $opt_o !~ /\//){
	$opt_o =~ s/\/$//;
	$outdir = $opt_o;
	$outpre = basename($opt_i); $outpre =~ s/\.xmap$//;
}
else{
	$outdir = dirname($opt_o);
	$outpre = basename($opt_o);
}
my ($cmd, $retCode);
$cmd = "mkdir -p $outdir";
$retCode = system($cmd);
die("**ERROR: Can not create directory \"$outdir\"\n") if($retCode != 0);

# get the other arguments
my $fasta_file = $opt_s;
my $hybrid_cmap_file = $opt_c;
my $NCBI_cut_site = "N";

my $ngs_namemap_file = $opt_m;
my $cut_coord_file = $opt_t;
my $padding_gap_len = ((!defined $opt_g) ? 13 : (($opt_g < 1) ? 13 : $opt_g));
my $padding_gap = ("N" x $padding_gap_len);

# set the output file names
my $agp_out_file = "$outdir/$outpre.agp";
my $begin_end_file = "$outdir/${outpre}_trimHeadTailGap.coord";
my $fasta_out_file = "$outdir/$outpre.fasta";
my $NCBI_fasta_out_file = "$outdir/${outpre}_NCBI.fasta";
my $unused_fasta_out = "$outdir/${outpre}_NOT_SCAFFOLDED.fasta";
my $gap_out_file = "$outdir/${outpre}.gap";

###########
our $sorted_xmap_file = &sortXMap($xmap_file, "$outdir/${outpre}_sorted.xmap");
our $xmap = readXMap($sorted_xmap_file);

if($cut_coord_file){
	print "Detect cut-coord file, updating fast and fasta name map file\n";
	($fasta_file, $ngs_namemap_file) = &updateNGSFiles($fasta_file, $ngs_namemap_file, $cut_coord_file, $outdir); 
}	

our $ngsMap = &getNGSMap($ngs_namemap_file);
our $seqMap = readFasta($fasta_file, $outdir);


our ($hybridCmap, $numcontig, $contigLength) = readCMap($hybrid_cmap_file);
my $cut_site = $hybridCmap->{channels}->{1};

processAlign($xmap, $gap_out_file, $agp_out_file, $begin_end_file, $padding_gap_len, $ngsMap);
printFasta($seqMap, $xmap, $fasta_out_file, $hybridCmap, $cut_site, $padding_gap);
printFasta($seqMap, $xmap, $NCBI_fasta_out_file, $hybridCmap, $NCBI_cut_site, $padding_gap);
&printUnUsedNGS($agp_out_file, $ngsMap, $begin_end_file, $unused_fasta_out, $seqMap);
&cleanUp($outdir);

#############Function for pre-processing inputs for the exporter ###################################
#sort xmap by refID and increasing refStartPos and decreasing refEndPos
#this simpifly checking if one contigs is completely embedded in another contigs
sub sortXMap
{
	my ($xmapFile, $sorted_xmap) = @_;
	my $headerCnt = `grep -c '^#' $xmapFile`; chomp $headerCnt;

	`head -n $headerCnt $xmapFile > $sorted_xmap`;  #export header to sorted xmap 
     $headerCnt = $headerCnt+1; #set to appropriate param for tail
	`tail -n+$headerCnt $xmapFile | sort -k3,3n -k6,6n -k7,7nr >> $sorted_xmap`; # RefContigID, RefStartPos, RefEndPos

	return($sorted_xmap);
}

#parse the ngs mapping file to create name mapping between contig id in cmap and xmap and its ngs name and other stats
sub getNGSMap
{
	my ($name_map_file, $cut_coord_file) = @_;
	my %name_map=();
	
	open(NAMEMAP, $name_map_file) || die("ERROR: cannot open sequence name map file $name_map_file");
	while(my $line = <NAMEMAP>){
		chomp $line;
		if($line !~ /^\d/){
			next;
		}
		my @tokens = split(/\t/, $line);
		#we store the ngs stat as 1) name; 2) length; 3) a flag whether it is used scaffolding; 4) the begin coord 
		#of this contig in case it is split by the pipline into multiple smaller contigs
		my $name = (split(/\s+/, $tokens[1]))[0];
		$name_map{$tokens[0]} = [$name, $tokens[2], 0, 1];
	}
	if(!defined $cut_coord_file || length($cut_coord_file) == 0){
		return \%name_map;
	}
	my %new_map=();
	print "cut_coord_file is: $cut_coord_file\n";
	open(CUT_COORD, $cut_coord_file) || die("ERROR: cannot open cut coordinate files $cut_coord_file\n");
	while(my $line = <CUT_COORD>){
		chomp $line;
		if($line =~ /^oldId/){ #skipping header
			next;
		}
		my @tokens = split(/\t/, $line);
		if(exists $name_map{$tokens[0]}){
			my $newname = $name_map{$tokens[0]}->[0];
			my $newlength = $tokens[2] - $tokens[1] + 1; 
			$new_map{$tokens[3]} = [$newname, $newlength, 0, $tokens[1]+1];
		}else{
			warn "Warning map with Id: $tokens[0] is not recognized to be a valide ID, check $cut_coord_file or $name_map_file\n";
		}
	}
	return \%new_map;
}

#when hybrid-scaffold cut ngs contigs to resolve conflicts
#we need to generate new fasta files and ngs name map files for 
#the cut ngs contigs
sub updateNGSFiles
{
	my ($fasta_file, $ngs_namemap_file, $cut_coord_file, $out_dir) = @_;
	
	my $cutted_fasta_file = "$out_dir/".basename($fasta_file).".cut.fasta";
	my $cutted_ngsNameMap_file = "$out_dir/".basename($ngs_namemap_file).".cut.txt";
	
	my %ngs_map = %{&getNGSMap($ngs_namemap_file, $cut_coord_file)};
	my $fasta_accessor = readFasta($fasta_file, $out_dir);

	open(my $cutted_fasta, ">".$cutted_fasta_file) || die("ERROR: cannot open file for writing: $cutted_fasta_file");
	open(my $cutted_ngs_nameMap, ">".$cutted_ngsNameMap_file) || die("ERROR: cannot open file for writing: $cutted_ngsNameMap_file");

	print $cutted_ngs_nameMap "CompntId\tCompntName\tCompntLength\n";
	my @ids = keys(%ngs_map);
	my @sort_ids = sort {$a <=> $b} @ids;
	#print(@sort_ids);
	foreach my $key(@sort_ids){
		my @ngscontig = @{$ngs_map{$key}};
		my $fasta_seq = $fasta_accessor->($ngscontig[0]);
		if($ngscontig[1] == length($fasta_seq)){
			print $cutted_fasta ">$ngscontig[0]\n";
			print $cutted_fasta $fasta_seq."\n";
			print $cutted_ngs_nameMap $key."\t".$ngscontig[0]."\t".$ngscontig[1]."\n";
			#print 			
		}else{
			my $new_ngs_name = $ngscontig[0]."_subseq_".$ngscontig[3].":".($ngscontig[3] + $ngscontig[1]-1);
			print $cutted_fasta ">$new_ngs_name\n";
			print $cutted_fasta substr($fasta_seq, $ngscontig[3] - 1, $ngscontig[1])."\n";
			print $cutted_ngs_nameMap $key."\t".$new_ngs_name."\t".$ngscontig[1]."\n";
		}	
	}
	print "DONE updating NGS files\n";
	close($cutted_fasta);
	close($cutted_ngs_nameMap);
	return ($cutted_fasta_file, $cutted_ngsNameMap_file);
}



#Read in a fasta file, return a function pointer which can be used 
#to access the fasta sequence by NGS name
sub readFasta
{
	my ($fasta_file, $out_dir) = @_;
	
	my ($fastaMap, $tmp_file) = processFastaFile($fasta_file, $out_dir);
	tie my @fasta_arry, 'Tie::File',  $tmp_file,  or die("ERROR: cannot process fasta file");
	my %fasta_map = %{$fastaMap};
	
	my $accesser = sub{
		my ($ID) = @_;
		if(exists $fasta_map{$ID}){
			my $ind = $fasta_map{$ID};
			return $fasta_arry[$ind];
		}else{
			print "Warning: cannot find sequence for ID: $ID\n";
			return "";
		}
	};	
	return $accesser;
}

#This function remove extra white space in a fasta file, so the sequence remains in one line
#It also creates an index table of seq name to the which line in the file
#corresponding to that sequence
sub processFastaFile
{
	my ($fasta_file, $out_dir) = @_;
	my $tmp_file = "$outdir/" . basename($fasta_file) . ".tmp.fa.tmp";
	open FASTA, $fasta_file || die("ERROR: cannot open fasta file");
	open TMP, ">".$tmp_file || die("ERROR: cannot access file system for temporary file");
	my %fastaTable=();
	my $curr="";
	my $currHeader="";
	my $numLines=1;
	my $line;
	while($line = <FASTA>){
		chomp $line;
		if($line =~ /^>/){
			if($curr ne ""){
				$fastaTable{substr($currHeader,1)}=$numLines;
				print TMP $currHeader."\n";
				print TMP $curr."\n";
				$numLines=$numLines+2;
			}
			my @tokens = split(/\s+/, $line);
			$currHeader = $tokens[0];		
			$curr="";
		}else{
			$curr=$curr.$line;
		}
	}
	#taking care of the last seq contig
	$fastaTable{substr($currHeader,1)}=$numLines;	
	print TMP $currHeader."\n";
	print TMP $curr."\n";
	
	close(FASTA);
	close(TMP);
	return (\%fastaTable, $tmp_file)
}

##########################Functions for processing the hybrid-scaffold to NGS contigs alignment (i.e. the xmap file) ##########################

#This function read-in the alignments btw hybrid-scaffold contigs and the NGS contigs and compute gap length
#bewteen each consecutive pair of scaffolded NGS contigs and output this information into a file

sub processAlign
{
	my ($xmap, $gapLen_out_file, $agp_out_file, $begin_end_file, $paddingGapLen, $ngs_map) = @_;
	$xmap = mapNGSStat($xmap, $ngs_map);
	$xmap = computeGapLengths($xmap);
	printGapFile($xmap, $gapLen_out_file);
	printAGPFile($xmap, $agp_out_file, $begin_end_file, $paddingGapLen, $ngs_map);
}

#This function maps the QryContigID in xmap to the ngs stats such as name and length 
#from the original ngs contigs 
sub mapNGSStat
{
	my ($xmap, $ngs_map) = @_;
	@{$xmap->{hits}->{NGSName}} = ("") x $xmap->{totalHits};
	@{$xmap->{hits}->{NGSLen}} = (0) x $xmap->{totalHits};
	@{$xmap->{hits}->{NGSStart}} = (1) x $xmap->{totalHits};
	for(my $i = 0; $i < $xmap->{totalHits}; $i++){
		my $ID = $xmap->{hits}->{QryContigID}[$i];
		if(exists $ngs_map->{$ID}){
			#print "NSG name: ".$ngs_map->{$ID}->[0]."\n"; 
			$xmap->{hits}->{NGSName}[$i] =$ngs_map->{$ID}->[0];
			$xmap->{hits}->{NGSLen}[$i] =$ngs_map->{$ID}->[1];
			$xmap->{hits}->{NGSStart}[$i] = $ngs_map->{$ID}->[3];
		}else{
			print "Warning: cannot find ngs name for QryContig ID: ".$ID.". Using stats from qry in Xmap instead.\n";
			#if we cannot find ngs for some reasons use stat from original xmap
			$xmap->{hits}->{NGSName}[$i] =$ngs_map->{hits}->{QryContigID}[$i];
			$xmap->{hits}->{NGSLen}[$i] =$ngs_map->{hits}->{QryLen}[$i];
		} 
	}
	return $xmap;
}

#adjust the ref start/end positions of alignments taking into account
#un-aligned sequence and direction of the query contig
sub adjustRefStart
{
	my($qryStart, $qryEnd, $refStart, $qryLen, $direction) = @_;
	return $refStart - beginGap($qryStart, $qryEnd, $qryLen, $direction);	
}

sub beginGap
{
	my($qryStart, $qryEnd, $qryLen, $direction) = @_;

	return ($direction eq '+') ? ($qryStart - 1) : ($qryLen - $qryStart);	
}

sub adjustRefEnd
{
	my($qryStart, $qryEnd, $refEnd, $qryLen, $direction) = @_;
	return $refEnd + EndGap($qryStart, $qryEnd, $qryLen, $direction);
}

sub EndGap
{
	my($qryStart, $qryEnd, $qryLen, $direction) = @_;

	return ($direction eq '+') ? ($qryLen - $qryEnd) : ($qryEnd - 1);
}

#compute gaps for a particular hybrid-scaffolded contig
#given the alignments of all ngs contigs to the hybrid
sub computeGapLengths
{
	my ($xmap) = @_;
	#the xmap data table will be used as the main table to store data
	#we added three-more columns to store the gap length information
	my $hits = $xmap->{hits};
	my $num_hits = $xmap->{totalHits};
	@{$hits->{GapLen}} = (0) x $num_hits;
	@{$hits->{AdjustedGapLen}} = (0) x $num_hits;
	@{$hits->{IsEmbedded}} = (0) x $num_hits;

	my ($currInd, $nextInd) = (0, 1);
	while($nextInd < $num_hits){
		if($hits->{RefContigID}[$nextInd] ne $hits->{RefContigID}[$currInd]){
			$currInd = $nextInd;
			$nextInd++;
			next;
		}
		my ($qryStartPos1, $qryEndPos1, $refStartPos1, $refEndPos1, $qryLen1, $dir1) = (
					$hits->{QryStartPos}[$currInd], $hits->{QryEndPos}[$currInd],
					$hits->{RefStartPos}[$currInd], $hits->{RefEndPos}[$currInd],
					$hits->{QryLen}[$currInd], $hits->{Orientation}[$currInd]
				);
		my ($qryStartPos2, $qryEndPos2, $refStartPos2, $refEndPos2, $qryLen2, $dir2) = (
					$hits->{QryStartPos}[$nextInd], $hits->{QryEndPos}[$nextInd],
					$hits->{RefStartPos}[$nextInd], $hits->{RefEndPos}[$nextInd],
					$hits->{QryLen}[$nextInd], $hits->{Orientation}[$nextInd]
				);
		
		my $alignGapLen = round($refStartPos2) - round($refEndPos1) - 1; 

		my $adjustedRefEnd1 = round(adjustRefEnd($qryStartPos1, $qryEndPos1, $refEndPos1, $qryLen1, $dir1));
		my $adjustedRefStart2 = round(adjustRefStart($qryStartPos2, $qryEndPos2, $refStartPos2, $qryLen2, $dir2));
		my $adjustedGapLen = $adjustedRefStart2 - $adjustedRefEnd1  - 1;
		
		#detecting contig that is completely covered by the previous contig
		#if this is the case we do not advance the current contig index until the next contig that is not embedded
		if(-1*$adjustedGapLen  > $qryLen2){
			$hits->{IsEmbedded}[$nextInd] = 1;
		}else{	
			$currInd = $nextInd;
		}
		# store the gap length information in the next contig rather than the first one to handle cases of multiple embedded contigs
		$hits->{GapLen}[$nextInd] = $alignGapLen;
		$hits->{AdjustedGapLen}[$nextInd] = $adjustedGapLen;
		$hits->{AdjustedGapBegin}[$nextInd] = $adjustedRefEnd1 + 1;
		$hits->{AdjustedGapEnd}[$nextInd] = $adjustedRefStart2 - 1;
		
		$nextInd++;
	}
	return $xmap;
}

###########################Functions for outputting agp/fasta or its auxilliary files###########################################

#print the gap informaton from alignments to a file
sub printGapFile
{
	my ($xmap, $out_put) = @_;
	open(my $fh, ">".$out_put);
	print $fh "#NGSId1\tNGSId2\tSuperScaffoldId\tXmapGapLength\tAdjustedGapLength\tNGSLength1\tNGSLength2\n";
	my $numAlign = $xmap->{totalHits};
	my $hits = $xmap->{hits};
	my ($prevInd, $currInd) = (0, 1);
	while($currInd < $numAlign){
		if( $hits->{RefContigID}[$prevInd] == $hits->{RefContigID}[$currInd]){
			print $fh $hits->{NGSName}[$prevInd]."\t".
				$hits->{NGSName}[$currInd]."\t".
				$hits->{RefContigID}[$currInd]."\t".
				$hits->{GapLen}[$currInd]."\t".
				$hits->{AdjustedGapLen}[$currInd]."\t".
				$hits->{QryLen}[$prevInd]."\t".
				$hits->{QryLen}[$currInd]."\n";
			if(!$hits->{IsEmbedded}[$currInd]){
				$prevInd = $currInd;
			}	
		}else{
			print $fh $hits->{NGSName}[$prevInd]."\t".
				$hits->{NGSName}[$prevInd]."\t".
				$hits->{RefContigID}[$prevInd]."\t".
				"0\t0\t".
				$hits->{QryLen}[$prevInd]."\t".
				$hits->{QryLen}[$prevInd]."\n";
				$prevInd = $currInd;
		}
		$currInd = $currInd + 1;
	}
}

#print header of AGP file
sub printAGPHeader
{
	my ($fh) = @_;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	my $agp_header = "##agp-version\t2.0\n".
				"# Organism:  \n".
				"# Platform:     \n".
				"# Model:        \n".
				"# Enzyme(s):    \n".
				"# BioSample:    \n".
				"# BioProject:   \n".
				"# Obj_Name\tObj_Start\tObj_End\tPartNum\tCompnt_Type\tCompntId_GapLength\tCompntStart_GapType\tCompntEnd_Linkage\tOrientation_LinkageEvidence\n";
	print $fh $agp_header;
}


sub printAGPAuxHeader
{
	my ($fh) = @_;
	my $header= "##agp-version\t2.0\n".
				"# Organism:   \n".
				"# Platform:     \n".
				"# Model:        \n".
				"# Enzyme(s):    \n".
				"# BioSample:    \n".
				"# BioProject:   \n".
				"Obj_Id\tHeadTrimmedLength\tTailTrimmedLength\n";
	print $fh $header;
}

#print the agp files
sub printAGPFile
{
	my ($xmap, $out_put, $aux_out_file, $dummyGapLen, $ngs_map)= @_;
	print "AGP output file: $out_put\n";
	open(my $fh, ">".$out_put);
	printAGPHeader($fh);
	open(my $aux_fh, ">".$aux_out_file);
	printAGPAuxHeader($aux_fh);
	my $numAlign = $xmap->{totalHits};
	my $currRefStart = 1;
	my $currRefEnd = 1;
	my $Ind = 1;
	
	#outputing ngs contigs
	my $currInd = 0;
	my $nextInd = 1;	
	while($currInd < $numAlign){
		#if contig is embbed inside another contig, skip this contig
		if($xmap->{hits}->{IsEmbedded}[$currInd]){
			#marking the ngs contig as embedded
			$ngs_map->{$xmap->{hits}->{QryContigID}[$currInd]}->[2] = -1;
			$currInd = $currInd + 1;
			next;
		}
		
		if($xmap->{hits}->{IsEmbedded}[$nextInd]){
			#marking the ngs contig as embeddedq
			$ngs_map->{$xmap->{hits}->{QryContigID}[$nextInd]}->[2] = -1;
			$nextInd = $nextInd + 1;
			next;
		}
					
		#printing begin gap for first contig
		#note first contig cannot be embedded so this condition is valid
		if($currInd == 0){
			my $begin_gap = $xmap->{hits}->{RefStartPos}[$currInd]-1;
			print $aux_fh "Super-Scaffold_".$xmap->{hits}->{RefContigID}[$currInd]."\t$begin_gap\t";
		}
		#printing ngs sequence
		#print "Index $i\n";
		$currRefEnd = $currRefStart + round($xmap->{hits}->{NGSLen}[$currInd]) - 1;
		print $fh "Super-Scaffold_".$xmap->{hits}->{RefContigID}[$currInd]."\t".
				"$currRefStart\t".	
				$currRefEnd."\t".
				"$Ind\tW\t".
				$xmap->{hits}->{NGSName}[$currInd]."\t".
				#"1\t".$xmap->{hits}->{QryLen}[$currInd]."\t".
				$xmap->{hits}->{NGSStart}[$currInd]."\t".($xmap->{hits}->{NGSStart}[$currInd]+$xmap->{hits}->{NGSLen}[$currInd]-1)."\t".
				$xmap->{hits}->{Orientation}[$currInd]."\n";
		$currRefStart = $currRefEnd + 1;
		#marking in the ngs table that this contig was scaffolded
		if(exists $ngs_map->{$xmap->{hits}->{QryContigID}[$currInd]}){
			$ngs_map->{$xmap->{hits}->{QryContigID}[$currInd]}->[2] = 1;
		}
		
		#printing gap		
		if($nextInd < $numAlign && $xmap->{hits}->{RefContigID}[$currInd] == $xmap->{hits}->{RefContigID}[$nextInd]){
			#my $gapLen = $xmap->{hits}->{GapLen}[$currInd];
			my $gapLen = $xmap->{hits}->{AdjustedGapLen}[$nextInd];
			
			#if two contigs begin/stop at the same label/position from alignment (by definition they have gap length of -1, contig right next to each other has gap length of zero), we assume it is not overlap
			#if($gapLen == -1){
			#	$gapLen = 0;
			#}
			
			if($gapLen < 0){
				$gapLen = $dummyGapLen;
			}elsif($gapLen < $dummyGapLen + 10){
				$gapLen = $dummyGapLen + 10
			}
			
			$currRefEnd = $currRefStart + $gapLen - 1;
			if($gapLen > 0){
				$Ind++;
				print $fh "Super-Scaffold_".$xmap->{hits}->{RefContigID}[$currInd]."\t".
					"$currRefStart\t".	
					"$currRefEnd\t".
					"$Ind\tN\t".$gapLen."\tscaffold\tyes\tmap\n";
			}
			$currRefStart = $currRefEnd + 1;		
		}else{								
			#handle beginning and end gap sequence
			#end gap for the previous hybrid contig
			my $end_gap = $xmap->{hits}->{RefLen}[$currInd] - $xmap->{hits}->{RefEndPos}[$currInd];
			#print "printing END gap: ".$xmap->{hits}->{RefLen}[$currInd]."\t". $xmap->{hits}->{RefEndPos}[$currInd]."\n";
			print $aux_fh round($end_gap)."\n"; 
			#begin gap for the current hybrid contig
			if($nextInd < $numAlign){ 
				my $begin_gap = $xmap->{hits}->{RefStartPos}[$nextInd]-1;
				print $aux_fh "Super-Scaffold_".$xmap->{hits}->{RefContigID}[$nextInd]."\t".round($begin_gap)."\t";
			}
			$currRefStart = 1;
			$currRefEnd = 0;
			$Ind = 0;
		}
		$currInd = $nextInd;
		$nextInd = $nextInd + 1;
		$Ind++;
	}
	close($fh);
	close($aux_fh);
}


#convert the interval objects to fasta file format
sub printFasta
{
	my ($fasta_map, $xmap, $out_file, $hybridCmap, $cut_site, $paddingGap) = @_;
	my $numAlign = $xmap->{totalHits};
	my $isBegin = 1;
	open(my $fh, ">".$out_file) || die "cannot open output fasta file $out_file\n";
	
	my $currInd = 0;
	my $nextInd = 1;
	
	my $minGap = $paddingGap; #min gap is padding + 10base
	for(my $i=0; $i < 10; $i++){
		$minGap = $minGap."N";
	} 
		
	while($currInd < $numAlign){
		#if contig is embbed inside another contig, skip this contig
		if($xmap->{hits}->{IsEmbedded}[$currInd]){
			$currInd = $currInd + 1;
			next;
		}
		
		if($xmap->{hits}->{IsEmbedded}[$nextInd]){
			$nextInd = $nextInd + 1;
			next;
		}

		my $numAlign = $xmap->{totalHits};
		#if contig is embbed inside another contig, skip this contig
		if($xmap->{hits}->{IsEmbedded}[$currInd]){
			next;
		}
		#handling the first line
		if($isBegin){
			print $fh ">Super-Scaffold_".$xmap->{hits}->{RefContigID}[$currInd]."\n";
			$isBegin=0;
		}
		
		#printing ngs sequence
#		print "NGSname $xmap->{hits}->{NGSName}[$currInd]\n";
		my $ngsseq = $fasta_map->($xmap->{hits}->{NGSName}[$currInd]);
		if($xmap->{hits}->{Orientation}[$currInd] eq '+'){
			print $fh $ngsseq;
		}else{
			print $fh reverseComplement($ngsseq);
		}

		#printing gap		
		if($nextInd < $numAlign && $xmap->{hits}->{RefContigID}[$currInd] == $xmap->{hits}->{RefContigID}[$nextInd]){
			#my $gapBegin = round($xmap->{hits}->{RefEndPos}[$currInd]);
			#my $gapEnd = round($xmap->{hits}->{RefStartPos}[$nextInd])-1; #gap end at one b.p.before next start
			my $gapBegin = round($xmap->{hits}->{AdjustedGapBegin}[$nextInd]);
			my $gapEnd = round($xmap->{hits}->{AdjustedGapEnd}[$nextInd]);
			my $gaplen = $xmap->{hits}->{AdjustedGapLen}[$nextInd];
			
			my $RefContigID = $xmap->{hits}->{RefContigID}[$currInd];
			my $gapSeq = $paddingGap;
			if($gaplen >=0){
				if($gaplen < length($minGap)){
					$gapSeq = $minGap
				}else{
			    	#if($gaplen >= 0){  
					$gapSeq = getSeqFromCmap([$gapBegin, $gapEnd], $hybridCmap, $RefContigID, $cut_site);
				}
			}
			#when neighboring contigs begin and end at same position, we assume they are not overlap (this by our current definition translate to gaplen of -1)
			#if($gaplen == -1){
			#    $gapSeq = "";
			#}
			print $fh $gapSeq;
		}else{								
			print $fh "\n";	
			if($nextInd < $numAlign){
				print $fh ">Super-Scaffold_".$xmap->{hits}->{RefContigID}[$nextInd]."\n";
			}
		}
		$currInd = $nextInd;
		$nextInd = $nextInd + 1;
		
	}
	close($fh);
}


#printing out a sequence accoriding to Cmap coordinate
#This is essentially a N-sequence with location of cut sites
#specified by the cut-site sequences of the enzyme/enzymes

sub getSeqFromCmap
{
	my ($interval, $cmap, $contigID, $cut_site) = @_;
	if($interval->[0] > $interval->[1]){
		print STDERR "Warning: Begin index is not smaller than end Index in interval: \[$interval->[0], $interval->[1]\]\n";
		return "";
	}
	if(!exists $cmap->{contigs}->{$contigID}){
		print "Cannot find hybrid scaffold Id ".$contigID." cmap file\n";
		return "";
	}
	my $contig =$cmap->{contigs}->{$contigID};
	my $seq="";
	my ($begin, $end) = getCoordinatesSubset($interval, $contig);

	#print "begin $begin end $end\n";
	my $coord = $interval->[0];
	for(my $i = $begin; $i < $end; $i++){
		#print "end pos: ".($contig->{Position}[$i])."\n"; 
		for(; $coord < $contig->{Position}[$i]; $coord++){
			$seq = $seq."N";
		}
		$seq=$seq.$cut_site;
		$coord = $coord+length($cut_site);
	}
	for(; $coord <= $interval->[1]; $coord++){
			$seq = $seq."N";
	}
	#adding nick site exceed gap length, need to trim the end
	if($coord > $interval->[1] + 1){
	    $seq = substr($seq, 0, length($seq) - ($coord-$interval->[1]-1));
	}
	return $seq;
}

#return begin and end index (right exclusive) in a cmap contig within 
#a specified coordinate range note that to avoid scanning 
#the whole cmap data all the time, one can "batch lookup" 
#the coordinates in successive calls of this function on
#many non-overlapping intervals in increasing coordinate order
#the method will remember the last coordinate left-off using
#a static/state variable

{
	my $prevInd = 0;
	my $prevContig;
	sub getCoordinatesSubset{
		my @interval = @{$_[0]};
		my $contig = $_[1];
		my $begin = 0;
		my $end = 0;
		#resetting the remembered index
		if(!(defined $prevInd) || !(defined $prevContig) || !($prevContig == $contig)){
			#print "resetting\n";
			$prevInd = 0;
		}
		#user can reset the index and decide to scan from a specified index
		if($#_ == 2){
			$prevInd = $_[2];
		}
		
		my @positions = @{$contig->{Position}};
		for(my $i = $prevInd; $i < $#positions; $i++){
			if($positions[$i] > $interval[0] && $begin <= 0){
				$begin = $i;
			}
			
			if($positions[$i] > $interval[1]){
				$end = $i;
				$prevInd = $end-1;
				last;
			}
			
			if($i == $#positions && $begin > 0){
				$end = $i+1;
				$prevInd = $end;
			}
		}
		$prevContig = $contig;
		if($begin > $end){
			$begin = 0;
			$end = 0;
		}
		return ($begin, $end);		
	}
}

sub reverseComplement
{
	my ($str) = @_;
	#print "forward: ".substr($str, length($str)-50, 50)."\n";
	my @match_forward = $str =~ /GCTCTTC/g;
	$str =~ tr/ATCGatcg/TAGCtagc/;
	$str = reverse($str);
	#my @match_reverse = $str =~ /GAAGAGC/g;
	#if($#match_forward != $#match_reverse){
	#	print "foward and reverse matches are not the same: $#match_forward\t$#match_reverse\n"
	#}
	#print "rcmp: ".substr(reverse($str), 0, 50)."\n";
	return $str;
}

#printing out the un-used ngs object as "singleton" to agp file as well as fasta
#file
#@TODO when only part of the ngs contig is not used only print out 
#part of the contigs
sub printUnUsedNGS
{
	my ($agp_out_file, $ngs_map, $aux_out_file, $fasta_out, $seqMap) = @_;
	open(my $fh, ">>".$agp_out_file) || die "Cannot open agp output file";
	open(my $aux_fh, ">>".$aux_out_file) || die "Cannot open agp auxilliary file";
	open(my $fasta_fh, ">".$fasta_out) || die "Cannot open fasta output file";
	
	my $unUsedCount=0;
	foreach my $key(keys %{$ngs_map}){
		my @ngs = @{$ngs_map->{$key}};
		my $ngsStart = $ngs[3];
		my $ngsEnd = $ngsStart + $ngs[1] - 1;
		if($ngs[2] <= 0){
			print $fh "$ngs[0]_obj\t1\t$ngs[1]\t1\tW\t$ngs[0]\t1\t$ngs[1]\t\+\n";
			print $aux_fh "$ngs[0]_obj\t0\t0\n";
			if(length($seqMap->($ngs[0])) > 0){				
				print $fasta_fh ">".$ngs[0]."_obj\n";
				print $fasta_fh $seqMap->($ngs[0])."\n";
			}else{
				warn("Sequence $ngs[0] has zero length, skipping it in output fasta\n");
			}
			$unUsedCount++;
		}
	}
	close($fh);
	close($aux_fh);
	close($fasta_fh);
}

sub cleanUp
{
	my ($dir) = @_;
	opendir(DIR, $dir);
	while(my $file = readdir(DIR)){
		if($file =~ /tmp\.fa\.tmp$/){
			print "Removing temp file: ". $file."\n";
			unlink $dir."/".$file;
		}
	}
	close(DIR);
}
