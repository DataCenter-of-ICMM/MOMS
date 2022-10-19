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

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script for assembling molecules in a BNX file to a CMAP (Consensus MAP) file
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
    -i, --input <str>    A BNX file for assembly (REQUIRED)
    -m, --model <str>    System model, i.e. "irys" or "saphyr" (default: saphyr)
    -g, --gsize <int>    The estimated genome size (M bp) (REQUIRED)
    -l, --len <int>      Minimum length of molecules (Kb) (default: 150)
    -s, --sites <int>    Minimum sites (labels) per molecule (default: 9)
    -p, --pvalue <float> The P-value threshold for assembly (default: 1e-10)
    -o, --output <str>   The output path and prefix (REQUIRED)
    -t, --threads <int>  The number of threads for parallel processing (default: 8)
    -c, --cluster        Using cluster for parallel computing (default: no)    
    --haplotype          From haplotype sample (default: no)
    --human              From human sample (default: no)
    -h, --help           Help

Options for invoking \$bionano/pipelineCL.py:
    -x, --xml <str>      Parameter file in XML format (default: \$bionano/xml/assembleArguments.xml)
    -e, --extent <int>   Number of extension and merge iterations N [0, 20] (default: 1)
    -b, --begin <int>    Start stage number (default: 0)
            1: sortBNX; 2: splitBNX; 3: Pairwise alignment; 4: assembly;
            5: refineA; 6: refinedB; 7: refineB1; 8: Merge0A;
            9+(I-1)*3: ext0(I); 10+(I-1)*3: ext1(I); 11+(I-1): mrg(I);
            9+N*3: refineFinal0; 10+N*3: refineFinal1; 11+N*3: alignmol

Example:
    $program -m saphyr --human --haplotype -i molecules_1.bnx -g 3000 -o assembly/enzyme1
USAGE

use vars qw($opt_i $opt_m $opt_hap $opt_hm $opt_g $opt_l $opt_s $opt_p $opt_o $opt_t $opt_c $opt_x $opt_e $opt_b $opt_h);
GetOptions( "i|input=s" => \$opt_i,
			"m|model=s" => \$opt_m,
			"haplotype" => \$opt_hap,
			"human" => \$opt_hm,
			"g|gsize=i" => \$opt_g,
			"l|len=i" => \$opt_l,
			"s|sites=i" => \$opt_s,
			"p|pvalue=f" => \$opt_p,
			"o|output=s" => \$opt_o,
			"t|threads=i" => \$opt_t,
			"c|cluster" => \$opt_c,
			"x|xml=s" => \$opt_x,
			"e|extent=i" => \$opt_e,
			"b|begin=i" => \$opt_b,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

BEGIN {
	select(STDERR);	$| = 1;
	select(STDOUT);	$| = 1;
}

# set the bionano tool path
my $progpath = abs_path(dirname($0));
my $bionano = abs_path(dirname($0) . "/../bionano");
my $toolpath = "$bionano/binary";
die("**ERROR: \"$toolpath\" directory is not found\n") unless(-d $toolpath);

my $revisor = "$progpath/reviseArguments.pl";
my $pipelineCL = "$bionano/pipelineCL.py";
my $assembler = "$toolpath/Assembler";
my $aligner  = "$toolpath/RefAligner";
die("**ERROR: can not find reviseArguments.pl at $progpath") unless(-f $revisor);
die("**ERROR: can not find pipelineCL.py at $bionano") unless(-f $pipelineCL);
die("**ERROR: Can not find Assembler at $toolpath\n") unless(-f $assembler);
die("**ERROR: Can not find RefAligner at $toolpath\n") unless(-f $aligner);

# check the required arguments
die("**ERROR: -i option must be specified\n") if(!defined $opt_i);
die("**ERROR: -g option must be specified\n") if(!defined $opt_g);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);

die("**ERROR: valid value of -m is either \"irys\" or \"saphyr\"\n") if((defined $opt_m) and !(grep{$_ eq lc($opt_m)} ("irys", "saphyr")));
my $bSaphyr = (!defined $opt_m || $opt_m ne "irys") ? 1 : 0;
my $bHaplotype = (defined $opt_hap) ? 1 : 0;
my $bHuman = (defined $opt_hm) ? 1 : 0;

# set the arguments
my $bnxfile = $opt_i;
die("**ERROR: \"$bnxfile\" is not a valid BNX file\n") if( ! -f $bnxfile);

my $genomesize = $opt_g;
my ($minlen, $minsites, $pvalue);
$minlen = (defined $opt_l) ? $opt_l : 150;
$minsites = (defined $opt_s) ? $opt_s : 9;
$pvalue = (defined $opt_p) ? $opt_p : 1e-10;
my $real_pvalue = $pvalue / $genomesize;
my $nthreads =  (defined $opt_t && ($opt_t =~ /^\d+$/)) ? (($opt_t < 1) ? 1 : $opt_t) : 8;
my $niter = (defined $opt_e && ($opt_e =~ /^\d+$/)) ? (($opt_e < 0) ? 0 : ($opt_e > 20) ? 20 : $opt_e) : 1;
my $bypass = (defined $opt_b && ($opt_b =~ /^\d+$/)) ? (($opt_b <= 1) ? 0 : ($opt_b >= 11 + $opt_e * 3) ? (10 + $opt_e * 3) : ($opt_b - 1)) : 0;

my $bZip = ($bnxfile =~ /\.gz$/);

my ($outdir, $outprefix);
if($opt_o =~ /\/$/){
	$opt_o =~ s/\/$//;
	$outdir = $opt_o;
	$outprefix = basename($bnxfile);
	$outprefix =~ s/\.gz$//;
	$outprefix =~ s/\.bnx$//;
}
else{
	$outdir = dirname($opt_o);
	$outprefix = basename($opt_o);
}
my $cmd = "mkdir -p $outdir";
my $retCode = system($cmd);
die("**ERROR: can not create directory $outdir") if($retCode != 0);
$outdir = abs_path($outdir);

my $optxml = (defined $opt_x) ? $opt_x : "$bionano/xml/assembleArguments.xml";
my $params="";
$params .= " -m " . ($bSaphyr ? "saphyr" : "irys");
$params .= " --haplotype" if($bHaplotype);
$params .= " --human" if($bHuman);
$cmd="$revisor $params -l $minlen -s $minsites -p $real_pvalue -i $optxml > $outdir/assembleArguments.xml";
$retCode = system($cmd);

my $clusterXml = "$bionano/xml/clusterArguments.xml";
$params = (defined $opt_c && -f $clusterXml) ? " -C $clusterXml" : "";

if($bZip){
	$cmd="zcat $bnxfile | python2 $pipelineCL -T $nthreads -N 6 -i $niter -t $toolpath -a $outdir/assembleArguments.xml$params -b /dev/stdin -V 0 -m -z -l $outdir -e $outprefix";
} else{
	$cmd="python2 $pipelineCL -T $nthreads -N 6 -i $niter -t $toolpath -a $outdir/assembleArguments.xml$params -b $bnxfile -V 0 -m -z -l $outdir -e $outprefix";
}
if($bypass == 0){
	$cmd .= " -w";
}
else{
	$cmd .= " -B $bypass";
}
print "$cmd\n";
system($cmd);

exit(0);
