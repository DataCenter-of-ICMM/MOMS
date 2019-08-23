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
$program: A perl script for aligning a pair of CMAP files
Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences

Usage: $program [options]
	-r, --reference <str> The reference CMAP file for alignment (REQUIRED)
	-q, --query <str>     The query CMAP file for alignment (REQUIRED)
	-o, --output <str>    The output path (REQUIRED)
	-t, --threads <int>   The number of threads for parallel processing (default: 8)
	-h, --help            Help

Options for invoking \$bionano/binary/RefAligner:
	-s, --stage           The stage of the alignment (first|BNG|final) (default: first)
	-c, --hashcolor <int> Set the hashcolor for the final stage (default: 1)
	-x, --xml <str>       Parameter file in XML format (default: \$bionano/xml/alignArguments.xml)

Example:
	$program -r ngs_assembly.cmap -q bng_assembly.cmap -o align
USAGE

use vars qw($opt_r $opt_q $opt_s $opt_c $opt_o $opt_t $opt_h $opt_x);
GetOptions( "r|reference=s" => \$opt_r,
			"q|query=s" => \$opt_q,
			"s|stage=s" => \$opt_s,
			"c|hashcolor=i" => \$opt_c,
			"o|output=s" => \$opt_o,
			"t|threads=i" => \$opt_t,
			"x|xml=s" => \$opt_x,
			"h|help" => \$opt_h);

die($usage) if($opt_h);

# set the bionano tool path
my $bionano = abs_path(dirname($0) . "/../bionano");
my $toolpath = abs_path("$bionano/binary");
die("**ERROR: \"$toolpath\" directory is not found\n") unless(-d $toolpath);

# set the BNG refaligner
my $refaligner  = "$toolpath/RefAligner";
die("**ERROR: Can not find RefAligner at $toolpath\n") unless(-f $refaligner);

# check the required arguments
die("**ERROR: -r option must be specified\n") if(!defined $opt_r);
die("**ERROR: -q option must be specified\n") if(!defined $opt_q);
die("**ERROR: -o option must be specified\n") if(!defined $opt_o);

# set the arguments
my $refcmap = $opt_r;
my $qrycmap = $opt_q;
my $optxml = (defined $opt_x) ? $opt_x : "$bionano/xml/alignArguments.xml";
my $stage = (defined $opt_s) ? $opt_s : "first";
my $color = (defined $opt_c) ? $opt_c : 1;

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
die("**ERROR: reference cmap \"$refcmap\" does not exist") unless(-f $refcmap);
die("**ERROR: query cmap \"$qrycmap\" does not exist") unless(-f $qrycmap);
die("**ERROR: parameter file \"$optxml\" does not exist") unless(-f $optxml);
# translate to absolute paths for later use
#$refcmap = abs_path($refcmap);
#$qrycmap = abs_path($qrycmap);
#$optxml = abs_path($optxml);

# prepare output directory
my $cmd = "mkdir -p $outdir";
my $retCode = system($cmd);
die("**ERROR: can not create directory $outdir") if($retCode != 0);
#$outdir = abs_path($outdir);

# process input files
my $XML = new XML::Simple(KeyAttr=>[]);
my $configs = $XML->XMLin($optxml);

my $paraDict = parseConfig($configs, "global");
my $maxmem = $paraDict->{'maxmem'}{val};
my $maxthreads = $paraDict->{'maxthreads'}{val};
my $maxvirtmem = (defined $paraDict->{'maxvirtmem'}{val}) ? $paraDict->{'maxvirtmem'}{val} : 0;
my $nthreads =  (defined $opt_t && ($opt_t =~ /^\d+$/)) ? (($opt_t < 1) ? 1 : $opt_t) : $maxthreads;

my $sharedParams = "-stdout -stderr -maxmem $maxmem -maxthreads $maxthreads -maxvirtmem $maxvirtmem";
if(defined $paraDict->{'RAmem'}{val}){
	$sharedParams .= " -RAmem " . $paraDict->{'RAmem'}{val};
}
$sharedParams .= " -hashcolor $color" if($stage eq "final");
my $paraStage = parseConfig($configs, "stage_${stage}");
my $params = makeParams($paraStage);

$cmd = "$refaligner -ref $refcmap -i $qrycmap -o $outdir/$outpre $sharedParams $params;";
print "$cmd\n";
system("$cmd");

#my $nhits = `grep -c '^[^#]' $outdir/$outpre.xmap 2>/dev/null | tr -d '\n'`;
#$nhits = 0 unless($nhits ne "");
#printf "%5s alignments found.\n", $nhits;

exit(0);
