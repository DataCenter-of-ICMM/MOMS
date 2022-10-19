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

# default values
my ($def_minlen, $def_minsites, $def_pvalue) = ("150", "9", "1e-5");
my ($def_fp, $def_fn, $def_sd, $def_sf) = ("1.0", "0.10", "0.0", "0.15");

my $program = basename($0);
my $usage = << "USAGE";
$program: A perl script used to revise arguments for running Bionano pipelineCL.py
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences
Usage: $program [options]
    -i, --input <str>     A XML file for revision (REQUIRED)
    -m, --model <str>     System model, i.e. "irys" or "saphyr" (default: saphyr)
    -l, --len <int>       Minimum length of molecules (Kb) (default: $def_minlen)
    -s, --sites <int>     Minimum sites (labels) per molecule (default: $def_minsites)
    -p, --pvalue <float>  The P-value threshold for assembly (default: $def_pvalue)
    --fp <int>            False Positives per 100kb (default: $def_fp)
    --fn <float>          False Negative Probability (default: $def_fn)
    --sd <float>          scaling error in root-kb (default: $def_sd)
    --sf <float>          fixed error in kb (default: $def_sf)
    --haplotype           Enable haplotype calling (default: no)
    --human               For human sapiens (default: no)
    -h, --help            Help
Example:
    $program -m irys -l 120 -s 8 -p 1e-6 -i oldArguments.xml > newArguments.xml
USAGE

use vars qw($opt_i $opt_m $opt_hap $opt_hm $opt_l $opt_s $opt_p $opt_fp $opt_fn $opt_sd $opt_sf $opt_h);
GetOptions( "i|input=s" => \$opt_i,
			"m|model=s" => \$opt_m,
			"haplotype" => \$opt_hap,
			"human" => \$opt_hm,
			"l|len=i" => \$opt_l,
			"s|sites=i" => \$opt_s,
			"p|pvalue=f" => \$opt_p,
			"fp=f" => \$opt_fp,
			"fn=f" => \$opt_fn,
			"sd=f" => \$opt_sd,
			"sf=f" => \$opt_sf,
			"h|help");

die($usage) if(defined $opt_h);

BEGIN {
	select(STDERR);	$| = 1;
	select(STDOUT);	$| = 1;
}

die("**ERROR: -i option must be specified\n") if(!defined $opt_i);
die("**ERROR: valid value of -m is either \"irys\" or \"saphyr\"\n") if((defined $opt_m) and !(grep{$_ eq lc($opt_m)} ("irys", "saphyr")));
my $bSaphyr = (!defined $opt_m || $opt_m ne "irys");
my $bHaplotype = (defined $opt_hap);
my $bHuman = (defined $opt_hm);

my $xmlfile = $opt_i;

my $minlen = (defined $opt_l) ? $opt_l : $def_minlen;
my $minsites = (defined $opt_s) ? $opt_s : $def_minsites;
my $pvalue = (defined $opt_p) ? $opt_p : $def_pvalue;

my $fp = (defined $opt_fp) ? $opt_fp : $def_fp;
my $fn = (defined $opt_fn) ? $opt_fn : $def_fn;
my $sd = (defined $opt_sd) ? $opt_sd : $def_sd;
my $sf = (defined $opt_sf) ? $opt_sf : $def_sf;

my ($res, $mergeres, $rres, $arres);
if($bSaphyr){
	$res = 3.1;
	$mergeres = 2.6;
	$rres = 0.9;
	$arres = 0.9;
}
else{
	$res  = 3.4;
	$mergeres = 2.9;
	$rres = 1.2;
	$arres = 1.05;
}

my ($bestrefwt, $addfragswt);
if($bHaplotype){
	$bestrefwt = 1;
	$addfragswt = 1;
}
else{
	$bestrefwt = 0;
	$addfragswt = 0;
}

my ($maxsites, @minSNRs);
if($bHuman){
	$maxsites = 500;
	@minSNRs = (2.0, 15.0, 0.5, 100, 0.8);
}
else{
	$maxsites = 200;
	@minSNRs = (2.0, 2.0, 0.5, 100, 0.8);
}

my $line;
my @sections = ();
my $sec;
my $orphan;
open(IN, "<$xmlfile") or die("Can not open $xmlfile for reading");
while ($line = <IN>){
	chomp($line);
	if($line =~ /^\s*<([A-Za-z_0-9]+)>/){
		push(@sections, $1);
		$sec = $1;
		print "$line\n";
		next;
	}
	if($line =~ /^\s*<\/([A-Za-z_0-9]+)>/){
		if( (scalar(@sections) <= 0) or ($1 ne $sections[-1]) ){
			$orphan = $1;
			last;
		}
		pop(@sections);
		if(scalar(@sections) > 0){
			$sec = $sections[-1];
		}
		else{
			undef $sec;
		}
		print "$line\n";
		next;
	}
	if(!defined $sec){
		print "$line\n";
		next;
	}
	if($sec eq "bnx_sort"){
		if($line =~ /<flag attr=\"-minlen\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-minlen\"\s+val0=\")[^\"]*(\")/$1$minlen$2/;
		}
		elsif($line =~ /<flag attr=\"-minsites\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-minsites\"\s+val0=\")[^\"]*(\")/$1$minsites$2/;
		}
		elsif($line =~ /<flag attr=\"-maxsites\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-maxsites\"\s+val0=\")[^\"]*(\")/$1$maxsites$2/;
		}
		elsif($line =~ /<flag attr=\"-MaxIntensity\"\s+/){
			if($bSaphyr){
				$line = "";
			}
		}
	}
	elsif($sec eq "noise0"){
		if($line =~ /<flag attr=\"-FP\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-FP\"\s+val0=\")[^\"]*(\")/$1$fp$2/;
		}
		elsif($line =~ /<flag attr=\"-FN\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-FN\"\s+val0=\")[^\"]*(\")/$1$fn$2/;
		}
		elsif($line =~ /<flag attr=\"-sd\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-sd\"\s+val0=\")[^\"]*(\")/$1$sd$2/;
		}
		elsif($line =~ /<flag attr=\"-sf\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-sf\"\s+val0=\")[^\"]*(\")/$1$sf$2/;
		}
		elsif($line =~ /<flag attr=\"-res\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-res\"\s+val0=\")[^\"]*(\")/$1$res$2/;
		}
	}
	elsif($sec eq "initialAssembly"){
		if($line =~ /<flag attr=\"-T\"\s+val0=\".*\"/){
			$line  =~ s/(<flag attr=\"-T\"\s+val0=\")[^\"]*(\")/$1$pvalue$2/;
		}
	}
	elsif($sec eq "extendRefine"){
		if($line =~ /<flag attr=\"-T\"\s+val0=\".*\"/){
			my $val = $pvalue / 10;
			$line  =~ s/(<flag attr=\"-T\"\s+val0=\")[^\"]*(\")/$1$val$2/;
		}
	}
	elsif($sec eq "autoNoise"){
		if($line =~ /<flag attr=\"-rres\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-rres\"\s+val0=\")[^\"]*(\")/$1$rres$2/;
		}
	}
	elsif($sec eq "autoNoise0"){
		if($line =~ /<flag attr=\"-minSNRestimate\"\s+/){
			if($bSaphyr){
				my $str = "";
				for(my $i=0; $i<scalar(@minSNRs); $i++){
					$str .= " val" . ($i+1) . "=\"$minSNRs[$i]\"";
				}
				$line =~ s/(<flag attr=\"-minSNRestimate\")\s+.*(\/>)/$1$str$2/;
			}
			else{
				$line = "";
			}
		}
		elsif($line =~ /<flag attr=\"-bpp\"\s+/){
			if(!$bSaphyr){
				$line = "";
			}
		}
	}
	elsif($sec eq "autoNoise1"){
		if($line =~ /<flag attr=\"-minSNRestimate\"\s+/){
			if(!$bSaphyr or !$bHuman){
				$line = "";
			}
		}
	}
	elsif($sec eq "refineBCommon"){
		if($line =~ /<flag attr=\"-rres\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-rres\"\s+val0=\")[^\"]*(\")/$1$rres$2/;
		}
	}
	elsif($sec eq "refineFinalCommon"){
		if($line =~ /<flag attr=\"-rres\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-rres\"\s+val0=\")[^\"]*(\")/$1$rres$2/;
		}
		elsif($line =~ /<flag attr=\"-BestRefWT\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-BestRefWT\"\s+val0=\")[^\"]*(\")/$1$bestrefwt$2/;
		}
		elsif($line =~ /<flag attr=\"-AddFragsWT\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-AddFragsWT\"\s+val0=\")[^\"]*(\")/$1$addfragswt$2/;
		}
		elsif($line =~ /<flag attr=\"-Haplotype\"\s+/){
			if(!$bHaplotype){
				$line = "";
			}
		}
	}
	elsif($sec eq "extensionCommon"){
		if($line =~ /<flag attr=\"-rres\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-rres\"\s+val0=\")[^\"]*(\")/$1$rres$2/;
		}
	}
	elsif($sec eq "merge"){
		if($line =~ /<flag attr=\"-T\"\s+val0=\".*\"/){
			my $val = $pvalue / 100;
			$line  =~ s/(<flag attr=\"-T\"\s+val0=\")[^\"]*(\")/$1$val$2/;
		}
		elsif($line =~ /<flag attr=\"-res\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-res\"\s+val0=\")[^\"]*(\")/$1$mergeres$2/;
		}
	}
	elsif($sec eq "lastMerge"){
		if($line =~ /<flag attr=\"-pairmergeHmap\"\s+/){
			if(!$bHaplotype){
				$line = "";
			}
		}
	}
	elsif($sec eq "characterizeDefault"){
		if($line =~ /<flag attr=\"-res\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-res\"\s+val0=\")[^\"]*(\")/$1$mergeres$2/;
		}
	}
	elsif($sec eq "alignmolCommon"){
		if($line =~ /<flag attr=\"-res\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-res\"\s+val0=\")[^\"]*(\")/$1$res$2/;
		}
		elsif($line =~ /<flag attr=\"-rres\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-rres\"\s+val0=\")[^\"]*(\")/$1$rres$2/;
		}
	}
	elsif($sec eq "alignmolvref"){
		if($line =~ /<flag attr=\"-rres\"\s+val0=\".*\"/){
			$line =~ s/(<flag attr=\"-rres\"\s+val0=\")[^\"]*(\")/$1$arres$2/;
		}
	}

	print "$line\n";
}
close IN;

die("Unpaired tag found: </$orphan>") if(defined $orphan);
