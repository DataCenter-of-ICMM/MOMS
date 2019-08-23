#!/bin/bash
# Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences
# MOMS is licensed under the Mulan PSL v1.
# You can use this software according to the terms and conditions of the Mulan PSL v1.
# You may obtain a copy of Mulan PSL v1 at:
#    http://license.coscl.org.cn/MulanPSL
# THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND, EITHER EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT, MERCHANTABILITY OR FIT FOR A PARTICULAR
# PURPOSE.
# See the Mulan PSL v1 for more details.

if [ $# -lt 6 ]; then
	echo "Usage: $0 ngspath bngpath inpath enpath outdir fillflag [nthreads [xmldir]]";
	exit 1;
fi

# for DEBUG
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$DIR/util.sh"
merger="$DIR/perl/cmap_merger.pl"
cmapfilter="$DIR/perl/cmap_filter.pl"
gapfiller="$DIR/perl/cmap_gap_filler.pl"
ngsaligner="$DIR/perl/cmap_aligner_two_passes.pl"
bngaligner="$DIR/perl/cmap_aligner.pl"
joiner="$DIR/perl/cmap_joiner.pl"
exporter="$DIR/perl/ExportAGP.pl"

ngspath=$1
bngpath=$2
inpath=$3
enpath=$4
outdir=$5
bfill=$(echo $6 | egrep '^[0-9]+$');
stfile="$outdir/status.txt"

nthreads=0
alignconf=""
mergeconf=""
if [ $# -ge 7 ]; then
	nthreads=$(echo $7 | egrep '^[0-9]+$');
	if [ $# -ge 8 ]; then
		xmldir=$8
		if [[ -d $xmldir ]]; then
			alignconf=" -x $xmldir/alignArguments.xml"
			mergeconf=" -x $xmldir/mergeArguments.xml"
		fi
	fi
fi

if [[ "$ngspath" != *\/*[^\/] ]]; then
	echo "ngspath must be in \"dir/prefix\" format";
	echo 'Error 1' > $stfile
	exit 1;
fi
ngsdir=${ngspath%/*}
ngspre=${ngspath##*/};

if [[ "$bngpath" != *\/*[^\/] ]]; then
	echo "bngpath must be in \"dir/prefix\" format";
	echo 'Error 2' > $stfile
	exit 1;
fi
bngdir=${bngpath%/*};
bngpre=${bngpath##*/};
alndir="$bngdir/align1";

ngssuf="cut";
bngsuf="adjusted_cut";

if [[ "$inpath" != *\/*[^\/] ]]; then
	echo "inpath must be in \"dir/prefix\" format";
	echo 'Error 3' > $stfile
	exit 1;
fi
indir=${inpath%/*}
infile=${inpath##*/};

if [[ "$enpath" != *\/*[^\/] ]]; then
	echo "enpath must be in \"dir/prefix\" format";
	echo 'Error 4' > $stfile
	exit 1;
fi
encdir=${enpath%/*}
encpre=${enpath##*/};

shareparams="";
if [ $nthreads -gt 0 ]; then
	shareparams+=" -t $nthreads";
fi

declare -A counts
for fst in `ls -1 ${ngsdir}/${ngspre}_*_$ngssuf.cmap`; do
	rname=${fst##*/${ngspre}_};
	rname=${rname%_$ngssuf.cmap};
	((counts[$rname]++));
done
for scd in `ls -1 ${bngdir}/${bngpre}_*_$bngsuf.cmap`; do
	qname=${scd##*/${bngpre}_};
	qname=${qname%_$bngsuf.cmap};
	((counts[$qname]++));
done

# force the enzymes in an order, which is very important for downstream analysis
eval "declare -a ordered=`orderedEnzymes`"
declare -A valid
for name in ${ordered[@]}; do
	((valid[$name]++));
done

for name in ${!counts[@]}; do
	if [ -z "${valid[$name]}" ]; then
		echo "$name is not a valid motif";
		echo 'Error 5' > $stfile
		exit 0;
	fi
done

# find the common enzyme names
for name in ${ordered[@]}; do
	if [[ -n "${counts[$name]}" ]]; then
		if [[ ${counts[$name]} -eq 2 ]]; then
			names+=($name);
		fi
	fi
done

mkdir -p $outdir
cmapdir="$outdir/CMAPs"
fastadir="$outdir/FASTAs"
ancdir="$outdir/anchor"

echo -e "Beginning to merge NGS cmaps and BN cmaps to construct single-enzyme scaffolds ...";
stime=$(ntime)
mkdir -p $cmapdir
for name in ${names[@]}; do
	errbin="$alndir/$name.errbin";
	ngscmap="$ngsdir/${ngspre}_${name}_$ngssuf.cmap";
	bngcmap="$bngdir/${bngpre}_${name}_$bngsuf.cmap";
	if [[ ! -f "$cmapdir/$name.hybrid.cmap" ]]; then
		mkdir -p $cmapdir/$name
		if [[ ! -f "$cmapdir/$name/step2.hybrid.cmap" ]]; then
			check "$merger -f $ngscmap -s $bngcmap -e $errbin -o $cmapdir/$name$mergeconf $shareparams > $cmapdir/$name/cmap_merge.stdout";
		fi
		if [[ -f "$cmapdir/$name/step2.hybrid.cmap" ]]; then
# get the used IDs
			ids="$cmapdir/$name/used_ids.txt"
			eval "(awk -F'\" \"' '{if(NR>1){printf \$2\"\\n\"; printf \$3\"\\n\"}}' $cmapdir/$name/step1.merge.pairs.txt | sort -nu > $ids)";

# categorize the NGS cmaps and BNG cmaps accordingly
			log="$cmapdir/$name/cmap_filter.stdout";
			:> $log
			cmd="$cmapfilter -c $ngscmap -i $ids -o $cmapdir/$name/${ngspre}_${name}.used";
			echo $cmd >> $log; eval "($cmd >> $log)"
			cmd="$cmapfilter -c $ngscmap -i $ids -exclude -o $cmapdir/$name/${ngspre}_${name}.non_used";
			echo $cmd >> $log; eval "($cmd >> $log)"
			bngcmap="$cmapdir/$name/${bngpre}_${name}_${bngsuf}_idshift.cmap";
			cmd="$cmapfilter -c $bngcmap -i $ids -o $cmapdir/$name/${bngpre}_${name}.used";
			echo $cmd >> $log; eval "($cmd >> $log)"
			cmd="$cmapfilter -c $bngcmap -i $ids -exclude -o $cmapdir/$name/${bngpre}_${name}.non_used";
			echo $cmd >> $log; eval "($cmd >> $log)"

# remove temporary files
			eval "rm -f $cmapdir/$name/Mrg_*.cmap";
# fill the gaps within cmaps
			if [ $bfill -eq 0 ]; then
				eval "(ln -sf $name/step2.hybrid.cmap $cmapdir/$name.hybrid.cmap)";
			else
				eval "$gapfiller -i $cmapdir/$name/step2.hybrid.cmap -a $cmapdir/$name/${bngpre}_${name}.non_used.cmap -o $cmapdir/$name.hybrid -t $nthreads > $cmapdir/$name.filled.log";
			fi
		fi
	fi
	ncmaps=$(count_cmaps $cmapdir/$name.hybrid.cmap);
	printf "%8s super-contigs found for enzyme $name.\n" $ncmaps;
done
etime=$(ntime)
echo -e "Single-enzyme scaffolds construction completed in $(duration $stime $etime) s.\n";

echo -e "Beginning alignment of NGS cmap to Hybrid CMAP ...";
stime=$(ntime)
mkdir -p $ancdir
for name in ${names[@]}; do
	scaffold="$cmapdir/$name.hybrid.cmap";
	ngscmap="$ngsdir/${ngspre}_${name}_$ngssuf.cmap";
	if [[ ! -f "$ancdir/$name/$name-NGS.xmap" ]]; then
		mkdir -p $ancdir/$name
		check "$ngsaligner -r $scaffold -q $ngscmap -o $ancdir/$name/$name-NGS $shareparams$alignconf > $ancdir/$name/$name-NGS-align.log";
	fi
	if [[ -f "$ancdir/$name/$name-NGS.xmap" ]]; then
		nhits=`grep -c '^[^#]' $ancdir/$name/$name-NGS.xmap 2>/dev/null | tr -d '\n'`;
		if [[ -n "$nhits" ]]; then
			refhits=$(count_xmaps $ancdir/$name/$name-NGS.xmap 3);
			qryhits=$(count_xmaps $ancdir/$name/$name-NGS.xmap 2);
			printf "%8s alignments found for enzyme $name between %6s NGS contigs and %6s super-contigs.\n" $nhits $qryhits $refhits
		fi
		if [[ ! -f "$ancdir/$name.xmap" ]]; then
			eval "ln -sf $name/$name-NGS.xmap $ancdir/$name.xmap"
		fi
	fi
done
etime=$(ntime)
echo -e "align sequences to hybrid scaffolds completed in $(duration $stime $etime) s.\n";

echo -e "Beginning alignment of BNG cmap to Hybrid CMAP ...";
stime=$(ntime)
mkdir -p $ancdir
for name in ${names[@]}; do
	scaffold="$cmapdir/$name.hybrid.cmap";
	bngcmap="$bngdir/${bngpre}_${name}_$bngsuf.cmap";
	if [[ ! -f "$ancdir/$name/$name-BNG.xmap" ]]; then
		mkdir -p $ancdir/$name
		check "$bngaligner -r $scaffold -q $bngcmap -o $ancdir/$name/$name-BNG $shareparams$alignconf > $ancdir/$name/$name-BNG-align.log";
	fi
	if [[ -f "$ancdir/$name/$name-BNG.xmap" ]]; then
		nhits=`grep -c '^[^#]' $ancdir/$name/$name-BNG.xmap 2>/dev/null | tr -d '\n'`;
		if [[ -n "$nhits" ]]; then
			refhits=$(count_xmaps $ancdir/$name/$name-BNG.xmap 3);
			qryhits=$(count_xmaps $ancdir/$name/$name-BNG.xmap 2);
			printf "%8s alignments found for enzyme $name between %6s BNG contigs and %6s super-contigs.\n" $nhits $qryhits $refhits
		fi
	fi
done
etime=$(ntime)
echo -e "align Bionano genome maps to hybrid scaffolds completed in $(duration $stime $etime) s.\n";

echo -e "Merging Hybrid CMAP with NGS not participated in the hybrid scaffold ...";
stime=$(ntime)
for name in ${names[@]}; do
	scaffold="$cmapdir/$name.hybrid.cmap";
	ngscmap="$cmapdir/$name/${ngspre}_${name}.non_used.cmap"
	check "$joiner -f $scaffold -s $ngscmap -o $cmapdir/$name/step2.hybrid-NGS-non_used"
done
etime=$(ntime)
echo -e "Merging Hybrid CMAP with naive NGS CMAP complete in $(duration $stime $etime) s.\n";

echo -e "Merging Hybrid CMAP with BNG CMAP not participated in the hybrid scaffold ...";
stime=$(ntime)
for name in ${names[@]}; do
	scaffold="$cmapdir/$name.hybrid.cmap";
	bngcmap="$cmapdir/$name/${bngpre}_${name}.non_used.cmap"
	check "$joiner -f $scaffold -s $bngcmap -o $cmapdir/$name/step2.hybrid-BNG-non_used"
done
etime=$(ntime)
echo -e "Merging Hybrid CMAP with naive NGS CMAP complete in $(duration $stime $etime) s.\n";

echo -e "Beginning construction of AGP and FASTA file of the scaffolded and unscaffolded sequences ...";
stime=$(ntime)
for name in ${names[@]}; do
	xmap="$ancdir/$name/$name-NGS.xmap";
	cmap="$ancdir/$name/$name-NGS_r.cmap";
	fasta="$indir/$infile";
	map="$encdir/${encpre}_${name}_key.txt";
	coord="$ngsdir/${name}_auto_cut_NGS_coord_translation.txt";
	if [[ ! -f "$fastadir/$name/$name.agp" || ! -f "$fastadir/$name/$name.fasta" ]]; then
		mkdir -p $fastadir/$name
		check "$exporter -i $xmap -c $cmap -s $fasta -m $map -t $coord -o $fastadir/$name/$name > $fastadir/$name/export.log";
	fi
	if [ -f $fastadir/$name/$name.agp ]; then
		eval "ln -sf $name/$name.agp $fastadir/$name.agp"
	fi
	if [ -f $fastadir/$name/$name.fasta ]; then
		eval "ln -sf $name/$name.fasta $fastadir/$name.fasta"
	fi
done
etime=$(ntime)
echo -e "AGP and FASTA generation complete in $(duration $stime $etime) s.\n";

for name in ${names[@]}; do
	if [ ! -f "$cmapdir/$name.hybrid.cmap" ]; then
		echo 'Error 6' > $stfile
		exit 1;
	fi
	eval "ln -sf CMAPs/$name.hybrid.cmap $outdir";
	if [ ! -f $fastadir/$name.agp ]; then
		echo 'Error 7' > $stfile
		exit 1;
	fi
	eval "ln -sf FASTAs/$name.agp $outdir";
	if [ ! -f $fastadir/$name.fasta ]; then
		echo 'Error 8' > $stfile
		exit 1;
	fi
	eval "ln -sf FASTAs/$name.fasta $outdir";
done

echo 'Done' > $stfile
exit 0;
