#!/bin/bash
# Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences
# MOMS is licensed under the Mulan PSL v1.
# You can use this software according to the terms and conditions of the Mulan PSL v1.
# You may obtain a copy of Mulan PSL v1 at:
#    http://license.coscl.org.cn/MulanPSL
# THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND, EITHER EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT, MERCHANTABILITY OR FIT FOR A PARTICULAR
# PURPOSE.
# See the Mulan PSL v1 for more details.

if [ $# -lt 5 ]; then
	echo "Usage: $0 refdir qrypath suffix outdir fillflag [nthreads [xmldir]]";
	exit 1;
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$DIR/util.sh"
aligner="$DIR/perl/cmap_aligner.pl"
cmapfilter="$DIR/perl/cmap_filter.pl"
gapfiller="$DIR/perl/cmap_gap_filler.pl"
fullsetter="$DIR/perl/cmap_fullset.pl"
unifier="$DIR/perl/cmap_uni.pl"

refdir=$1
qrypath=$2
suffix=$3
outdir=$4
bfill=$(echo $5 | egrep '^[0-9]+$');
stfile="$outdir/status.txt"

nthreads=0
alignconf=""
if [ $# -ge 6 ]; then
	nthreads=$(echo ${6} | egrep '^[0-9]+$');
	if [ $# -ge 7 ]; then
		xmldir=$7
		if [[ -d $xmldir ]]; then
			alignconf=" -x $xmldir/alignArguments.xml";
		fi
	fi
fi

if [[ "$qrypath" != *\/*[^\/] ]]; then
	echo "qrypath must be in \"dir/prefix\" format";
	echo 'Error 1' > $stfile
	exit 1;
fi
qrydir=${qrypath%/*};
qrypre=${qrypath##*/};

shareparams="";
if [ $nthreads -gt 0 ]; then
	shareparams+=" -t $nthreads";
fi

alndir="$outdir/single"
echo -e "Beginning alignment of BNG contigs to the final scaffold ... ";
stime=$(ntime)
i=0;
for refmap in `ls -1 $refdir/mono/*.cmap`; do
	name=${refmap##*/};
	name=${name%.cmap};
	names+=($name);
	mkdir -p  "$alndir/$name";
	if [ ! -f $alndir/$name/$name.xmap ]; then
		qrymap="$qrydir/${qrypre}_${name}_$suffix.cmap";
		check "$aligner -s BNG -r $refmap -q $qrymap -o $alndir/$name/$name $shareparams$alignconf > $alndir/$name/$name-align.log";
	fi
	if [ -f $alndir/$name/$name.xmap ]; then
		nhits=`grep -c '^[^#]' $alndir/$name/$name.xmap 2>/dev/null | tr -d '\n'`;
		if [[ -n "$nhits" ]]; then
			refhits=$(count_xmaps $alndir/$name/$name.xmap 3);
			qryhits=$(count_xmaps $alndir/$name/$name.xmap 2);
			printf "%5s alignments found for enzyme $name between %6s BNG contigs and %6s super-contigs.\n" $nhits $qryhits $refhits
		fi
		eval "ln -sf $name/$name.xmap $alndir";
		((i++));
	fi
done
etime=$(ntime)
if [ $i -lt ${#names[@]} ]; then
	echo 'Error 2' > $stfile
	exit 1;
fi
echo -e "Alignment completed in $(duration $stime $etime) s.\n";

if [ $bfill -ne 0 ]; then
	filldir="$outdir/gapfilled"

	mkdir -p $filldir
	echo  -e "Beginning to fill gaps using non-aligned BNG contigs ... ";
	stime=$(ntime)
	i=0;
	for name in ${names[@]}; do
		if [ ! -f $filldir/$name.filled.cmap ]; then
			check "grep -v \"^#\" $alndir/$name.xmap | cut -f2 | sort | uniq | $cmapfilter -c $qrydir/${qrypre}_${name}_$suffix.cmap -i - -exclude -o $filldir/$name.non_used > $filldir/$name.nonused.log";
			check "$gapfiller -i $refdir/mono/$name.cmap -a $filldir/$name.non_used.cmap --addonly -o $filldir/$name.filled $shareparams > $filldir/$name.filled.log";
		fi
		if [ -f $filldir/$name.filled.cmap ]; then
			((i++));
		fi
	done
	etime=$(ntime)
	if [ $i -lt ${#names[@]} ]; then
		echo 'Error 3' > $stfile
		exit 1;
	fi
	echo -e "Gap-filling completed in $(duration $stime $etime) s.\n";

	eval "rm -rf $outdir/mono";
	mkdir -p "$outdir/mono";
	for name in ${names[@]}; do
		eval "ln -sf ../gapfilled/$name.filled.cmap $outdir/mono/$name.cmap";
	done

	cmd="$fullsetter"
	for name in ${names[@]}; do
		cmd="$cmd -i $outdir/mono/$name.cmap";
	done
	cmd="$cmd -o $outdir/multicolors";
	check "$cmd";
else
	eval "rm -rf $outdir/mono; cp -rp $refdir/mono/ $outdir";
	check "cp -p $refdir/multicolors.cmap $outdir"
fi

check "$unifier -i $outdir/multicolors.cmap -o $outdir/unified";

echo 'Done' > $stfile
exit 0;
