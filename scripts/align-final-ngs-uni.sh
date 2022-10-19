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

if [ $# -lt 7 ]; then
	echo "Usage: $0 refdir qrypath suffix infasta enzymes endir outdir [nthreads [xmldir]]";
	exit 1;
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$DIR/util.sh"
aligner="$DIR/perl/cmap_aligner_two_passes.pl"
collector="$DIR/perl/collectUnusedContigs.pl"
fa2cmap="$DIR/perl/fa2cmap_multi_color.pl"
unifier="$DIR/perl/cmap_uni.pl"
xmapfilter="$DIR/perl/xmap_filter_by_cosine.pl"
idmap="$DIR/perl/xmap_idmap.pl"

refdir=$1
qrypath=$2
suffix=$3
infasta=$4
IFS=',' read -r -a array <<< "$5"
enzymenames=`for i in ${!array[@]}; do str=$(echo ${array[$i]} | cut -d: -f2); echo -n "_$str"; done; echo`
enzymenames=${enzymenames#_}
enzymeseqs=`for i in ${!array[@]}; do str=$(echo ${array[$i]} | cut -d: -f2,3); echo -n "$str $(( i + 1 )) "; done; echo`
endir=$6
outdir=$7
stfile="$outdir/status.txt"
nthreads=0
alignconf=""
if [ $# -ge 8 ]; then
	nthreads=$(echo ${8} | egrep '^[0-9]+$');
	if [ $# -ge 9 ]; then
		xmldir=$9
		if [[ -d $xmldir ]]; then
			alignconf=" -x $xmldir/alignArguments.xml"
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

# force the enzymes in an order, which is very important for downstream analysis
eval "declare -a ordered=`orderedEnzymes`"
declare -A valid
for name in ${ordered[@]}; do
	((valid[$name]++));
done
declare -A used
for cmap in `ls -1 $refdir/mono/*.cmap`; do
	name=${cmap##*/};
	name=${name%.cmap};
	if [ -z "${valid[$name]}" ]; then
		echo "$name is not a valid motif";
		echo "Error 1" > $stfile;
		exit 1;
	fi
	used[$name]=1;
done
for name in ${ordered[@]}; do
	if [[ -n "${used[$name]}" ]]; then
		names+=($name);
	fi
done

sgldir="$outdir/single"

echo -e "Beginning alignment of NGS contigs to the final scaffold using single channel ... ";
stime=$(ntime)
mkdir -p $sgldir
j=0;
for ((i=0; i<${#names[@]}; i++)); do
	name=${names[$i]};
	mkdir -p  "$sgldir/$name";
	if [ ! -f $sgldir/$name/$name.xmap ]; then
		refmap="$refdir/mono/$name.cmap";
		qrymap="$qrydir/${qrypre}_${name}_$suffix.cmap";
		check "$aligner -r $refmap -q $qrymap -o $sgldir/$name/$name $shareparams$alignconf > $sgldir/$name/$name-align.log";
	fi
	if [ -f $sgldir/$name/$name.xmap ]; then
		nhits=`grep -c '^[^#]' $sgldir/$name/$name.xmap 2>/dev/null | tr -d '\n'`;
		if [[ -n "$nhits" ]]; then
			refhits=$(count_xmaps $sgldir/$name/$name.xmap 3);
			qryhits=$(count_xmaps $sgldir/$name/$name.xmap 2);
			printf "%8s alignments are found for enzyme $name between %6s NGS contigs and %6s super-contigs.\n" $nhits $qryhits $refhits
		fi
		eval "ln -sf $name/$name.xmap $sgldir";
		((j++));
	fi
done
etime=$(ntime)
if [ $j -lt ${#names[@]} ]; then
	echo 'Error 2' > $stfile
	exit 1;
fi
eval "cat $sgldir/*.xmap | grep -v '^#' | cut -f2 | sort -n | uniq > $sgldir/used_NGS_id.txt";
echo -e "Single-channel aligment complete in $(duration $stime $etime) s.\n";

# prepare the combined coordinate translation file
combinedtrans="$outdir/combined_NGS_coord_translation.txt";
eval "(cat $qrydir/*_auto_cut_NGS_coord_translation.txt | sort -k1,1n -k2,2n -k3,3nr | uniq > $combinedtrans)";

echo -e "Beginning alignment of unused NGS contigs to the final scaffold using unified channel ... ";
stime=$(ntime)
unidir="$outdir/unified"
mkdir -p $unidir
if [ ! -f $unidir/final.xmap ]; then
	check "$collector -i $infasta -m $endir/${qrypre}_${enzymenames}_key.txt -t $combinedtrans -u $sgldir/used_NGS_id.txt -min 20000 -max 200000 -o $unidir/unused_fit_contigs > $unidir/contig-filtering.log"; # $enzymenames
	check "$fa2cmap -i $unidir/unused_fit_contigs.fa -e ${enzymeseqs} -o $unidir/encoding" # $enzymeseqs
	eval "(cd $unidir; ln -sf  encoding/unused_fit_contigs_${enzymenames}.cmap unused_fit_contigs_multicolors.cmap)"; # enzymenames
	check "$unifier -i $unidir/unused_fit_contigs_multicolors.cmap -o $unidir/unused_fit_contigs_unified"
	if [ ! -f $unidir/aligned.xmap ]; then
		check "$aligner -r $refdir/unified.cmap -q $unidir/unused_fit_contigs_unified.cmap -o $unidir/aligned --skip $shareparams$alignconf > $unidir/aligned.log"
	fi
	check "$xmapfilter -x $unidir/aligned.xmap -r $refdir/multicolors.cmap -q $unidir/unused_fit_contigs_multicolors.cmap > $unidir/filtered.xmap";
	check "$idmap -i $unidir/filtered.xmap -k $unidir/encoding/unused_fit_contigs_${enzymenames}_key.txt -o $unidir/final.xmap";
fi
if [ -f $unidir/final.xmap ]; then
	nhits=`grep -c '^[^#]' $unidir/final.xmap 2>/dev/null | tr -d '\n'`;
	if [[ -n "$nhits" ]]; then
		refhits=$(count_xmaps $unidir/final.xmap 3);
		qryhits=$(count_xmaps $unidir/final.xmap 2);
		printf "%8s alignments are found for between %6s NGS contigs and %6s super-contigs.\n" $nhits $qryhits $refhits
	fi
fi
etime=$(ntime)
echo -e "Unified-channel aligment complete in $(duration $stime $etime) s.\n";

# prepare the combined XMAP file
finalxmap="$outdir/combined.xmap"
header=("# XMAP File Version:\t0.2"
		"# Label Channels:\t$j"
		"# Reference Maps From:\t$refdir/mono/*.cmap"
		"# Query Maps From:\t$qrydir/${qrypre}_*_$siffix.cmap"
		"#h XmapEntryID\tQryContigID\tRefContigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\tQryLen\tRefLen\tLabelChannel\tAlignment"
		"#f int        \tint        \tint        \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat     \tstring \tfloat \tfloat \tint         \tstring");
# print the header
`:>$finalxmap`;
for ((i=0; i<${#header[@]}; i++)); do
	line=${header[$i]};
	echo -e $line >> $finalxmap;
done
# print the content
k=0;
for ((i=0; i<${#names[@]}; i++)); do
	name=${names[$i]};
	((j=$i+1))
	eval "grep -v '^#' $sgldir/$name.xmap | sort -k3,3n -k6,6n -k7,7nr | awk -v b=\"$k\" -v no=\"$j\" '{printf b+NR; for(i=2; i<=NF; i++){printf \"\\t\"; if(i!=13)printf \$i; else printf no;} printf \"\\n\"}' >> $finalxmap";
	nr=`grep -vc '^#' $sgldir/$name.xmap`;
	((k+=$nr))
done
eval "grep -v '^#' $unidir/final.xmap | sort -k3,3n -k6,6n -k7,7nr | awk -v b=\"$k\" -v no=0 '{printf b+NR; for(i=2; i<=NF; i++){printf \"\\t\"; if(i!=13)printf \$i; else printf no;} printf \"\\n\"}' >> $finalxmap";

echo 'Done' > $stfile
exit 0;
