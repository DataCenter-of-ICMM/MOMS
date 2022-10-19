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

if [ $# -lt 3 ]; then
	echo "Usage: $0 refpath qrypath outdir [nthreads [xmldir]]";
	exit 1;
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$DIR/util.sh"
aligner="$DIR/perl/cmap_aligner.pl"
identifier="$DIR/perl/conflicts_identifier.pl"
cutter="$DIR/perl/conflicts_cutter.pl"

refpath=$1
qrypath=$2
outdir=$3
stfile="$outdir/status.txt"

nthreads=0
alignconf=""
xmlconf=""
if [ $# -ge 4 ]; then
	nthreads=$(echo ${4} | egrep '^[0-9]+$');
	if [ $# -ge 5 ]; then
		xmldir=$5
		if [[ -d $xmldir ]]; then
			alignconf=" -x $xmldir/alignArguments.xml";
			xmlconf=" -x $xmldir/conflictsArguments.xml"
		fi
	fi
fi
if [[ "$refpath" != *\/*[^\/] ]]; then
	echo "refpath must be in \"dir/prefix\" format";
	echo 'Error 1' > $stfile
	exit 1;
fi
refdir=${refpath%/*};
refpre=${refpath##*/};
if [[ "$qrypath" != *\/*[^\/] ]]; then
	echo "qrypath must be in \"dir/prefix\" format";
	echo 'Error 2' > $stfile
	exit 1;
fi
qrydir=${qrypath%/*};
qrypre=${qrypath##*/};

shareparams="";
if [ $nthreads -gt 0 ]; then
	shareparams+=" -t $nthreads";
fi

# alignment
declare -A counts
for ref in `ls -1 ${refdir}/${refpre}_*.cmap`; do
	rname=${ref##*/${refpre}_};
	rname=${rname%.cmap};
	((counts[$rname]++));
done
for qry in `ls -1 ${qrydir}/${qrypre}_*_adjusted.cmap`; do
	qname=${qry##*/${qrypre}_};
	qname=${qname%_adjusted.cmap};
	((counts[$qname]++));
done

# force the enzymes in an order, which is very important for downstream analysis
eval "declare -a ordered=`orderedEnzymes`"

# find the common enzyme names
for name in ${ordered[@]}; do
	if [[ -n "${counts[$name]}" ]]; then
		if [ ${counts[$name]} -eq 2 ]; then
			names+=($name);
		fi
	fi
done

# the 2nd round alignments
alndir="$outdir/align1";
mkdir -p $alndir
echo -e "Beginning initial NGS CMAP to rescaled BioNano CMAP alignment ... ";
stime=$(ntime)
for name in ${names[@]}; do
	refmap="$refdir/${refpre}_$name.cmap"
	qrymap="$qrydir/${qrypre}_${name}_adjusted.cmap"
	if [ ! -f "$alndir/$name.xmap" ]; then
		check "$aligner -r $refmap -q $qrymap -o $alndir/$name $shareparams$alignconf > /dev/null";
	fi
	nhits=`grep -c '^[^#]' $alndir/$name.xmap 2>/dev/null | tr -d '\n'`;
	if [[ -n "$nhits" ]]; then
		refhits=$(count_xmaps $alndir/$name.xmap 3);
		qryhits=$(count_xmaps $alndir/$name.xmap 2);
		printf "%8s alignments are found for enzyme $name between %6s BNG contigs and %6s NGS contigs.\n" $nhits $qryhits $refhits
	fi
done
etime=$(ntime)
echo -e "Initial rescaled alignment complete in $(duration $stime $etime) s.\n";

cfldir="$outdir/conflicts"
mkdir -p $cfldir

declare -A used
# identification of conflicts
echo -e "Beginning conflicts identification ...";
stime=$(ntime)
for xmap in `ls -1 $alndir/*.xmap`; do
	name=${xmap##*/};
	name=${name%.xmap};
	((used[$name]++));
	refmap="$alndir/${name}_r.cmap";
	qrymap="$alndir/${name}_q.cmap";
	ref0map="$refdir/${refpre}_$name.cmap";
	qry0map="$qrydir/${qrypre}_${name}_adjusted.cmap";
	ref0maps=$(count_cmaps $ref0map);
	qry0maps=$(count_cmaps $qry0map);
	check "$identifier -i $xmap -r $refmap -q $qrymap -r0 $ref0map -q0 $qry0map -o $cfldir/$name$xmlconf >/dev/null";
	refhits=$(count_xmaps $cfldir/$name.xmap 3); [ $refhits -eq 0 ] && refhits="NONE"
	qryhits=$(count_xmaps $cfldir/$name.xmap 2); [ $qryhits -eq 0 ] && qryhits="NONE"
	printf "%8s of %6s BNG contigs have been flagged as conflicting for enzyme $name.\n" $qryhits $qry0maps
	printf "%8s of %6s NGS contigs have been flagged as conflicting for enzyme $name.\n" $refhits $ref0maps
done
etime=$(ntime)
echo -e "Conflicts identification complete in $(duration $stime $etime) s.\n";

echo -e "Beginning conflicts cutting ...";
stime=$(ntime)
i=0;
for name in ${names[@]}; do
	xmap="$alndir/$name.xmap";
	refmap="$alndir/${name}_r.cmap";
	qrymap="$alndir/${name}_q.cmap";
	ref0map="$refdir/${refpre}_$name.cmap";
	qry0map="$qrydir/${qrypre}_${name}_adjusted.cmap";
	conflict="$cfldir/${name}_conflicts.txt";
	((i++));
	check "$cutter -i $xmap -r $refmap -q $qrymap -r0 $ref0map -q0 $qry0map -c $conflict -s $i,${#names[@]} -o $outdir/$name $xmlconf >/dev/null";
	refmaps=$(count_cmaps $outdir/${refpre}_${name}_cut.cmap);
	qrymaps=$(count_cmaps $outdir/${qrypre}_${name}_adjusted_cut.cmap);
	printf "%8s BNG contigs found for enzyme $name after conflicts cutting.\n" $qrymaps
	printf "%8s NGS contigs found for enzyme $name after conflicts cutting.\n" $refmaps
done
etime=$(ntime)
echo -e "Conflicts cutting complete in $(duration $stime $etime) s.\n";

allenzymes=$( IFS="_"; echo "${names[*]}" );
eval "(cp $refdir/${refpre}_${allenzymes}.cmap $outdir/${refpre}_${allenzymes}.cmap)";

echo 'Done' > $stfile
exit 0;
