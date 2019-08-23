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

if [ $# -lt 3 ]; then
	echo "Usage: $0 refpath qrypath outdir [nthreads [xmldir]]";
	exit 1;
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$DIR/util.sh"
aligner="$DIR/perl/cmap_aligner.pl"
refaligner="$DIR/bionano/binary/RefAligner"

refpath=$1
qrypath=$2
outdir=$3
stfile="$outdir/status.txt"

nthreads=0
alignconf=""
if [ $# -ge 4 ]; then
	nthreads=$(echo ${4} | egrep '^[0-9]+$');
	if [ $# -ge 5 ]; then
		xmldir=$5
		if [[ -d $xmldir ]]; then
			alignconf=" -x $xmldir/alignArguments.xml";
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

declare -A counts
for ref in `ls -1 ${refdir}/${refpre}_*.cmap`; do
	rname=${ref##*/${refpre}_};
	rname=${rname%.cmap};
	((counts[$rname]++));
done
for qry in `ls -1 ${qrydir}/${qrypre}_*.cmap`; do
	qname=${qry##*/${qrypre}_};
	qname=${qname%.cmap};
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
# alignment
echo -e "Beginning initial NGS CMAP to BioNano CMAP alignment ... ";
stime=$(ntime)
for name in ${names[@]}; do
	refmap="$refdir/${refpre}_$name.cmap"
	qrymap="$qrydir/${qrypre}_$name.cmap"
	if [ ! -f "$outdir/align0/$name.xmap" ]; then
		check "$aligner -s first -r $refmap -q $qrymap -o $outdir/align0/$name $shareparams$alignconf > /dev/null";
	fi
done
etime=$(ntime)
echo -e "Initial alignment complete in $(duration $stime $etime) s.\n";

# rescale the query CMAPs
echo -e "Rescaling BioNano CMAP ... ";
stime=$(ntime)
for name in ${names[@]}; do
	qrymap="$qrydir/${qrypre}_$name.cmap"
	errbin="$outdir/align0/$name.errbin"
	adjusted="$outdir/align0/${qrypre}_${name}_adjusted"
	if [ ! -f "$adjusted.cmap" ]; then
		check "$refaligner -merge -i $qrymap -readparameters $errbin -o $adjusted -stdout -stderr > /dev/null";
	fi
	if [ -f "$adjusted.cmap" ]; then
		eval "ln -sf align0/${qrypre}_${name}_adjusted.cmap $outdir";
	fi
done
etime=$(ntime)
echo -e "Rescaling complete in $(duration $stime $etime) s.\n";

echo 'Done' > $stfile
exit 0;

:<<"END"
# the 2nd round alignments
echo -e "Beginning initial NGS CMAP to rescaled BioNano CMAP alignment ... ";
stime=$(ntime)
for name in ${names[@]}; do
	refmap="$refdir/${refpre}_$name.cmap"
	qrymap="$outdir/align0/${qrypre}_${name}_adjusted.cmap"
	if [ ! -f "$outdir/align1/$name.xmap" ]; then
		check "$aligner -r $refmap -q $qrymap -o $outdir/align1/$name $shareparams$alignconf > /dev/null";
	fi
	nhits=`grep -c '^[^#]' $outdir/align1/$name.xmap 2>/dev/null | tr -d '\n'`;
	if [[ -n "$nhits" ]]; then
		refhits=$(count_xmaps $outdir/align1/$name.xmap 3);
		qryhits=$(count_xmaps $outdir/align1/$name.xmap 2);
		printf "%8s alignments are found for enzyme $name between %6s BNG contigs and %6s NGS contigs.\n" $nhits $qryhits $refhits
	fi
	eval "ln -sf align1/$name.xmap $outdir";
done
etime=$(ntime)
echo -e "Initial rescaled alignment complete in $(duration $stime $etime) s.\n";

echo 'Done' > $stfile
exit 0;
END
