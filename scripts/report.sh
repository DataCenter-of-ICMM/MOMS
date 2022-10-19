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
	echo "Usage: $0 xmappath fastapath cmappath encpath outpath";
	exit 1;
fi

DIR="$( cd  "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$DIR/util.sh"
exporter="$DIR/perl/ExportAGP_TwoEnzyme.pl"

xmappath=$1
fastapath=$2
cmappath=$3
encpath=$4
outpath=$5
xmapdir=${xmappath%/*}
xmapfile=${xmappath##*/}
fastadir=${fastapath%/*}
fastafile=${fastapath##*/}
cmapdir=${cmappath%/*}
cmapfile=${cmappath##*/}
encdir=${encpath%/*}
encpre=${encpath##*/}
outdir=${outpath%/*}
outpre=${outpath##*/}
rptdir=${outdir%/*}
stfile="$outdir/status.txt"

if [[ "$xmappath" != *\/*[^\/] ]]; then
	echo "xmappath must be in \"dir/prefix\" format";
	echo 'Error 1' > $stfile
	exit 1;
fi
if [[ "$fastapath" != *\/*[^\/] ]]; then
	echo "fastapath must be in \"dir/prefix\" format";
	echo 'Error 2' > $stfile
	exit 1;
fi
if [[ "$cmappath" != *\/*[^\/] ]]; then
	echo "cmappath must be in \"dir/prefix\" format";
	echo 'Error 3' > $stfile
	exit 1;
fi
if [[ "$encpath" != *\/*[^\/] ]]; then
	echo "encpath must be in \"dir/prefix\" format";
	echo 'Error 4' > $stfile
	exit 1;
fi
if [[ "$outpath" != *\/*[^\/] ]]; then
	echo "outpath must be in \"dir/prefix\" format";
	echo 'Error 5' > $stfile
	exit 1;
fi

declare -A used
for map in `ls -1 $fastadir/assembly_*.cmap`; do
	name=${map##*/assembly_};
	name=${name%.cmap};
	((used[$name]++));
done

# force the enzymes in an order, which is very important for downstream analysis
eval "declare -a ordered=`orderedEnzymes`"
for name in ${ordered[@]}; do
	if [ -n "${used[$name]}" ]; then
		names+=($name);
	fi
done

allenzymes=$( IFS="_"; echo "${names[*]}" );
keyfile="${encpre}_${allenzymes}_key.txt";

mkdir -p $outdir
echo  -e "Beginning to export AGP and FASTA files for multi-channel CMAPs ..."
stime=$(ntime)
check "$exporter -i $xmapdir/$xmapfile -c $cmapdir/$cmapfile -s $fastadir/$fastafile -o $outdir/$outpre -m $encdir/$keyfile -t $xmapdir/combined_NGS_coord_translation.txt 2>&1 > $outdir/export.log";
etime=$(ntime)
echo -e "Exporting AGP and FASTA files complete in $(duration $stime $etime) s.\n";

check "ls $rptdir/Step*/*.cmap | xargs $DIR/perl/calc_cmap_stats.pl -f -o $outdir/cmap_file_stats.txt";
check "ls -1 $outdir/$outpre*.fa* | xargs $DIR/perl/calc_fasta_stats.pl -f > $outdir/fasta_file_stats.txt";
check "$DIR/perl/calc_xmap_stats.pl -q $xmapdir/$xmapfile -r $cmapdir/$cmapfile > $outdir/xmap_file_stats.txt";

echo 'Done' > $stfile
exit 0;
