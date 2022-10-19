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

if [ $# -lt 4 ]; then
	echo "Usage: $0 scfdir scfsuf ngspath outdir [nthreads]";
	exit 1;
fi

# for DEBUG
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$DIR/util.sh"
xmapfilter="$DIR/perl/xmap_filter.pl"
associator="$DIR/perl/xmap_associator.pl"
lmfitter="$DIR/perl/fitsLinearModel.pl"
allocator="$DIR/perl/edges_allocator.pl"
layout="$DIR/perl/clusterAndLayout.pl"
multimerger="$DIR/perl/cmap_merger_multi_color.pl"
subset="$DIR/perl/cmap_subset.pl"

scfdir=$1
scfsuf=$2
ngspath=$3
outdir=$4
stfile="$outdir/status.txt"

nthreads=0
if [ $# -ge 5 ]; then
	nthreads=$(echo $5 | egrep '^[0-9]+$');
fi

shareparams="";
if [ $nthreads -gt 0 ]; then
	shareparams+=" -t $nthreads";
fi

if [[ "$ngspath" != *\/*[^\/] ]]; then
	echo "ngspath must be in \"dir/prefix\" format";
	echo 'Error 1' > $stfile
	exit 1;
fi
ngsdir=${ngspath%/*}
ngspre=${ngspath##*/};

# force the enzymes in an order, which is very important for downstream analysis
eval "declare -a ordered=`orderedEnzymes`"
declare -A valid
for name in ${ordered[@]}; do
	((valid[$name]++));
done
declare -A defined
for cmap in `ls -1 ${scfdir}/*.${scfsuf}.cmap`; do
	name=${cmap##*/};
	name=${name%.${scfsuf}.cmap};
	if [ -z "${valid[$name]}" ]; then
		echo "$name is not a valid motif";
		echo 'Error 1' > $stfile
		exit 1;
	fi
	defined[$name]=1;
done
for name in ${ordered[@]}; do
	if [[ -n "${defined[$name]}" ]]; then
		names+=($name);
	fi
done

mkdir -p $outdir

alndir="$scfdir/anchor";
swdir="$outdir/sandwich";
mono="$outdir/mono";

if [[ ${#names[@]} -lt 2 ]]; then
	echo 'Error 2' > $stfile
	exit 0;
fi

echo -e "Beginning to build relationship between the BN contigs from different enzyme assemblies ...";
stime=$(ntime)
mkdir -p $swdir
for ((i=0; i<${#names[@]}-1; i++)); do
	name1=${names[$i]};
	coord1="$ngsdir/${name1}_auto_cut_NGS_coord_translation.txt";
	for ((j=i+1; j<${#names[@]}; j++)); do
		name2=${names[$j]};
		coord2="$ngsdir/${name2}_auto_cut_NGS_coord_translation.txt";
		pairname="${name1}_${name2}";

		assprefix="glues_$pairname";
		fitprefix="edges_$pairname";
		check "$associator -f $alndir/$name1.xmap -s $alndir/$name2.xmap -t1 $coord1 -t2 $coord2 -o $swdir/$pairname/$assprefix";
		check "$lmfitter -i $swdir/$pairname/$assprefix.tsv -o $swdir/$pairname/$fitprefix";
	done
done
etime=$(ntime)
echo -e "Scaffolds/Contigs relationship building complete in $(duration $stime $etime) s.\n";

echo -e "Beginning to make the layout of the final assembly ...";
stime=$(ntime)
check "$allocator $swdir/*/edges_*.tsv $swdir/*/glues_*.tsv $swdir";
references="";
for ((i=0; i<${#names[@]}; i++)); do
	name=${names[$i]};
	references+=" -r $alndir/$name/${name}-NGS_r.cmap";
done
check "$layout $references -e $swdir/edges.tsv -g $swdir/glues.tsv -o $swdir/paths";
check "$multimerger $references -l $swdir/paths.tsv -o $outdir/multicolors";

mkdir -p $mono
for ((i=0; i<${#names[@]}; i++)); do
	name=${names[$i]};
	color=$(($i + 1))
	check "$subset -i $outdir/multicolors.cmap -c $color -o $mono/$name";
done
etime=$(ntime)
echo -e "Cmap assembly complete in $(duration $stime $etime) s.\n";

echo 'Done' > $stfile
exit 0;
