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
	echo "Usage: $0 indir outdir genome-size molecular-len nsites pvalue [bHaplotype [bHuman [bCluster [bVerbose [nthreads]]]]]";
	exit;
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$DIR/util.sh"
assembler="$DIR/perl/bnx_assembler.pl"

indir=$1
outdir=$2
gsize=$3
mlen=$4
nsites=$5
pvalue=$6
bHaplotype=0
bHuman=0
bCluster=0
bVerbose=0
nthreads=0
if [ $# -ge 7 ]; then
	bHaplotype=$(echo $7 | egrep '^[0-9]+$')
	if [ $# -ge 8 ]; then
		bHuman=$(echo $8 | egrep '^[0-9]+$')
		if [ $# -ge 9 ]; then
			bCluster=$(echo $9 | egrep '^[0-9]+$')
			if [ $# -ge 10 ]; then
				bVerbose=$(echo ${10} | egrep '^[0-9]+$')
				if [ $# -ge 11 ]; then
					nthreads=$(echo ${11} | egrep '^[0-9]+$')
				fi
			fi
		fi
	fi
fi

share_params="-g $gsize -l $mlen -s $nsites -p $pvalue"
[ $bHaplotype -ne 0 ] && share_params+=" --haplotype";
[ $bHuman -ne 0 ] && share_params+=" --human";
for bnx in `ls -1 $indir/molecules_*.bnx`; do
	fname=${bnx##*/molecules_};
	enzyme=${fname%.bnx};

	echo "assembling for enzyme $enzyme ... ";
	stime=$(ntime)

	BNXV=$(grep -m 1 '^# BNX File Version' $bnx | sed 's/.*\.//');
	[ "${BNXV:-0}" -ge 3 ] && model="saphyr" || model="irys";

	params="-m $model -i $bnx $share_params -o $outdir/enzyme_$enzyme/";
	if [ $nthreads -gt 0 ]; then
		params+=" -t $nthreads";
	fi
	if [ $bCluster -gt 0 ]; then
		params+=" --cluster";
	fi
	check "mkdir -p $outdir/enzyme_$enzyme";
	if [ $bVerbose -eq 0 ]; then
		params+=" > $outdir/enzyme_$enzyme/bnx_assembler.stdout";
	fi
	check "$assembler $params";
	etime=$(ntime)
	cmap="$outdir/enzyme_$enzyme/contigs/molecules_${enzyme}_refineFinal1/MOLECULES_${enzyme}_REFINEFINAL1.cmap";
	if [ -f $cmap ]; then
		check "cp $cmap $outdir/assembly_${enzyme}.cmap";
		echo "done ($(duration $stime $etime)s)";
	else
		echo "FAILED ($(duration $stime $etime)s)";
	fi
done
