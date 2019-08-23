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
	echo "Usage: $0 input.fa enzyme1.cmap[,enzyme2.cmap[,...]] outdir";
	exit 1;
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";
source "$DIR/util.sh"

infasta=$(abs_path $1)
IFS=',' read -r -a bngs <<< "$2"
outdir=$3
stfile="$outdir/status.txt"

mkdir -p $outdir
`(cd $outdir; ln -sf $infasta ${infasta##*/})`
`(cd $outdir; awk 'BEGIN{len=0}{if($0 ~ /^[^>]/){len+=length($0);}}END{print len}' $infasta > contig-length.txt)`

eval "declare -A enzymes=`enzymeNames`"
declare -A defined
cnt=0;
for i in ${!bngs[@]}; do
	bng=$(abs_path ${bngs[$i]});
	MOTIF=$(grep -m 1 '^# Nickase Recognition Site' $bng | sed 's/.*:[[:blank:]]\+//' | sed 's/;.*$//');
	enzyme=${enzymes[${MOTIF^^}]}
	if [[ -z "$enzyme" ]]; then
		echo "Unrecognizable nickase motif \"$MOTIF\" found in \"$bng\". Skip processing.";
		continue;
	fi
	enzymeuc=${enzyme^^};
	if [[ -f "$outdir/assembly_$enzymeuc.cmap" ]]; then
		echo "Duplicate nickase motif \"$MOTIF\" found in \"$bng\". Skip processing.";
		continue;
	fi
	eval "(cd $outdir; ln -sf $bng assembly_$enzymeuc.cmap)"
	if [ $? -eq 0 ]; then
		((cnt++));
		defined[$enzymeuc]=${MOTIF^^}
	fi
done

if [ $cnt -ne ${#bngs[@]} ]; then
	echo 'Error 1' > $stfile
	exit 1;
fi

# force the enzymes in an order, which is very important for downstream analysis
eval "declare -a ordered=`orderedEnzymes`"
enzymefile="$outdir/enzymes.txt";
:> $enzymefile;
cnt=0;
for enzyme in ${ordered[@]}; do
	if [[ -n "${defined[$enzyme]}" ]]; then
		((cnt++));
		echo "$cnt:$enzyme:${defined[$enzyme]}" >> $enzymefile;
	fi
done

if [ $cnt -ne ${#bngs[@]} ]; then
	echo 'Error 2' > $stfile
	exit 1;
fi

echo 'Done' > $stfile
exit 0;
