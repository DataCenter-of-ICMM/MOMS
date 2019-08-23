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
	echo "Usage: $0 input.fa enzymes outdir [minNicks [minKlen]]";
	exit;
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
fa2cmap="$DIR/perl/fa2cmap_multi_color.pl"

infasta=$1
IFS=',' read -r -a array <<< "$2"
enzymes=`for i in ${!array[@]}; do str=$(echo ${array[$i]} | cut -d: -f2,3); echo -n "$str $(( i + 1 )) "; done; echo`
outdir=$3
stfile="$outdir/status.txt"

minnicks=0
minklen=0
if [ $# -ge 4 ]; then
	minnicks=$(echo $4 | egrep '^[0-9]+$');
	if [ $# -ge 5 ]; then
		minklen=$(echo $5 | egrep '^[0-9]+$');
	fi
fi

convert(){
	params="-i $1 -e $2"
	if [ $3 -gt 0 ]; then
		params+=" -m $3"
	fi
	if [ $4 -gt 0 ]; then
		params+=" -M $4"
	fi
	params+=" -o $5"

	cmd="$fa2cmap $params";
	echo "$cmd"
	eval "($cmd)"
	result=$?
	if [ "$result" -ne 0 ]; then
		echo 'Error 1' > $stfile
		exit $result
	fi
}

echo "Encoding with all enzyme channels ..."
convert $infasta "$enzymes" $minnicks $minklen $outdir

if [ ${#array[@]} -gt 1 ]; then
	echo "Encoding with single enzyme channels ..."
	for i in ${!array[@]}; do
		str=$(echo ${array[$i]} | cut -d: -f2,3);
		enzyme="$str 1";
		convert $infasta "$enzyme" $minnicks $minklen $outdir
	done
fi

echo "100% done" > $stfile

exit 0
