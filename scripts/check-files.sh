#!/usr/bin/bash
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
	echo "Usage: $0 list.txt search-dir suffix";
	exit 1;
fi

listfile=$1
searchdir=$2
suffix=$3

if [ ! -f $listfile ]; then
	echo "Error: $listfile does not exist";
	exit 1
fi

declare -A defined
cnt1=0;
for obj in `cut -d':' -f1 $listfile`; do
	name=${obj^^};
	defined[$name]=1;
	((cnt1++))
done

cnt2=0
for fnd in `ls -1 $searchdir/*.$suffix`; do
	name=${fnd##*/};
	name=${name%.$suffix}
	if [[ ${defined[$name]} != 0 ]]; then
		((cnt2++))
	fi
done

if [[ $cnt1 == $cnt2 ]]; then
	exit 0;
fi

exit 1
