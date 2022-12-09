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

if [ $# -lt 1 ]; then
	echo "USAGE: $0 result_dir_of_MOMS"
	echo "   OR  $0 encpath respath"
	exit 1
fi

if [ $# -ge 2 ]; then
	d=$1;
	if [ -z "$d" ]; then
		>&2 echo "ERROR: \"$d\" does not exist."
		exit 255;
	elif [ ! -d "$d" ]; then
		>&2 echo "ERROR: \"$d\" is not a valid directory."
		exit 255;
	fi
	qrydir=$(ls -1d $2/conflicts 2>/dev/null | sed -n '1p');
	if [ -z "$qrydir" ]; then
		>&2 echo "ERROR: \"$2/conflicts\" does not exist."
		exit 255;
	elif [ ! -d $qrydir ]; then
		>&2 echo "ERROR: \"$2/conflicts\" is not a valid directory.";
		exit 255;
	fi
else
	dir=$1;
	d=$(ls -1d $dir/Step-*_NGS_CMAP_encoding/ 2>/dev/null | sed -n '1p');
	if [ -z "$d" ]; then
		>&2 echo "ERROR: \"$dir/Step-*_NGS_CMAP_encoding\" does not exist."
		exit 255;
	elif [ ! -d $d ]; then
		>&2 echo "ERROR: \"$dir/Step-*_NGS_CMAP_encoding\" is not a valid directory."
		exit 255;
	fi
	qrydir=$(ls -1d $dir/Step-*_Chimeras_resolution/conflicts 2>/dev/null | sed -n '1p');
	if [ -z "$qrydir" ]; then
		>&2 echo "ERROR: \"$dir/Step-*_Chimeras_resolution/conflicts\" does not exist."
		exit 255;
	elif [ ! -d $qrydir ]; then
		>&2 echo "ERROR: \"$dir/Step-*_Chimeras_resolution/conflicts\" is not a valid directory.";
		exit 255;
	fi
fi
f=$(ls -1 $d/*_*.cmap 2>/dev/null | sed -n '1p');
if [ -z "$f" ]; then
	>&2 echo "ERROR: \"$d\" does not contain valid CMAP files.";
	exit 255;
fi

refdir=${f%/*};
prefix=${f##*/}; prefix=${prefix%%_*};

nWidth=80
offset='  '
function formatLine(){
	str=$*
	len=$( expr $nWidth - ${#str} - 2 )
	printf "$offset|$str%${len}s|\n"
}

function count_cmaps(){
    local result=$(grep '^[^#]' $1 2>/dev/null| cut -f1 | uniq | wc -l)
    echo "$result"
}

function count_xmaps(){
	infile=${1/\\*/*}
	local result=$(grep '^[^#]' $infile 2>/dev/null| cut -f$2 | sort -n | uniq | wc -l)
	echo "$result"
}

sEnzymes=( $( ls -1 $qrydir/*.xmap | sed 's/^.*\///' | sed 's/\..*$//' | tr '\n' ' ' ) )
for e in ${sEnzymes[@]}
do
	if [ -z $nContigs ]; then
		ref="${refdir}/${prefix}_$e.cmap"
		if [ -s "$ref" ]; then
			nContigs=$(count_cmaps $ref)
		fi
	fi
done

horizontal=$(printf %${nWidth}s | tr ' ' '-' | sed "s/^/$offset/" )
title='Conflicts between Bionano maps and NGS contigs'
len=$( expr $nWidth - ${#title} - 2 )
left=$( expr $len / 2 )
right=$( expr $len - $left )
title=$( printf "$offset|%${left}s\033[33m$title\033[0m%${right}s|" )

echo -e "\n$horizontal\n$title\n$horizontal"
nConflictContigs=$(count_xmaps "$qrydir/*.xmap" 3)
nConflicts=$( for f in `ls -1 $qrydir/*.xmap`; do grep -v '^#' $f | cut -f3 | sort -n; echo; done | awk 'BEGIN{id="";total=0}{if($1!=id){if(id!="")printf("%s\t%d\n", id, cnt); id=$1; cnt=1} else{cnt++}}' | sort -k1,1n | awk 'BEGIN{id="";total=0}{if($1!=id){if(id!="")total+=max; id=$1; max=$2}else{if($2>max)max=$2}}END{if(id!="")total+=max; print total}' )
formatLine " $nConflicts conflicts in $nConflictContigs / $nContigs contigs were flagged as conflicting; of these:"
for e in ${sEnzymes[@]}
do
	count=$( grep -v '^#' $qrydir/$e.xmap | wc -l )
	count2=$(count_xmaps $qrydir/$e.xmap 3)
	formatLine "   $count in $count2 contigs were supported by enzyme $e"
done

nResolvedContigs=$( cat  ${qrydir%/*}/*_auto_cut_NGS_coord_translation.txt | awk '{if($1~/^[0-9]/) print}' | cut -f1 | uniq -c | sed 's/^ \+//' | awk '{if($1>1)print $2}' | sort -n | uniq | wc -l );
nResolved=$( cut -f1 ${qrydir%/*}/*_auto_cut_NGS_coord_translation.txt | awk 'BEGIN{id="";total=0}{if($1!=id){if(cnt>0)printf("%s\t%d\n", id, cnt); id=$1; cnt=0} else{cnt++}}' | sort -k1,1n | awk 'BEGIN{id="";total=0}{if($1!=id){if(id!="")total+=max; id=$1; max=$2}else{if($2>max)max=$2}}END{if(id!="")total+=max; print total}' )
formatLine " $nResolved conflicts in $nResolvedContigs contigs were confirmed and resolved; of these:"
for e in ${sEnzymes[@]}
do
	count=$( cut -f1 ${qrydir%/*}/${e}_auto_cut_NGS_coord_translation.txt | awk 'BEGIN{id="";total=0}{if($1!=id){if(cnt>0)total+=cnt; id=$1; cnt=0} else{cnt++}}END{print total}' );
	count2=$( cut -f1 ${qrydir%/*}/${e}_auto_cut_NGS_coord_translation.txt | uniq -c | sed 's/^ \+//' | awk '{if($1>1)print $2}' | wc -l );
	formatLine "   $count in $count2 contigs were identified by enzyme $e"
done
echo -e "$horizontal\n"

echo -e "Detailed files are stored in $qrydir\n";

exit 0
