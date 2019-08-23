ntime() {
	echo "$(date +%s%N | cut -b1-13)";
}

duration(){
	elapsed=`echo "scale=3;($2-$1)/1000" | bc -l`;
	echo "$elapsed";
}

check() {
	echo "Running command: ${1}"
	stime=$(ntime)
	eval "(${1})"
	RESULT=$?
	if [ "${RESULT}" -ne 0 ]; then
		etime=$(ntime)
		echo "FAILED ($(duration $stime $etime)s)";
		exit ${RESULT}
	fi
}

abs_path(){
	echo "$(cd $(dirname $1); pwd)/$(basename $1)";
}

count_cmaps() {
	local result=$(grep '^[^#]' $1 2>/dev/null| cut -f1 | uniq | wc -l)
	echo "$result"
}

count_xmaps() {
	local result=$(grep '^[^#]' $1 2>/dev/null| cut -f$2 | sort -n | uniq | wc -l)
	echo "$result"
}

## for enzyme ordering
orderedEnzymes()
{
	local -a arr=(BSPQI BSSSI BBVCI BSMI BSRDI BSECI HAMHI DLE1)
	declare -p arr | sed -e 's/^declare -a [^=]*=//'
}

## for retrieving enzyme name
enzymeNames()
{
	local -A names=(["GCTCTTC"]="BspQI" ["CACGAG"]="BssSI" ["CCTCAGC"]="BbvCI" ["GAATGC"]="BsmI" ["GCAATG"]="BsrDI" ["ATCGAT"]="bseCI" ["GGATCC"]="BamHI" ["CTTAAG"]="DLE1")
	declare -p names | sed -e 's/^declare -A [^=]*=//'
}

