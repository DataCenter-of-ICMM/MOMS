#!/bin/bash

USAGE="USAGE: ./compressDeNovo_forIrysView.sh --targetFolder <path/to/assembly/output> --outputFolder <path/to/compressFileOutput> --prefix <name_for_output_file> [--gzip <force gzip output>] [--verbose <enable verbose output>]"
echo ""

PBZIP2=NO
VERBOSE=NO
GZIPFORCED=NO
HELP=NO
start=`date +%s`

echo "COMMAND: ${0##*/} ${@}"
echo ""

if hash pbzip2 2>/dev/null; 
then
	PBZIP2=YES
fi

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -t|-tf|--targetFolder)
    TARGET="$2"
    shift # past argument
    ;;
    -p|--prefix)
    PREFIX="$2"
    shift # past argument
    ;;
    -o|-of|--outputFolder)
    OUTPUTFOLDER="$2"
    shift # past argument
    ;;
    -h|--help)
    HELP=YES
    ;;
    -g|-gz|--gz|--gzip)
    PBZIP2=NO
    GZIPFORCED=YES
    ;;
    -v|--verbose)
    VERBOSE=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

if [ $HELP == YES ];
then
	echo $USAGE
	echo ""
	exit 0
fi

if [ -d "$TARGET" ];
then
	TARGETFOLDER=$(readlink -m -n $TARGET)
	echo "Compressing assembly folder: $TARGETFOLDER"
	#echo ""
else
	echo $USAGE
	echo ""
	exit -1
fi

if [ -d "$OUTPUTFOLDER" ];
then
	OUTPUTFOLDER=$(readlink -m -n $OUTPUTFOLDER)
	echo "Using output folder: $OUTPUTFOLDER"
	#echo ""
else
	mkdir -p $OUTPUTFOLDER
	if [ $? -ne 0 ]; 
	then
		echo "ERROR: Could not create directory $OUTPUTFOLDER"
		#echo ""
		echo $USAGE
		echo ""
		exit -1
	fi
	OUTPUTFOLDER=$(readlink -m -n $OUTPUTFOLDER)
	echo "Using output folder: $OUTPUTFOLDER"
	#echo ""
fi

if [ -z "$PREFIX" ];
then
	echo "ERROR: --prefix not defined"
	echo $USAGE
	echo ""
	exit -1
else
	if [ $PBZIP2 == YES ];
	then
		echo "Creating output file: $OUTPUTFOLDER/$PREFIX.tar.bz2"
	else
		echo "Creating output file: $OUTPUTFOLDER/$PREFIX.tar.gz"
		if [ $GZIPFORCED == YES ];
		then
			echo -e "\tGzip outout forced"
		fi
	fi
	echo ""
fi

if [ $VERBOSE == YES ];
then
	echo "Verbose output enabled"
	echo ""
fi

#MAIN COMPRESSION COMMANDS

CURRPATH=$(pwd)
cd $TARGETFOLDER
cd ..
#tar -cfv --no-recursion $OUTPUTFOLDER/$PREFIX.tar.gz $TARGETFOLDER/*
#CWD=$(pwd)
#echo "CWD= $CWD"

BASENAME=$(basename $TARGETFOLDER)

if [ $VERBOSE == NO ];
then
	echo -n "Collecting files..."
	find $BASENAME -maxdepth 1 -type f -print0 | xargs -0 tar --exclude='*.tar.gz' --exclude='*_of_*.bnx' --exclude='*.map' --exclude='*_refined_*' --exclude='*_group*' --exclude='all_sorted.bnx' -cf $OUTPUTFOLDER/$PREFIX.tar
	tar --exclude='*.tar.gz' --exclude='*.map' --exclude='*_refined_*' --exclude='*_group*' --append --file=$OUTPUTFOLDER/$PREFIX.tar $BASENAME/ref $BASENAME/contigs/alignmolvref/merge $BASENAME/contigs/exp_refineFinal1/alignref_final $BASENAME/contigs/alignmolvref/copynumber $BASENAME/contigs/exp_refineFinal1_sv
	tar --append --file=$OUTPUTFOLDER/$PREFIX.tar $BASENAME/contigs/exp_refineFinal1/alignmol/merge/EXP_REFINEFINAL1_merge.map
	tar --append --file=$OUTPUTFOLDER/$PREFIX.tar $BASENAME/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap
	cd $BASENAME/contigs/exp_refineFinal1/alignmol/merge 
	tar --exclude='*.tar.gz' --exclude='*.map' --exclude='*_refined_*' --exclude='*_group*' -zcf $OUTPUTFOLDER/molecules.tar.gz *
        cd $OUTPUTFOLDER
	tar --append --file=$OUTPUTFOLDER/$PREFIX.tar molecules.tar.gz
	echo "done!"
	echo ""

	echo -n "Compressing..."
	if [ $PBZIP2 == NO ];
	then
		gzip -f $OUTPUTFOLDER/$PREFIX.tar
		gzip -t $OUTPUTFOLDER/$PREFIX.tar.gz
	else
		pbzip2 -f $OUTPUTFOLDER/$PREFIX.tar
		pbzip2 -t $OUTPUTFOLDER/$PREFIX.tar.bz2
	fi
	echo "done!"
	echo ""
else
	echo "Collecting files..."
	find $BASENAME -maxdepth 1 -type f -print0 | xargs -0 tar --exclude='*.tar.gz' --exclude='*_of_*.bnx' --exclude='*.map' --exclude='*_refined_*' --exclude='*_group*' --exclude='all_sorted.bnx' -vcf $OUTPUTFOLDER/$PREFIX.tar
	tar --exclude='*.tar.gz' --exclude='*.map' --exclude='*_group*' --verbose --append --file=$OUTPUTFOLDER/$PREFIX.tar $BASENAME/ref $BASENAME/contigs/alignmolvref/merge $BASENAME/contigs/*_refineFinal1/alignref_final $BASENAME/contigs/alignmolvref/copynumber $BASENAME/contigs/*_refineFinal1_sv $BASENAME/contigs/auto_noise
	tar --append --verbose --file=$OUTPUTFOLDER/$PREFIX.tar $BASENAME/contigs/*_refineFinal1/alignmol/merge/*_REFINEFINAL1_merge.map
	tar --append --verbose --file=$OUTPUTFOLDER/$PREFIX.tar $BASENAME/contigs/*_refineFinal1/*_REFINEFINAL1.cmap
	alignmol=$BASENAME/contigs/*_refineFinal1/alignmol/merge
	if [ -d $alignmol ]; then
	    echo "Collecting" $alignmol
	    cd $alignmol
	    tar --exclude='*.tar.gz' --exclude='*.map' --exclude='*_refined_*' --exclude='*_group*' -zcvf $OUTPUTFOLDER/molecules.tar.gz *
            #cd $OUTPUTFOLDER
	    #tar --append --verbose --file=$OUTPUTFOLDER/$PREFIX.tar molecules.tar.gz
            cd $OUTPUTFOLDER/..
	    tar --append --verbose --file=$OUTPUTFOLDER/$PREFIX.tar $BASENAME/molecules.tar.gz
	else
	    echo "Alignmol dir missing:" $alignmol
	fi
	echo "...done!"
	echo ""

	echo "Compressing..."
	if [ $PBZIP2 == NO ];
	then
		gzip -fv $OUTPUTFOLDER/$PREFIX.tar
		gzip -t $OUTPUTFOLDER/$PREFIX.tar.gz
	else
		pbzip2 -fv $OUTPUTFOLDER/$PREFIX.tar
		pbzip2 -t $OUTPUTFOLDER/$PREFIX.tar.bz2
	fi
	echo "...done!"
	echo ""
fi
rm $OUTPUTFOLDER/molecules.tar.gz
cd $CURRPATH
end=`date +%s`
runtime=$((end-start))

echo "COMPLETE! in $runtime seconds"
echo ""
