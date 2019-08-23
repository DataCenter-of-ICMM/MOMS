#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if  [[ ! -f "$DIR/bionano/binary/RefAligner" ]]; then
	LOCATION="~/tools/pipeline/1.0/RefAligner/1.0/";
	echo -e "\033[31mPlease specify the location of RefAligner:\033[0m [eg. $LOCATION]";
	read LOCATION
	while [[ ! -d "$LOCATION" || ! -f "$LOCATION/RefAligner" ]]; do
		if [ -z $LOCATION ]; then
			break
		fi
		echo -e "Can not find RefAligner @ the specified location";
		read LOCATION
	done
	if [[ -d "$LOCATION" && -f "$LOCATION/RefAligner" ]]; then
		mkdir -p $DIR/bionano/binary
		ln -sf $LOCATION/RefAligner $DIR/bionano/binary
		echo -e " \033[32mOK\033[0m."
	fi
else
	echo -e "RefAligner exists";
fi
