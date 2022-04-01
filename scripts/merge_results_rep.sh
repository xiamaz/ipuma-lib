#!/bin/bash

if [ -z $2 ]; then
	echo "$0 [dataset prefix] [destination csv]"
	exit
fi

EXPORTSCRIPT=./scripts/extract.py
MERGED_OUTPUT=$2

if [ -f $MERGED_OUTPUT ]; then
	echo "${MERGED_OUTPUT} already exists."
	exit
fi

for repdir in $1*; do
	repnum=${repdir##*_}
	# $EXPORTSCRIPT --move_failed 1 "$repdir" "/tmp/${repnum}.csv"
	$EXPORTSCRIPT "$repdir" "/tmp/${repnum}.csv"
	if [ $repnum = "rep1" ]; then
		cat "/tmp/${repnum}.csv" >> $MERGED_OUTPUT
	else
		tail -n+1 "/tmp/${repnum}.csv" >> $MERGED_OUTPUT
	fi
done
