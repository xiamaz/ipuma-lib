#!/bin/bash

EXPORTSCRIPT=./scripts/extract.py
MERGED_OUTPUT=output.csv

if [ -f $MERGED_OUTPUT ]; then
	echo "${MERGED_OUTPUT} already exists."
	exit
fi

for repdir in $1*; do
	repnum=${repdir##*_}
	$EXPORTSCRIPT --move_failed 1 "$repdir" "${repnum}.csv"
	cat "${repnum}.csv" >> $MERGED_OUTPUT
done