#!/bin/bash

CMD="./build/bin/ipusw --numDevices 32 --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 2000 --bufsize 400000 --fillAlgo greedy -- /global/D1/projects/ipumer/datasets/compare/dna/DNA_2_200_ref.fasta /global/D1/projects/ipumer/datasets/compare/dna/DNA_2_200_qer.fasta"

COUNT=100

ERROR_EXECS=0
for i in `seq 1 ${COUNT}`; do 
	$CMD
	ret=$?
	if [ $ret -ne 0 ]; then
		ERROR_EXECS=$((${ERROR_EXECS}+1))
	fi
done

echo "Failed ${ERROR_EXECS}/${COUNT}"