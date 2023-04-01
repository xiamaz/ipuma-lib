#!/bin/bash
mkdir -p /tmp/simplecorrectness
/home/zhaom/ipuma-lib/build/bin/ipusw \
 --numVertices 1472  --fillAlgo greedy --seedLength 17 \
 --maxComparisonsPerVertex 400 \
 --bandPercentageXDrop 0.45 \
 --maxSequenceLength 20000 \
 --output /tmp/simplecorrectness/scores_ipu.json \
 --vtype xdroprestrictedseedextend \
 --useMulticomparison \
 --comparisons /global/D1/projects/ipumer/inputs_ab/cmps_celegans_multi.json \
 --sequences /global/D1/projects/ipumer/inputs_ab/seqs_celegans_multi.json