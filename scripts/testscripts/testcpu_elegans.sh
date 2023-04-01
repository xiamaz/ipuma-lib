#!/bin/bash
mkdir -p /tmp/simplecorrectness
numactl -N 1 -m 1 /home/zhaom/ipuma-lib/build/bin/cpusw \
 --output /tmp/simplecorrectness/scores_ksw_elegans.json \
 --algo ksw2 \
 --threads 48 \
 --xDrop 400 \
 --comparisons /global/D1/projects/ipumer/inputs_ab/cmps_celegans_multi.json \
 --sequences /global/D1/projects/ipumer/inputs_ab/seqs_celegans_multi.json