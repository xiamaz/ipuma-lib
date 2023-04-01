#!/bin/bash
/home/zhaom/ipuma-lib/build/bin/ipusw \
 --numVertices 1472  --fillAlgo fillfirst --seedLength 17 \
 --maxComparisonsPerVertex 10 \
 --bandPercentageXDrop 0.45 \
 --maxSequenceLength 20000 \
 --output scores_ipu.json \
 --vtype xdroprestrictedseedextend \
 --generateSequences \
 --generatorCount 10000 \
 --generatorLength 20000 \
 --generatorSimilarity 0
 # /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seqs_seed_H.txt \
 # /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seqs_seed_V.txt \
 # /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seeds_H1.txt \
 # /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seeds_V1.txt
 # /global/D1/projects/ipumer/inputs_ab/seqs_seed_H.txt /global/D1/projects/ipumer/inputs_ab/seqs_seed_V.txt /global/D1/projects/ipumer/inputs_ab/seeds_H.txt /global/D1/projects/ipumer/inputs_ab/seeds_V.txt