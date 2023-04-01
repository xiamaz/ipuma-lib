#!/bin/bash
#gdb --args /home/zhaom/ipuma-lib/build/bin/cpusw \
/home/zhaom/ipuma-lib/build/bin/cpusw \
 --output /tmp/simplecorrectness/scores_cpu_ipuma.json \
 --algo libgaba \
 --threads 96 \
 --seqsH /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seqs_seed_H.txt \
 --seqsV /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seqs_seed_V.txt \
 --seedsH1 /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seeds_H1.txt \
 --seedsV1 /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seeds_V1.txt \
 --seedsH2 /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seeds_H2.txt \
 --seedsV2 /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seeds_V2.txt
 # --generateSequences \
 # --generatorCount 100 \
 # --generatorLength 100 \
 # --generatorSimilarity 1