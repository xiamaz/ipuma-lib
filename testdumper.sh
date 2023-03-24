#!/bin/bash
DSDUMPER=./build/bin/dsdumper

LOGAN_ARGS="--seqsH /home/zhaom/ipuma-lib/output/logan_100k/seqH_100k.txt \
 --seqsV /home/zhaom/ipuma-lib/output/logan_100k/seqV_100k.txt \
 --seedsH1 /home/zhaom/ipuma-lib/output/logan_100k/seedH_100k.txt \
 --seedsV1 /home/zhaom/ipuma-lib/output/logan_100k/seedV_100k.txt"

GENERATOR1_ARGS="--generatorCount 20000 \
 --generatorSeqLen 20000 \
 --generatorSimilarity 1"
GENERATOR0_ARGS="--generatorCount 20000 \
 --generatorSeqLen 20000 \
 --generatorSimilarity 0"
GENERATOR_ARGS="--generatorCount 20000 \
 --generatorSeqLen 20000 \
 --generatorSimilarity 0.85"

TARGS="--generatorCount 2 \
 --generatorSeqLen 20000 \
 --generatorSimilarity 0.85"

dsname=logan
dsargs="$LOGAN_ARGS"
$DSDUMPER $dsargs \
 --outputCmps "/global/D1/projects/ipumer/datasets_giulia/cmps_${dsname}_multi.json" \
 --outputSequences "/global/D1/projects/ipumer/datasets_giulia/seqs_${dsname}_multi.json" \

dsname=simulated0
dsargs="$GENERATOR0_ARGS"
$DSDUMPER $dsargs \
 --outputCmps "/global/D1/projects/ipumer/datasets_giulia/cmps_${dsname}_multi.json" \
 --outputSequences "/global/D1/projects/ipumer/datasets_giulia/seqs_${dsname}_multi.json" \

dsname=simulated1
dsargs="$GENERATOR1_ARGS"
$DSDUMPER $dsargs \
 --outputCmps "/global/D1/projects/ipumer/datasets_giulia/cmps_${dsname}_multi.json" \
 --outputSequences "/global/D1/projects/ipumer/datasets_giulia/seqs_${dsname}_multi.json" \

dsname=simulated85
dsargs="$GENERATOR_ARGS"
$DSDUMPER $dsargs \
 --outputCmps "/global/D1/projects/ipumer/datasets_giulia/cmps_${dsname}_multi.json" \
 --outputSequences "/global/D1/projects/ipumer/datasets_giulia/seqs_${dsname}_multi.json" \