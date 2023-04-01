#!/bin/bash
NUMSEQ=10000
LENSEQ=10000
OUTDIR=/tmp/parity
mkdir -p $OUTDIR
/home/zhaom/ipuma-lib/build/bin/ipusw \
 --numVertices 1472  --fillAlgo fillfirst --seedLength 17 \
 --maxComparisonsPerVertex 10 \
 --bandPercentageXDrop 0.45 \
 --maxSequenceLength 20000 \
 --output $OUTDIR/scores_ipuint_sim05.json \
 --vtype xdroprestrictedseedextend \
 --generateSequences \
 --generatorCount $NUMSEQ \
 --generatorLength $LENSEQ \
 --generatorSimilarity 0.5

ALGO=seqan
/home/zhaom/ipuma-lib/build/bin/cpusw \
 --output $OUTDIR/scores_cpu_${ALGO}_sim05.json \
 --algo seqan \
 --threads 96 \
 --generateSequences \
 --generatorCount $NUMSEQ \
 --generatorLength $LENSEQ \
 --generatorSimilarity 0.5

ALGO=ipumacpu
/home/zhaom/ipuma-lib/build/bin/cpusw \
 --output $OUTDIR/scores_cpu_${ALGO}_sim05_orderflipped.json \
 --algo $ALGO \
 --threads 96 \
 --generateSequences \
 --generatorCount $NUMSEQ \
 --generatorLength $LENSEQ \
 --generatorSimilarity 0.5