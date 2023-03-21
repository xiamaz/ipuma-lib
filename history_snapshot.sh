#!/bin/bash

set -euo pipefail

/home/lukb/miniconda3/bin/cmake --build `pwd`/build --config Release --target ipusw
out=`pwd`/single_tile
mkdir -p $out

globals="--numVertices 1 \
 --fillAlgo fillfirst \
 --complexityAlgo xdrop \
 --seedLength 17 \
 --bandPercentageXDrop 0.45 \
 --maxSequenceLength 19295 \
 --vertexBufferSize 160000 \
 --vtype xdroprestrictedseedextend"

synth="--generateSequences \
 --generatorCount 30000 \
 --generatorLength 15000 "

set -x
for X in 15 20; do

echo "X=$X"

build/bin/ipusw \
  $globals \
  $synth \
 --xDrop $X \
 --generatorSimilarity 0.85 |& tee "$out"/15_${X}.log

# build/bin/ipusw \
#   $globals \
#   $synth \
#  --xDrop $X \
#  --generatorSimilarity 0.98 |& tee "$out"/2_${X}.log

# build/bin/ipusw \
#   $globals \
#   $synth \
#  --xDrop $X \
#  --generatorSimilarity 0.99 |& tee "$out"/1_${X}.log

# build/bin/ipusw \
#   $globals \
#   $synth \
#  --xDrop $X \
#  --generatorSimilarity 0.995 |& tee "$out"/05_${X}.log


# build/bin/ipusw \
#   $globals \
#  --xDrop $X \
#  --seqsH /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seqs_seed_H.txt \
#  --seqsV /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seqs_seed_V.txt \
#  --seedsH1 /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seeds_H1.txt \
#  --seedsV1 /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seeds_V1.txt \
#  --seedsH2 /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seeds_H2.txt \
#  --seedsV2 /global/D1/projects/ipumer/xdrop_inputs_elba_ecoli/seeds_V2.txt |& tee "$out"/elba_${X}.log

done
