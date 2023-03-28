#!/bin/bash
set -uo pipefail
PROJECT_DIR=$(git rev-parse --show-toplevel)
BUILD_DIR="$PROJECT_DIR/build_$(hostname)"
# --- Build the CPUSW Binary on the hostname system
mkdir -p "$BUILD_DIR"
pushd "$BUILD_DIR"
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --target ipusw
popd
BIN="$PROJECT_DIR/build_$(hostname)/bin/ipusw"

if [[ ! -f $BIN ]]; then
	echo "Binary not found build failed."
	exit 1
fi

export POPLAR_LOG_LEVEL=ERR

ECOLI_ARGS="--comparisons /global/D1/projects/ipumer/inputs_ab/cmps_ecoli_multi_single.json \
--sequences /global/D1/projects/ipumer/inputs_ab/seqs_ecoli_multi_single.json"

ECOLI100_ARGS="--comparisons /global/D1/projects/ipumer/inputs_ab/cmps_elba100_multi.json \
 --sequences /global/D1/projects/ipumer/inputs_ab/seqs_elba100_multi.json"

ELEGANS_ARGS="--comparisons /global/D1/projects/ipumer/inputs_ab/cmps_celegans_multi.json \
 --sequences /global/D1/projects/ipumer/inputs_ab/seqs_celegans_multi.json"

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

 declare -A dsmap
dsmap[logan]="$LOGAN_ARGS"
dsmap[ecoli]="$ECOLI_ARGS"
dsmap[simulated1]="$GENERATOR1_ARGS"
dsmap[simulated0]="$GENERATOR0_ARGS"
dsmap[simulated85]="$GENERATOR_ARGS"
dsmap[ecoli100]="$ECOLI100_ARGS"
dsmap[elegans]="$ELEGANS_ARGS"

OUTDIR="$PROJECT_DIR/output/ipusw_benchmark_final/$(hostname)"
mkdir -p "$OUTDIR"

run() {
	RUNNAME="${DSNAME}_${ALGO}_x${XDROP}_decom${DECOMPOSE}"
  if [[ $DEVICES > 1 ]]; then
    RUNNAME="${RUNNAME}_ipu${DEVICES}"
  fi
	echo "Running $RUNNAME "
  if [[ $DECOMPOSE = "y" ]]; then
    EXTRAARGS="--decomposeMulticomparisons"
  else
    EXTRAARGS=""
  fi
  scorefile="$OUTDIR/${RUNNAME}_scores.json"
  if [[ -f $scorefile ]]; then
    echo "Scorefile $RUNNAME already computed"
    return
  fi
	if [[ $NUMACTL = "y" ]]; then
   numactl -N 1 -m 1 "$BIN" \
     --numVertices 1472  --fillAlgo greedy --seedLength 17 \
     --maxComparisonsPerVertex 400 \
     --bandPercentageXDrop 0.45 \
     --maxSequenceLength 20000 \
     --numDevices "$DEVICES" \
     --xDrop "$XDROP" \
     --vtype xdroprestrictedseedextend \
   	  --output "$scorefile" \
   	$DSARGS $EXTRAARGS |& tee "$OUTDIR/${RUNNAME}.log"
	else
   "$BIN" \
     --numVertices 1472  --fillAlgo greedy --seedLength 17 \
     --maxComparisonsPerVertex 400 \
     --bandPercentageXDrop 0.45 \
     --maxSequenceLength 20000 \
     --xDrop "$XDROP" \
     --vtype xdroprestrictedseedextend \
   	 --output "$OUTDIR/${RUNNAME}_scores.json" \
   	 $DSARGS $EXTRAARGS |& tee "$OUTDIR/${RUNNAME}.log"
	fi
}

NUMACTL=y
ALGO=ipuma
for DSNAME in "${!dsmap[@]}"; do
for DECOMPOSE in n y; do
for XDROP in 5 10 15 20 50; do
for DEVICES in 1 2 4 8 16 32; do
  DSARGS="${dsmap[$DSNAME]}"
  run
done
done
done
done