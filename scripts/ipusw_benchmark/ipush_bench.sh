#!/bin/bash
# DATASET_DIR="/global/D1/projects/ipumer/inputs_ab"
DATASET_DIR="/home/zhaom/alextit/datasets"
set -uo pipefail
PROJECT_DIR=$(git rev-parse --show-toplevel)
BUILD_DIR="$PROJECT_DIR/build_$(hostname)"
# --- Build the CPUSW Binary on the hostname system
mkdir -p "$BUILD_DIR"
pushd "$BUILD_DIR"
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --target ipusw
popd

BIN="$BUILD_DIR/bin/ipusw"
OUTDIR="$PROJECT_DIR/output/ipusw_benchmark_gc/$(hostname)"

if [[ ! -f $BIN ]]; then
	echo "Binary not found build failed."
	exit 1
fi

ECOLI_ARGS="--comparisons ${DATASET_DIR}/cmps_ecoli_multi_single.json \
--sequences ${DATASET_DIR}/seqs_ecoli_multi_single.json"

ECOLI100_ARGS="--comparisons ${DATASET_DIR}/cmps_ecoli100_multi.json \
 --sequences ${DATASET_DIR}/seqs_ecoli100_multi.json"

ELEGANS_ARGS="--comparisons ${DATASET_DIR}/cmps_celegans_multi.json \
 --sequences ${DATASET_DIR}/seqs_celegans_multi.json"

LOGAN_ARGS="--comparisons ${DATASET_DIR}/cmps_logan_multi.json \
 --sequences ${DATASET_DIR}/seqs_logan_multi.json"

SIMULATED85_ARGS="--comparisons ${DATASET_DIR}/cmps_simulated85_multi.json \
 --sequences ${DATASET_DIR}/seqs_simulated85_multi.json"

# GENERATOR1_ARGS="--generatorCount 20000 \
#  --generatorSeqLen 20000 \
#  --generatorSimilarity 1"
# GENERATOR0_ARGS="--generatorCount 20000 \
#  --generatorSeqLen 20000 \
#  --generatorSimilarity 0"
# GENERATOR_ARGS="--generatorCount 20000 \
#  --generatorSeqLen 20000 \
#  --generatorSimilarity 0.85"

 declare -A dsmap
# dsmap[logan]="$LOGAN_ARGS"
dsmap[ecoli]="$ECOLI_ARGS"
# dsmap[simulated1]="$GENERATOR1_ARGS"
# dsmap[simulated0]="$GENERATOR0_ARGS"
dsmap[simulated85]="$SIMULATED85_ARGS"
# dsmap[elegans]="$ELEGANS_ARGS"
# dsmap[ecoli100]="$ECOLI100_ARGS"
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
  DSARGS="${dsmap[$DSNAME]}"
  DECOMPOSE=n
  XDROP=5
  DEVICES=1
  run
done

# for DSNAME in "${!dsmap[@]}"; do
# for DECOMPOSE in n y; do
# for XDROP in 5 10 15 20 50; do
# for DEVICES in 1 2 4 8 16 32; do
#   DSARGS="${dsmap[$DSNAME]}"
#   run
# done
# done
# done
# done