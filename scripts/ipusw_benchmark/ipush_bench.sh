#!/bin/bash
set -euo pipefail
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

ECOLI_ARGS="--comparisons /global/D1/projects/ipumer/inputs_ab/cmps_ecoli_multi_single.json \
--sequences /global/D1/projects/ipumer/inputs_ab/seqs_ecoli_multi_single.json"

ECOLI100_ARGS="--comparisons /global/D1/projects/ipumer/inputs_ab/cmps_elba100_multi.json \
 --sequences /global/D1/projects/ipumer/inputs_ab/seqs_elba100_multi.json"

ELEGANS_ARGS="--comparisons /global/D1/projects/ipumer/inputs_ab/cmps_celegans_multi.json \
 --sequences /global/D1/projects/ipumer/inputs_ab/seqs_celegans_multi.json"

GENERATOR_ARGS="--generatorCount 10000 \
 --generatorSeqLen 20000 \
 --generatorSimilarity 0.85"

OUTDIR="$PROJECT_DIR/output/ipusw_benchmark/$(hostname)"
mkdir -p "$OUTDIR"

run() {
	RUNNAME="${DSNAME}_${ALGO}_x${XDROP}_decom${DECOMPOSE}"
	echo "Running $RUNNAME "
  if [[ $DECOMPOSE = "y" ]]; then
    EXTRAARGS="--decomposeMulticomparisons"
  else
    EXTRAARGS=""
  fi
	if [[ $NUMACTL = "y" ]]; then
   numactl -N 1 -m 1 "$BIN" \
     --numVertices 1472  --fillAlgo greedy --seedLength 17 \
     --maxComparisonsPerVertex 400 \
     --bandPercentageXDrop 0.45 \
     --maxSequenceLength 20000 \
     --xDrop "$XDROP" \
     --vtype xdroprestrictedseedextend \
   	  --output "$OUTDIR/${RUNNAME}_scores.json" \
   	$DSARGS $EXTRAARGS |& tee "$OUTDIR/${RUNNAME}.log"
	else
   numactl -N 1 -m 1 "$BIN" \
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
DECOMPOSE=y
DSNAME=ecoli
DSARGS=$ECOLI_ARGS
XDROP=20
run

XDROP=10
run

XDROP=15
run

XDROP=5
run

DSNAME=ecoli100
DSARGS=$ECOLI100_ARGS
XDROP=20
run

XDROP=10
run

XDROP=15
run

XDROP=5
run

DSNAME=elegans
DSARGS=$ELEGANS_ARGS
XDROP=20
run

XDROP=10
run

XDROP=15
run

XDROP=5
run

DSNAME=simulated85
DSARGS=$GENERATOR_ARGS
XDROP=20
run
XDROP=15
run
XDROP=10
run
XDROP=5
run

DECOMPOSE=n
DSNAME=ecoli
DSARGS=$ECOLI_ARGS
XDROP=20
run

XDROP=10
run

XDROP=15
run

XDROP=5
run

DSNAME=ecoli100
DSARGS=$ECOLI100_ARGS
XDROP=20
run
XDROP=10
run

XDROP=15
run

XDROP=5
run

DSNAME=elegans
DSARGS=$ELEGANS_ARGS
XDROP=20
run
XDROP=10
run

XDROP=15
run

XDROP=5
run

DSNAME=simulated85
DSARGS=$GENERATOR_ARGS
XDROP=20
run
XDROP=15
run
XDROP=10
run
XDROP=5
run