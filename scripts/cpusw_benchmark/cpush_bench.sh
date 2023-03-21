#!/bin/bash
set -euo pipefail
PROJECT_DIR=$(git rev-parse --show-toplevel)
BUILD_DIR="$PROJECT_DIR/build_$(hostname)"
# --- Build the CPUSW Binary on the hostname system
mkdir -p "$BUILD_DIR"

pushd "$PROJECT_DIR/extern/genometools"
make cleanup
popd

pushd "$PROJECT_DIR/extern/libgaba"
make clean
popd

pushd "$BUILD_DIR"
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --target cpusw
popd
BIN="$PROJECT_DIR/build_$(hostname)/bin/cpusw"

THREADS=$1
NUMACTL=y

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

OUTDIR="$PROJECT_DIR/output/cpusw_benchmark/$(hostname)"
mkdir -p "$OUTDIR"

run() {
	RUNNAME="${DSNAME}_${ALGO}_x${XDROP}_t$THREADS"
	echo "Running $RUNNAME using $THREADS threads"
	if [[ $NUMACTL = "y" ]]; then
		numactl -N 1 -m 1 "$BIN" \
			--algo "$ALGO" \
			--threads "$THREADS" \
			--xDrop "$XDROP" \
			--output "$OUTDIR/${RUNNAME}_scores.json" \
			$DSARGS |& tee "$OUTDIR/${RUNNAME}.log"
	else
		"$BIN" \
			--algo "$ALGO" \
			--threads "$THREADS" \
			--xDrop "$XDROP" \
			--output "$OUTDIR/${RUNNAME}_scores.json" \
			$DSARGS |& tee "$OUTDIR/${RUNNAME}.log"
	fi
}

NUMACTL=y
DSNAME=ecoli
DSARGS=$ECOLI_ARGS

ALGO=seqan
XDROP=20
run
XDROP=15
run
XDROP=10
run
XDROP=5
run

ALGO=ksw2
XDROP=400
run
ALGO=ksw2
XDROP=20
run
ALGO=ksw2
XDROP=15
run
ALGO=ksw2
XDROP=10
run
ALGO=ksw2
XDROP=5
run

DSNAME=ecoli100
DSARGS=$ECOLI100_ARGS

ALGO=seqan
XDROP=15
run

ALGO=seqan
XDROP=5
run

ALGO=seqan
XDROP=20
run

ALGO=seqan
XDROP=10
run

ALGO=ksw2
XDROP=400
run

ALGO=ksw2
XDROP=20
run

ALGO=ksw2
XDROP=15
run

ALGO=ksw2
XDROP=10
run

ALGO=ksw2
XDROP=5
run

DSNAME=elegans
DSARGS=$ELEGANS_ARGS

ALGO=seqan
XDROP=20
run

ALGO=seqan
XDROP=10
run

ALGO=seqan
XDROP=15
run

ALGO=seqan
XDROP=5
run

ALGO=ksw2
XDROP=400
run

ALGO=ksw2
XDROP=15
run

ALGO=ksw2
XDROP=5
run

ALGO=ksw2
XDROP=10
run

ALGO=ksw2
XDROP=20
run

DSNAME=simulated85
DSARGS=$GENERATOR_ARGS

ALGO=seqan
XDROP=20
run

ALGO=seqan
XDROP=10
run

ALGO=seqan
XDROP=15
run

ALGO=seqan
XDROP=5
run

ALGO=ksw2
XDROP=400
run

ALGO=ksw2
XDROP=15
run

ALGO=ksw2
XDROP=5
run

ALGO=ksw2
XDROP=10
run

ALGO=ksw2
XDROP=20
run

# ALGO=ipumacpu
# XDROP=15
# run
# 
# ALGO=libgaba
# XDROP=15
# run
# 
# ALGO=genometools
# XDROP=15
# run