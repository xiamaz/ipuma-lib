#!/bin/bash
OUTDIR=output/logan_correctness_kswscoreonly_defaultvals
mkdir -p $OUTDIR

dataset=logan
dsargs="--seqsH /home/zhaom/ipuma-lib/output/logan_100k/seqH_100k.txt \
 --seqsV /home/zhaom/ipuma-lib/output/logan_100k/seqV_100k.txt \
 --seedsH1 /home/zhaom/ipuma-lib/output/logan_100k/seedH_100k.txt \
 --seedsV1 /home/zhaom/ipuma-lib/output/logan_100k/seedV_100k.txt"

run() {
	runname=cpu_${algo}_${dataset}_x${xdrop}
	scorefile="$OUTDIR/scores_$runname.json"
	if [[ ! -f $scorefile ]]; then
		echo "Running $runname"
		/home/zhaom/ipuma-lib/build/bin/cpusw \
		 --output $OUTDIR/scores_$runname.json \
		 --algo $algo \
		 --xDrop $xdrop \
		 --matchValue 2 \
		 --mismatchValue -4 \
		 --gapInit -4 \
		 --gapExtend -2 \
		 --threads 32 \
		 $dsargs |& tee $OUTDIR/$runname.log
	fi
}

# algo=ipumacpu
# xdrop=10
# run
# xdrop=50
# run
# xdrop=500
# run

algo=ksw2
xdrop=10
run
# xdrop=50
# run
# xdrop=500
# run

 dataset=sim85
 dsargs="--generatorCount 100 \
 --generatorSeqLen 2500 \
 --generatorSimilarity 0.85"

algo=ipumacpu
xdrop=10
run
xdrop=50
run
xdrop=500
run

# algo=ksw2
# xdrop=10
# run
# xdrop=50
# run
# xdrop=500
# run