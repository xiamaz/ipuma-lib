#!/bin/bash

set -euo pipefail

X=20

#--fillAlgo roundrobin greedy fillfirst
pushd build && ninja; popd

# ./scripts/lib.sh logs/ipu_elba$X

CMD="build/bin/ipusw --numDevices 1 --numThreads 2 --tilesUsed 1472"
$CMD    --vtype multibandxdrop \
        --forwardOnly \
        --maxBatches $((6*10)) \
        --bufsize 173000 \
        --fillAlgo fillfirst \
        --xDrop 100 \
        --bandPercentageXDrop 0.5 \
        --maxAB 17920  \
        /global/D1/projects/ipumer/xdrop_inputs/elbaAs.txt /global/D1/projects/ipumer/xdrop_inputs/elbaBs.txt \
        | tee elbaX${X}.log
        # /global/D1/projects/ipumer/xdrop_inputs/elbaAs.txt /global/D1/projects/ipumer/xdrop_inputs/elbaBs.txt \