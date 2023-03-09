#!/bin/bash

set -euo pipefail

/home/lukb/miniconda3/bin/cmake --build `pwd`/build --config Release --target ipuma-core


pushd `pwd`/build/test/
numactl --cpunodebind 0 ./ipuma-core |& tee run.log
gcups_inner=$(cat run.log | grep JOBLOG | grep -o '{.*' | jq .gcups_inner) 
cycles_inner=$(cat run.log | grep JOBLOG | grep -o '{.*' | jq .cycles_inner) 



total_cycles=$(echo "$cycles_inner" | paste -sd+ | bc)
echo "Hist:"
echo "$gcups_inner" | hist -x
echo "MCycles" $(( "$total_cycles" / 10 ** 7))

