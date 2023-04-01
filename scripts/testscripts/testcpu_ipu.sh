#!/bin/bash
#gdb --args /home/zhaom/ipuma-lib/build/bin/cpusw \
/home/zhaom/ipuma-lib/build/bin/cpusw \
 --output /tmp/simplecorrectness/scores_ksw.json \
 --algo ipumacpu \
 --threads 1 \
 --generatorCount 100 \
 --generatorSeqLen 2500 \
 --generatorSimilarity 0.85
