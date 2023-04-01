#!/bin/bash
#gdb --args /home/zhaom/ipuma-lib/build/bin/cpusw \
/home/zhaom/ipuma-lib/build/bin/cpusw \
 --output /tmp/simplecorrectness/scores_ksw.json \
 --algo genometools \
 --threads 1 \
 --xDrop 5 \
 --generatorCount 100 \
 --generatorSeqLen 20000 \
 --generatorSimilarity 0.6