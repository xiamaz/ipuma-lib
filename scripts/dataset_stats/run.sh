#!/bin/bash

DSNAME=simulated85
python dataset_stats.py ../data/cmps_${DSNAME}_multi.json.gz ../data/seqs_${DSNAME}_multi.json.gz > dataset_stats/${DSNAME}_stats.json

DSNAME=ecoli
python dataset_stats.py ../data/cmps_${DSNAME}_multi_single.json.gz ../data/seqs_${DSNAME}_multi_single.json.gz > dataset_stats/${DSNAME}_stats.json

DSNAME=elba100
python dataset_stats.py ../data/cmps_${DSNAME}_multi.json.gz ../data/seqs_${DSNAME}_multi.json.gz > dataset_stats/${DSNAME}_stats.json

DSNAME=celegans
python dataset_stats.py ../data/cmps_${DSNAME}_multi.json.gz ../data/seqs_${DSNAME}_multi.json.gz > dataset_stats/${DSNAME}_stats.json
