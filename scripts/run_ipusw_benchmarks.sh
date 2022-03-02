#!/bin/bash
BIN=./build/bin/ipusw
OUTPUT=/global/D1/projects/ipumer/datasets/results/ipu_synthetic_benchmarks
OVERWRITE=

INPUT1=/global/D1/projects/ipumer/datasets/DNA-big-As.fasta
INPUT2=/global/D1/projects/ipumer/datasets/DNA-big-Bs.fasta
name=dna_large
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo roundrobin -- ${INPUT1} ${INPUT2} | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/DNA-big-As.fasta
INPUT2=/global/D1/projects/ipumer/datasets/DNA-big-Bs.fasta
name=dna_large_greedy
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo greedy -- ${INPUT1} ${INPUT2} | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_1_200_ref.fasta
INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_1_200_qer.fasta
name=dna_1_200
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo roundrobin -- ${INPUT1} ${INPUT2} | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/compare/protein/PROTEIN_200_ref.fasta
INPUT2=/global/D1/projects/ipumer/datasets/compare/protein/PROTEIN_200_que.fasta
name=protein_200
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo roundrobin -- ${INPUT1} ${INPUT2} | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/compare/protein/PROTEIN_200_ref.fasta
INPUT2=/global/D1/projects/ipumer/datasets/compare/protein/PROTEIN_200_que.fasta
name=protein_200_greedy
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo greedy -- ${INPUT1} ${INPUT2} | tee ${output_log}
fi