#!/bin/bash
BIN=./build/bin/ipusw
OUTPUT=/global/D1/projects/ipumer/datasets/results/ipu_synthetic_benchmarks
OVERWRITE=

for NUM_IPU in 1 2 4 8 16 32 64; do
# INPUT1=/global/D1/projects/ipumer/datasets/DNA-big-As.fasta
# INPUT2=/global/D1/projects/ipumer/datasets/DNA-big-Bs.fasta
# name=dna_large
# output_log=${OUTPUT}/${name}.log
# if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
# 	echo "Running ${name}"
# 	${BIN} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo roundrobin -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
# fi
# 
# INPUT1=/global/D1/projects/ipumer/datasets/DNA-big-As.fasta
# INPUT2=/global/D1/projects/ipumer/datasets/DNA-big-Bs.fasta
# name=dna_large_greedy
# output_log=${OUTPUT}/${name}.log
# if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
# 	echo "Running ${name}"
# 	${BIN} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo greedy -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
# fi

INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_1_200_ref.fasta
INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_1_200_qer.fasta
name=ipu${NUM_IPU}_dna_1_200
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo roundrobin -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_150_ref.fasta
INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_150_qer.fasta
name=ipu${NUM_IPU}_dna_2_150
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo roundrobin -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_200_ref.fasta
INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_200_qer.fasta
name=ipu${NUM_IPU}_dna_2_200
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo roundrobin -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_250_ref.fasta
INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_250_qer.fasta
name=ipu${NUM_IPU}_dna_2_250
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo roundrobin -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_200_ref.txt
INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_200_que.txt
name=ipu${NUM_IPU}_protein_200
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo roundrobin -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_200_ref.txt
INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_200_que.txt
name=ipu${NUM_IPU}_protein_200_greedy
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo greedy -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_400_ref.txt
INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_400_que.txt
name=ipu${NUM_IPU}_protein_400
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo roundrobin -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_400_ref.txt
INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_400_que.txt
name=ipu${NUM_IPU}_protein_400_greedy
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo greedy -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_600_ref.txt
INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_600_que.txt
name=ipu${NUM_IPU}_protein_600
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo roundrobin -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_600_ref.txt
INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_600_que.txt
name=ipu${NUM_IPU}_protein_600_greedy
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo greedy -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/compare/base/PROTEIN-As.txt
INPUT2=/global/D1/projects/ipumer/datasets/compare/base/PROTEIN-Bs.txt
name=ipu${NUM_IPU}_protein_unfiltered
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo roundrobin -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi

INPUT1=/global/D1/projects/ipumer/datasets/compare/base/PROTEIN-As.txt
INPUT2=/global/D1/projects/ipumer/datasets/compare/base/PROTEIN-Bs.txt
name=ipu${NUM_IPU}_protein_unfiltered_greedy
output_log=${OUTPUT}/${name}.log
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} --numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo greedy -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi
done