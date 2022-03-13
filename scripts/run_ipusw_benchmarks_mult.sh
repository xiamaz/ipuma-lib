#!/bin/bash
BIN=./build/bin/ipusw
OUTPUT=/global/D1/projects/ipumer/datasets/results/ipu_synthetic_benchmarks
OVERWRITE=
PRINTOUT=

mkdir -p ${OUTPUT}

DEVNUM=(1 2 8 32 64)
TH_FACTORS=(1 5 10)

run() {
name=ipu${NUM_IPU}_th${NUM_THREADS}_${dsname}_${fillAlgo}_${DDUP}
output_log=${OUTPUT}/${name}.log
if [ ${PRINTOUT} ]; then
	echo "Print run command for ${name}: ${BIN} ${config} -- ${INPUT1} ${INPUT2}"
else if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	echo "CMD: ${BIN} ${config} -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}"
	${BIN} ${config} -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi
fi
}
remoteMemory=stream
for NUM_IPU in ${DEVNUM[@]}; do
for TH_FACTOR in ${TH_FACTORS[@]}; do
for DDUP in yesdup nodup; do
if [ $DDUP = "yesdup" ]; then
	DUPLICATE_DS="--duplicateDatasets"
	if [ $NUM_IPU -le 2 ]; then
		continue
	fi
else
	DUPLICATE_DS=""
fi
NUM_THREADS=$((NUM_IPU * TH_FACTOR))
# DNA experiments
fillAlgo=roundrobin
BUFSIZE=34000
MAX_BATCHES=100
config="--numDevices ${NUM_IPU} --numThreads ${NUM_THREADS} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches ${MAX_BATCHES} --bufsize ${BUFSIZE} --fillAlgo ${fillAlgo} ${DUPLICATE_DS}"
	INPUT1=/global/D1/projects/ipumer/datasets/compare/base/DNA-big-As.txt
	INPUT2=/global/D1/projects/ipumer/datasets/compare/base/DNA-big-Bs.txt
	dsname=dna_large
	run

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_1_200_ref.fasta
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_1_200_qer.fasta
	dsname=dna_1_200
	run

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_150_ref.fasta
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_150_qer.fasta
	dsname=dna_2_150
	run

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_200_ref.fasta
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_200_qer.fasta
	dsname=dna_2_200
	run

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_250_ref.fasta
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_250_qer.fasta
	dsname=dna_2_250
	run

# Protein experiments
fillAlgo=greedy
BUFSIZE=140000
MAX_BATCHES=220
config="--numDevices ${NUM_IPU} --numThreads ${NUM_THREADS} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches ${MAX_BATCHES} --bufsize ${BUFSIZE} --fillAlgo ${fillAlgo} ${DUPLICATE_DS}"

	INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_200_ref.txt
	INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_200_que.txt
	dsname=protein_200
	run

	INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_400_ref.txt
	INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_400_que.txt
	dsname=protein_400
	run

	INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_600_ref.txt
	INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_600_que.txt
	dsname=protein_600
	run

	INPUT1=/global/D1/projects/ipumer/datasets/compare/base/PROTEIN-As.txt
	INPUT2=/global/D1/projects/ipumer/datasets/compare/base/PROTEIN-Bs.txt
	dsname=protein_unfiltered
	run
done
done
done