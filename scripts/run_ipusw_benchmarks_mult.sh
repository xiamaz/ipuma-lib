#!/bin/bash
BIN=./build/bin/ipusw
OUTPUT=/global/D1/projects/ipumer/datasets/results/ipu_synthetic_benchmarks
OVERWRITE=
PRINTOUT=

mkdir -p ${OUTPUT}

THREADS=(1 2 4 8 16 32 64)
# THREADS=(1 2 4 8 16 32)

run() {
name=ipu${NUM_IPU}_${dsname}_${fillAlgo}
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
BUFSIZE=400000
MAX_BATCHES=2000

for NUM_IPU in ${THREADS[@]}; do

for fillAlgo in roundrobin greedy; do
# DNA experiments
config="--numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches ${MAX_BATCHES} --bufsize ${BUFSIZE} --fillAlgo ${fillAlgo}"
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
config="--numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches ${MAX_BATCHES} --bufsize ${BUFSIZE} --fillAlgo ${fillAlgo}"

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
