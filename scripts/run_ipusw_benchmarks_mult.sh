#!/bin/bash
BIN=./build/bin/ipusw
OUTPUT_STEM=/global/D1/projects/ipumer/datasets/results/ipu_synthetic_benchmarks_paper_final
OVERWRITE=
PRINTOUT=

# DEVNUM=(1 2 4 8 16 32 64)
DEVNUM=(1 2 4 8 16 32 64)
REPMAX=5

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

for REPNUM in `seq $REPMAX`; do
OUTPUT=${OUTPUT_STEM}_rep${REPNUM}
remoteMemory=stream
mkdir -p ${OUTPUT}
for NUM_IPU in ${DEVNUM[@]}; do
for DDUP in yesdup nodup; do
if [ $DDUP = "yesdup" ]; then
	DUPLICATE_DS="--duplicateDatasets"
	if [ ${NUM_IPU} -lt 2 ]; then
		continue
	fi
else
	DUPLICATE_DS=""
fi
# NUM_THREADS=$((NUM_IPU * TH_FACTOR))
NUM_THREADS=`nproc`
# DNA experiments
fillAlgo=roundrobin
BUFSIZE=34000
MAX_BATCHES=180
config="--numDevices ${NUM_IPU} --numThreads ${NUM_THREADS} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches ${MAX_BATCHES} --bufsize ${BUFSIZE} --fillAlgo ${fillAlgo} ${DUPLICATE_DS}"
	INPUT1=/global/D1/projects/ipumer/datasets/compare/base/DNA-big-As.txt
	INPUT2=/global/D1/projects/ipumer/datasets/compare/base/DNA-big-Bs.txt
	dsname=dna_large
	run

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/8x/DNA_2_150_ref.fasta.2x.4x.8x.txt
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/8x/DNA_2_150_qer.fasta.2x.4x.8x.txt
	dsname=dna_2_150
	run

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/128x/DNA_2_200_ref.fasta.2x.4x.8x.16x.32x.64x.128x.txt
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/128x/DNA_2_200_qer.fasta.2x.4x.8x.16x.32x.64x.128x.txt
	dsname=dna_2_200
	run

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/32x/DNA_2_250_ref.fasta.2x.4x.8x.16x.32x.txt
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/32x/DNA_2_250_qer.fasta.2x.4x.8x.16x.32x.txt
	dsname=dna_2_250
	run

	# INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_1_200_ref.fasta
	# INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_1_200_qer.fasta
	# dsname=dna_1_200
	# run

	# INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_150_ref.fasta
	# INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_150_qer.fasta
	# dsname=dna_2_150
	# run

	# INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_200_ref.fasta
	# INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_200_qer.fasta
	# dsname=dna_2_200
	# run

	# INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_250_ref.fasta
	# INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_250_qer.fasta
	# dsname=dna_2_250
	# run


# Protein experiments
fillAlgo=greedy
BUFSIZE=170000
MAX_BATCHES=300
config="--numDevices ${NUM_IPU} --numThreads ${NUM_THREADS} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 1500 --datatype aa --similarity blosum50 --maxBatches ${MAX_BATCHES} --bufsize ${BUFSIZE} --fillAlgo ${fillAlgo} ${DUPLICATE_DS}"

	INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_200_ref.txt
	INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_200_que.txt
	dsname=protein_200
	run

	INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_400_ref.txt
	INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_400_que.txt
	dsname=protein_400
	run

	INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_600_ref.txt
	INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_600_que.txt
	dsname=protein_600
	run

	INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN-longer_que.txt
	INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN-longer_ref.txt
	dsname=protein_unfiltered
	run

	INPUT1=/global/D1/projects/ipumer/As.txt
	INPUT2=/global/D1/projects/ipumer/Bs.txt
	dsname=protein_full
	run
done
done
done
