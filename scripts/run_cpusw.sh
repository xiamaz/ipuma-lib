#!/bin/bash
BIN=./build/bin/cpusw
OUTPUT=/global/D1/projects/ipumer/datasets/results/cpu_synthetic_benchmarks_paper_final_withoutpreprocess_outer_includeprofile
OVERWRITE=
PRINTOUT=

mkdir -p ${OUTPUT}

# THREADS=(1 2 4 8 16 32 64 128)
THREADS=(48 96)
NUMACTL=(0 0-1)

run() {
	name=ssw${NUM_THREADS}_${dsname}_numa${NUMA}
	output_log=${OUTPUT}/${name}.log
	if [ ${PRINTOUT} ]; then
		echo "Print run command for ${name}"
		echo "CMD: numactl -N ${NUMA} -m ${NUMA} -- ${BIN} ${config} -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}"
	elif [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
		echo "Running ${name}"
		echo "CMD: numactl -N ${NUMA} -m ${NUMA} -- ${BIN} ${config} -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}"
		numactl -N ${NUMA} -m ${NUMA} -- ${BIN} ${config} -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
	fi
}

for NUMA in ${NUMACTL[@]}; do
	NUM_THREADS=48
	if [ $NUMA = "0-1" ]; then
		NUM_THREADS=96
	fi

	config="--threads ${NUM_THREADS}"
# DNA experiments
	INPUT1=/global/D1/projects/ipumer/datasets/compare/base/DNA-big-As.txt
	INPUT2=/global/D1/projects/ipumer/datasets/compare/base/DNA-big-Bs.txt
	dsname=dna_large
	run

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/8x/DNA_2_150_ref.fasta.2x.4x.8x.txt
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/8x/DNA_2_150_qer.fasta.2x.4x.8x.txt
	dsname=dna_2_150_8x
	run

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/128x/DNA_2_200_ref.fasta.2x.4x.8x.16x.32x.64x.128x.txt
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/128x/DNA_2_200_qer.fasta.2x.4x.8x.16x.32x.64x.128x.txt
	dsname=dna_2_200_128x
	run

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/32x/DNA_2_250_ref.fasta.2x.4x.8x.16x.32x.txt
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/32x/DNA_2_250_qer.fasta.2x.4x.8x.16x.32x.txt
	dsname=dna_2_250_32x
	run

# Protein experiments
	config="--threads ${NUM_THREADS} --datatype aa"

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

	# INPUT1=/global/D1/projects/ipumer/datasets/compare/base/PROTEIN-As.txt
	# INPUT2=/global/D1/projects/ipumer/datasets/compare/base/PROTEIN-Bs.txt
	# dsname=protein_unfiltered
	# run

 	INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN-longer_que.txt
 	INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN-longer_ref.txt
 	dsname=protein_unfiltered
 	run

 	INPUT1=/global/D1/projects/ipumer/As.txt
 	INPUT2=/global/D1/projects/ipumer/Bs.txt
 	dsname=protein_full
 	run
done
