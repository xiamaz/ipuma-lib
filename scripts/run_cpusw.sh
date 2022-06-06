#!/bin/bash
BIN=./build/bin/cpusw
RUNSET_NAME="${1:-cpuswrun}"
OUTPUT="./output/${RUNSET_NAME}"
# OUTPUT=/global/D1/projects/ipumer/datasets/results/cpu_synthetic_benchmarks_paper_final_withoutpreprocess_outer_includeprofile
OVERWRITE=
PRINTOUT=

mkdir -p ${OUTPUT}

# setting for two socket systems
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
	INPUT1=./download/DNA-big-As.txt
	INPUT2=./download/DNA-big-Bs.txt
	dsname=dna_large
	run

	INPUT1=./download/DNA_2_150_ref.txt
	INPUT2=./download/DNA_2_150_qer.txt
	dsname=dna_2_150
	run

	INPUT1=./download/DNA_2_200_ref.txt
	INPUT2=./download/DNA_2_200_qer.txt
	dsname=dna_2_200
	run

	INPUT1=./download/DNA_2_250_ref.txt
	INPUT2=./download/DNA_2_250_qer.txt
	dsname=dna_2_250
	run

# Protein experiments
	config="--threads ${NUM_THREADS} --datatype aa"

	INPUT2=./download/PROTEIN_200_ref.txt
	INPUT1=./download/PROTEIN_200_que.txt
	dsname=protein_200
	run

	INPUT2=./download/PROTEIN_400_ref.txt
	INPUT1=./download/PROTEIN_400_que.txt
	dsname=protein_400
	run

	INPUT2=./download/PROTEIN_600_ref.txt
	INPUT1=./download/PROTEIN_600_que.txt
	dsname=protein_600
	run

	INPUT1=./download/PROTEIN-longer_que.txt
	INPUT2=./download/PROTEIN-longer_ref.txt
	dsname=protein_unfiltered
	run

	INPUT1=./download/As_new.txt
	INPUT2=./download/Bs_new.txt
	dsname=protein_full
	run
done
