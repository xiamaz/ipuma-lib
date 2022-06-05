#!/bin/bash
BIN=./build/bin/ipusw
OUTPUT_STEM=/global/D1/projects/ipumer/datasets/results/ipu_synthetic_benchmarks_paper_final_rebalance_container
OVERWRITE=
PRINTOUT=
CONTAINER=y

# DEVNUM=(1 2 4 8 16 32 64)
DEVNUM=(1 2 4 8 16 32 64)
REPMAX=5


gen_bin() {
	if [ ! ${CONTAINER} ]; then
		echo "${BIN}"
		return
	fi
	in1=$(realpath ${INPUT1})
	in2=$(realpath ${INPUT2})
	echo "gc-docker -- -e IPUOF_CONFIG_PATH=${IPUOF_CONFIG_PATH} --rm -it -v ${IPUOF_CONFIG_PATH}:${IPUOF_CONFIG_PATH} -v ${in1}:${in1} -v ${in2}:${in2} lukburchard/ipuma-lib:latest  /build/bin/ipusw"
}

run() {
	in1=$(realpath ${INPUT1})
	in2=$(realpath ${INPUT2})
	name=ipu${NUM_IPU}_th${NUM_THREADS}_${dsname}_${fillAlgo}_${DDUP}
	output_log=${OUTPUT}/${name}.log
	if [ ${PRINTOUT} ]; then
		echo "Print run command for ${name}: ${BIN} ${config} -- ${in1} ${in2}"
	elif [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
		echo "Running ${name}"
		echo "CMD: $(gen_bin) ${config} -- ${in1} ${in2} 2>&1 | tee ${output_log}"
		$(gen_bin) ${config} -- ${in1} ${in2} 2>&1 | tee ${output_log}
	fi
}

for REPNUM in `seq $REPMAX`; do
OUTPUT=${OUTPUT_STEM}_rep${REPNUM}
mkdir -p ${OUTPUT}
for NUM_IPU in "${DEVNUM[@]}"; do
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
# config="--numDevices ${NUM_IPU} --numThreads ${NUM_THREADS} --tilesUsed 1472 --vtype multiasm --maxBatches ${MAX_BATCHES} --bufsize ${BUFSIZE} --fillAlgo ${fillAlgo} ${DUPLICATE_DS}"
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
fillAlgo=greedy
BUFSIZE=170000
MAX_BATCHES=300
config="--numDevices ${NUM_IPU} --numThreads ${NUM_THREADS} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 1500 --datatype aa --similarity blosum50 --maxBatches ${MAX_BATCHES} --bufsize ${BUFSIZE} --fillAlgo ${fillAlgo} ${DUPLICATE_DS}"

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
done
done