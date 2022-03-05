#!/bin/bash
BIN=./build/bin/ipusw
OUTPUT=/global/D1/projects/ipumer/datasets/results/ipu_synthetic_benchmarks
OVERWRITE=

mkdir -p ${OUTPUT}

THREADS=(1 2 3 4 8 16 32 64)

runIpu() {
if [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
	echo "Running ${name}"
	${BIN} ${config} -- ${INPUT1} ${INPUT2} 2>&1 | tee ${output_log}
fi
}

for NUM_IPU in ${THREADS[@]}; do

for fillAlgo in roundrobin greedy; do
	INPUT1=/global/D1/projects/ipumer/datasets/compare/base/DNA-big-As.txt
	INPUT2=/global/D1/projects/ipumer/datasets/compare/base/DNA-big-Bs.txt
	name=ipu${NUM_IPU}_dna_large_${fillAlgo}
	config="--numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo ${fillAlgo}"
	output_log=${OUTPUT}/${name}.log
	runIpu

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_1_200_ref.fasta
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_1_200_qer.fasta
	config="--numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo ${fillAlgo}"
	name=ipu${NUM_IPU}_dna_1_200_${fillAlgo}
	output_log=${OUTPUT}/${name}.log
	runIpu

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_150_ref.fasta
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_150_qer.fasta
	config="--numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo ${fillAlgo}"
	name=ipu${NUM_IPU}_dna_2_150_${fillAlgo}
	output_log=${OUTPUT}/${name}.log
	runIpu

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_200_ref.fasta
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_200_qer.fasta
	config="--numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo ${fillAlgo}"
	name=ipu${NUM_IPU}_dna_2_200_${fillAlgo}
	output_log=${OUTPUT}/${name}.log
	runIpu

	INPUT1=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_250_ref.fasta
	INPUT2=/global/D1/projects/ipumer/datasets/compare/dna/DNA_2_250_qer.fasta
	config="--numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches 1000 --bufsize 300000 --fillAlgo ${fillAlgo}"
	name=ipu${NUM_IPU}_dna_2_250_${fillAlgo}
	output_log=${OUTPUT}/${name}.log
	runIpu

	INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_200_ref.txt
	INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_200_que.txt
	config="--numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo ${fillAlgo}"
	name=ipu${NUM_IPU}_protein_200_${fillAlgo}
	output_log=${OUTPUT}/${name}.log
	runIpu

	INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_400_ref.txt
	INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_400_que.txt
	config="--numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo ${fillAlgo}"
	name=ipu${NUM_IPU}_protein_400_${fillAlgo}
	output_log=${OUTPUT}/${name}.log
	runIpu

	INPUT1=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_600_ref.txt
	INPUT2=/global/D1/projects/ipumer/datasets/protein-txt/PROTEIN_600_que.txt
	config="--numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo ${fillAlgo}"
	name=ipu${NUM_IPU}_protein_600_${fillAlgo}
	output_log=${OUTPUT}/${name}.log
	runIpu

	INPUT1=/global/D1/projects/ipumer/datasets/compare/base/PROTEIN-As.txt
	INPUT2=/global/D1/projects/ipumer/datasets/compare/base/PROTEIN-Bs.txt
	config="--numDevices ${NUM_IPU} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxAB 2000 --datatype aa --similarity blosum50 --maxBatches 1000 --bufsize 300000 --fillAlgo ${fillAlgo}"
	name=ipu${NUM_IPU}_protein_unfiltered_${fillAlgo}
	output_log=${OUTPUT}/${name}.log
	runIpu
done

done