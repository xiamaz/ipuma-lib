#!/bin/bash
RUNDIR="scripts/seqan"
COMMIT="eed206e113f2aa1a2ef7a2441f19e56ecafa23fc"
ALIGN_BENCH_URL="https://github.com/rrahn/align_bench/archive/$COMMIT.zip"
ALIGN_BENCH_DIR="$RUNDIR/align_bench"
BIN="${ALIGN_BENCH_DIR}/build/bin/align_bench_par"
OUTPUT_PATH="${RUNDIR}/results"

mkdir -p $OUTPUT_PATH

# THREADS=(1 2 4 8 16 32 64 128)
THREADS=(48 96)
NUMACTL=(0 0-1)

if [ ! -f $BIN ]; then
	mkdir -p $RUNDIR

	if [ ! -d $ALIGN_BENCH_DIR ]; then
		git clone --recurse-submodules https://github.com/rrahn/align_bench.git $ALIGN_BENCH_DIR
	fi

	cd $ALIGN_BENCH_DIR
	git checkout $COMMIT
	mkdir -p build
	cd build
	cmake -DCMAKE_BUILD_TYPE=Release -DSEQAN_ARCH_SSE4=ON -DSEQAN_ARCH_AVX2=ON -DUSE_UME_SIMD=OFF -DSEQAN_ARCH_AVX512_CNL=ON ../
	make -j `nproc`
fi

run() {
	name=seqan${NUM_THREADS}_${dsname}_numa${NUMA}
	output_log=${OUTPUT_PATH}/${name}.log
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
	config="--threads ${NUM_THREADS} -v -d local -i 32 -a dna -r 5"
# DNA experiments
	INPUT1=./download/DNA-big-As.fasta
	INPUT2=./download/DNA-big-Bs.fasta
	dsname=dna_large
	run

	INPUT1=./download/DNA_2_150_ref.fasta
	INPUT2=./download/DNA_2_150_qer.fasta
	dsname=dna_2_150
	run

	INPUT1=./download/DNA_2_200_ref.fasta
	INPUT2=./download/DNA_2_200_qer.fasta
	dsname=dna_2_200
	run

	INPUT1=./download/DNA_2_250_ref.fasta
	INPUT2=./download/DNA_2_250_qer.fasta
	dsname=dna_2_250
	run

# Protein experiments
	config="--threads ${NUM_THREADS} -v -d local -i 32 -a aa -r 5"

	INPUT2=./download/PROTEIN_200_ref.fasta
	INPUT1=./download/PROTEIN_200_que.fasta
	dsname=protein_200
	run

	INPUT2=./download/PROTEIN_400_ref.fasta
	INPUT1=./download/PROTEIN_400_que.fasta
	dsname=protein_400
	run

	INPUT2=./download/PROTEIN_600_ref.fasta
	INPUT1=./download/PROTEIN_600_que.fasta
	dsname=protein_600
	run

	INPUT1=./download/PROTEIN-longer_que.fasta
	INPUT2=./download/PROTEIN-longer_ref.fasta
	dsname=protein_unfiltered
	run

	INPUT1=./download/As_new.fasta
	INPUT2=./download/Bs_new.fasta
	dsname=protein_full
	run
done

# export the generated data
logfiles=( $OUTPUT_PATH/*.log )
sed -n 's/Policy/numa,dataset,&/p' ${logfiles[1]}
for logfile in $OUTPUT_PATH/*.log; do
	logname=`basename $logfile`
	logstem="${logname%.*}"
	numa=${logstem##*_}
	dataset=${logstem#*_}
	dataset=${dataset%_*}
	if [ $numa = "numa0" ]; then
		numa=1
	elif [ $numa = "numa0-1" ]; then
		numa=2
	fi
	sed -n "s/parallel/$numa,$dataset,&/p" $logfile
done
