#!/bin/bash
BIN=./build/bin/ipusw

OUTPUT_BASEDIR="${1:-./logs/testrun}"

# overwrite existing output directories
OVERWRITE=
# dry run without execution
PRINTOUT=
# unset CONTAINER if script should NOT be run in docker
CONTAINER=


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
        local repnum="${1}"
	OUTPUT=${OUTPUT_BASEDIR}_rep${repnum}
	mkdir -p ${OUTPUT}
	in1=$(realpath ${INPUT1})
	in2=$(realpath ${INPUT2})
	name=`date -Ins`
	output_log=${OUTPUT}/${name}"${3}".log
	if [ ${PRINTOUT} ]; then
		echo "Print run command for ${name}: ${BIN} ${config} -- ${in1} ${in2}"
	elif [ $OVERWRITE ] || [ ! -f ${output_log} ]; then
		echo "Running ${name}"
		echo "CMD: $(gen_bin) ${config} -- ${in1} ${in2} 2>&1 | tee ${output_log}"
		nice -n18 $(gen_bin) ${2} -- ${in1} ${in2} 2>&1 | tee ${output_log}
	fi
}

REPMAX=5

for REPNUM in `seq $REPMAX`; do
	# DNA experiments
        config="
        --numDevices 1 \
        --numThreads 2 \
        --tilesUsed 1472
        --vtype multibandxdrop \
        --forwardOnly \
        --maxBatches $((6*2)) \
        --bufsize 173000 \
        --fillAlgo fillfirst \
        --maxAB 17920"
        echo $config
	# config="--numDevices ${NUM_IPU} --numThreads ${NUM_THREADS} --tilesUsed 1472 --vtype multiasm --forwardOnly --maxBatches ${MAX_BATCHES} --bufsize ${BUFSIZE} --fillAlgo ${fillAlgo} ${DUPLICATE_DS}"
	# config="--numDevices ${NUM_IPU} --numThreads ${NUM_THREADS} --tilesUsed 1472 --vtype multiasm --maxBatches ${MAX_BATCHES} --bufsize ${BUFSIZE} --fillAlgo ${fillAlgo} ${DUPLICATE_DS}"
	INPUT1=/global/D1/projects/ipumer/xdrop_inputs/elbaAs.txt
	INPUT2=/global/D1/projects/ipumer/xdrop_inputs/elbaBs.txt
	run $REPNUM "${config}" elba 
done
