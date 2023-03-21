#!/bin/bash
#SBATCH -p habanaq    # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH --sockets-per-node 1 # number of cores
#SBATCH --cores-per-socket 36
#SBATCH --ntasks-per-node=128
#gid  #SBATCH -o slurm.%N.%j.out # STDOUT
# #SBATCH -e slurm.%N.%j.err # STDERR
PROJECT_DIR=$(git rev-parse --show-toplevel)
bash "$PROJECT_DIR/scripts/cpusw_benchmark/cpush_bench.sh" 72