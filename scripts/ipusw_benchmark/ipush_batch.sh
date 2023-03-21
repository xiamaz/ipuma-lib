#!/bin/bash
#SBATCH -p ipuq    # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH --sockets-per-node 2 # number of cores
#SBATCH --cores-per-socket 24
#SBATCH --ntasks-per-node=96
#gid  #SBATCH -o slurm.%N.%j.out # STDOUT
# #SBATCH -e slurm.%N.%j.err # STDERR
PROJECT_DIR=$(git rev-parse --show-toplevel)
bash "$PROJECT_DIR/scripts/ipusw_benchmark/ipush_bench.sh"