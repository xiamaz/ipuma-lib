#!/bin/bash
#SBATCH -p hgx2q    # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH --sockets-per-node 1 # number of cores
#SBATCH --cores-per-socket 1
#SBATCH --ntasks-per-node=128
#gid  #SBATCH -o slurm.%N.%j.out # STDOUT
# #SBATCH -e slurm.%N.%j.err # STDERR

cmake --build `pwd`/build --config Release --target ipusw