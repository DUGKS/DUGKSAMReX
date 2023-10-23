#!/bin/bash
#SBATCH -J DUGKSAMReX
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -p amd_512
#. ./loadmod.sh
#srun ./DUGKS.3d.intel.MPI.ex ./Benchmark/benchmark.ZalesakDisk
# case=FlatInterface
# case=QuiescentDrop
# case=SpinodalDecomposition
# case=DropletCoalescence
# case=BasicCase
case=UnsteadyTGV
DIM=2
N=4
nohup mpiexec -n $N ./DUGKS.${DIM}d.gnu.MPI.ex ./Benchmark/${case}.benchmark &