#!/bin/bash

# when using GPUs
#module load nano CUDA GCC OpenMPI MPI-settings/CUDA Doxygen texlive

# CPU-only
#module load nano GCC OpenMPI Doxygen texlive
module load user_spack/22.2.1
module unload intel-mpi/2019-intel
module load mpi.intel/2019.12_gcc
module load gcc/11

