#!/bin/bash
#SBATCH --partition=qot
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=5-00:00:00
#SBATCH --job-name=correl2e
#SBATCH --mail-user=ayounis@ur.rochester.edu
#SBATCH --mail-type=end

module unload gcc lapack openblas hdf5
module load gcc/11.2.0/b1 lapack/3.9.0/b2 openblas/0.3.10/b1 hdf5/1.12.1/b1

export OMP_NUM_THREADS=12
export PWD=/scratch/ayounis/QuantaRay/bin

$PWD/correl2e part qm_data.h5
