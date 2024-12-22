#!/bin/bash
#SBATCH --partition=qot
#SBATCH --cpus-per-task=10
#SBATCH --mem=15G
#SBATCH --time=5-00:00:00
#SBATCH --job-name=qray
#SBATCH --mail-user=ayounis@ur.rochester.edu
#SBATCH --mail-type=end

module unload gcc openmpi hdf5
module load gcc/11.2.0/b1 openmpi/2.1.1/b1 hdf5/1.12.1/b1

export OMP_NUM_THREADS=10
export PWD=/scratch/ayounis/QuantaRay/bin

$PWD/qrayQ input.deck
