#!/bin/bash
#SBATCH --partition=qot
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=pysc

module unload anaconda3
module load anaconda3/2021.05

python3 pdb.py
