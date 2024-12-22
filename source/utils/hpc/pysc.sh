#!/bin/bash
#SBATCH --partition=qot
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=pysc

module unload anaconda3
module load anaconda3/2021.11

python3 pdb.py
