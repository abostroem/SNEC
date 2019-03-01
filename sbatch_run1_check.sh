#! /bin/bash -l 
#SBATCH -t 8:00:00
#SBATCH --mail-type=begin 
#SBATCH --mail-type=end 
#SBATCH --mail-user=kabostroem@ucdavis.edu 
#SBATCH --job-name='Check model output'
#SBATCH -o slurm_%A.out # Standard output
#SBATCH -e slurm_%A.err # Standard error
srun ./run1_check.py