#! /bin/bash -l 
#SBATCH -t 8:00:00
#SBATCH --mail-type=begin 
#SBATCH --mail-type=end 
#SBATCH --mail-user=kabostroem@ucdavis.edu 
#SBATCH --job-name='Check model output'
#SBATCH -o name_%A_%a.out # Standard output
#SBATCH -e name_%A_%a.err # Standard error
SRUN python run1_check.py