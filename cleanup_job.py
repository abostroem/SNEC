import glob
import os
import shutil
'''
Clean directory after a snec run by moving files into directories outside of the main one
'''
def clean_slurm_output():
    os.makedirs('slurm_out', exist_ok=True)
    flist = glob.glob('slurm*.out')
    for ifile in flist:
        shutil.move(ifile, os.path.join('slurm_out', ifile))
        
        
def clean_sbatch_scripts():
    dirlist = glob.glob('sbatch_out/*')
    dirnum = len(dirlist)+1

    os.makedirs(os.path.join('sbatch_out', 'run{}'.format(dirnum)))
    flist = glob.glob('snec*.sh')
    for ifile in flist:
        shutil.move(ifile, os.path.join('sbatch_out', 'run{}'.format(dirnum), ifile))
        
if __name__ == "__main__":
    clean_slurm_output()
    clean_sbatch_scripts()
    