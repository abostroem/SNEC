#!/usr/bin/env python

import os
from sys import argv
from string import Template
import stat

BASEPATH = '/Users/bostroem/Desktop/research/not_my_code/SNEC-1.01'

script,SNname,MixNi,MM,Evalues = argv
ofile = open('snec_master.cmd', 'w')
ofile.write(' \
#! /bin/bash -l \n \
#SBATCH -N 1   # node count \n \
#SBATCH -t 12:00:00 \n \
#SBATCH --mail-type=begin \n \
#SBATCH --mail-type=end \n \
#SBATCH --mail-user=kabostroem@ucdavis.edu \n \
#SBATCH --job-name=${SN}_${Nimixing}_${Mvalue} \n \
srun ./snec${SLURM_ARRAY_TASK_ID}.sh \n'.format(SNname, MixNi, MM))
ofile.close()
os.chmod('snec_master.cmd', stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | stat.S_IWUSR | stat.S_IRUSR)



energies = Evalues.strip('(').strip(')').strip("'").split(' ')
for indx, ienergy in enumerate(energies):
    ofile = open('snec{}.sh'.format(indx+1), 'w')
    energy_dir =  'E_{}r'.format(ienergy)
    mixing_dir = 'mixing_{}'.format(MixNi)
    mass_dir = 'M{}'.format(MM)
    ofile.write('cd {} \n'.format(os.path.join(BASEPATH, SNname, mixing_dir, mass_dir, energy_dir)))
    ofile.write('./snec &>snec.out & \n')
    ofile.close()
    os.chmod('snec{}.sh'.format(indx+1), stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | stat.S_IWUSR | stat.S_IRUSR)



