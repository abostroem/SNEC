#!/usr/bin/env python

import os
from sys import argv
from string import Template
import stat

if len(argv) == 4:
    script,SNname,KK,Rvalues = argv
    start_indx = 0
elif len(argv) == 5:
    script,SNname,KK,Rvalues,start_indx = argv
    start_indx = int(start_indx)
else:
    print('Incorrect number of arguments, expecting SNname, KK, Rvalues, and possibly a starting index')
    sys.exit()

BASEPATH = '/home/bostroem/SNEC/'
#BASEPATH = '/Users/bostroem/Desktop/research/not_my_code/SNEC-1.01'
ofile = open('snec_master_wind.cmd', 'w')
ofile.write('\
#!/bin/bash -l \n\
#SBATCH -N 1   # node count \n\
#SBATCH -t 12:00:00 \n\
#SBATCH --mail-type=begin \n\
#SBATCH --mail-type=end \n\
#SBATCH --mail-user=kabostroem@ucdavis.edu \n\
#SBATCH --job-name={}_{} \n\
srun snec_wind${}SLURM_ARRAY_TASK_ID{}.sh \n'.format(SNname, KK,'{', '}'))
ofile.close()
os.chmod('snec_master_wind.cmd', stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | stat.S_IWUSR | stat.S_IRUSR)

r_values = Rvalues.strip('(').strip(')').split(' ')

for indx, iradius in enumerate(r_values):
    ofile = open('snec_wind{}.sh'.format(indx+1+start_indx), 'w')
    radius_dir =  'R_{}'.format(iradius.strip("'"))
    k_dir = 'K{}'.format(KK)
    ofile.write("#!/bin/bash -l \n")
    sndir = '{}_wind'.format(SNname)
    ofile.write('cd {} \n'.format(os.path.join(BASEPATH, sndir, k_dir, radius_dir)))
    ofile.write('./snec &>snec.out \n')
    ofile.close()
    os.chmod('snec_wind{}.sh'.format(indx+1), stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | stat.S_IWUSR | stat.S_IRUSR)



