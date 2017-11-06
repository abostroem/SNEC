#!/usr/bin/env python

import os
from sys import argv
from string import Template

script,SNname,MixNi,MM,Evalues = argv

PARFILE = Template("""\
#! /bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=16  # core count
#SBATCH -t 12:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=vsg@princeton.edu
#SBATCH --job-name=${SN}_${Nimixing}_${Mvalue}

BASEPATH=/scratch/gpfs/vsg/${SN}/mixing_${Nimixing}/M${Mvalue}
RUNS=${energies}

for r in $${RUNS[@]}; do
    cd $$BASEPATH/E_$$r
    time ./snec &> snec.out &
done

wait
""")

open("/scratch/gpfs/vsg/snec.cmd", "w").write(
    PARFILE.substitute(
         SN            = str(SNname),
         Nimixing      = str(MixNi),
         Mvalue        = str(MM),
         energies      = str(Evalues))
    )
