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
#SBATCH --mail-user=kabostroem@ucdavis.edu
#SBATCH --job-name=${SN}_${Nimixing}_${Mvalue}

BASEPATH=/home/bostroem/SNEC/${SN}/mixing_${Nimixing}/M${Mvalue}
RUNS=${energies}

for r in $${RUNS[@]}; do
    cd $$BASEPATH/E_$$r
    time ./snec &> snec.out &
done

wait
""")

open("/home/bostroem/SNEC/snec.cmd", "w").write(
    PARFILE.substitute(
         SN            = str(SNname),
         Nimixing      = str(MixNi),
         Mvalue        = str(MM),
         energies      = str(Evalues))
    )
