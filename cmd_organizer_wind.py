#!/usr/bin/env python

import os
from sys import argv
from string import Template

script,SNname,KK,Rvalues = argv

PARFILE = Template("""\
#! /bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=16  # core count
#SBATCH -t 12:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=vsg@princeton.edu
#SBATCH --job-name=${SN}_${Kvalue}

BASEPATH=/scratch/gpfs/vsg/${SN}_wind/K${Kvalue}
RUNS=${radii}

for r in $${RUNS[@]}; do
    cd $$BASEPATH/R_$$r
    time ./snec &> snec.out &
done

wait
""")

open("/scratch/gpfs/vsg/snec_wind.cmd", "w").write(
    PARFILE.substitute(
         SN            = str(SNname),
         Kvalue        = str(KK),
         radii      = str(Rvalues))
    )
