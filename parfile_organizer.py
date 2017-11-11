#!/usr/bin/env python

import os
from sys import argv
from string import Template

script,basename,EE,mixing,Ni_mass = argv

M_ZAMS = ["9.0","9.5","10.0","10.5","11.0","11.5","12.0","12.5","13.0","13.5","14.0","14.5","15.0","15.5","16.0","16.5","17.0","17.5","18.0","18.5","19.0","19.5","20.0","20.5","21.0","21.5","22.0","22.5","23.0","23.5","24.0","24.5","25.0","25.5","26.0","26.5","27.0","27.5","28.0","28.5","29.0","29.5","30.0"]

M_ex = [1.4,1.4,1.423,1.483,1.404,1.551,1.478,1.568,1.614,1.615,1.652,1.688,1.886,1.916,1.504,1.518,1.530,1.911,1.905,1.810,1.646,1.773,1.804,1.536,1.480,1.569,1.817,2.061,2.121,2.134,2.093,2.019,1.930,1.930,1.930,1.930,1.930,1.930,1.930,1.930,1.930,1.930,1.930]

#M_ex = [1.4,1.4,1.423,1.483,1.404,1.551,1.478,1.568,1.614,1.615,1.652,1.688,1.886,1.916,1.504,1.518,1.530,1.911,1.905,1.810,1.646,1.773,1.804,1.536,1.480,1.569,2.3,2.4,2.4,2.4,2.4,2.4,2.4]

PARFILE = Template("""\
#____________LAUNCH_____________

outdir              = "Data"

#___________PROFILE_____________

profile_name 		= $profile_short

comp_profile_name	= $profile_iso_dat

#__________EXPLOSION_____________

initial_data 		= "Thermal_Bomb"

#Options:
#"Piston_Explosion"
#"Thermal_Bomb"

piston_vel          = 5.0d9
piston_tstart       = 0.0d0
piston_tend         = 1.0d-2

final_energy        = $fin_energy
bomb_tstart         = 0.0d0
bomb_tend           = $bomb_time
bomb_mass_spread    = $bomb_mass_size #(in solar mass)
bomb_start_point    = 1

#_____________GRID_______________

imax         = $resolution

gridding = "from_file_by_mass"

#Options:
#"uniform_by_mass"
#"from_file_by_mass"

mass_excision = 1
mass_excised = $excised_mass #in solar mass, provided mass_excision = 1

#___________EVOLUTION_____________

radiation = 1
eoskey = $eos

#Options:
#1 - ideal eos
#2 - Paczynski

helm_table_name = "src/helmholtz_eos/helm_table.dat"

Ni_switch = $Ni_switch_par
Ni_mass = $Ni_total_mass 			#(in solar mass)
Ni_boundary_mass = $Ni_mixing          #(in solar mass, here carefull with the excised mass)
                                       #(attention - smoothing is going to change it, if applied)
Ni_period = 1.0d4

saha_ncomps = 3

boxcar_smoothing = $boxcar

opacity_floor_envelope = 0.01d0
opacity_floor_core     = 0.24d0

#___________TIMING_______________

ntmax               = 10000000000000

tend                = $endtime

dtout               = 1.7d5
dtout_scalar        = 1.7d4
dtout_check         = 1.7d5

ntout               = -1
ntout_scalar        = -1
ntout_check         = -1

ntinfo              = 1000

dtmin               = 1.0d-10
dtmax               = 3.0d2

#____________TEST_________________

sedov 		    = 0
""")

shortfile = [f for f in os.listdir(basename + "/profiles") if (f.startswith("s") and f.endswith(".short"))][0]
isofile = [f for f in os.listdir(basename + "/profiles") if (f.startswith("s") and f.endswith(".iso.dat"))][0]
M_ZAMS_string = shortfile[1:-6]
open(basename + "/parameters", "w").write(
    PARFILE.substitute(
        profile_short      = "\""+"profiles/"+shortfile+"\"",
        profile_iso_dat      = "\""+"profiles/"+isofile+"\"",
        fin_energy      = str(EE)+'d51',
        bomb_time       = '1.0d0',
        bomb_mass_size      = '0.02d0',
        resolution      = '1000',
        excised_mass  = str(M_ex[M_ZAMS.index(M_ZAMS_string)]),
        Ni_switch_par = '1',
        Ni_total_mass = str(Ni_mass),
        Ni_mixing  = str(mixing),
        boxcar = '1',
        eos = '2',
        endtime = '2.592d7')
        )
