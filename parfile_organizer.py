#!/usr/bin/env python

import os
import glob
import sys
from string import Template

def write_parfile(i_ni_mass, i_ni_mix, imass, ienergy, parameters, basepath):

    PARFILE = Template(
    """\
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

#helm_table_name = "src/helmholtz_eos/helm_table.dat"

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
    profile_flist = glob.glob(os.path.join(basepath, 'profiles', 's*'))
    shortfile = None
    isofile = None
    for ifile in profile_flist:
        if ifile.endswith('short'): 
            shortfile=os.path.basename(ifile)
        if ifile.endswith('dat'):
            isofile=os.path.basename(ifile)
    if (isofile is None) or (shortfile is None):
        print(profile_flist)
        import pdb; pdb.set_trace()

    with open(os.path.join(basepath,"parameters"), "w") as ofile:
        if imass not in parameters['excised_mass'].keys():
            print('Please add an entry to the excised_mass dictionary defined in user_def.param for mass {}'.format(imass))
            sys.exit()
        ofile.write(
            PARFILE.substitute(
                profile_short      =  '"{}"'.format(os.path.join("profiles",shortfile)),
                profile_iso_dat      = '"{}"'.format(os.path.join("profiles",isofile)),
                fin_energy      = str(ienergy)+'d51',
                bomb_time       = '{}d0'.format(parameters['bomb_time']),
                bomb_mass_size      = '{}d0'.format(parameters['bomb_mass_size']),
                resolution      = str(parameters['resolution']),
                excised_mass = str(parameters['excised_mass'][imass]),
                Ni_switch_par = str(parameters['Ni_switch_par']),
                Ni_total_mass = '{}'.format(i_ni_mass),
                Ni_mixing  = '{}'.format(i_ni_mix),
                boxcar = str(parameters['boxcar']),
                eos = str(parameters['eos']),
                endtime = '{:1.3E}'.format(parameters['endtime']).replace('E+0', 'd')
            )
        )
