import glob
import os
import sys
import shutil

import numpy as np
import yaml

import profile_manipulation


def get_input_parameters(filename='user_def.param'):
    '''
    Read user supplied file to get initial setup parameters
    Input:
        filename: name of user supplied file with inputs
    Output:
        parameters: dictionary of parameters
    '''
    ofile = open('user_def.param', 'r')
    parameters = yaml.load(ofile)

    for iparam in parameters['expand_params']:
        parameters[iparam] = np.arange(*parameters[iparam])
    return parameters
    
def create_directory(i_ni_mass, i_ni_mix, imass, ienergy, idensity, iradius, overwrite=False):
    '''
    Create directory structure, one directory for each model run
    '''
    path = os.path.join('snec_models', 
                        'Ni_mass_{:1.2f}'.format(i_ni_mass),
                        'Ni_mixing_{:1.1f}'.format(i_ni_mix),
                        'M{:2.1f}'.format(imass), #TODO - figure out formating so if string is single digit, prepended with a 0
                        'E_{:1.3f}'.format(ienergy),
                        'K_{:2.1f}'.format(idensity), #TODO - figure out formating so if string is single digit, prepended with a 0
                        'R_{}'.format(int(iradius)))
    if os.path.exists(path) and overwrite is False:
        overwrite=input('overwrite existing path (y, n): {} '.format(path))
        if overwrite == 'y':
            overwrite = True
        else:
            sys.exit()
    if (os.path.exists(path) and overwrite is True):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)
    os.makedirs(os.path.join(path, 'profiles'), exist_ok=True)
    os.makedirs(os.path.join(path, 'Data'), exist_ok=True)
    return path


    
def copy_files(basepath, imass, idensity, iradius):
    profile_basename = os.path.join('sukhbold_profiles_wind',
                    'M{:2.1f}'.format(imass), #TODO - figure out formating so if string is single digit, prepended with a 0
                    'K_{:2.1f}'.format(idensity), #TODO - figure out formating so if string is single digit, prepended with a 0
                    'R_{}'.format(int(iradius)))
    profile_list = glob.glob(os.path.join(profile_basename, 's*')) #get the .iso and .dat files
    for ifile in profile_list:
        shutil.copyfile(ifile, os.path.join(basepath, 'profiles', os.path.basename(ifile)))
    shutil.copytree('tables', os.path.join(basepath, 'tables'))


                
parameters = get_input_parameters()
profile_manipulation.add_wind(parameters)
for i_ni_mass in parameters['ni_mass']:
    for i_ni_mix in parameters['ni_mass']:
        for imass in parameters['mass']:
            for ienergy in parameters['explosion_energy']:
                for idensity in parameters['density_1D']:
                    for iradius in parameters['wind_extent']:
                        basepath=create_directory(i_ni_mass, i_ni_mix, imass, ienergy, idensity, iradius, overwrite=True)
                        copy_files(basepath, imass, idensity, iradius)
                    
                        

                                            
    
