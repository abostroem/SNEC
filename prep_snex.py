import argparse
import glob
import os
import stat
import sys
import shutil

import numpy as np
import yaml

import profile_manipulation
import parfile_organizer

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
                        'Ni_mass_{:1.4f}'.format(i_ni_mass),
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
    '''
    Copy relevant profiles, executable, and opacity tables to each directory
    '''
    profile_basename = os.path.join('sukhbold_profiles_wind',
                    'M{:2.1f}'.format(imass), #TODO - figure out formating so if string is single digit, prepended with a 0
                    'K_{:2.1f}'.format(idensity), #TODO - figure out formating so if string is single digit, prepended with a 0
                    'R_{}'.format(int(iradius)))
    profile_list = glob.glob(os.path.join(profile_basename, 's*')) #get the .iso and .dat files
    for ifile in profile_list:
        shutil.copyfile(ifile, os.path.join(basepath, 'profiles', os.path.basename(ifile)))
    shutil.copytree('tables', os.path.join(basepath, 'tables'))
    shutil.copyfile('snec', os.path.join(basepath, 'snec'))

def write_sbatch_executable(basepath, array_num, snec_dir):
    '''
    Create a file to execute snec for each array job
    '''
    ofile = open('snec{}.sh'.format(array_num), 'w')
    ofile.write('cd {}\n'.format(os.path.join(snec_dir, basepath)))
    ofile.write('.snec &>snec.out')
    ofile.close()
    os.chmod('snec{}.sh'.format(array_num), stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | stat.S_IWUSR | stat.S_IRUSR)

def write_sbatch_job(timeout):
    '''
    Creates script to be called by sbatch
    '''
    ofile = open('snec_master.sh', 'w')
    ofile.write(
    '''#! /bin/bash -l 
#SBATCH -t {}
#SBATCH --mail-type=begin 
#SBATCH --mail-type=end 
#SBATCH --mail-user=kabostroem@ucdavis.edu 
#SBATCH --job-name=M={}-{},E={}-{},K={}-{}, R={}-{}
srun ./snec${}SLURM_ARRAY_TASK_ID{}.sh 
   '''.format(timeout, parameters['mass'][0], parameters['mass'][-1],
                                            parameters['explosion_energy'][0], parameters['explosion_energy'][-1],
                                            parameters['density_1D'][0], parameters['density_1D'][-1],
                                            parameters['wind_extent'][0], parameters['wind_extent'][-1],
                                            '{', '}')
        )
    ofile.close()
    os.chmod('snec_master.sh', stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | stat.S_IWUSR | stat.S_IRUSR)
    
    
if __name__ == "__main__":      

    parser = argparse.ArgumentParser()
    parser.add_argument('timeout', type=int, nargs=1,
                        help='the maximum time for a job to run, this should be twice what you expect it to take')
    args = parser.parse_args()
    snec_dir = os.getcwd()     
    parameters = get_input_parameters()
    profile_manipulation.add_wind(parameters)
    array_num=0
    for indx_ni_mass, i_ni_mass in enumerate(parameters['ni_mass']):
        for indx_ni_mix, i_ni_mix in enumerate(parameters['ni_mix']):
            for indx_mass, imass in enumerate(parameters['mass']):
                for indx_energy, ienergy in enumerate(parameters['explosion_energy']):
                    for indx_density, idensity in enumerate(parameters['density_1D']):
                        for indx_radius, iradius in enumerate(parameters['wind_extent']):
                            basepath=create_directory(i_ni_mass, i_ni_mix, imass, ienergy, idensity, iradius, overwrite=True)
                            copy_files(basepath, imass, idensity, iradius)
                            parfile_organizer.write_parfile(i_ni_mass, i_ni_mix, imass, ienergy, parameters, basepath)
                            array_num +=1
                            write_sbatch_executable(basepath, array_num, snec_dir)
    write_sbatch_job(args.timeout)
    print('Ready to run, start your jobs with the command:')
    print('sbatch --partition high --array=1-{}'.format(array_num))
                    
                        

                                            
    
