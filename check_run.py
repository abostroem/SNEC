import os
import subprocess
import tarfile
import shutil
import time

# I use tar and zip interchangably

def get_last_line(filename, handle=False):
    if handle is False:
        f = open(filename, 'rb')
    else:
        f = filename
    
    first = f.readline()        # Read the first line.
    f.seek(-2, os.SEEK_END)     # Jump to the second last byte.
    while f.read(1) != b"\n":   # Until EOL is found...
        f.seek(-2, os.SEEK_CUR) # ...jump back the read byte plus one more.
    last = f.readline().decode('ascii')        # Read last line. and convert to string
    return float(last.split()[0])
        
        
def check_run(parameters, base_model_dir):
    start_time = time.time()
    missing_list = []
    incomplete_list = []
    change_list = []
    error_list = []
    cur_dir = os.getcwd()
    for indx_ni_mass, i_ni_mass in enumerate(parameters['ni_mass']):
        for indx_ni_mix, i_ni_mix in enumerate(parameters['ni_mix']):
            for indx_mass, imass in enumerate(parameters['mass']):
                print('finished mass {} for Ni mass {}'.format(imass, i_ni_mass))
                for indx_energy, ienergy in enumerate(parameters['explosion_energy']):
                    for indx_density, idensity in enumerate(parameters['density_1D']):
                        for indx_radius, iradius in enumerate(parameters['wind_extent']):
                            loop_start = time.time()
                            path = os.path.join(base_model_dir,
                                'Ni_mass_{:1.4f}'.format(i_ni_mass),
                                'Ni_mixing_{:1.1f}'.format(i_ni_mix),
                                'M{:2.1f}'.format(imass), #TODO - figure out formating so if string is single digit, prepended with a 0
                                'E_{:1.3f}'.format(ienergy),
                                'K_{:2.1f}'.format(idensity), #TODO - figure out formating so if string is single digit, prepended with a 0
                                'R_{}'.format(int(iradius)))
                            tarfilename = '{}.tar.gz'.format(path)
                            unzip_dir_exists = os.path.exists(path)
                            zip_dir_exists = os.path.isfile(tarfilename)
                            print(path)
                            #Neither directory exists
                            # --> record directory in missing_list
                            if (unzip_dir_exists is False) and (zip_dir_exists is False):
                                #print('Neither')
                                missing_list.append(path)
                            #Only an unzipped directory exists
                            # --> Look at last time step
                            # --> if last time step is the end time
                            # --> if last time is not the end time
                            #        add file to incomplete list
                            # --> zip the directory
                            # --> add a message to the log file
                            elif (unzip_dir_exists is True) and (zip_dir_exists is False):
                                #print('only unzip')
                                last_time_step = get_last_line(os.path.join(path, 'Data', 'lum_observed.dat'))
                                if last_time_step != parameters['endtime']:
                                    incomplete_list.append('last entry: {} for {}'.format(last_time_step, path))
                                #tar the directory
                                log_message = 'Zipping {}'.format(path)
                                change_list.append(log_message)
                                os.chdir(os.path.split(path)[0])
                                tar = tarfile.open(tarfilename, 'w:gz')
                                tar.add(os.path.split(path)[1])
                                tar.close()
                                os.chdir(cur_dir)
                                shutil.rmtree(path)
                            #Only zip file
                            # --> Extract the lum_observed file
                            # --> Look at the last time recorded
                            # --> If last time recorded is not the end time, record in incomplete list
                            elif (unzip_dir_exists is False) and (zip_dir_exists is True):
                                #print('only zip')
                                ##try:
                                ##    tar = tarfile.open(tarfilename, 'r')
                                ##except (tarfile.ReadError,tarfile.CompressionError,tarfile.StreamError,tarfile.ExtractError):
                                ##    error_list.append(tarfile)
                                ##    continue
                                ##lum_file = os.path.join('R_{}'.format(int(iradius)), 'Data', 'lum_observed.dat')
                                ##ofile = tar.extractfile(lum_file)
                                ##last_time_step = get_last_line(ofile, handle=True)
                                ##tar.close()
                                ##if last_time_step != parameters['endtime']:
                                ##    incomplete_list.append('{}; last time step: {}'.format(path, last_time_step))
                            #If both exist   
                            # --> Extract the lum_observed file
                            # --> Look at the last time recorded
                            # --> Look at last time recorded in the unzipped lum_observed file
                            # --> If zip file has more data than unzipped file
                            #        remove unzipped file
                            #        record removal in log file
                            #        check that kept data is complete, if not, add to incomplete list
                            # --> If unzipped file has more data than the zipped file
                            #        remove the zipped file
                            #        zip the unzipped file
                            #        record remove and zipping in log file
                            #        check that kept data is complete, if not, add to incomplete list
                        
                            elif (unzip_dir_exists is True) and (zip_dir_exists is True):
                                ##print('both')
                                untar_file = os.path.join(path, 'Data', 'lum_observed.dat')
                                lum_file = os.path.join('R_{}'.format(int(iradius)), 'Data', 'lum_observed.dat')
                                last_time_step_untar = get_last_line(untar_file)
                                try:
                                    tar = tarfile.open(tarfilename, 'r')
                                except (tarfile.ReadError,tarfile.CompressionError,tarfile.StreamError,tarfile.ExtractError):
                                    error_list.append(tarfilename)
                                    continue
                                ofile = tar.extractfile(lum_file)
                                last_time_step_tar = get_last_line(ofile, handle=True)
                                tar.close()
                                if last_time_step_tar >= last_time_step_untar:
                                    shutil.rmtree(path)
                                    change_list.append('Deleted {}'.format(path))
                                    if last_time_step_tar < parameters['endtime']:
                                        incomplete_list.append('{}; last time step: {}'.format(path, last_time_step_tar))

                                else:
                                    shutil.rmtree(tarfilename)
                                    os.chdir(os.path.split(path)[0])
                                    tar = tarfile.open(tarfilename, 'w:gz')
                                    tar.add(os.path.split(path)[1])
                                    tar.close()
                                    os.chdir(cur_dir)
                                    shutil.rmtree(path)
                                    log_message = 'Zipping {}'.format(path)
                                    change_list.append(log_message)
                                    change_list.append('Deleted {}'.format(tarfilename))
                                    if last_time_step_untar < parameters['endtime']:
                                        incomplete_list.append('{}; last time step: {}'.format(path, last_time_step_untar))
                            loop_end = time.time()
                            #print('loop time = {}'.format(loop_end - loop_start))
    ofile_change = open('change_log.txt', 'w')
    for iline in change_list:
        ofile_change.write('{}\n'.format(iline))
    ofile_change.close()

    ofile_missing = open('missing_files.txt', 'w')
    for iline in missing_list:
        ofile_missing.write('{}\n'.format(iline))
    ofile_missing.close()

    ofile_incomplete = open('incomplete_files.txt', 'w')
    for iline in incomplete_list:
        ofile_incomplete.write('{}\n'.format(iline))
    ofile_incomplete.close()
    
    ofile_error = open('error_files.txt', 'w')
    for iline in error_list:
        ofile_error.write('{}\n'.format(iline))
    ofile_error.close()
    end_time = time.time()
    print('total time = {}'.format(end_time - start_time))                    


#Run locally
#Catch error for improperly tarred file
                            