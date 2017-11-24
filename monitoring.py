import os
import glob

flist = sorted(glob.glob('snec_models/*/*/*/*/*/*/Data/lum_observed.dat'))

for ifile in flist:
    with open(ifile, 'r') as ofile:
        all_lines = ofile.readlines()
        with open(ifile.replace('Data/lum_observed.dat', 'parameters'), 'r') as  param_file:
            all_lines = param_file.readlines()
            for iline in all_lines:
                if iline.startswith('tend'):
                    tend = float(iline.split('=')[1])
            print('{:2.2f} for {}'.format(float(all_lines[-1].split()[0])/tend, ifile))

