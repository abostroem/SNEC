import os
import glob
import shutil

base_dir = '/dark/SNEC/snec_models'
from_dir = '/home/bostroem/SNEC_bare_15oz/SNEC/snec_models'

input('THIS CODE ISNT TESTED - VERIFY THAT IT WORKS BEFORE RUNNING (enter to continue, ctl+c to exit ')
for ini_mass in ni_masses:
    new_path = os.path.join(base_dir, 'Ni_mass_{:1.4f}'.format(i_ni_mass))
    if not os.path.exists(new_path):
        os.mkdir(new_path)
    for ini_mix in ni_mix:
        new_path = os.path.join(new_path, 'Ni_mixing_{:1.1f}'.format(i_ni_mix))
        if not os.path.exists(new_path):
            os.mkdir(new_path)
        for imass in masses:
            new_path = os.path.join(new_path, 'M{:2.1f}'.format(imass))
            if not os.path.exists(new_path):
                os.mkdir(new_path)
            for k in kvalues:
                new_path = os.path.join(new_path, 'K_{:2.1f}'.format(idensity))
                if not os.path.exists(new_path):
                    os.mkdir(new_path)
                for r in radii:
                    new_path = os.path.join(new_path, 'R_{}'.format(int(iradius)))
                    if not os.path.exists(new_path):
                        os.mkdir(new_path)
                    new_files = glob.glob(os.path.join(new_path.replace(base_dir, from_dir), '*.gz'))
                    for ifile in new_files:
                        if not os.path.exists(ifile.replace(from_dir, base_dir)):
                            shutil.copyfile(ifile.replace(from_dir, base_dir))
                        else:
                            print('file already exists, not copying: {}'.format(ifile))
                        