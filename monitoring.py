import os
import glob

flist = sorted(glob.glob('asassn15oz/mixing_5.0/M??/E*/Data/lum_observed.dat'))
tend=1.728E7

for ifile in flist:
    ofile = open(ifile, 'r')
    all_lines = ofile.readlines()
    ofile.close()
    print('{:2.2f} for {}'.format(float(all_lines[-1].split()[0])/tend, ifile))

flist = sorted(glob.glob('asassn15oz2/mixing_5.0/M??/E*/Data/lum_observed.dat'))
tend=2.592E7
for ifile in flist:
    ofile = open(ifile, 'r')
    all_lines = ofile.readlines()
    ofile.close()
    print('{:2.2f} for {}'.format(float(all_lines[-1].split()[0])/tend, ifile))

flist = sorted(glob.glob('asassn15oz_wind/*/*/Data/lum_observed.dat'))
tend=2.592E7
for ifile in flist:
    ofile = open(ifile, 'r')
    all_lines = ofile.readlines()
    ofile.close()
    print('{:2.2f} for {}'.format(float(all_lines[-1].split()[0])/tend, ifile))
