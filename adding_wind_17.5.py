import os, sys
from sys import argv
import subprocess

#import numpy as np
#from pylab import *
from matplotlib import rcParams
from scipy import interpolate


#Run as many times as you want with different KK values
script,KK = argv

msol = 1.98e33
rsol = 6.96e10
ggrav = 6.6742e-8
velocity_of_wind = 1.0e6

### --------------------- Parameters --------------------

delta_r = 1 #Radius step in solar radii

R_extent = ['700','800','900','1000','1100','1200','1300','1400','1500','1600','1700','1800','1900','2000','2100','2200']
#R_extent = ['2300','2400','2500','2600','2700','2800','2900','3000','3100','3200','3300','3400','3500','3600','3700','3800']

for k in range(len(R_extent)):

    mainfolder = os.path.join('profiles_wind_17','K{}_R{}'.format(KK,R_extent[k])) #Change

    os.makedirs(mainfolder, exist_ok=True)

    M_ZAMS = [17] #Change
    string = ["17"] #Change

    rho_attach_gl = []
    vel_esc_gl = []
    mdot_gl = []
    mass_in_wind_gl = []
    tau_wind_gl = []


    for i, imass in enumerate(M_ZAMS):

        fname = "sukhbold_profiles/s{0}/profiles/s{0:2.1f}.short".format(imass) #change
        fname_iso = "sukhbold_profiles/s{0}/profiles/s{0:2.1f}.iso.dat".format(imass) #change

        os.makedirs(os.path.join(mainfolder,'s{}'.format(imass)), exist_ok=True)
        os.makedirs(os.path.join(mainfolder,'s{}'.format(imass),'profiles'), exist_ok=True)

        ### --------------------- Read the profile of the core --------------------

        mass = []
        radius = []
        temp = []
        rho = []
        rho_log = []
        vel = []

        for l in open(fname, 'r').readlines():
            s = l.split()
            if len(s) != 1:
                mass.append(float(s[1])/msol)
                radius.append(float(s[2])/rsol)
                temp.append(float(s[3]))
                rho.append(float(s[4]))
                rho_log.append(log10(float(s[4])))
                vel.append(float(s[5]))
        #Attach to profile when densities are the same
        rho_wind_attach = float(KK)/((radius[-1]*rsol)**2)

        vel_esc = sqrt(2*ggrav*mass[-1]*msol/(radius[-1]*rsol))
        mdot_msol_yr_max = rho_wind_attach*4*math.pi*(radius[-1]*rsol)**2*velocity_of_wind*365*24*3600/msol

        rho_attach_gl.append(rho_wind_attach)
        vel_esc_gl.append(vel_esc)
        mdot_gl.append(mdot_msol_yr_max)

        ### -------------------- Building profiles with the wind ----------------------

        radius_wind = [] #in solar radii
        rho_wind = []
        rho_log_wind = []
        mass_wind = [] #in solar masses
        temp_wind = []
        vel_wind = []

        for l in range(10000000): #(just big number to cover all the profile)
            if l < len(radius) and rho[l] > rho_wind_attach:
                radius_wind.append(radius[l])
                rho_wind.append(rho[l])
                rho_log_wind.append(rho_log[l])
                mass_wind.append(mass[l])
                temp_wind.append(temp[l])
                vel_wind.append(vel[l])
            else:
                if (radius_wind[l-1]+delta_r) > float(R_extent[k]):
                    break
                radius_wind.append(radius_wind[l-1]+delta_r)
                rho_wind.append(rho_wind_attach*(radius[-1]/radius_wind[l])**2)
                rho_log_wind.append(log10(rho_wind[l]))
                mass_wind.append(mass_wind[l-1]+4.0*math.pi/3.0*(radius_wind[l]**3 -
                            radius_wind[l-1]**3)*rho_wind[l]*rsol**3/msol)
                temp_wind.append(temp[-1])
                vel_wind.append(velocity_of_wind)

        mass_in_wind = mass_wind[-1] - mass[-1]
        mass_in_wind_gl.append(mass_in_wind)

        ### ------------ calculating the optical depth of the wind --------------
        #This is informational - SNEC will recalculate more carefully
        kappa = 1.0e-4
        tau_wind = 0
        for l in range(len(radius_wind)-1):
            if radius_wind[len(radius_wind)-1-l-1] < radius[-1]:
                break
            tau_wind = tau_wind + kappa*rsol*(radius_wind[len(radius_wind)-1-l]
                -radius_wind[len(radius_wind)-1-l-1])*rho_wind[len(radius_wind)-1-l-1]

        tau_wind_gl.append(tau_wind)

        ### -------------------- reading a composition file ------------------------

        iso_lines = []

        n_line = 0

        for l in open(fname_iso, 'r').readlines():
            if n_line == 0:
                s = l.split()
                ncells = int(s[0])
                niso = int(s[1])
            elif n_line == 1:
                A_line = l
            elif n_line == 2:
                Z_line = l
            else:
                iso_lines.append(l)
                s = l.split()
                last_line = s

            n_line = n_line + 1

        ### --------------------- Output the profiles -----------------------------
        # note, that here we convert mass and radius in cgs units, used by the code.
        # in the script they are in the solar masses and radii for convenience of plotting
            
        outfile = open(os.path.join(mainfolder,'s{}'.format(M_ZAMS[i]),'profiles','s{}_{}.short'.format(string[i],R_extent[k])),"w")

        outfile.write(str(len(radius_wind)) + '\n')

        for l in range(len(radius_wind)):
            if l == len(radius_wind) - 1:
                outfile.write(str(l+1) + '   ' + str(mass_wind[l]*msol) + '   ' + str(radius_wind[l]*rsol)
                          + '   ' + str(temp_wind[l]) + '   ' + str(rho_wind[l]) + '   '
                          + str(vel_wind[l]))
            else:
                outfile.write(str(l+1) + '   ' + str(mass_wind[l]*msol) + '   ' + str(radius_wind[l]*rsol)
                              + '   ' + str(temp_wind[l]) + '   ' + str(rho_wind[l]) + '   '
                              + str(vel_wind[l]) + '\n')

        outfile.write('\n')

        outfile.close()

        outfile_iso = open(mainfolder+'/s'+str(M_ZAMS[i])+'/profiles'+'/s'+string[i]+'_'+str(R_extent[k])+'.iso.dat',"w")

        outfile_iso.write(str(ncells+1) + '   ' + str(niso) + '\n')
        outfile_iso.write(A_line)
        outfile_iso.write(Z_line)

        for l in range(len(iso_lines)):
            outfile_iso.write(iso_lines[l])
        #outfile_iso.write('\n') # for some reason it was used before...

        outfile_iso.write(str(mass_wind[-1]*msol) + '        ' + str(radius_wind[-1]*rsol) + '          ')

        for l in range(2,len(last_line)):
            outfile_iso.write(last_line[l] + '   ')

        outfile_iso.write('\n')

        outfile_iso.close()

    ######### ------------ table with the wind properties

    outfile_info = open(mainfolder+'/info.dat',"w")

    for l in range(len(string)):
        outfile_info.write(string[l] + ": rho attach=" + str(rho_attach_gl[l]) + ", vel esc=" +
                      str(vel_esc_gl[l]) + ", mdot=" + str(mdot_gl[l]) + ", mass in wind=" +
                      str(mass_in_wind_gl[l]) + ", tau wind=" + str(tau_wind_gl[l]) + '\n')

    outfile_info.close()
