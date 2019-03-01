import sys
import os
import math

import astropy.constants as c
import astropy.units as u
import numpy as np


def generate_profile_directory_name(imass):
    if int(imass) == imass: #imass is an integer value
        str_mass = '{}'.format(int(imass))
        mass_format='single'
    elif np.isclose(imass%.1, 0.1): #imass has a single decimal
        str_mass = '{:2.1f}'.format(imass)
        mass_format='single'
    elif np.isclose(imass%0.05, 0.05): #imass has two decimals
        str_mass = '{:2.2f}'.format(imass)
        mass_format='double'
    if mass_format == 'single':
        profile_basename = os.path.join('sukhbold_profiles', 's{}'.format(str_mass), 'profiles', 's{:2.1f}'.format(imass)) 
    elif mass_format == 'double':
        profile_basename = os.path.join('sukhbold_profiles','s{}'.format(str_mass), 'profiles', 's{:2.2f}'.format(imass))
    return profile_basename

def add_wind(parameters):
    msol = c.M_sun.to(u.gram).value
    rsol = c.R_sun.to(u.cm).value
    ggrav = c.G.to(u.cm**3/u.gram/u.second**2).value
    velocity_of_wind = parameters['wind_velocity']
    #Verified that any differences in the output file for M=17, K=10, R=2200 were from the
    #change in constants to more exact values; with old constants files were identical

    ### --------------------- Parameters --------------------



    delta_r = parameters['delta_r']

    R_extent = parameters['wind_extent']

    for imass in parameters['mass']:
        for idensity in parameters['density_1D']:
            if idensity == 0:  #Copy profiles without adding wind
                mainfolder = os.path.join('sukhbold_profiles_wind',
                    'M{:2.1f}'.format(imass), #TODO - figure out formating so if string is single digit, prepended with a 0
                    'K_0.0', #TODO - figure out formating so if string is single digit, prepended with a 0
                    'R_0000',
                    'profiles')

                os.makedirs(mainfolder, exist_ok=True)
                base_profile_dir = generate_profile_directory_name(imass)

                fname = base_profile_dir+'.short'
                fname_iso = base_profile_dir+'.iso.dat'
                shutil.copyfile(fname, mainfolder)
                shutil.copyfile(fname_iso, mainfolder)
            else:
                idensity = idensity*10**17
                for iradius in parameters['wind_extent']:
                    mainfolder = os.path.join('sukhbold_profiles_wind',
                        'M{:2.1f}'.format(imass), #TODO - figure out formating so if string is single digit, prepended with a 0
                        'K_{:2.1f}'.format(idensity/10**17), #TODO - figure out formating so if string is single digit, prepended with a 0
                        'R_{}'.format(int(iradius)))
    
                    os.makedirs(mainfolder, exist_ok=True)
    
                    rho_attach_gl = []
                    vel_esc_gl = []
                    mdot_gl = []
                    mass_in_wind_gl = []
                    tau_wind_gl = []
    
                    base_profile_dir = generate_profile_directory_name(imass)
    
    
                    fname = base_profile_dir+'.short'
                    fname_iso = base_profile_dir+'.iso.dat'
    
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
                            rho_log.append(math.log10(float(s[4])))
                            vel.append(float(s[5]))
                    #Attach to profile when densities are the same
                    rho_wind_attach = float(idensity)/((radius[-1]*rsol)**2)
    
                    vel_esc = math.sqrt(2*ggrav*mass[-1]*msol/(radius[-1]*rsol))
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
                            if (radius_wind[l-1]+delta_r) > float(iradius):
                                break
                            radius_wind.append(radius_wind[l-1]+delta_r)
                            rho_wind.append(rho_wind_attach*(radius[-1]/radius_wind[l])**2)
                            rho_log_wind.append(math.log10(rho_wind[l]))
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
            
                    outfile = open(os.path.join(mainfolder,'{}.short'.format(os.path.basename(base_profile_dir))),"w")
    
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
    
                    outfile_iso = open(os.path.join(mainfolder,'{}.iso.dat'.format(os.path.basename(base_profile_dir))),"w")
    
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
    
                    outfile_info = open(os.path.join(mainfolder,'info.dat'),"w")
    
                    outfile_info.write("{}: rho attach={}, vel esc={}, mdot={}, mass in wind={}, tau wind={}\n".format(
                                        imass, rho_attach_gl[0], vel_esc_gl[0], mdot_gl[0], mass_in_wind_gl[0], tau_wind_gl[0]))
    
                    outfile_info.close()
    