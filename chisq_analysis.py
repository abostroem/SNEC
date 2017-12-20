import os

from astropy.io import ascii as asc
import astropy.units as u

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from visualization import zscale
import define_filters


class SnecAnalysis(object):
    def __init__(self, snname, base_dir, S2_start, S2_end, ni_mass,
                 ni_mixing, masses, energies, time_offsets, 
                 Kvalues, radii, fig_dir=None):
        self.name = snname
        self.base_dir = base_dir
        self.S2_start = S2_start
        self.S2_end = S2_end
        self.ni_mass = ni_mass
        self.ni_mixing = ni_mixing
        self.masses = masses
        self.energies = energies
        self.time_offsets = time_offsets
        self.Kvalues = Kvalues
        self.radii = radii
        self.fig_dir = fig_dir

    def get_breakout_time(self, model_dir):
        ofile = open(os.path.join(model_dir, 'info.dat'), 'r')
        all_lines = ofile.readlines()
        if len(all_lines)>6:
            time_breakout = float((all_lines[5].split('=')[1]).strip('seconds\n'))
        else: #SN never got to breakout
            time_breakout = None
        return time_breakout
    
    def prepare_model_data(self, model_dir):
        model_mag_tbdata = asc.read(os.path.join(model_dir,'magnitudes.dat'),
                                    names=['time', 'temp', 'weird_mag', 
                                           'u', 'g', 'r', 'i', 'z', 'U', 
                                           'B', 'V', 'R', 'I'])
        #Observers call t_explosion the time when the explosion is first visible, therefore
        #the models should use this same conventions
        time_breakout = self.get_breakout_time(model_dir)
        if time_breakout is not None:
            model_mag_tbdata['time'] = ((model_mag_tbdata['time']-time_breakout)*u.second).to(u.day).value
        else:
            model_mag_tbdata=None
        return model_mag_tbdata


    def calc_chisq(self, sn_lc):
        skip_filters = []
        start_now = True
        if not os.path.exists('chisq_table.txt'):
            ofile = open('chisq_table.txt', 'w')
        else:
            ofile = open('chisq_table.txt', 'r')
            all_lines = ofile.readlines()
            if len(all_lines) > 0:
                last_line = all_lines[-1]
                split_last_line = last_line.split(',')
                last_complete_dir = os.path.join(self.base_dir, 
                                 'Ni_mass_{:1.4f}'.format(float(split_last_line[0])),
                                 'Ni_mixing_{:1.1f}'.format(float(split_last_line[1])),
                                 'M{:2.1f}'.format(float(split_last_line[2])),
                                 'E_{:1.3f}'.format(float(split_last_line[3])),
                                 'K_{:2.1f}'.format(float(split_last_line[4])), 
                                 'R_{}'.format(int(split_last_line[5])),
                                 'Data')
                start_now = False
        with open('chisq_table.txt', 'a') as output_ofile:
            with open('missing_mag_files.txt', 'a') as missing_ofile:
                chisq = np.ones((len(self.ni_mass), 
                                 len(self.ni_mixing), 
                                 len(self.masses), 
                                 len(self.energies),
                                 len(self.Kvalues),
                                 len(self.radii),
                                 len(self.time_offsets)))*np.nan
                for ni_mindx, i_ni_mass in enumerate(self.ni_mass):
                    for ni_indx, i_ni_mix in enumerate(self.ni_mixing):
                        for mindx, imass in enumerate(self.masses):
                            for eindx, ienergy in enumerate(self.energies):
                                for kindx, idensity in enumerate(self.Kvalues):
                                    for rindx, iradius in enumerate(self.radii):
                                        #read in model data 
                                        #TODO: make this and prep_snec look at the same piece of code
                                        model_dir = os.path.join(self.base_dir, 
                                                                 'Ni_mass_{:1.4f}'.format(i_ni_mass),
                                                                 'Ni_mixing_{:1.1f}'.format(i_ni_mix),
                                                                 'M{:2.1f}'.format(imass),
                                                                 'E_{:1.3f}'.format(ienergy),
                                                                 'K_{:2.1f}'.format(idensity), 
                                                                 'R_{}'.format(int(iradius)),
                                                                 'Data')
                                        if start_now is False:
                                            if model_dir != last_complete_dir:
                                                continue
                                            else:
                                                start_now = True
                                                continue
                                        if not os.path.exists(os.path.join(model_dir, 'magnitudes.dat')):
                                            missing_ofile.write("Missing Mag File: NiMass={},NiMix={},M={},E={},K={}, R={}\n".format(i_ni_mass, i_ni_mix, imass, ienergy, idensity, iradius))
                                            continue
                                        model_mag_tbdata = self.prepare_model_data(model_dir)
                                        if model_mag_tbdata is not None: #Successful explosion model
                                            chisq_filters = []
                                            #Loop over time shifts
                                            for tindx, toffset in enumerate(self.time_offsets):
                                                #Loop over filters
                                                for ifilter in sn_lc.abs_mag.keys():
                                                    s2_indx = (sn_lc.phase[ifilter] >= self.S2_start)& (sn_lc.phase[ifilter] <=self.S2_end)
                                                    if len(sn_lc.phase[ifilter][s2_indx]) > 5:
                                                        if model_mag_tbdata['time'][-1]+toffset > self.S2_end:
                                                            interp_mod_mag = np.interp(sn_lc.phase[ifilter][s2_indx]+toffset, 
                                                                                       model_mag_tbdata['time'], 
                                                                                       model_mag_tbdata[ifilter])
                                                            chisq_tmp = np.sum(((sn_lc.abs_mag[ifilter][s2_indx]-interp_mod_mag)/sn_lc.abs_mag_err[ifilter][s2_indx])**2)
                                                            chisq_filters.append(chisq_tmp)
                                                        else:
                                                            missing_ofile.write("Failed (LC too short) Model: NiMass={},NiMix={},M={},E={},K={}, R={}\n".format(i_ni_mass, i_ni_mix, imass, ienergy, idensity, iradius))
                                                            chisq_filters.append(1E10) #if a model doesn't explode it should never be the best model
                                                
                                                    else:
                                                        skip_filters.append(ifilter)  
                                                chisq[ni_mindx, ni_indx, mindx, eindx, kindx, rindx, tindx] = np.sum(np.array(chisq_filters))
                                                output_ofile.write('{},{},{},{},{},{},{},{}\n'.format(i_ni_mass, i_ni_mix, imass, ienergy, idensity, iradius, toffset, chisq[ni_mindx, ni_indx, mindx, eindx, kindx, rindx, tindx]))
                                                output_ofile.flush()          
                                        else:
                                            missing_ofile.write("Failed (unexploded) Model: NiMass={},NiMix={},M={},E={},K={}, R={}\n".format(i_ni_mass, i_ni_mix, imass, ienergy, idensity, iradius))
        print('Skipped Filters {} b/c <5 points in fit region'.format(set(skip_filters)))
        self.min_indx_base_mod = np.where(chisq == np.nanmin(chisq))
        self.best_ni_mass = self.ni_mass[self.min_indx_base_mod[0][0]]
        self.best_ni_mix = self.ni_mixing[self.min_indx_base_mod[1][0]]
        self.best_mass = self.masses[self.min_indx_base_mod[2][0]]
        self.best_energy = self.energies[self.min_indx_base_mod[3][0]]
        self.best_Kvalue = self.Kvalues[self.min_indx_base_mod[4][0]]
        self.best_radius = self.radii[self.min_indx_base_mod[5][0]]
        self.best_time_offset = self.time_offsets[self.min_indx_base_mod[6][0]]
        print('Best chi square for Ni mixing={}, Ni mass = {}, mass={}, energy={}, density = {}, radius = {}, time offset={}'.format(
                    self.best_ni_mass,
                    self.best_ni_mix,
                    self.best_mass,
                    self.best_energy,
                    self.best_Kvalue,
                    self.best_radius,
                    self.best_time_offset))

        self.best_model_dir = os.path.join(self.base_dir, 
                                         'Ni_mass_{:1.4f}'.format(self.best_ni_mass),
                                         'Ni_mixing_{:1.1f}'.format(self.best_ni_mix),
                                         'M{:2.1f}'.format(self.best_mass),
                                         'E_{:1.3f}'.format(self.best_energy),
                                         'K_{:2.1f}'.format(self.best_Kvalue), 
                                         'R_{}'.format(int(self.best_radius)),
                                         'Data')

        self.best_model_tbdata = self.prepare_model_data(self.best_model_dir)
        self.model_chisq = chisq
        


    def plot_lightcurve(self, sn_lc, band='all'):
        if band == 'all':
            bands = sn_lc.abs_mag.keys()
        elif len(bands) == 1:
            bands = [band]
        else:
            bands = band
        offset = 0
        filter_dict = define_filters.define_filters()
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for iband in filter_dict.keys():
            if iband in bands:
                l = ax.errorbar(sn_lc.phase[iband], sn_lc.abs_mag[iband]-offset, sn_lc.abs_mag_err[iband], fmt='o', linestyle='none', label='{} data+{}'.format(iband, offset))
                ax.plot(self.best_model_tbdata['time'], self.best_model_tbdata[iband]-offset, color=l[0].get_color(), ls='--', label='{}+{}'.format(iband, offset))

        ax.set_xlabel('Phase (days)')
        ax.set_ylabel('Absolute Magnitude + Offset')
        ax.set_title('Light Curve and SNEC Models for {}'.format(self.name))
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.legend(loc=3, ncol=2)
        plt.savefig(os.path.join(self.fig_dir, '{}_snec_lc.pdf'.format(self.name)))
        
 
# Write every step
# Check that file exists
# Write issues to log file   
# Restart from last check    
                
        
            
        
                         

                    
                
