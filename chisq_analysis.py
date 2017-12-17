import os

from astropy.io import ascii as asc
from astropy.io import fits
import astropy.units as u

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import supernova as sn
from visualization import zscale
import define_filters


class SnecAnalysis(object):
    def __init__(self, snname, base_dir, S2_start, S2_end, 
                 ni_mixing, masses, energies, time_offsets, 
                 Kvalues=None, Kvalues_str=None, radii=None, fig_dir=None):
        self.name = snname
        self.base_dir = base_dir
        self.S2_start = S2_start
        self.S2_end = S2_end
        self.ni_mixing = ni_mixing
        self.masses = masses
        self.energies = energies
        self.time_offsets = time_offsets
        self.Kvalues = Kvalues
        self.Kvalues_str = Kvalues_str
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


    def calc_chisq_base_model(self, sn_lc):
        chisq = np.ones((len(self.masses), len(self.energies), len(self.time_offsets), len(self.ni_mixing)))*np.nan
        skip_filters = []
        for ni_indx, ni_mix in enumerate(self.ni_mixing):
            for mindx, imass in enumerate(self.masses):
                for eindx, ienergy in enumerate(self.energies):
                    #read in model data 
                    if ienergy == 1.0:
                        ienergy = '1.00'
                    model_dir = os.path.join(self.base_dir, 
                                                     sn_lc.name, 
                                                     'mixing_{}'.format(ni_mix),
                                                     'M{}'.format(imass),
                                                     'E_{:}'.format(ienergy),
                                                     'Data')
                    if not os.path.exists(os.path.join(model_dir, 'magnitudes.dat')):
                        print("Model doesn't exist: M={},E={}".format(imass, ienergy))
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
                                        print("Failed (LC too short) Model: mass={},energy={}, ni_mix={}".format(imass, ienergy,ni_mix))
                                        chisq_filters.append(1E10)
                                    chisq[mindx, eindx, tindx, ni_indx] = np.sum(np.array(chisq_filters))    
                                else:
                                    skip_filters.append(ifilter)            
                    else:
                            print("Failed (unexploded) Model: mass={},energy={}, ni_mix={}".format(imass, ienergy, ni_mix))
        print('Skipped Filters {} b/c <5 points in fit region'.format(set(skip_filters)))
        self.min_indx_base_mod = np.where(chisq == np.nanmin(chisq))
        self.best_mass = self.masses[self.min_indx_base_mod[0][0]]
        self.best_energy = self.energies[self.min_indx_base_mod[1][0]]
        self.best_time_offset = self.time_offsets[self.min_indx_base_mod[2][0]]
        self.best_ni_mixing = self.ni_mixing[self.min_indx_base_mod[3][0]]
        print('Best chi square for mass={}, energy={}, time offset={}, Ni mixing={}'.format(
                    self.best_mass,
                    self.best_energy,
                    self.best_time_offset,
                    self.best_ni_mixing))

        self.best_base_model_dir = os.path.join(self.base_dir, 
                                                 sn_lc.name, 
                                                 'mixing_{}'.format(self.best_ni_mixing),
                                                 'M{}'.format(self.best_mass),
                                                 'E_{:}'.format(self.best_energy),
                                                 'Data')
        self.best_base_model_tbdata = self.prepare_model_data(self.best_base_model_dir)
        self.base_model_chisq = chisq


    def calc_chisq_wind(self, sn_lc):
        chisq = np.ones((len(self.Kvalues), len(self.radii)))*np.nan
        for kindx, iK in enumerate(self.Kvalues):
            for rindx, iradius in enumerate(self.radii):
                #read in model data 
                model_dir = os.path.join(self.base_dir, 
                                                 '{}_wind'.format(sn_lc.name), 
                                                 'K{}'.format(self.Kvalues_str[kindx]),
                                                 'R_{:}'.format(iradius),
                                                 'Data')
                if not os.path.exists(os.path.join(model_dir, 'magnitudes.dat')):
                    print("Model doesn't exist: radius={}, K={}".format(iradius, iK))
                    continue
                model_mag_tbdata = self.prepare_model_data(model_dir)
                if model_mag_tbdata is not None: #Successful explosion model
                    chisq_filters = []
                    #Loop over filters
                    for ifilter in sn_lc.abs_mag.keys():
                        if ifilter not in ['U', 'u', 'B', 'b']:
                            S2_indx = (sn_lc.phase[ifilter] <=self.S2_end)
                            if model_mag_tbdata['time'][-1] > self.S2_end:
                                interp_mod_mag = np.interp(sn_lc.phase[ifilter][S2_indx], 
                                                           model_mag_tbdata['time'], 
                                                           model_mag_tbdata[ifilter])
                                chisq_tmp = np.sum(((sn_lc.abs_mag[ifilter][S2_indx]-interp_mod_mag)/sn_lc.abs_mag_err[ifilter][S2_indx])**2)
                                chisq_filters.append(chisq_tmp)
                            else:
                                print("Failed (LC too short) Model: K={},radius={}".format(iK, iradius))
                                chisq_filters.append(1E10)
                            chisq[kindx, rindx] = np.sum(np.array(chisq_filters))
                else:
                    print("Failed (unexploded) Model: K={},radius={}".format(iK, iradius))
        self.min_indx_csm_mod = np.where(chisq == np.nanmin(chisq))
        
        self.best_Kvalue = self.Kvalues[self.min_indx_csm_mod[0][0]]
        self.best_radius = self.radii[self.min_indx_csm_mod[1][0]] 
        print('Best chi square for K={}, radius={}'.format(
                    self.best_Kvalue,self.best_radius))
        self.best_csm_model_dir = os.path.join(self.base_dir, 
                                                 '{}_wind'.format(sn_lc.name), 
                                                 'K{}'.format(self.Kvalues_str[self.min_indx_csm_mod[0][0]]),
                                                 'R_{:}'.format(self.best_radius),
                                                 'Data')
        self.csm_model_chisq = chisq
        self.best_csm_model_tbdata = self.prepare_model_data(self.best_csm_model_dir)
        

    def plot_chisq_base_mod(self):
        num_xticks=6
        num_yticks=8
        #TODO - update to use different Ni Masses, this about how to display this
        with PdfPages(os.path.join(self.fig_dir, '{}_chisq_base_mod.pdf').format(self.name)) as pp:
            vmin, vmax = zscale(self.base_model_chisq[:,:, self.min_indx_base_mod[3][0],0])
            for tindx, toffset in enumerate(self.time_offsets):
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1) 
                im = ax.imshow(self.base_model_chisq[:,:,tindx,0], 
                            interpolation='nearest', 
                            #vmin=vmin, vmax=vmax,
                            extent=[self.energies[0]-.2, self.energies[-1]+.2, self.masses[0]-.2, self.masses[-1]+.2],
                            aspect='auto')   
                if toffset == self.best_time_offset:
                    ax.axvline(self.best_energy, color='r')
                    ax.axhline(self.best_mass, color='r')
                ax.set_xlabel('Energy')
                ax.set_ylabel('Mass')
                ax.set_title('Chi Square Grid for Mass/Energy of Base Model; Offset Time = {}'.format(toffset))
                fig.colorbar(im)
                xticks = ax.get_xticks()
                yticks = ax.get_yticks()
                pp.savefig(fig)
                plt.close()
            
    def plot_chisq_csm_mod(self):
        vmin, vmax = zscale(self.csm_model_chisq[:,:])
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1) 
        im = ax.imshow(self.csm_model_chisq[:,:], 
                    interpolation='nearest', 
                    #vmin=vmin, vmax=vmax,
                    extent=[self.radii[0], self.radii[-1], float(self.Kvalues[0]), float(self.Kvalues[-1])],
                    aspect='auto')   
        plt.axvline(self.best_radius, color='r')
        plt.axhline(self.best_Kvalue, color='r')         
        ax.set_xlabel('Radius (Rsun)')
        ax.set_ylabel('K Value')
        title = 'Chi Square Grid for K/Radius of CSM model, M={}, E={}, M_Ni={}'
        ax.set_title(title.format(self.best_mass, self.best_energy, self.best_ni_mixing))
        fig.colorbar(im)
        plt.savefig((os.path.join(self.fig_dir, '{}_chisq_csm_mod.pdf').format(self.name)))
        plt.close()
            
            
    def plot_lightcurve(self, sn_lc, base_model=True, csm_model=True, band='all'):
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
                if base_model is True:
                    ax.plot(self.best_base_model_tbdata['time'], self.best_base_model_tbdata[iband]-offset, color=l[0].get_color(), ls='--', label='{} Base+{}'.format(iband, offset))
                if csm_model is True:
                    ax.plot(self.best_csm_model_tbdata['time'], self.best_csm_model_tbdata[iband]-offset, color=l[0].get_color(), label='{} CSM+{}'.format(iband, offset))
                offset+=2
        ax.set_xlabel('Phase (days)')
        ax.set_ylabel('Absolute Magnitude + Offset')
        ax.set_title('Light Curve and SNEC Models for {}'.format(self.name))
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.legend(loc=3, ncol=2)
        plt.savefig(os.path.join(self.fig_dir, '{}_snec_lc.pdf'.format(self.name)))
        
        
                
        
            
        
                         

                    
                
