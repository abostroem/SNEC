# Install
* In make.inc comment/uncomment line for specific machine
* Run make clean
* Run make
* Create a Data folder
This produces an executable *snec* that can be run

# Batch processing
TBD


# Output
In the .xg files for each time frame the first column of the output is mass coordinate in grams (apart from 'mass.xg', where the first column is radius and the second is mass). The second column is the following (in the order the files are listed in output.F90):

* vel.xg 			- velocity
* rho.xg 			- density
* ye.xg	    	    	- electron fraction (ratio of electrons to baryons)
* press.xg   	     	- pressure
* cs2.xg  			- speed of sound squared
* Q.xg 			- artificial viscosity
* eps.xg 			- specific internal energy
* temp.xg 			- temperature
* lum.xg 			- luminosity at each point (as given by Eq.4 of the notes)
* tau.xg 			- optical depth
* delta_time.xg 	- timestep as computed at each point in the timestep.F90 (Eq.34 of the notes)
* radius.xg 		- radial coordinate
* kappa.xg 		- opacity (after the opacity floor is applied) (see Sec.3.3 of the notes)
* kappa_table.xg 	- opacity (before the opacity floor is applied)
* logR_op.xg 		- logarithm of R, defined in the Sec.3.3, see calculations in opacity.F90
* logT.xg 			- logarithm of temperature (sometimes more convenient than just temperature)
* p_rad.xg 			- radiation pressure
* Ni_deposit_function.xg  - deposition function of the energy of gamma rays, as given by Eq.46 of the notes
* He_1.xg He_2.xg He_3.xg - ionization fractions of helium (_1 is neutral)
* H_1.xg, H_2.xg 		- ionization fractions of hydrogen (_1 is neutral)
* free_electron_frac.xg   - fraction of free electrons
* E_shell.xg , time_diff.xg , time_exp.xg - some quantities, used for analysis, as given in the Sec.7 of the notes (in particular, see references [20] and [21])
* photosphere_tracer.xg 	- =1 at the position of the photosphere, =0 everywhere else (just to keep track of the motion of the photoshere with time)

For all .dat files the first column is time in seconds, the second is the following (in the order they appear in output.F90):

* T_eff.dat 			- effective temperature (as in Eq.86 of the notes)
* Ni_total_luminosity.dat 	- total amount of energy per second deposited in the system by gamma rays
* lum_observed.dat        	- observed luminosity (sum of the photospheric luminosity and Ni contribution, Eq.83 of the notes)
* index_photo.dat 		- grid index of the photosphere location
* lum_photo.dat 		- luminosity at the photosphere
* mass_photo.dat 		- mass coordinate of the photosphere
* vel_photo.dat 		- velocity at the photosphere
* rad_photo.dat 		- radius of the photosphere
* opacity_corrupted.dat      - =1 if the photosphere passes through the region where the tabulated data for the opacity is absent (calculated in analysis.F90, see discussion in Sec.3.3 of the notes), =0 otherwise
* index_lumshell.dat 	- grid index of the luminosity shell (calculated in analysis.F90,for definition see http://adsabs.harvard.edu/abs/2010ApJ...725..904N)
* mass_lumshell.dat 	- mass coordinate of the luminosity shell

In addition, there are few files with the initial data, generated in problem.F90 (rho_initial.dat, rad_initial.dat, mass_initial.dat, delta_mass_initial.dat, press_initial.dat). First column in these files is the number of the gridpoint, the second column gives density, radius, mass coordinate, mass resolution and pressure, respectively. The file conservation.dat gives the energy balance, as described in the Sec.2.2.4 of the notes.