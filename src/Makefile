include ../make.inc


SOURCES=problem.F90\
	allocate_vars.F90\
	eos_switch.F90\
	eos_content.F90\
	grid.F90\
	read_profile.F90\
	read_profile_compositions.F90\
	output.F90\
	timestep.F90\
	blstep.F90\
	artificial_viscosity.F90\
	opacity.F90\
	bicubic_interpolator.F90\
	hydro.F90\
	hydro_rad.F90\
	outinfo.F90\
	input_parser.F90\
	conservation.F90\
	analysis.F90\
	read_BolCorr.F90\
	opal_opacity.F90\
	nickel.F90\
	simple_saha.F90\
	saha_tables.F90\
	luminosity.F90\
	bomb_profile.F90\
	arrays.F90\
	boxcar.F90\
	shock_capture.F90\
	snec.F90


CLEANSTUFF=rm -rf *.o *.mod

OBJECTS=$(SOURCES:.F90=.o)

all: blmod.o outinfo.o snec

blmod.o: blmod.F90
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) -c blmod.F90	

snec: $(OBJECTS) 
	$(F90) $(LDFLAGS) -o ../snec  blmod.o $(OBJECTS) $(EXTRAOBJECTS) $(EXTRALIBS) $(LIBS) 

$(OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@

snec.o: blmod.o  

input_parser.o: blmod.o 


clean: 
	$(CLEANSTUFF)
