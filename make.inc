# some definitions controlling compilation
export F90=gfortran
export F90FLAGS=-g -O3 -Warray-bounds  -fbounds-check
#export F90FLAGS=-g -O0 #for debugging run with lldb
export LDFLAGS=$(F90FLAGS) 
#use the following line when compiling on unix
#export LAPACKLIBS=-llapack  #-lg2c
#use the following line when compiling on mac
export LAPACKLIBS=-framework Accelerate
export LIBS+=$(LAPACKLIBS)

#don't use lapack for the moment
#export LAPACKLIBS=-llapack  -lg2c
#export LIBS+=$(LAPACKLIBS)

