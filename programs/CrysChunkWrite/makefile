# makefile for building on OS X with gfortran 
#   gfortran installed with homebrewed gcc
# the dependencies (ARPACK crysfml and lapack95 are installed in ~/lib)
HLIB = /Users/dmr3/lib
DIR = -L$(HLIB)/CrysENM -lcrysenm -L$(HLIB)/crysfml/GFortran64/LibC -lcrysfml -L$(HLIB)/ARPACK -larpack_osx 

# os x links into lapack and blas
LIBS = -framework Accelerate

# includes for fortran modules
MODSCRYSFML =  -I$(HLIB)/crysfml/GFortran64/LibC 
MODSCRYSENM =  -I$(HLIB)/CrysENM

F90 = gfortran -O3 #-finstrument-functions #-i-static  

chunkwrite: chunkwrite.f90; $(F90) -o chunkwrite chunkwrite.f90  $(MODSCRYSFML) $(MODSCRYSENM) $(DIR) $(LIBS)  
crystal_coords: crystal_coords.f90; $(F90) -o crystal_coords crystal_coords.f90  $(MODSCRYSFML) $(DIR) $(LIBS)  

# some other compiler options used in the past
# -CB   pretty sweet.  Checks the bounds of the arrays at runtime!
#-openmp # -lm -i-dynamic 
#shared-intel #-parallel  -par_threshold0 -par_report2 -ipo  

clean: rm -f *.o *.mod
