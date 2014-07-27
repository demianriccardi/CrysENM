# the directories containing the OpenGL libraries, f90gl libraries, GLUT
HLIB = /Users/demianriccardi/lib
DIR = -L$(HLIB)/CrysENM -lcrysenm -L$(HLIB)/crysfml/GFortran/LibC -lcrysfml -L$(HLIB)/ARPACK -larpack_osx -L$(HLIB)/LAPACK95 -llapack95

LIBS = -framework Accelerate
MODSCRYSFML =  -I$(HLIB)/crysfml/GFortran/LibC 
MODSCRYSENM =  -I$(HLIB)/CrysENM
MODLAPACK95 =  -I$(HLIB)/LAPACK95

F90 = gfortran -O3 #-finstrument-functions #-i-static  

#bfactcorr: bfactcorr.f90; $(F90) -o bfactcorr bfactcorr.f90 $(MODSCRYSFML) $(MODSCRYSENM) $(MODLAPACK95) $(DIR) $(LIBS)  

chunkwrite: chunkwrite.f90; $(F90) -o chunkwrite chunkwrite.f90  $(MODSCRYSFML) $(MODSCRYSENM) $(MODLAPACK95) $(DIR) $(LIBS)  

#all: bfactcorr

# -CB   pretty sweet.  Checks the bounds of the arrays at runtime!
#-openmp # -lm -i-dynamic 
#shared-intel #-parallel  -par_threshold0 -par_report2 -ipo  

clean: rm -f *.o *.mod
