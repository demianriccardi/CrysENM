gfortran -c mod_constants.f90    -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_math.f90         -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_types.f90        -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_linalg.f90       -I/Users/demianriccardi/lib/crysfml/GFortran/LibC -I/Users/demianriccardi/lib/LAPACK95/lapack95_modules 
gfortran -c mod_arpack.f90       -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_inout.f90        -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_crysbuild.f90    -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_correl_dist.f90  -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_tempfact.f90     -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_kirchoff.f90     -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_hessblock.f90    -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_hessian.f90      -I/Users/demianriccardi/lib/crysfml/GFortran/LibC 
gfortran -c mod_curvefit.f90     -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_symop.f90        -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_bnm.f90          -I/Users/demianriccardi/lib/crysfml/GFortran/LibC 
gfortran -c mod_tls.f90          -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_lattdyn.f90      -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_valvec_store.f90 -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_vcov_store.f90   -I/Users/demianriccardi/lib/crysfml/GFortran/LibC
gfortran -c mod_intensity.f90    -I/Users/demianriccardi/lib/crysfml/GFortran/LibC

ar cr libcrysenm.a *.o
