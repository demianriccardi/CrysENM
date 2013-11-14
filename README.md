CrysENM
=======

Fortran 90 toolkit for application of elastic network models to proteins in the crystalline state. 
I wrote this as a NLM postdoctoral trainee in the laboratory of George N. Phillips, Jr from 2007-2010.  
It is in need of some TLC (not TLS, wink).  TODO: Refactor, reduce, write testing framework, remove compiler dependencies, 
expand, DOCUMENT, and improve.  

Dependencies: 

1. intel fortran compiler: uses MKL libraries for hefty matrix diagonalization

2. ARPACK: http://www.caam.rice.edu/software/ARPACK/

3. CrysFML

It will work, but you have to know what you're doing!  

I will push programs that use the library as I clean them up.  Contact me if you are interested!  
A contact got me this far!

if used productively, please cite:

Riccardi D., Cui Q., Phillips Jr. G. N., 
“Application of elastic network models to proteins in the crystalline state”, Biophys. J., 96, 464 (2009)

Riccardi D., Cui Q., Phillips Jr. G. N., 
“Evaluating elastic network models of crystalline biological molecules with temperature factors, correlated motions, and diffuse X-ray scattering”, Biophys. J., 99, 2616 (2010)
