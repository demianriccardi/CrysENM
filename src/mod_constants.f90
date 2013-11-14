module mod_constants     
! DMR: 11-11-2008
!   just a bunch of constants
! kind=dp from crysfml 
!  INTEGER, PARAMETER :: DP = KIND (1.d0) ! Alternatively
use cfml_globaldeps,                 only: dp,sp

Implicit none

real(DP), parameter :: Deg_Per_Rad     = 57.295779513082320876798155_DP
real(DP), parameter :: Rad_Per_Deg     =  0.017453292519943295769237_DP

real(DP), parameter :: pi              =  3.141592653589793238462643_DP

real(DP), parameter :: zero      = 0.0_DP
real(DP), parameter :: one       = 1.0_DP
real(DP), parameter :: two       = 2.0_DP
real(DP), parameter :: three     = 3.0_DP
real(DP), parameter :: four      = 4.0_DP
real(dp), parameter :: mp        = 1.660538782D-27    ! Kg
real(dp), parameter :: clight    = 2.99792458D08     ! m/s
real(dp), parameter :: Nav       = 6.02214179D23     ! mol-1
!real(dp), parameter :: CONVFRQ   = 2045.4828288_DP/c*2pi
real(dp), parameter :: CONVFRQ   = 108.5913586_DP ! sqrt(mol amu kcal-1)*(angs rad-1 cm-1)
real(DP), parameter :: twopi              =  two*pi
!real(dp), parameter :: CONVFRQ2   = 2045.4828288_DP/c
!real(dp), parameter :: CONVFRQ2   = 682.2996288_DP ! sqrt(mol amu kcal-1)*(angs rad-1 cm-1)
! this factor of two pi is still an issue
real(dp), parameter :: hplanck   =   6.62606896D-34     !J s       wikipedia
real(dp), parameter :: hbar      =   1.05457163D-34     !J s rad-1 wikipedia
real(dp), parameter :: hkcal     =   9.53707619D-14 ! kcal s mol-1      
real(dp), parameter :: hbarkcal  =   1.51787282D-14 ! kcal s mol-1 rad-1 
real(dp), parameter :: boltzc    =   1.3806503D-23 ! 
real(dp), parameter :: Rkcal     =   1.9872D-03 !kcal mol-1 K-1 
real(dp), parameter :: Rjoule    =   8.31447215D00 !J mol-1 K-1 
real(DP), parameter :: u2b             =  8.0d0*pi*pi
real(DP), parameter :: b2u             =  one/u2b
real(DP), parameter :: cal2joule       =  4.184d0
real(DP), parameter :: joule2cal       =  one/cal2joule
end module 
