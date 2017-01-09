module mod_types
use cfml_atom_typedef,              only: atom_type, atom_list_type
!init_atom_type
use cfml_globaldeps,                 only: sp,dp

!cdyn_atom_type
type, extends(atom_type) :: Protein_Atom_Type
  character(len=25)  :: resn=" "
  integer            :: ires,iblock,iaunit,ichain
end type Protein_Atom_Type

Type, public :: Protein_Atom_List_Type
    integer                                  :: natoms
    type(Protein_Atom_Type),dimension(:),allocatable :: atom
End type Protein_Atom_List_Type

!b=sum_r/ndat - sum_c/ndat
type :: Neigh_List_Type
  integer                                  :: nneighs
  type(protein_atom_list_Type),dimension(:),allocatable :: neigh
  integer,dimension(2)                     :: aneigh,bneigh,cneigh ! range of neighbors 
                                              ! ie. -1 to 1 etc.. unit cell start is always
                                              ! zero so zero must be in series
end type Neigh_list_type

type :: uc_Type ! need to remove this!
  integer                       :: naunits,natoms
  type(Protein_Atom_Type), dimension(:), allocatable  :: atom 
  real(dp), dimension(3)        :: avec,bvec,cvec ! cartesian lattice vectors
  integer,dimension(2)          :: aneigh,bneigh,cneigh ! range of neighbors 
  integer                       :: nneighs ! total neighbors within cutoff
end type uc_type

type :: asym_List_Type
  integer                                  :: naunits,natoms
  integer                                  :: nblocks ! contained in all asyms 
  type(Protein_Atom_list_Type),dimension(:),allocatable :: aunit 
  real(dp), dimension(3) :: avec,bvec,cvec ! cartesian lattice vectors
  integer,dimension(2)   :: aneigh,bneigh,cneigh ! range of neighbors 
  integer                :: nneighs ! total neighbors within cutoff
  real(sp), allocatable,dimension(:,:)   :: block_cents
end type asym_list_type

type :: uc_list_type
  integer                                  :: nneighs
  type(asym_list_Type),dimension(:),allocatable :: neigh
  integer,dimension(2)                     :: aneigh,bneigh,cneigh ! range of neighbors 
end type uc_list_type

type :: isoaniso
  integer  :: nviab ! number of atoms with occupancy of 1.0
  integer  :: naniso ! number of atoms with aniso entries
  integer  :: ndotcc   ! number of atoms with aniso entries < 0.5
  real(dp) :: isocorr,avgdotab,anisocorr  ! correlation between exper and theor Bs
  real(dp) :: avganisoexp,avganisothr,avgccmod,avgsuij ! see merrit ACTA CRys D55,1997 
  real(dp), dimension(:),  allocatable   :: exp_flucu, thr_flucu, exp_aniso,thr_aniso  
  integer,  dimension(:),  allocatable   :: resid ! number of the residue with same indexing
  real(dp), dimension(:,:),allocatable   :: atcorr
       ! flucuations in angs^2 and aniso unitless
end type isoaniso

type :: csr
!CSR Sparse Matrix Storage type see mkl_ug
  integer                                  :: nnzero
  real(dp),dimension(:),allocatable        :: values
  integer ,dimension(:),allocatable        :: columns
  integer ,dimension(:),allocatable        :: pointerB
  integer ,dimension(:),allocatable        :: pointerE
end type csr 

type :: sparse
! generc sparse Matrix Storage Format see mkl_ug
! typical use is for hessian and kirchoff, so added in distance array to make
! life faster and easier
  integer                                  :: nnzero,ndim,nmulti,nbonds,rdim,cdim ! ndim for square 
     ! nmulti is for multiple ints
  real(dp),dimension(:),       allocatable :: values
  complex(kind=8),dimension(:),allocatable :: zvalues
  real(dp),dimension(:),       allocatable :: distance 
  real(dp),dimension(:,:),     allocatable :: tab_dxyz 
  integer ,dimension(:),       allocatable :: columns
  integer ,dimension(:),       allocatable :: rows
  real(sp)                                 :: sparsity 
end type sparse

type :: inp_par
  character(len=25)  :: gnm="gnm",enm="enm",eigtyp
  character(len=25)  :: coord="cart"
  character(len=25)  :: tpbc,runtype,genm,weigh,fileroot,bztyp,ctyp
  character(len=25)  :: fctyp ! heavy, hinsen, power
  character(len=25)  :: atom_anl = "CA", isovcov_scale = "lsq"! atoms used for analysis
  integer            :: fcpower
  real(dp)           :: rcut_start, rcut_end, drcut,perc,fcnst
  real(dp)           :: distcorr_max
  integer            :: nviab       ! number of viable iso
  integer            :: nfrqs       ! number of frqs used, arpack
  integer            :: maxfrqs=1200 ! maximum number of frqs, arpack
  integer            :: print_unit  
  integer            :: print_level 
  integer            :: nnzeros,natoms
  integer            :: first_mode,last_mode
  integer            :: multiplicity
  integer            :: dumbo
  integer            :: qdiv
  integer            :: maxatinblock = 999999
  integer            :: bnm
  integer            :: nblocks = 1 ! number of blocks per asym unit
  logical            :: block_read = .false.
  integer            :: nbz ! number of brilluoin zone vectors sampled
  logical            :: aniso,dosgo,multint
  real(sp)           :: time 
  integer            :: aunit=1  ! for tracking which asym unit I'm dealing with
  real(dp)           :: sparsity,gnm_sparsity, enm_sparsity 
  real(dp)           :: qmax,qvec(3),fqvec(3),animate_mag
  integer            :: animate_nframes
  real(dp)           :: isocorr_cut,isocorr_gamm,isocorr_c0
  real(sp), allocatable,dimension(:,:) :: block_cents
  real(dp)           :: temperature
  integer            :: hihfkikflilf(6)
  real(dp)           :: dh=1.0d0,dk=1.0d0,dl=1.0d0
  !real(dp)           :: dh=0.0d0,dk=0.0d0,dl=0.0d0
 
  integer            :: qdir(3)
end type inp_par

type :: dos
  integer                                  :: nfrqs,nbins,total_frqs 
          !nfrqs: number of frequencies in histo, total_frqs: total number frqs in system
  real(dp),dimension(:),allocatable        :: coor(:),freqs(:),fgw(:),gw(:),bigw(:)
  real(dp)                                 :: maxfrq,delta_frq,cell_vol
  integer, dimension(:),allocatable        :: histo(:)
end type dos

contains

subroutine init_protein_atom_type (A)
    type (Protein_Atom_Type), intent(in out)   :: A
    A%resn = " "
    A%ires = 0
    A%iblock = 0
    A%Lab      =" "
    A%ChemSymb =" "
    A%SfacSymb =" "
    A%Wyck     ="."
    A%Active   =.true.
    A%Z        =0
    A%Mult     =0
    A%X        =0.0
    A%X_Std    =0.0
    A%MX       =0.0
    A%LX       =0
    A%Occ      =0.0
    A%Occ_Std  =0.0
    A%MOcc     =0.0
    A%LOcc     =0
    A%Biso     =0.0
    A%Biso_std =0.0
    A%MBiso    =0.0
    A%LBiso    =0
    A%Utype    ="none"
    A%ThType   ="isotr"
    A%U        =0.0
    A%U_std    =0.0
    A%Ueq      =0.0
    A%MU       =0.0
    A%LU       =0
    A%Charge   =0.0
    A%Moment   =0.0
    A%Ind      =0
    A%NVar     =0
    A%VarF     =0.0
    A%LVarF    =0
    A%mVarF    =0.0
    A%AtmInfo  ="None"
    return
end subroutine init_protein_atom_type

Subroutine Allocate_protein_Atom_List(N,A)
   !---- Arguments ----!
   integer,               intent(in)       :: n  
   type (protein_atom_list_type), intent(in out)   :: A  

   !---- Local Variables ----!
   integer :: i,ier

   A%natoms = n
   if (allocated(A%Atom)) deallocate(A%Atom)
   allocate (A%atom(n),stat=ier)

   do i=1,n
      call init_protein_atom_type(A%atom(i))
   end do

   return
End Subroutine Allocate_protein_atom_list

subroutine inflate_atom_list (a,p)
! there must be a better way than this...
    type (atom_list_type), intent(in) :: a
    type (protein_atom_list_type), intent(out) :: p
    integer :: i

    call Allocate_protein_Atom_List(a%natoms, p)

    do i = 1, p%natoms
        p%atom(i)%lab      = a%atom(i)%Lab 
        p%atom(i)%ChemSymb = a%atom(i)%ChemSymb 
        p%atom(i)%SfacSymb = a%atom(i)%SfacSymb
        p%atom(i)%Wyck     = a%atom(i)%Wyck
        p%atom(i)%Active   = a%atom(i)%Active 
        p%atom(i)%Z        = a%atom(i)%Z       
        p%atom(i)%Mult     = a%atom(i)%Mult    
        p%atom(i)%X        = a%atom(i)%X       
        p%atom(i)%X_Std    = a%atom(i)%X_Std   
        p%atom(i)%MX       = a%atom(i)%MX      
        p%atom(i)%LX       = a%atom(i)%LX      
        p%atom(i)%Occ      = a%atom(i)%Occ     
        p%atom(i)%Occ_Std  = a%atom(i)%Occ_Std 
        p%atom(i)%MOcc     = a%atom(i)%MOcc    
        p%atom(i)%LOcc     = a%atom(i)%LOcc    
        p%atom(i)%Biso     = a%atom(i)%Biso    
        p%atom(i)%Biso_std = a%atom(i)%Biso_std
        p%atom(i)%MBiso    = a%atom(i)%MBiso   
        p%atom(i)%LBiso    = a%atom(i)%LBiso   
        p%atom(i)%Utype    = a%atom(i)%Utype   
        p%atom(i)%ThType   = a%atom(i)%ThType  
        p%atom(i)%U        = a%atom(i)%U       
        p%atom(i)%U_std    = a%atom(i)%U_std   
        p%atom(i)%Ueq      = a%atom(i)%Ueq     
        p%atom(i)%MU       = a%atom(i)%MU      
        p%atom(i)%LU       = a%atom(i)%LU      
        p%atom(i)%Charge   = a%atom(i)%Charge  
        p%atom(i)%Moment   = a%atom(i)%Moment  
        p%atom(i)%Ind      = a%atom(i)%Ind     
        p%atom(i)%NVar     = a%atom(i)%NVar    
        p%atom(i)%VarF     = a%atom(i)%VarF    
        p%atom(i)%LVarF    = a%atom(i)%LVarF   
        p%atom(i)%mVarF    = a%atom(i)%mVarF   
        p%atom(i)%AtmInfo  = a%atom(i)%AtmInfo 
    end do               

end subroutine inflate_atom_list

subroutine dos_append(freqs,states)
! grows the array of frequencies
real(dp) , allocatable, intent(in)     :: freqs(:)
type(dos),              intent(in out) :: states
real(dp) , allocatable                 :: tmparray(:)
integer :: ialloc


!first, if allocated, copy states%freqs to tmparray
if (allocated(states%freqs)) then

  allocate(tmparray(states%nfrqs), stat=ialloc)
  if (ialloc /= 0) STOP "*** not enough memory ***"
  tmparray=states%freqs
  deallocate(states%freqs)

end if

! increment number of frqs by number to be added
states%nfrqs = states%nfrqs+size(freqs)

! allocate states%freqs; at this point states%freqs is not allocated 
allocate(states%freqs(states%nfrqs),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

if (allocated(tmparray)) then

  states%freqs(1:size(tmparray))              = tmparray(:)
  states%freqs(size(tmparray)+1:states%nfrqs) = freqs(:)

else

  states%freqs(:) = freqs(:)

end if

end subroutine

subroutine dos_statalloc(states)
! grows the array of frequencies
type(dos),              intent(in) :: states

if (allocated(states%freqs)) print *, "nfrqs", size(states%freqs), states%nfrqs 
if (allocated(states%histo)) print *, "nhisto", size(states%histo)
if (allocated(states%gw)) print *, "ngw", size(states%gw)
if (allocated(states%fgw)) print *, "nfgw", size(states%fgw)
if (allocated(states%bigw)) print *, "nbigw", size(states%bigw)
if (allocated(states%coor)) print *, "ncoor", size(states%coor)

end subroutine

subroutine dos_histgen(states)
! generate histogram, and calculate density of states
use mod_math, only: trapezoidal_int
type(dos), intent(in out) :: states
real(dp) :: tot
integer :: i,j,val,nbins

states%maxfrq    = maxval(states%freqs)
states%delta_frq = states%maxfrq/real(states%nbins-1,kind=dp) ! there're nbin - 1 nonzeroes
nbins = states%nbins

if (allocated(states%histo)) deallocate(states%histo)
if (allocated(states%gw))    deallocate(states%gw)
if (allocated(states%fgw))    deallocate(states%fgw)
if (allocated(states%bigw))    deallocate(states%bigw)
if (allocated(states%coor))    deallocate(states%coor)
allocate(states%histo(nbins),states%gw(nbins),states%fgw(nbins), &
         states%coor(nbins),states%bigw(nbins), stat=ialloc)
if (ialloc /= 0) STOP "dos_histgen> malloc"

states%histo =0
states%gw = zero 
states%fgw = zero
states%bigw = zero

do i=1,states%nfrqs
  val               = NINT((states%freqs(i))/states%delta_frq)+1 ! the 1 is to avoid zero
  !print *, states%freqs(i),states%delta_frq,NINT((states%freqs(i))/states%delta_frq)+1,val, size(states%histo)
  states%histo(val) = states%histo(val)+1
end do

tot=0

do i=1,states%nbins
  states%coor(i) = real(i-1)*states%delta_frq
  tot            = tot + real(states%histo(i))/real(states%nfrqs)
! add up the fractions of bigG
  states%bigw(i) = tot
  states%gw(i)   = real(states%histo(i))/(real(states%nfrqs)*states%delta_frq)
  states%fgw(i)  = real(states%histo(i))/(states%delta_frq)
end do

! don't need to integrate-integrate function gw to get BiGw because we're treating gw
! via histograms which gives us the total number of freqs within a deltafreq bin
! multiply by deltafreq and you have the number of freqs!  add them up and you have
! biGw
!call trapezoidal_int(states%coor,states%gw,states%bigw)

end subroutine

subroutine dos_hist_regen(constant,states)
! multiply freqs by constant and regenerate the histogram and dos 
use mod_math, only: trapezoidal_int
type(dos), intent(in out) :: states
real(dp) , intent(in) :: constant
real(dp) :: tot
integer :: i,j,val

states%freqs = states%freqs*constant
states%maxfrq    = maxval(states%freqs)
states%delta_frq = states%maxfrq/real(states%nbins-1,kind=dp) ! there're nbin - 1 nonzeroes

states%histo = 0

do i=1,states%nfrqs
  val               = NINT((states%freqs(i))/states%delta_frq)+1 ! the 1 is to avoid zero
  !print *, states%freqs(i),states%delta_frq,NINT((states%freqs(i))/states%delta_frq)+1,val, size(states%histo)
  states%histo(val) = states%histo(val)+1
end do

tot=0

do i=1,states%nbins
  states%coor(i) = real(i-1)*states%delta_frq
  tot            = tot + real(states%histo(i))/real(states%nfrqs)
! add up the fractions of bigG
  states%bigw(i) = tot
  states%gw(i)   = real(states%histo(i))/(real(states%nfrqs)*states%delta_frq)
  states%fgw(i)  = real(states%histo(i))/(states%delta_frq)
end do

! don't need to integrate-integrate function gw to get BiGw because we're treating gw
! via histograms which gives us the total number of freqs within a deltafreq bin
! multiply by deltafreq and you have the number of freqs!  add them up and you have
! biGw
!call trapezoidal_int(states%coor,states%gw,states%bigw)

end subroutine

subroutine dos_init(nbins,states)
integer  , intent(in)  :: nbins
type(dos), intent(out) :: states
integer :: ialloc

states%nfrqs     = 0
states%nbins     = nbins
states%maxfrq    = 0.0d0
states%delta_frq = 0.0d0

if (allocated(states%histo)) deallocate(states%histo)
if (allocated(states%gw))    deallocate(states%gw)
if (allocated(states%fgw))    deallocate(states%fgw)
allocate(states%histo(nbins),states%gw(nbins),states%fgw(nbins), &
         states%coor(nbins),states%bigw(nbins), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

states%histo = 0
states%gw    = zero

end subroutine

subroutine dos_analysis(input,states)
!DMR add in some curve fitting stuff
use mod_curvefit
use mod_constants , only: rkcal, hkcal, convfrq,clight
integer, parameter :: max_temp=400
type(inp_par), intent(in) :: input
type(dos),     intent(in) :: states
real(dp), allocatable, dimension(:) :: temps, cv
integer :: i,jtemp,ier
real(dp) :: hfreq,kt,bestpart,hvkt,maxperc,firstperc,deltaperc
type(datafunc) :: adfunc,gdfunc


allocate(temps(max_temp),cv(max_temp), stat=ier)

if (input%print_level .ge. 2) then
      open(unit=22,file=trim(input%fileroot)//"-"//trim(input%genm)//"-"//trim(input%tpbc)//"-"// &
                        trim(input%fctyp) //"-"//trim(input%bztyp)//"-dos.txt", &
      status='unknown',action='write', iostat=ier)

      open(unit=23,file=trim(input%fileroot)//"-"//trim(input%genm)//"-"//trim(input%tpbc)//"-"// &
                        trim(input%fctyp) //"-"//trim(input%bztyp)//"-cv.txt", &
      status='unknown',action='write', iostat=ier)
end if
! write out histogram
if (input%print_level .ge. 2) then
  do i=1,states%nbins
    write (22,'(I4,I7,F10.5,F10.5,F10.5,F15.5)') i,states%histo(i),  &
         states%coor(i),states%gw(i),states%biGw(i) ,states%fgw(i)
  end do
end if

! calculate heat capacity
! something is just not right
    cv=zero
  do jtemp = 1, 400
    temps(jtemp) = real(jtemp,kind=dp)
! start from 2 to avoid the zero
    do i=2,states%nbins !DMR
      hfreq    = hkcal * states%coor(i) * 1.0D12
      !hfreq    = hkcal * real(i)*states%delta_frq * 1.0D12
      kt       = real(jtemp)*rkcal
      bestpart = one/(two*cosh(hfreq/kt) - two )
      cv(jtemp)    = cv(jtemp) + &
                     (one/(kt*real(jtemp)))*states%gw(i)*(hfreq**2)*bestpart*states%delta_frq
    end do
    if (cv(jtemp) .le. 1.0D-16) cv(jtemp)=zero
if (input%print_level .ge. 2) then
    write(23,'(F10.2,E15.5)') real(jtemp), Cv(jtemp)
end if

  end do

close(22)
close(23)

print '(A12)', "BEGINFIT CV>"
print '(A3,2A6,A10,A6,2A10,A6,A10,A6,A10,A6)', "Cv>","frst%","last%","a","b","msd", &
                             "a_log","b_log","msd_log", "xcent","yshft", "#pnts"

maxperc = 0.1d0
firstperc = 1.0d-02
deltaperc = 5.0d-02

call data_adjust(temps,cv,rkcal,firstperc,maxperc,adfunc)
call minsearch(adfunc)
call loglogfit(adfunc)

print '(A3,2F6.2,E10.3,F6.2,2E10.3,F6.2,E10.3,F6.2,E10.3,I6)', "Cv>", firstperc,maxperc, &
                                        adfunc%a, adfunc%b,adfunc%msd, &
                                        adfunc%a_log,adfunc%b_log, adfunc%msd_log, &
                                        adfunc%acoor_offset,adfunc%func_offset, &
                                        size(adfunc%func)


print '(A10)', "ENDFIT CV>"
print *, " "
print '(A12)', "BEGINFIT GW>"
print '(A3,2A6,A10,A6,2A10,A6,A10,A6,A10,A6)', "Gw>","frst%","last%","a","b","msd", &
                             "a_log","b_log","msd_log", "xcent","yshft", "#pnts"

maxperc = 0.2d0
firstperc = 3.0d-02

call data_adjust(states%coor,states%biGw,one,firstperc,maxperc,gdfunc)
call minsearch(gdfunc)
call loglogfit(gdfunc)

print '(A3,2F6.2,E10.3,F6.2,2E10.3,F6.2,E10.3,F6.2,E10.3,I6)', "Gw>",firstperc,maxperc, &
                                          gdfunc%a, gdfunc%b,gdfunc%msd, &
                                          gdfunc%a_log,gdfunc%b_log, gdfunc%msd_log, &
                                          gdfunc%acoor_offset,gdfunc%func_offset, &
                                          size(gdfunc%func)


print '(A10)', "ENDFIT GW>"


end subroutine

end module

