module mod_crysbuild
! DMR 01-11-2008
! DMR update: 04-30-2009 ...  now images has recentered unit cell construction
!                             while asyms_im has no recentering!  more molecular in nature
! subroutines that construct the contents of the unit cell and a chunk of crystal that,
! for now, it the unit cell and its nearest neighbors.  Mainly driver routines for the 
! crysFML library.
!
!Acta Cryst. (2005). A61, C22
!CrysFML: a Library to Develop Crystallographic Programs in
!Fortran 95. Powder diffraction examples
!Juan Rodríguez-Carvajala, Javier González-Platasb, aLaboratoire Léon
!Brillouin (CEA-CNRS), CEA/Saclay France. bDepto. Física
!Fundamental, Universidad de La Laguna, Tenerife, Spain. E-mail:
!Juan.Rodriguez-Carvajal@cea.fr
!CrysFML is a set of Fortran 95 modules to be used in
!crystallographic and diffraction computing programs [1]. Modern
!array syntax and new features of Fortran 95 are used through the
!modules. We take advantage of all object oriented programming
!techniques already available in Fortran (user-defined types,
!encapsulation, overload of procedures and functions). The lacking
!features (e.g. inheritance and class methods) will be easily
!implemented as soon as compilers of the new standard become
!available. We aim to preserve the efficiency, the simplicity and the
!adequacy of modern Fortran for numerical calculations. All aspects of
!symmetry and handling of reflections and structure factor calculations
!are treated in dedicated modules. Main programs using the adequate
!modules may perform more or less complicated calculations with only
!few lines of code. The documentation is written in the source code. A
!document, in HTML format can be generated using a program.
!We shall present an overview of the present status of the library
!and a series of examples useful for powder diffraction: simple
!crystallographic calculations, bond-valence sums, aids to space group
!determination, profile fitting, powder diffraction simulations, kernel of
!the Rietveld method, etc.
![1] Rodríguez-Carvajal J., González-Platas J., Compcomm Newsletter 2003, 1,50.
!
!Subroutines
!  unit_cell_gen(cell,asym_un,SpG,asym_rcnt,unit_cell)
!  image_gen(cell,unit_cell,images)
!  cart_to_fract (cell, atoms_in)
!  fract_to_cart (cell, atoms_in)
!

! crysFML modules
use cfml_crystallographic_symmetry,  only: space_group_type, get_orbit,applyso
use cfml_Scattering_Chemical_Tables, only: Get_Atomic_Mass
!use cfml_Atom_typedef,                only: Atom_List_Type, init_atom_type
use cfml_crystal_metrics,              only: Crystal_Cell_Type
use cfml_globaldeps,                   only: sp,dp

! external modules
use mod_types,                  only: uc_list_type,Neigh_List_Type,asym_List_Type, init_protein_atom_type, protein_atom_list_type 
use mod_math,                   only: dist_sqr
use mod_inout,                  only: fileroot

implicit none

private :: unit_cell_gen, image_gen
public ::  crys_build_new,asymlist_to_atomlist,chunk_gen,&
           mass_vector_build,mass_matrix_build,crys_build, &
           asym_build, cart_to_fract, fract_to_cart,zmass_weightit

Interface     cart_to_fract
   Module Procedure asymlist_cart_to_fract, atomlist_cart_to_fract
End Interface cart_to_fract 

interface     block_assign
  Module Procedure block_assign_atasy, block_assign_asyp1
end interface

Interface     fract_to_cart
   Module Procedure asymlist_fract_to_cart, atomlist_fract_to_cart
End Interface fract_to_cart 

interface mass_vector_build
  module procedure mass_vector_build_uc ,mass_vector_build_atom
end interface mass_vector_build


contains 

subroutine mass_weightit(genm,atoms,matin,matout,zd)
use mod_types,     only: sparse
use mod_inout,     only: sparse_refinit
character(len=3),    intent(in)      :: genm
type(protein_atom_list_type), intent(in)     :: atoms
type(sparse),intent(in)              :: matin
type(sparse),intent(out)             :: matout
type(sparse)                         :: rtmasses
character(len=1)                     :: zd
integer                              :: i,j,nnzero,ialloc,ndim

call mass_matrix_build(atoms,rtmasses,genm)
call sparse_refinit(matin,matout,zd)

do i=1,matin%nnzero

if (zd .eq. "d") then
    matout%values(i) = rtmasses%values(matin%columns(i)) &
                    *matin%values(i)* &
                    rtmasses%values(matin%rows(i))
else if (zd .eq. "z") then
    matout%zvalues(i) = rtmasses%values(matin%columns(i)) &
                    *matin%zvalues(i)* &
                    rtmasses%values(matin%rows(i))

end if

end do

end subroutine

subroutine zmass_weightit(genm,atoms,matin)
use mod_types,     only: sparse
use mod_inout,     only: sparse_refinit
character(len=3),         intent(in)     :: genm
type(protein_atom_list_type),     intent(in)     :: atoms
complex(dp) ,allocatable, intent(inout)  :: matin(:,:)
type(sparse)                         :: rtmasses
integer                              :: i,j,k,l,ndim,nat

ndim = size(matin,1)

call mass_matrix_build(atoms,rtmasses,genm)

do i=1,ndim
  do j=i,ndim
    matin(i,j) = rtmasses%values(i)*matin(i,j)*rtmasses%values(j)
    if (i.ne.j) matin(j,i) = rtmasses%values(j)*matin(j,i)*rtmasses%values(i)
  end do
end do

end subroutine

subroutine mass_vector_build_atom (atom_list,masses)
use mod_types,     only : sparse
use mod_constants, only : one,zero
! returns U sparse-matrix with masses along diag
type(protein_atom_list_type),  intent(in)  :: atom_list
real(dp), allocatable, intent(out) :: masses(:)
integer                            :: nrepeat
real(sp)                           :: mass
integer                            :: i,j,k,ndim, ialloc

ndim=atom_list%natoms

allocate(masses(ndim),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

masses = zero

do i=1,atom_list%natoms
call Get_Atomic_Mass(atom_list%atom(i)%chemsymb,mass)
masses(i)  = real(mass,kind=dp)
end do

end subroutine

subroutine mass_vector_build_uc (atom_list,masses)
use mod_types,     only : sparse,uc_type
use mod_constants, only : one,zero
! returns U sparse-matrix with masses along diag
type(uc_type),  intent(in)  :: atom_list
real(dp), allocatable, intent(out) :: masses(:)
integer                            :: nrepeat
real(sp)                           :: mass
integer                            :: i,j,k,ndim, ialloc

ndim=atom_list%natoms

allocate(masses(ndim),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

masses = zero

do i=1,atom_list%natoms
call Get_Atomic_Mass(atom_list%atom(i)%chemsymb,mass)
masses(i)  = real(mass,kind=dp)
end do

end subroutine

subroutine invrt_rtmass_vector (atom_list,inv_rtmass)
use mod_types,     only : sparse
use mod_constants, only : one,zero
! returns U sparse-matrix with masses along diag
type(protein_atom_list_type),  intent(in)  :: atom_list
real(dp), allocatable, intent(out) :: inv_rtmass(:)
real(sp)                           :: mass
integer                            :: i,j,k,ndim, ialloc

ndim=atom_list%natoms

call mass_vector_build (atom_list,inv_rtmass)

do i=1,atom_list%natoms
inv_rtmass(i)  = one/dsqrt(inv_rtmass(i))
end do

end subroutine

subroutine atvec_COM(atom_list,com)
! center of mass calc
use mod_constants, only : one,zero
type(protein_atom_list_type),  intent(in)   :: atom_list
real(dp),              intent(out)  :: com(3)
real(dp), allocatable               :: masses(:)
real(dp) :: sum_mass
integer   :: i

com = zero
call mass_vector_build (atom_list,masses)
sum_mass = sum(masses) 

do i =1, atom_list%natoms
com = com + atom_list%atom(i)%x*masses(i)
end do

com = com/sum_mass

end subroutine

subroutine mass_matrix_build (atom_list,masses,genm)
use mod_types,     only : sparse
use mod_constants, only : one
! returns U sparse-matrix with masses along diag
type(protein_atom_list_type), intent(in)  :: atom_list
type(sparse),         intent(out) :: masses
character(len=3),     intent(in)  :: genm
integer                           :: nrepeat
real(sp)                          :: mass
integer                           :: i,j,k,nnzero, ialloc

if (genm .eq. "enm") then
nrepeat=3
else 
nrepeat=1
end if

nnzero=nrepeat*atom_list%natoms

masses%nnzero=nnzero

allocate(masses%values(nnzero), &
        masses%columns(nnzero), masses%rows(nnzero),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

k=1
do i=1,atom_list%natoms
call Get_Atomic_Mass(atom_list%atom(i)%chemsymb,mass)
do j=0,nrepeat-1
    masses%values(k+j)  = one/dsqrt(real(mass,kind=dp))
    masses%rows(k+j)    = k+j
    masses%columns(k+j) = k+j
end do
k=k+nrepeat
end do

end subroutine

subroutine crys_build_new(input,cell,spg,asym_un,asyms,unit_cell)
use mod_types , only: inp_par, asym_list_type, uc_type 
! build p1 unit cell of crystal from the asymetric unit 
! takes in cartesian coordinates 
! operates in fractional coords
! returns cartesian coords
!   asyms:     list of all atoms p1 unit cell indexed by asymmetric unit
!   unit_cell: list of atoms in  p1 unit cell without indexing
!
!   both of these data structures contain information for cartesian lattice displacement
!   for use with periodic boundary conditions and visualization
!
!   this subroutine will work best for biological units that don't contain the more minimal
!   assymetric units found in small molecule crystal like urea:
!         asym u: OCNH2   molecule: OCN2H4
!   another subroutine for this case should be fixed up (there's one that works already, but
!     it needs to be explicitly written as such) 
type(inp_par)          , intent(in out)     :: input
type (crystal_cell_type),intent(in)     :: Cell
type (space_group_type), intent(in)     :: SpG
type (protein_atom_list_type),   intent(in out) :: asym_un ! inout tofrom fractional coords  
type (asym_list_type),   intent(out)    :: asyms 
type (asym_list_type),   intent(out)    :: unit_cell
real(dp)                                :: rcut
integer                                 :: i,j

! convert asym unit from cartesian to fractional space
call cart_to_fract(cell,asym_un)

! generate recentered asym unit and full unit cell
call asym_gen(input,cell,asym_un,spg,asyms)
call asyms_to_p1asym(asyms,unit_cell)

! convert pieces back to cartesian coords
call fract_to_cart(cell,asym_un)
call fract_to_cart(cell,asyms)
call fract_to_cart(cell,unit_cell)

end subroutine

subroutine crys_build(input,cell,spg,asym_un,images)
use mod_types , only: inp_par, neigh_list_type
! USES RECENTERING
! build hunk of crystal from the asymetric unit 
! takes in cartesian coordinates 
! operates in fractional coords
! returns cartesian coords
! 
! todo: there is a bit of redundancy and memory wasting here.  If it proved
!       necessary, the images datastructure could be defined solely in memory
!       and then the other pieces could be defined via pointers to the relavent
!       pieces of images.  not needed for now
type(inp_par)          , intent(in)     :: input
type (crystal_cell_type),intent(in)     :: Cell
type (space_group_type), intent(in)     :: SpG
type (protein_atom_list_type),   intent(in out) :: asym_un ! inout tofrom fractional coords  
type (neigh_list_type),  intent(out)    :: images 
type (protein_atom_list_type)                   :: asym_rcnt, unit_cell
real(dp)                                :: rcut
integer                                 :: i,j

! convert asym unit from cartesian to fractional space
call cart_to_fract(cell,asym_un)

! generate recentered asym unit and full unit cell
call unit_cell_gen(input,cell,asym_un,SpG,asym_rcnt,unit_cell)
call image_gen(input,cell,unit_cell,images)

! convert pieces back to cartesian coords
call fract_to_cart(cell,asym_un)
do i=0,images%nneighs
call fract_to_cart(cell,images%neigh(i))
end do

deallocate(asym_rcnt%atom)
deallocate(unit_cell%atom)

end subroutine

subroutine asym_build(input,cell,spg,asym_un,asyms_im)
use mod_types,  only: inp_par, uc_list_type
! build hunk of crystal from the asymetric unit 
! takes in cartesian coordinates 
! operates in fractional coords
! returns cartesian coords

type(inp_par)          , intent(in out)     :: input !passed to asym_gen
type (Crystal_Cell_Type),intent(in)     :: Cell
type (space_group_type), intent(in)     :: SpG
type (protein_atom_list_type),   intent(in out) :: asym_un ! needs to be convert to and from fractional coords  
type (protein_atom_list_type)                   :: asym_rcnt,unit_cell
type (uc_list_Type),  intent(out)       :: asyms_im 
type (asym_list_Type)                  :: asyms 
real(dp)                                :: rcut  !passed to asym_gen
integer                                 :: i,j


rcut = input%rcut_end

! convert asym unit from cartesian to fractional space
call cart_to_fract(cell,asym_un)

! generate recentered asym unit and full unit cell
call piecewise_uc_gen(input,cell,asym_un,SpG,unit_cell)
deallocate(unit_cell%atom)
! generate hunk of crystal
if( (trim(input%tpbc) .eq. "asy") .or. & 
    (trim(input%tpbc) .eq. "pbc")  .or. &
    (trim(input%tpbc) .eq. "bvk") ) then
call asym_gen(input,cell,asym_un,SpG,asyms)
else
call nonasym_gen(rcut,cell,asym_un,SpG,asyms)
end if
call asym_im_gen(input,cell,asyms,asyms_im)

! convert pieces back to cartesian coords
call fract_to_cart(cell,asym_un)

do i=0,asyms_im%nneighs
do j=1,asyms_im%neigh(i)%naunits
    call fract_to_cart(cell,asyms_im%neigh(i)%aunit(j))
end do
end do

do j=1,asyms%naunits
call fract_to_cart(cell,asyms%aunit(j))
end do

end subroutine

subroutine unit_cell_gen(input,cell,asym_un,SpG,asym_rcnt,unit_cell)
use mod_types , only: inp_par
use cfml_math_general, only: modulo_lat
! 
! build unit cell with RECENTERING from the asymetric unit with fractional coordinates
! returns recentered asym_unit and unit cell with fractional coordinates
! DMR: 4-13-09  fixed for special sym positions
type (inp_par),          intent(in)     :: input
type (Crystal_Cell_Type),intent(in)     :: Cell
type (protein_atom_list_type),   intent(in out)     :: asym_un 
type (space_group_type), intent(in)     :: SpG
type (protein_atom_list_type),   intent(out)    :: asym_rcnt, unit_cell
integer                                 :: i,j,jx,jy,jz,nat
integer                                 :: nmultip,ier,natoms
real,     allocatable                   :: orb(:,:) ! for get_orbit subroutine 
real(dp)                                :: dist,dxyz(3),xyzi(3),xyzj(3),sumfrc(3)

!initiate recentered asym unit
asym_rcnt%natoms=asym_un%natoms
if (allocated(asym_rcnt%atom)) deallocate(asym_rcnt%atom)
allocate (asym_rcnt%atom(asym_rcnt%natoms),stat=ier)

do i=1,asym_rcnt%natoms
asym_rcnt%atom(i)=asym_un%atom(i)
asym_rcnt%atom(i)%x = modulo_lat(asym_un%atom(i)%x)
end do

call number_of_p1atoms(asym_rcnt,spg,natoms)

! initialize unit cell 
unit_cell%natoms=natoms

if (allocated(unit_cell%atom)) deallocate(unit_cell%atom)
allocate (unit_cell%atom(unit_cell%natoms),stat=ier)
do i=1,unit_cell%natoms
call init_protein_atom_type(unit_cell%atom(i))
end do

allocate(orb(3,SpG%multip), stat=ier)
! what follow works with small molecules with special spots
! problem is that is switches up the ordering
! fill up unit cell
!nat=0
!do i=1, asym_rcnt%natoms
!  call Get_Orbit(asym_rcnt%atom(i)%x(:),SpG,nmultip,orb)
!call Get_Orbit_nocent(asym_rcnt%atom(i)%x(:),SpG,nmultip,orb)
!  do j=1,nmultip
!    nat=nat+1
!    unit_cell%atom(nat)        = asym_rcnt%atom(i)
!    unit_cell%atom(nat)%x      = orb(:,j)
!  end do
!end do

nat=0
do j=1,SpG%multip
do i=1, asym_rcnt%natoms
    call Get_Orbit(asym_rcnt%atom(i)%x(:),SpG,nmultip,orb)
    nat=nat+1
    unit_cell%atom(nat)        = asym_rcnt%atom(i)
    unit_cell%atom(nat)%x      = orb(:,j)
end do
end do



deallocate(orb)
!test if any atoms are overlapping.
call fract_to_cart(cell,unit_cell)
do j=1,unit_cell%natoms-1
do i=j+1,unit_cell%natoms
    dxyz(:) = unit_cell%atom(i)%x(:) - unit_cell%atom(j)%x(:)    
    dist=dist_sqr(dxyz)
    if (dist .le. 0.5d0) then
            print '(A16,1x,I5,1x,I5,1x,A1)','overlap atoms', i,j, ':'
            print '(A2,1x,I5,1x,A3,1x,3F12.5)', 'i:',i, &
                unit_cell%atom(i)%lab,unit_cell%atom(i)%x(:)
            print '(A2,1x,I5,1x,A3,1x,3F12.5)', 'j:',j, &
                unit_cell%atom(j)%lab,unit_cell%atom(j)%x(:)
    end if
end do  
end do    

call cart_to_fract(cell,unit_cell)

end subroutine

subroutine piecewise_uc_gen(input,cell,asym_un,SpG,unit_cell)
use mod_types , only: inp_par
use cfml_math_general, only: modulo_lat
! build unit cell from the asymetric unit with fractional coordinates
! returns recentered asym_unit and unit cell with fractional coordinates
! DMR: 4-13-09  fixed for special sym positions
type (inp_par),          intent(in)     :: input
type (Crystal_Cell_Type),intent(in)     :: Cell
type (protein_atom_list_type),   intent(in)     :: asym_un 
type (space_group_type), intent(in)     :: SpG
type (protein_atom_list_type),   intent(out)    :: unit_cell
integer                                 :: i,j,jx,jy,jz,nat
integer                                 :: nmultip,ier,natoms
real,     allocatable                   :: orb(:,:) ! for get_orbit subroutine 
real(dp)                                :: dist,dxyz(3),xyzi(3),xyzj(3),sumfrc(3)

!initiate recentered asym unit

! determine number of atoms in primitive unit cell
call number_of_p1atoms(asym_un,spg,natoms)
! initialize unit cell 
unit_cell%natoms=natoms

if (allocated(unit_cell%atom)) deallocate(unit_cell%atom)
allocate (unit_cell%atom(unit_cell%natoms),stat=ier)
do i=1,unit_cell%natoms
call init_protein_atom_type(unit_cell%atom(i))
end do

! fill up unit cell
allocate(orb(3,SpG%multip), stat=ier)
nat=0
do i=1, asym_un%natoms
!call Get_Orbit(asym_rcnt%atom(i)%x(:),SpG,nmultip,orb)
call Get_Orbit(asym_un%atom(i)%x(:),SpG,nmultip,orb)
!call Get_Orbit_nocent(asym_un%atom(i)%x(:),SpG,nmultip,orb)
do j=1,nmultip
    nat=nat+1
    unit_cell%atom(nat)        = asym_un%atom(i)
    unit_cell%atom(nat)%x      = orb(:,j)
end do
end do

deallocate(orb)
!test if any atoms are overlapping.
call fract_to_cart(cell,unit_cell)
do j=1,unit_cell%natoms-1
do i=j+1,unit_cell%natoms
    dxyz(:) = unit_cell%atom(i)%x(:) - unit_cell%atom(j)%x(:)    
    dist=dist_sqr(dxyz)
    if (dist .le. 0.5d0) then
            print '(A16,1x,I5,1x,I5,1x,A1)','overlap atoms', i,j, ':'
            print '(A2,1x,I5,1x,A3,1x,3F12.5)', 'i:',i, &
                unit_cell%atom(i)%lab,unit_cell%atom(i)%x(:)
            print '(A2,1x,I5,1x,A3,1x,3F12.5)', 'j:',j, &
                unit_cell%atom(j)%lab,unit_cell%atom(j)%x(:)
    end if
end do  
end do    

call cart_to_fract(cell,unit_cell)

end subroutine

subroutine number_of_p1atoms(asym,spg,natoms)
! determine how many atoms are in the unit cell
! this is NOT simply the multiplicity of space group times natoms: 
!             unit_cell%natoms=asym_rcnt%natoms*Spg%multip
! because there may be special sym positions (especially in systems with higher symm)           
type (protein_atom_list_type),   intent(in)    :: asym
type (space_group_type), intent(in)    :: SpG
integer,                 intent(out)   :: natoms
real,     allocatable                  :: orb(:,:) ! for get_orbit subroutine 
integer :: i,nmultip,ier

allocate(orb(3,SpG%multip), stat=ier)
natoms=0
do i=1, asym%natoms
!  print *, 'atom ', i, ' ', asym%atom(i)%biso
  call Get_Orbit(asym%atom(i)%x(:),SpG,nmultip,orb)
  !call Get_Orbit_nocent(asym%atom(i)%x(:),SpG,nmultip,orb)
!  print *, orb 
  natoms = natoms + nmultip
end do

!print *, 'shitnatoms ', natoms, ' asym natoms ', asym%natoms 

end subroutine

subroutine asymlist_to_uc(asyms,unit_cell)
use mod_types, only : asym_list_Type,uc_type
type (asym_list_Type),   intent(in)   :: asyms
type (uc_Type),          intent(out)  :: unit_cell
integer :: i,j,k,nat,natoms,ier

unit_cell%natoms  = asyms%natoms
unit_cell%avec    = asyms%avec
unit_cell%bvec    = asyms%bvec
unit_cell%cvec    = asyms%cvec
unit_cell%aneigh  = asyms%aneigh
unit_cell%bneigh  = asyms%bneigh
unit_cell%cneigh  = asyms%cneigh
unit_cell%nneighs = asyms%nneighs


if (allocated(unit_cell%atom)) deallocate(unit_cell%atom)
allocate (unit_cell%atom(unit_cell%natoms),stat=ier)
do i=1,unit_cell%natoms
  call init_protein_atom_type(unit_cell%atom(i))
end do

nat = 0
do i = 1,asyms%naunits
  do j = 1,asyms%aunit(1)%natoms
    nat = nat + 1
    unit_cell%atom(nat)=asyms%aunit(i)%atom(j)
    unit_cell%atom(nat)%iaunit = i
  end do
end do

end subroutine

subroutine asymlist_to_atomlist(asyms,unit_cell)
use mod_types, only : asym_list_Type
use mod_inout, only : pdbwriter_new
type (asym_list_Type),   intent(in)   :: asyms
type (protein_atom_list_type),   intent(out)  :: unit_cell
integer :: i,j,k,nat,natoms,ier

unit_cell%natoms = asyms%natoms

if (allocated(unit_cell%atom)) deallocate(unit_cell%atom)
allocate (unit_cell%atom(unit_cell%natoms),stat=ier)
do i=1,unit_cell%natoms
  call init_protein_atom_type(unit_cell%atom(i))
end do


nat = 0
do i = 1,asyms%naunits
  do j = 1,asyms%aunit(1)%natoms
    nat = nat + 1
    unit_cell%atom(nat)=asyms%aunit(i)%atom(j)
    unit_cell%atom(nat)%iaunit = i
  end do
end do

end subroutine

subroutine asymlist_to_uctype(asyms,unit_cell)
use mod_types, only : asym_list_Type
use mod_inout, only : pdbwriter_new
type (asym_list_Type),   intent(in)   :: asyms
type (protein_atom_list_type),   intent(out)  :: unit_cell
integer :: i,j,k,nat,natoms,ier

unit_cell%natoms = asyms%natoms

if (allocated(unit_cell%atom)) deallocate(unit_cell%atom)
allocate (unit_cell%atom(unit_cell%natoms),stat=ier)
do i=1,unit_cell%natoms
  call init_protein_atom_type(unit_cell%atom(i))
end do


nat = 0
do i = 1,asyms%naunits
  do j = 1,asyms%aunit(1)%natoms
    nat = nat + 1
    unit_cell%atom(nat)=asyms%aunit(i)%atom(j)
    unit_cell%atom(nat)%iaunit = i
  end do
end do

end subroutine

subroutine asyms_to_p1asym(asyms,unit_cell)
use mod_types, only : asym_list_Type
use mod_inout, only : pdbwriter_uc_crys,pdbwriter
! this subroutine copies the contents of all asym units
! into one asym unit...  ie.  reducing the symm to p1
type (asym_list_Type),   intent(in)    :: asyms
type (asym_list_Type),   intent(out)   :: unit_cell
integer                                :: i,j,k,ier

! determine number of atoms in primitive unit cell
unit_cell = asyms
if (allocated(unit_cell%aunit)) deallocate(unit_cell%aunit)
unit_cell%naunits = 1
unit_cell%nblocks = asyms%nblocks
allocate(unit_cell%aunit(1),stat = ier)
if (ier /= 0) stop "asyms_to_p1asym> aunit malloc"
unit_cell%aunit(1)%natoms = unit_cell%natoms

if (allocated(unit_cell%aunit(1)%atom)) deallocate(unit_cell%aunit(1)%atom)
allocate (unit_cell%aunit(1)%atom(unit_cell%natoms),stat=ier)
if (ier /= 0) STOP "*** not enough memory ***"

do j=1,unit_cell%natoms
  call init_protein_atom_type(unit_cell%aunit(1)%atom(j))
end do

! copy all asym units into one unit cell
i = 0
do k = 1, asyms%naunits
  do j = 1,asyms%aunit(1)%natoms
    i = i + 1
    unit_cell%aunit(1)%atom(i) = asyms%aunit(k)%atom(j)
  end do
end do


end subroutine

integer function iposinteq(i)
integer, intent(in)  :: i

if (i .lt. 1) then
  iposinteq = 1
else
  iposinteq = i
end if

end function

subroutine p1asym_to_chunk(na,nb,nc,unit_cell,chunk)
use mod_types, only : asym_list_Type
use mod_inout, only : pdbwriter_uc_crys,pdbwriter
! this subroutine copies the contents of all asym units
! into one asym unit...  ie.  reducing the symm to p1
integer, intent(in)                    :: na,nb,nc
type (asym_list_Type),   intent(in)    :: unit_cell
type (asym_list_Type),   intent(out)   :: chunk 
integer                                :: i,j,k,l,ii,jj,ier,ia,ib,ic,nun,maxblock

ia = iposinteq(na)
ib = iposinteq(nb)
ic = iposinteq(nc)

nun = (ia-1+1)*(ib-1+1)*(ic-1+1)

! determine number of atoms in primitive unit cell
chunk = unit_cell 
if (allocated(chunk%aunit)) deallocate(chunk%aunit)

chunk%naunits = nun
chunk%natoms  = unit_cell%natoms*nun 
chunk%nblocks = unit_cell%nblocks*nun
allocate(chunk%aunit(nun), stat = ier)
if (ier /= 0) stop "p1asym_to_chunk> aunit malloc"

maxblock = 1
do i = 1, unit_cell%aunit(1)%natoms
  if (unit_cell%aunit(1)%atom(i)%iblock .gt. maxblock) maxblock = unit_cell%aunit(1)%atom(i)%iblock
end do

if (unit_cell%nblocks .ne. maxblock) then
  print *, "p1asym_to_chunk> something wrong with the nblock param of unit_cell"
  stop
end if

do i = 1,nun
  chunk%aunit(i)%natoms = unit_cell%aunit(1)%natoms
  if (allocated(chunk%aunit(i)%atom)) deallocate(chunk%aunit(i)%atom)
  allocate (chunk%aunit(i)%atom(chunk%aunit(i)%natoms),stat=ier)
  if (ier /= 0) STOP "*** not enough memory ***"
  do j=1,chunk%aunit(i)%natoms
    call init_protein_atom_type(chunk%aunit(i)%atom(j))
  end do
end do
! next copy in the stuff
! copy all asym units into one unit cell
jj=0
l =0
do i = 1,ia
  do j = 1,ib
    do k = 1,ic
      jj = jj+1
      do ii=1,unit_cell%aunit(1)%natoms
        l = l+1
        chunk%aunit(jj)%atom(ii) = unit_cell%aunit(1)%atom(ii)
        chunk%aunit(jj)%atom(ii)%x = unit_cell%aunit(1)%atom(ii)%x + &
                                     unit_cell%avec*real(i-1,kind=dp) + &
                                     unit_cell%bvec*real(j-1,kind=dp) + &
                                     unit_cell%cvec*real(k-1,kind=dp)
        chunk%aunit(jj)%atom(ii)%iblock = unit_cell%aunit(1)%atom(ii)%iblock+(jj-1)*maxblock 

!        write(666,*) chunk%aunit(jj)%atom(ii)%resn, l,chunk%aunit(jj)%atom(ii)%ires,&
!                     chunk%aunit(jj)%atom(ii)%iblock
      end do
    end do
  end do
end do




end subroutine

subroutine asym_shrinker(input,asyms_in,asyms_out)
use mod_types, only : asym_list_Type,inp_par
type (inp_par),          intent(in)    :: input
type (asym_list_Type),   intent(in)    :: asyms_in
type (asym_list_Type),   intent(out)   :: asyms_out
integer                                :: i,j,k,ier

asyms_out = asyms_in
if (allocated(asyms_out%aunit)) deallocate(asyms_out%aunit)
allocate(asyms_out%aunit(asyms_in%naunits),stat = ier)
if (ier /= 0) stop "asyms_to_p1asym> aunit malloc"

do j = 1, asyms_in%naunits
  call atom_shrinker(input,asyms_in%aunit(j),asyms_out%aunit(j))
end do

asyms_out%natoms = asyms_out%aunit(1)%natoms*asyms_out%naunits

end subroutine

subroutine atom_shrinker(input,atoms_in,atoms_out)
! this subroutine copies the contents of all asym units
! into one asym unit...  ie.  reducing the symm to p1
use mod_types, only : inp_par
type (inp_par),          intent(in)    :: input
type (protein_atom_list_type),   intent(in)    :: atoms_in
type (protein_atom_list_type),   intent(out)   :: atoms_out
integer                                :: i,j,k,ier,nat

nat = 0
! count up the atom type
if (trim(input%atom_anl) .ne. "all") then
  do i = 1, atoms_in%natoms
    if (trim(input%atom_anl) .eq. atoms_in%atom(i)%lab ) nat = nat + 1
  end do
else
  nat = atoms_in%natoms
end if

atoms_out%natoms = nat

if (allocated(atoms_out%atom)) deallocate(atoms_out%atom)
allocate (atoms_out%atom(nat),stat=ier)
if (ier /= 0) STOP "*** not enough memory ***"

do j=1,nat
  call init_protein_atom_type(atoms_out%atom(j))
end do

j = 0
do i = 1, atoms_in%natoms
  if (trim(input%atom_anl) .ne. "all") then
    if (trim(input%atom_anl) .eq. atoms_in%atom(i)%lab ) then
      j = j+1
      atoms_out%atom(j)= atoms_in%atom(i)
    end if
  else
    atoms_out%atom(i)= atoms_in%atom(i)
  end if
end do

end subroutine

subroutine asym_gen(input,cell,asym,SpG,unit_cell)
use mod_types, only : inp_par,asym_list_Type
use mod_inout, only : pdbwriter_new
use mod_symop, only : vecmat3_rotate
! DMR:  WORKING ON THIS NOW... get the data structure for unit cell squared away
!       so we can test and try the faster neighbor approach
type(inp_par)          , intent(in out)    :: input
type (Crystal_Cell_Type),intent(in)    :: Cell
type (protein_atom_list_type),   intent(in)    :: asym
type (space_group_type), intent(in)    :: SpG
type (asym_list_Type),   intent(out)  :: unit_cell
integer                                :: inc,i,j,jj,h,k,jx,jy,jz,neigh
integer                                :: nmultip,nmultip1,nmultip2,ier
real,     allocatable                  :: orb(:,:) ! for get_orbit subroutine 
real(dp)                               :: dist,dxyz(3),xyzi(3),xyzj(3),dcut
logical                                :: selfinteraction=.false.
integer                                :: nself,natoms,nat
real(dp), allocatable :: bleh(:,:)
real(dp),dimension(3) :: avec,bvec,cvec ! cartesian vecs
integer :: nneighs, na(2),nb(2),nc(2)

!dcut = rcut*rcut

! determine number of atoms in primitive unit cell
call number_of_p1atoms(asym,spg,natoms)
unit_cell%naunits = natoms/asym%natoms
unit_cell%natoms  = natoms
!unit_cell%naunits=SpG%multip
if(unit_cell%naunits .ne. spg%multip) then
  print *, 'WARNING: there are special positions. were you expecting this?'
  print *, 'calculated number of asym units:', unit_cell%naunits
  print *, 'space group multiplicity:',spg%multip
  print *, 'wyckoff num_orbit:', spg%wyckoff%num_orbit
endif

! DMR gotta get them centers in here!!
! input%block_cents

if (allocated(unit_cell%aunit)) deallocate(unit_cell%aunit)
allocate (unit_cell%aunit(1:unit_cell%naunits),stat=ier)
if (ier /= 0) STOP "*** not enough memory ***"


do i=1,unit_cell%naunits
  unit_cell%aunit(i)%natoms=asym%natoms
end do

!allocate atoms
do i=1,unit_cell%naunits
  if (allocated(unit_cell%aunit(i)%atom)) deallocate(unit_cell%aunit(i)%atom)
  allocate (unit_cell%aunit(i)%atom(asym%natoms),stat=ier)
  if (ier /= 0) STOP "*** not enough memory ***"
  do j=1,asym%natoms
    call init_protein_atom_type(unit_cell%aunit(i)%atom(j))
  end do
end do

if (allocated(input%block_cents)) then
  allocate(unit_cell%block_cents(unit_cell%naunits*size(input%block_cents(:,1)),3),stat = ier)
  if (ier /= 0) STOP "*** not enough memory ***"
  unit_cell%block_cents = 0.0d0
  do i =1, size(input%block_cents(:,1))
    input%block_cents(i,:) = matmul( input%block_cents(i,:), &
                                               transpose(Cell%orth_Cr_cel) )
  end do
  unit_cell%block_cents(1:size(input%block_cents(:,1)),:) = input%block_cents(:,:)
else
  allocate(unit_cell%block_cents(unit_cell%naunits,3),stat = ier)
  if (ier /= 0) STOP "*** not enough memory ***" 
  unit_cell%block_cents = 0.0d0
end if

allocate(orb(3,SpG%multip), stat=ier)

! fill up unit cell
do j=1,spg%multip
  do i=1, asym%natoms
    call Get_Orbit(asym%atom(i)%x(:),SpG,nmultip,orb)
    !call Get_Orbit_nocent(asym%atom(i)%x(:),SpG,nmultip,orb)
    unit_cell%aunit(j)%atom(i)   = asym%atom(i)
    unit_cell%aunit(j)%atom(i)%x = orb(:,j)
    unit_cell%aunit(j)%atom(i)%iaunit = j
    ! rotate anisotropic temp factor
    call vecmat3_rotate(cell,spg,j,unit_cell%aunit(j)%atom(i)%u)
    ! increment block id appropriately for bnm 
    unit_cell%aunit(j)%atom(i)%iblock = asym%atom(i)%iblock + input%nblocks*(j-1) 
  end do
end do

! temp
!call fract_to_cart(cell,unit_cell)
!  do i = 1, size(unit_cell%block_cents(:,1))
!    unit_cell%block_cents(i,:) = matmul(unit_cell%block_cents(i,:),transpose(Cell%Cr_orth_cel))
!    input%block_cents(i,:)     = unit_cell%block_cents(i,:)
!  end do

!call pdbwriter_new(unit_cell,"shitbird")

!stop

! fix up the block centers
if (allocated(input%block_cents)) then
  do j = 1,spg%multip
    do i = 1, size(input%block_cents(:,1))
      call Get_Orbit(unit_Cell%block_cents(i,:),SpG,nmultip,orb) 
      !call Get_Orbit_nocent(unit_Cell%block_cents(i,:),SpG,nmultip,orb) 
      unit_cell%block_cents(i+input%nblocks*(j-1),:) = orb(:,j) 
    end do
  end do
end if

deallocate(orb)

! now for this unit cell, lets use lattice translations to maximize the amount of each molecule
!   between 0 and 1
!   just use a quick avg...  should be close enough
allocate(bleh(spg%multip,3))
bleh = 0.0d0 

do j = 1,spg%multip
  do jj = 1, asym%natoms 
    xyzj = unit_cell%aunit(j)%atom(jj)%x
    bleh(j,:) = bleh(j,:) + xyzj 
  end do
end do

! comp avg fract coor for each asym
bleh = bleh/asym%natoms

do j = 1, spg%multip
  do jj = 1,asym%natoms
    unit_cell%aunit(j)%atom(jj)%x(1) = unit_cell%aunit(j)%atom(jj)%x(1)-floor(bleh(j,1))
    unit_cell%aunit(j)%atom(jj)%x(2) = unit_cell%aunit(j)%atom(jj)%x(2)-floor(bleh(j,2))
    unit_cell%aunit(j)%atom(jj)%x(3) = unit_cell%aunit(j)%atom(jj)%x(3)-floor(bleh(j,3))
  end do
  if (allocated(input%block_cents)) then
    do i = 1,size(input%block_cents(:,1))
      unit_cell%block_cents(i+input%nblocks*(j-1),1) = unit_cell%block_cents(i+input%nblocks*(j-1),1) - floor(bleh(j,1))
      unit_cell%block_cents(i+input%nblocks*(j-1),2) = unit_cell%block_cents(i+input%nblocks*(j-1),2) - floor(bleh(j,2))
      unit_cell%block_cents(i+input%nblocks*(j-1),3) = unit_cell%block_cents(i+input%nblocks*(j-1),3) - floor(bleh(j,3))
    end do
  end if
end do

if (allocated(input%block_cents)) then
  deallocate(input%block_cents)
  allocate(input%block_cents(size(unit_cell%block_cents(:,1)),3),stat=ier)
  do i = 1, size(unit_cell%block_cents(:,1))
    unit_cell%block_cents(i,:) = matmul(unit_cell%block_cents(i,:),transpose(Cell%Cr_orth_cel))
    input%block_cents(i,:)     = unit_cell%block_cents(i,:)
  end do
end if

! add some info to the unit cell so that cryst can be generated on the fly
avec = (/1.0d0, 0.0d0, 0.0d0/)
bvec = (/0.0d0, 1.0d0, 0.0d0/)
cvec = (/0.0d0, 0.0d0, 1.0d0/)

avec = matmul(CELL%Cr_Orth_cel,avec)
bvec = matmul(CELL%Cr_Orth_cel,bvec)
cvec = matmul(CELL%Cr_Orth_cel,cvec)

unit_cell%avec = avec
unit_cell%bvec = bvec
unit_cell%cvec = cvec

call count_neighs(input,cell,nneighs,na,nb,nc)

unit_cell%nneighs = nneighs
unit_cell%aneigh  = na
unit_cell%bneigh  = nb
unit_cell%cneigh  = nc

unit_cell%nblocks = input%nblocks*unit_cell%naunits

if( (trim(input%tpbc) .eq. "bvk") .or. (trim(input%tpbc) .eq. "pbc") ) input%nblocks = unit_cell%nblocks

end subroutine

subroutine block_assign_atasy (asym,asyms,nblocks)
use mod_types, only : asym_list_Type
type (protein_atom_list_type),   intent(in)     :: asym
type (asym_list_Type),   intent(inout)  :: asyms
integer,                 intent(in)     :: nblocks
integer :: i,j

do j=1,asyms%naunits
  do i=1, asym%natoms
    ! increment block id appropriately for bnm 
    asyms%aunit(j)%atom(i)%iblock = asym%atom(i)%iblock + nblocks*(j-1)
  end do
end do
asyms%nblocks = nblocks*asyms%naunits

end subroutine

subroutine block_assign_asyp1 (asyms,unit_cell)
use mod_types, only : asym_list_Type
type (asym_list_Type),   intent(in)     :: asyms
type (asym_list_Type),   intent(inout)  :: unit_cell
integer :: i,j,k

unit_cell%nblocks = asyms%nblocks

k = 0
do j=1,asyms%naunits
  do i=1, asyms%aunit(j)%natoms
    ! increment block id appropriately for bnm 
    k = k+1
    unit_cell%aunit(1)%atom(k)%iblock = asyms%aunit(j)%atom(i)%iblock
  end do
end do

end subroutine


subroutine nonasym_gen(rcut,cell,asym,SpG,unit_cell)
real(dp)               , intent(in)    :: rcut
type (Crystal_Cell_Type),intent(in)    :: Cell
type (protein_atom_list_type),   intent(inout)    :: asym !DMR test inout 12-2-08
type (space_group_type), intent(in)    :: SpG
type (asym_list_Type),   intent(out)  :: unit_cell
integer                                :: inc,i,j,h,k,jx,jy,jz,neigh
integer                                :: nmultip,nmultip1,nmultip2,ier
real,     allocatable                  :: orb(:,:) ! for get_orbit subroutine 
real(dp)                               :: dist,xyzi(3),xyzj(3),dcut,com(3), mag
logical                                :: selfinteraction=.false.
integer                                :: nself

dcut = rcut*rcut
print *, 'nonasym_gen called.'

! allocate the asym units
unit_cell%naunits=SpG%multip
if (allocated(unit_cell%aunit)) deallocate(unit_cell%aunit)
allocate (unit_cell%aunit(1:unit_cell%naunits),stat=ier)
if (ier /= 0) STOP "*** not enough memory ***"

do i=1,unit_cell%naunits
  unit_cell%aunit(i)%natoms=asym%natoms
end do

!allocate atoms
do i=1,unit_cell%naunits
  if (allocated(unit_cell%aunit(i)%atom)) deallocate(unit_cell%aunit(i)%atom)
  allocate (unit_cell%aunit(i)%atom(asym%natoms),stat=ier)
  if (ier /= 0) STOP "*** not enough memory ***"
  do j=1,asym%natoms
    call init_protein_atom_type(unit_cell%aunit(i)%atom(j))
  end do
end do

! test the min distance
! fill up asym unit arrays
allocate(orb(3,SpG%multip), stat=ier)

do j=1,SpG%multip
!  print *, 'shit', spg%symop(j)%tr
  do i=1, asym%natoms
    unit_cell%aunit(j)%atom(i)        = asym%atom(i)
    unit_cell%aunit(j)%atom(i)%x      = matmul(spg%symop(j)%rot,asym%atom(i)%x(:))  !+ &
!                                        spg%symop(j)%tr
  end do
end do
deallocate(orb)

end subroutine

subroutine asym_im_gen(input,cell,unit_cell,images)
use mod_types, only : inp_par,uc_list_type,asym_list_Type
! build up the crystal chunk: unit cell plus nearest neighbors  
! not sure if this is going overboard
!   build crystal as array of unit cells
!      unitcell is an array of asymunits
!        asymunit is and arry of atoms
type (inp_par),               intent(in)   :: input
type (Crystal_Cell_Type),     intent(in)   :: Cell
type (asym_list_Type),       intent(in)   :: unit_cell
type (uc_list_type),   intent(out)         :: images
integer                                    :: nat_asym 
integer                                    :: inc,i,j,k,kk,jx,jy,jz,neigh,ier
integer, dimension(2) :: na,nb,nc
integer               :: nneighs

call count_neighs(input,cell,nneighs,na,nb,nc)
images%aneigh=na ; images%bneigh=nb ; images%cneigh=nc
images%nneighs = nneighs

nat_asym = unit_cell%aunit(1)%natoms



!first allocate the neighs
if (allocated(images%neigh)) deallocate(images%neigh)
allocate (images%neigh(0:images%nneighs),stat=ier)
if (ier /= 0) STOP "*** not enough memory ***"

!next allocate the asyms

do i=0,images%nneighs
    images%neigh(i)%naunits=unit_cell%naunits
    if (allocated(images%neigh(i)%aunit)) deallocate(images%neigh(i)%aunit)
    allocate (images%neigh(i)%aunit(unit_cell%naunits),stat=ier)
end do

!next allocate the atoms
do i=0,images%nneighs
  do j=1,unit_cell%naunits
    images%neigh(i)%aunit(j)%natoms=unit_cell%aunit(j)%natoms
    if (allocated(images%neigh(i)%aunit(j)%atom)) deallocate(images%neigh(i)%aunit(j)%atom)
    allocate (images%neigh(i)%aunit(j)%atom(nat_asym),stat=ier)
    if (ier /= 0) STOP "*** not enough memory ***"
    do k=1,unit_cell%aunit(j)%natoms
      call init_protein_atom_type(images%neigh(i)%aunit(j)%atom(k))
    end do
  end do
end do

! copy in the unit cell
inc=0
do j=1,unit_cell%naunits
  images%neigh(inc)  = unit_cell
  do k=1,unit_cell%aunit(j)%natoms
    images%neigh(inc)%aunit(j)%atom(k)  = unit_cell%aunit(j)%atom(k)
  end do
end do

do jx=-1,1
  do jy=-1,1
    do jz=-1,1

      if ((jx .eq. 0) .and. (jy .eq. 0) .and.(jz .eq. 0)) cycle
        inc=inc+1
        do j=1,unit_cell%naunits
          do k=1,unit_cell%aunit(j)%natoms
            images%neigh(inc)%aunit(j)%atom(k)      = unit_cell%aunit(j)%atom(k)
            images%neigh(inc)%aunit(j)%atom(k)%x(1) = unit_cell%aunit(j)%atom(k)%x(1)+real(jx)
            images%neigh(inc)%aunit(j)%atom(k)%x(2) = unit_cell%aunit(j)%atom(k)%x(2)+real(jy)
            images%neigh(inc)%aunit(j)%atom(k)%x(3) = unit_cell%aunit(j)%atom(k)%x(3)+real(jz)
          end do
        end do

    end do
  end do
end do

end subroutine

subroutine count_neighs ( input,cell,nneighs,na,nb,nc)
use mod_types, only: inp_par
type (inp_par)            , intent(in)   :: input
type (Crystal_Cell_Type),intent(in)    :: Cell
integer, intent(out) :: nneighs,na(2),nb(2),nc(2)
integer :: nadir,nbdir,ncdir

nadir=int((input%rcut_end-1.0D-05)/cell%cell(1))+1
nbdir=int((input%rcut_end-1.0D-05)/cell%cell(2))+1
ncdir=int((input%rcut_end-1.0D-05)/cell%cell(3))+1

na=(/-nadir, nadir/)
nb=(/-nbdir, nbdir/)
nc=(/-ncdir, ncdir/)

nneighs = (2*nadir+1)*(2*nbdir + 1)*(2*ncdir + 1) - 1 ! subtract one for zerozerozero

end subroutine

subroutine image_gen(input,cell,unit_cell,images)
use mod_types, only: inp_par,neigh_list_Type
! build up the crystal chunk: unit cell plus nearest neighbors  
! all are included even though there are probably some tricks that can be used
! as in charmm where a half plus the corners are included.
type (inp_par)            , intent(in)   :: input
type (Crystal_Cell_Type),intent(in)    :: Cell
type (protein_atom_list_type),   intent(in)    :: unit_cell
type (neigh_list_Type),   intent(out)  :: images
integer                                :: inc,i,j,jx,jy,jz,neigh,ier,nadir,nbdir,ncdir
integer :: nneighs,na(2),nb(2),nc(2)
! determine how big of a chunk we need to build
! I subtract off a smidgeon, in case rcut_end=cell(:)
call count_neighs(input,cell,nneighs,na,nb,nc)
images%aneigh=na ; images%bneigh=nb ; images%cneigh=nc
images%nneighs = nneighs

if (allocated(images%neigh)) deallocate(images%neigh)
allocate (images%neigh(0:images%nneighs),stat=ier)
if (ier /= 0) STOP "*** not enough memory ***"

do i=0,images%nneighs
  images%neigh(i)%natoms=unit_cell%natoms
end do

do i=0,images%nneighs
  if (allocated(images%neigh(i)%atom)) deallocate(images%neigh(i)%atom)
  allocate (images%neigh(i)%atom(unit_cell%natoms),stat=ier)
  if (ier /= 0) STOP "*** not enough memory ***"
  do j=1,unit_cell%natoms
    call init_protein_atom_type(images%neigh(i)%atom(j))
  end do
end do

! copy in the unit cell
inc=0

do i=1, unit_cell%natoms
  images%neigh(inc)%atom(i)  = unit_cell%atom(i)
end do

do jx=images%aneigh(1),images%aneigh(2)
  do jy=images%bneigh(1),images%bneigh(2)
    do jz=images%cneigh(1),images%cneigh(2)
      if ((jx .eq. 0) .and. (jy .eq. 0) .and.(jz .eq. 0)) cycle
      inc=inc+1
      do i=1, unit_cell%natoms
        images%neigh(inc)%atom(i)      = unit_cell%atom(i)
        images%neigh(inc)%atom(i)%x(1) = unit_cell%atom(i)%x(1)+real(jx)
        images%neigh(inc)%atom(i)%x(2) = unit_cell%atom(i)%x(2)+real(jy)
        images%neigh(inc)%atom(i)%x(3) = unit_cell%atom(i)%x(3)+real(jz)
      end do
    end do
  end do
end do

if (input%print_level .gt. 1) print *, images%nneighs, ' neighbors'

end subroutine

subroutine chunk_gen(input,cell,unit_cell,chunk)
use mod_types, only: inp_par, neigh_list_Type
use mod_inout, only:pdbwriter_cryst 
! build up the crystal chunk: unit cell plus neighbors (defined within chunk datastructure)  
type (inp_par)          , intent(in)   :: input
type (Crystal_Cell_Type),intent(in)    :: Cell
type (protein_atom_list_type),   intent(in out)  :: unit_cell
type (neigh_list_Type),   intent(in out)  :: chunk
integer                                :: inc,i,j,jx,jy,jz,neigh,ier

!convert to fractional coords
call cart_to_fract (cell, unit_cell)

!compute number of neighs
inc=0
do jx=chunk%aneigh(1),chunk%aneigh(2)
  do jy=chunk%bneigh(1),chunk%bneigh(2)
    do jz=chunk%cneigh(1),chunk%cneigh(2)
      if ((jx .eq. 0) .and. (jy .eq. 0) .and.(jz .eq. 0)) cycle
      inc=inc+1
    end do
  end do
end do

chunk%nneighs = inc

if (allocated(chunk%neigh)) deallocate(chunk%neigh)
allocate (chunk%neigh(0:chunk%nneighs),stat=ier)
if (ier /= 0) STOP "*** not enough memory ***"

do i=0,chunk%nneighs
  chunk%neigh(i)%natoms=unit_cell%natoms
end do

do i=0,chunk%nneighs
  if (allocated(chunk%neigh(i)%atom)) deallocate(chunk%neigh(i)%atom)
  allocate (chunk%neigh(i)%atom(unit_cell%natoms),stat=ier)
  if (ier /= 0) STOP "*** not enough memory ***"
  do j=1,unit_cell%natoms
    call init_protein_atom_type(chunk%neigh(i)%atom(j))
  end do
end do

! copy in the unit cell
inc=0

do i=1, unit_cell%natoms
  chunk%neigh(inc)%atom(i)  = unit_cell%atom(i)
end do

!add in neighbors

do jx=chunk%aneigh(1),chunk%aneigh(2)
  do jy=chunk%bneigh(1),chunk%bneigh(2)
    do jz=chunk%cneigh(1),chunk%cneigh(2)
      if ((jx .eq. 0) .and. (jy .eq. 0) .and.(jz .eq. 0)) cycle
      inc=inc+1
      do i=1, unit_cell%natoms
        chunk%neigh(inc)%atom(i)      = unit_cell%atom(i)
        chunk%neigh(inc)%atom(i)%x(1) = unit_cell%atom(i)%x(1)+real(jx)
        chunk%neigh(inc)%atom(i)%x(2) = unit_cell%atom(i)%x(2)+real(jy)
        chunk%neigh(inc)%atom(i)%x(3) = unit_cell%atom(i)%x(3)+real(jz)
      end do
    end do
  end do
end do

!convert to cartesian coords
call fract_to_cart (cell, unit_cell)

do i=0,chunk%nneighs
  call fract_to_cart (cell, chunk%neigh(i))
end do

call pdbwriter_cryst(chunk,trim(input%fileroot)//"-chunky")

end subroutine

! this cart to fract should be able to be overloaded with protein_atom_list_type
! or uc_type as they are similarly called
subroutine atomlist_cart_to_fract (cell, atoms_in)
type (Crystal_Cell_Type),intent(in)     :: Cell
type (protein_atom_list_type),   intent(in out) :: atoms_in
integer                                 :: i

    do i=1,atoms_in%natoms
      atoms_in%atom(i)%x(:)=matmul(atoms_in%atom(i)%x(:),transpose(Cell%orth_Cr_cel))
    end do

end subroutine

subroutine atomlist_fract_to_cart (cell, atoms_in)
type (Crystal_Cell_Type),intent(in)     :: Cell
type (protein_atom_list_type),   intent(in out) :: atoms_in
integer                                 :: i

    do i=1,atoms_in%natoms
      atoms_in%atom(i)%x(:)=matmul(atoms_in%atom(i)%x(:),transpose(Cell%Cr_orth_cel))
    end do

end subroutine

subroutine asymlist_fract_to_cart (cell, atoms_in)
type (Crystal_Cell_Type),intent(in)     :: Cell
type (Asym_list_Type),   intent(in out) :: atoms_in
integer                                 :: i,j

  do j = 1,atoms_in%naunits
    do i = 1,atoms_in%aunit(j)%natoms
      atoms_in%aunit(j)%atom(i)%x(:) = matmul( atoms_in%aunit(j)%atom(i)%x(:), &
                                               transpose(Cell%Cr_orth_cel) )
    end do
  end do

end subroutine

subroutine asymlist_cart_to_fract (cell, atoms_in)
type (Crystal_Cell_Type),intent(in)     :: Cell
type (Asym_list_Type),   intent(in out) :: atoms_in
integer                                 :: i,j

  do j = 1,atoms_in%naunits
    do i = 1,atoms_in%aunit(j)%natoms
      atoms_in%aunit(j)%atom(i)%x(:) = matmul( atoms_in%aunit(j)%atom(i)%x(:), &
                                               transpose(Cell%orth_Cr_cel) )
    end do
  end do

end subroutine

subroutine vec_c2f (cell, vec)
! converts vector (x1,y1,z1,.... xn,yn,zn) from cart to fractional coordinates
type (Crystal_Cell_Type),intent(in)     :: Cell
real(dp), allocatable,   intent(in out) :: vec(:)
integer                                 :: i,j,k,natoms

natoms = size(vec)/3

k=1
  do j = 1,natoms
    vec(k:k+2) = matmul(vec(k:k+2),transpose(Cell%orth_Cr_cel))
    k = k + 3
  end do

end subroutine

subroutine vec_f2c (cell, vec)
! converts vector (x1,y1,z1,.... xn,yn,zn) from fractional to cart coordinates
type (Crystal_Cell_Type),intent(in)     :: Cell
real(dp), allocatable,   intent(in out) :: vec(:)
integer                                 :: i,j,k,natoms

natoms = size(vec)/3
k=1
  do j = 1,natoms
    vec(k:k+2) = matmul(vec(k:k+2),transpose(Cell%Cr_orth_cel))
    k = k + 3
  end do

end subroutine

subroutine mass_by_residue(chemsymb,resn,mass)
! table of residue masses
! input: 
!   chemsymb: chemical symbol C,N,O... 
!   resn: name of residue     LYS,LEU ...
!   output: mass
!   default condition, when residue is not found, returns atomic mass

character(len=*) , intent(in)  :: chemsymb,resn
real(dp),          intent(out) :: mass
real(sp) :: mass_default

call Get_Atomic_Mass(chemsymb,mass_default)

select case(resn)
  case("ASP")
    mass = 114.08078d0
  case("PRO")
    mass = 98.12486d0
  case("ILE")
    mass = 113.15977d0
  case("LYS")
    mass = 129.18241d0
  case("TRP")
    mass = 183.18959d0
  case("GLY")
    mass = 57.05201d0
  case("CYS")
    mass = 103.13895d0
  case("PHE")
    mass = 147.17683d0
  case("GLN")
    mass = 128.13096d0
  case("SER")
    mass = 87.07835d0
  case("ASN")
    mass = 114.10402d0
  case("VAL")
    mass = 99.13283d0
  case("LEU")
    mass = 113.15977d0
  case("TYR")
    mass = 163.17623d0
  case("GLU")
    mass = 128.10772d0
  case("ARG")
    mass = 157.19581d0
  case("THR")
    mass = 101.10529d0
  case("ALA")
    mass = 71.07895d0
  case("MET")
    mass = 131.19283d0
  case("HIS")
    mass = 137.14129d0
  case("PHIS")
    mass = 138.14926d0
  case default
    mass = mass_default 
end select

end subroutine

subroutine residue_masses (atoms_in, masses)
! this sub returns an array of residue masses summed over atoms in res
! for this sort of mass weighting, the protein should be course grained
! to one atom per residue, we'll accomplish this with the requirement that all labels be the same
!       
type (protein_atom_list_type),   intent(in)  :: atoms_in
real(dp), allocatable,   intent(out) :: masses(:)
integer                              :: i,ier
character(len=6) :: label1

label1 = atoms_in%atom(1)%lab

allocate(masses(atoms_in%natoms), stat=ier)
if(ier /=0) stop "residue_masses> malloc"

do i = 1, atoms_in%natoms
  if(trim(atoms_in%atom(i)%lab) .eq. trim(label1)) then
    call mass_by_residue(atoms_in%atom(i)%chemsymb,trim(atoms_in%atom(i)%resn),masses(i))
  else
    stop "residue_masses> mass weighting by residue called, but more than 1 atper res?"
  end if
end do

end subroutine

subroutine atom_mass_constr(input,atoms,masses)
use mod_types,     only: inp_par
use mod_constants, only: one
type (inp_par),          intent(in)  :: input 
type (protein_atom_list_type),   intent(in)  :: atoms
real(dp),allocatable,    intent(out) :: masses(:)
integer :: ier

select case(trim(input%weigh))
  case("mass")
    call mass_vector_build (atoms,masses)
  case("byatom")
    call mass_vector_build (atoms,masses)
  case("byres") 
    call residue_masses (atoms,masses)
  case default
    if (allocated(masses)) deallocate(masses)
    allocate(masses(atoms%natoms),stat=ier)
    if(ier/=0) stop "atom_mass_constr> malloc"
    masses = one
end select

end subroutine

subroutine spherical_carve(iat,max_coor,dmax,Cell,Spg,A,coord)
! DMR: 05-24-2012
! get spherical carve of crystal about single atom (iat)
! derived from  Set_TDist_Coordination(max_coor,Dmax, Cell, Spg, A) in Geom_calc
       !---- Arguments ----!
    integer,                  intent(in)   :: iat,max_coor
    real(kind=dp),            intent(in)   :: dmax
    type (Crystal_cell_Type), intent(in)   :: Cell
    type (Space_Group_Type),  intent(in)   :: SpG
    type (protein_atom_list_type),    intent(inout):: A  ! inout to be able to convert to fract if needed
    character (len=25)   ,    intent(in)   :: coord

    !---- Local Variables ----!
    integer                              :: i,j,k,lk,i1,i2,i3,nn,L
    integer,       dimension(3)          :: ic1,ic2
    real(kind=dp), dimension(3)          :: xx,x1,xo,Tn,xr, QD
    real(kind=dp)                        :: T,dd,epsi
    real(kind=dp), dimension(3,max_coor) :: uu
    
    epsi = 0.001

    if (trim(coord) .eq. "cart") then
      call cart_to_fract(cell,A)
    end if
    ! we work with fractional, but report cartesian

    qd(:)=1.0/cell%rcell(:)
    ic2(:)= int(dmax/cell%cell(:))+1
    ic1(:)=-ic2(:)

    xo(:)=a%atom(iat)%x(:)
    write(6, '(A3,1x,3F10.3,1x,A5)') a%atom(iat)%ChemSymb,matmul(cell%cr_orth_cel,xo),a%atom(iat)%lab
    do k=1,a%natoms
       lk=1
       uu(:,lk)=xo(:)
       do j=1,Spg%Multip
          xx=ApplySO(Spg%SymOp(j),a%atom(k)%x)
          do i1=ic1(1),ic2(1)
             do i2=ic1(2),ic2(2)
                do_i3:do i3=ic1(3),ic2(3)
                      Tn(1)=real(i1); Tn(2)=real(i2); Tn(3)=real(i3)
                      x1(:)=xx(:)+tn(:)
                      do l=1,3
                         t=abs(x1(l)-xo(l))*qd(l)
                         if (t > dmax) cycle  do_i3
                      end do
                      do nn=1,lk
                         if (sum(abs(uu(:,nn)-x1(:)))  <= epsi) cycle  do_i3
                      end do
                      xr = matmul(cell%cr_orth_cel,x1-xo)
                      dd=sqrt(dot_product(xr,xr))
                      if (dd > dmax .or. dd < 0.001) cycle
                      lk=lk+1
                      uu(:,lk)=x1(:)
                      write(6, '(A3,1x,3F10.3,1x,A5)') a%atom(k)%ChemSymb,matmul(cell%cr_orth_cel,x1),a%atom(k)%lab
                      !Coord_Info%Dist(i,Coord_Info%Coord_Num(i))=dd
                end do do_i3 !i3
             end do !i2
          end do !i1
       end do !j
    end do !k

    if (trim(coord) .eq. "cart") then
      call fract_to_cart(cell,A)
    endif
    return
  end subroutine

end module

