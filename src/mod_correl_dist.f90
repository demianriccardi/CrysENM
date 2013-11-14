module mod_correl_dist
! Computation of correlation from ad hoc type stuff to GNM and ENM 
! DMR 10-30-2007
!
!
! TO DO:
! 1. distance based correlations   
!   a. absolute corr
!   b. distance scaled (distance squared?)
! 2. GNM
! 3. ENM
! 4. Interface for eigenvects from CHARMM
! 5. Crystal dyn interface
!        

use cfml_globaldeps,                 only: dp,sp
use cfml_Atom_typedef,               only: Atom_List_Type
use mod_types,                       only: inp_par,asym_List_Type
use mod_math
use mod_constants

implicit none

interface corrtab_dexp
  module procedure corrtabdexp_iso,corrtabdexp_pbc
end interface corrtab_dexp

contains

subroutine corrtabdexp_iso(input,atoms, &
                                  correlation)
! DMR 11/14/08
! isotropic correlation matrix constructor
! C_kk' = C0 x exp(-r_kk' / gamm)
! rcut is the cutoff
type(inp_par), intent(in)         :: input
type(atom_list_type),intent(in)   :: atoms
real(dp),intent(out), allocatable   :: correlation(:,:)
integer                         :: i,j,natom,ialloc
real(dp)                        :: dist,dxyz(3),rcut,gamm,c0

rcut=input%isocorr_cut
gamm=input%isocorr_gamm
c0=input%isocorr_c0

natom=atoms%natoms

if (allocated(correlation)) deallocate(correlation)
allocate(correlation(natom,natom),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"


correlation=zero

do i=1,natom
  do j=i,natom
    dxyz(:) = atoms%atom(i)%x(:)-atoms%atom(j)%x(:)
    dist=distance(dxyz)
    if(dist .le. rcut) then
      correlation(i,j)= c0 * dexp(-dist/gamm)
    end if
      correlation(j,i)=correlation(i,j)
  end do
end do

end subroutine

subroutine corrtabdexp_pbc(input,unit_cell, &
                                  correlation)
! DMR 5/09
! this subroutine needs to be checked again... what about multiple interactions 
! isotropic correlation matrix constructor
! C_kk' = C0 x exp(-r_kk' / gamm)
! rcut is the cutoff
type(inp_par), intent(in)           :: input
type(asym_list_type),intent(in)     :: unit_cell
real(dp),intent(out), allocatable   :: correlation(:,:)
integer                             :: i,ii,j,jj,k,natoms,naunits,ialloc
real(dp)                            :: dist,dxyz(3),rcut,gamm,c0
real(dp)                            :: dist2,dxyz2(3)
real(dp), dimension(3)              :: vec,avec,bvec,cvec,xyzi,xyzj
integer :: ja,jb,jc

print *, 'corrtab_pbc> allows multiple interactions...  double check this conceptually'

rcut=input%isocorr_cut
gamm=input%isocorr_gamm
c0=input%isocorr_c0

natoms  = unit_cell%aunit(1)%natoms
naunits = unit_cell%naunits
print *, "number of asy units", naunits
print *, "number of atoms", natoms

if (allocated(correlation)) deallocate(correlation)
allocate(correlation(natoms,natoms),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

correlation=zero

do ii=1,natoms

  xyzi = unit_cell%aunit(1)%atom(ii)%x
  do k = 1, naunits
  do jj=ii,natoms

    ! DMR do for relevant lattice displacements 
    do ja=unit_cell%aneigh(1),unit_cell%aneigh(2)
    avec = real(ja,kind=dp)*unit_cell%avec
    do jb=unit_cell%bneigh(1),unit_cell%bneigh(2)
    bvec = real(jb,kind=dp)*unit_cell%bvec
    do jc=unit_cell%cneigh(1),unit_cell%cneigh(2)
    cvec = real(jc,kind=dp)*unit_cell%cvec
    vec = avec + bvec + cvec

    xyzj = unit_cell%aunit(k)%atom(jj)%x+vec

    dxyz = xyzi-xyzj
    dist=distance(dxyz)

! use something along these lines to avoid multiple interactions
!    do ii=0,images%nneighs
!      do jj=ii,images%nneighs
!        dxyz2(:) = images%neigh(ii)%atom(i)%x(:)-images%neigh(jj)%atom(j)%x(:)
!        dist2=distance(dxyz2)
!        if (dist2 .lt. dist) then
!          dist = dist2
!        end if
!      end do
!    end do

    if(dist .le. rcut) then
      correlation(ii,jj)= correlation(ii,jj) + c0 * dexp(-dist/gamm)
    end if
      correlation(jj,ii)=correlation(ii,jj)
    end do
    end do
    end do
  end do
  end do
end do

print *, "exiting isocorr setup"

end subroutine

subroutine bfact_corrtab(biso,correlation)
! DMR 11/14/08
! multiply the isotropic correlation matrix by corresponding bfactors
! dsqrt[biso(k)biso(k')] * correlation(k,k')
real(dp),intent(in), allocatable     :: biso(:)
real(dp),intent(inout), allocatable  :: correlation(:,:)
integer                              :: i,j,natom,ialloc

natom=size(biso)

do i=1,natom
  do j=i,natom
    correlation(i,j) = dsqrt(biso(i)*biso(j))*correlation(i,j)  
    correlation(j,i)=correlation(i,j)
  end do
end do

end subroutine

subroutine correl_distscl(rcut,atoms, &
                                  correlation)
! absolute correlation (ie: 1.0) based on a cutoff of distance
real(dp),intent(in)             :: rcut
type(atom_list_type),intent(in)   :: atoms
real(dp),intent(out), allocatable   :: correlation(:,:)
integer                         :: i,j,natom,ialloc
real(dp)                        :: dist,dst_sqr,dcut,dxyz(3)

natom=atoms%natoms

if (allocated(correlation)) deallocate(correlation)
allocate(correlation(natom,natom),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

dcut=rcut**2

correlation=zero

do i=1,natom
  do j=i,natom
    dxyz = atoms%atom(i)%x(:) - atoms%atom(j)%x(:)
    dst_sqr=dist_sqr(dxyz)
    if(dst_sqr .le. dcut) then
      correlation(i,j)=one/dst_sqr
    end if
      correlation(j,i)=correlation(i,j)
  end do
end do

end subroutine

subroutine correl_distabs_pbc(rcut,unit_cell, &
                                  correlation)
! absolute correlation (ie: 1.0) based on a cutoff of distance
real(dp),intent(in)               :: rcut
type(asym_List_Type),intent(in)   :: unit_cell
real(dp),intent(out), allocatable     :: correlation(:,:)
integer                           :: i,ii,j,jj,natoms,ialloc
real(dp)                          :: dist,dst_sqr,dcut,dxyz(3)
real(dp), dimension(3)          :: vec,avec,bvec,cvec,xyzi,xyzj
integer :: ja,jb,jc

natoms   = unit_cell%natoms

if (allocated(correlation)) deallocate(correlation)
allocate(correlation(natoms,natoms),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"


correlation=zero

do ii=1,natoms

  xyzi = unit_cell%aunit(1)%atom(ii)%x

  do jj=ii,natoms

    ! DMR do for relavent lattice displacements 
    do ja=unit_cell%aneigh(1),unit_cell%aneigh(2)
    avec = real(ja,kind=dp)*unit_cell%avec
    do jb=unit_cell%bneigh(1),unit_cell%bneigh(2)
    bvec = real(jb,kind=dp)*unit_cell%bvec
    do jc=unit_cell%cneigh(1),unit_cell%cneigh(2)
    cvec = real(jc,kind=dp)*unit_cell%cvec
    vec = avec + bvec + cvec

    xyzj = unit_cell%aunit(1)%atom(jj)%x+vec

    dxyz = xyzi-xyzj
    dist=distance(dxyz)
    if(dist .le. rcut) then
        correlation(ii,jj)=one
    end if
    correlation(jj,ii)=correlation(ii,jj)
    end do
    end do
    end do
  end do
end do

end subroutine

subroutine correl_distscl_pbc(rcut,unit_cell, &
                                  correlation)
! absolute correlation (ie: 1.0) based on a cutoff of distance
real(dp),intent(in)               :: rcut
type(asym_List_Type),intent(in)   :: unit_cell
real(dp),intent(out), allocatable     :: correlation(:,:)
integer                           :: i,ii,j,jj,natoms,ialloc
real(dp)                          :: dist,dst_sqr,dcut,dxyz(3)
real(dp), dimension(3)          :: vec,avec,bvec,cvec,xyzi,xyzj
integer :: ja,jb,jc

natoms   = unit_cell%natoms

if (allocated(correlation)) deallocate(correlation)
allocate(correlation(natoms,natoms),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"


correlation=zero

do ii=1,natoms

  xyzi = unit_cell%aunit(1)%atom(ii)%x

  do jj=ii,natoms

    ! DMR do for relavent lattice displacements 
    do ja=unit_cell%aneigh(1),unit_cell%aneigh(2)
    avec = real(ja,kind=dp)*unit_cell%avec
    do jb=unit_cell%bneigh(1),unit_cell%bneigh(2)
    bvec = real(jb,kind=dp)*unit_cell%bvec
    do jc=unit_cell%cneigh(1),unit_cell%cneigh(2)
    cvec = real(jc,kind=dp)*unit_cell%cvec
    vec = avec + bvec + cvec

    xyzj = unit_cell%aunit(1)%atom(jj)%x+vec

    dxyz = xyzi-xyzj
    dist=distance(dxyz)
    if(dist .le. rcut) then
        correlation(ii,jj)=one/(dist*dist)
    end if
    correlation(jj,ii)=correlation(ii,jj)
    end do
    end do
    end do
  end do
end do

end subroutine

end module
