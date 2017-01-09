module mod_bnm
! DMR 04-27-2009
!
! for blocks of atoms computes:
!   moment of inertia
!   inertia tensor
!   center of mass
!   transformation matrix (P) to/fro cartesian and T/R   
!       
use mod_types,                  only: inp_par,sparse,protein_atom_list_type
use mod_constants

implicit none

!type block
!  integer :: natoms
!  integer, allocatable :: atinblock(:) ! this is contain "pointers" of each atom in the block
                                   ! useful if blocks contain atoms that are not in order
!  integer, allocatable :: ixyz(:)  ! give the indices corresponding to ats, ie. the location
                                   ! of block
!  real(dp) :: com(3),rotc(3),cent(3),momI(3),tnsI(3,3),eigI(3,3),tmass
!  type(sparse),allocatable :: pmat(:,:)
!  logical :: cent_read
!end type 
type block
  integer :: natoms
  integer, allocatable :: atinblock(:) ! this is contain "pointers" of each atom in the block
                                   ! useful if blocks contain atoms that are not in order
  integer, allocatable :: ixyz(:)  ! give the indices corresponding to ats, ie. the location
                                   ! of block
  real(dp) :: com(3),rotc(3),cent(3),momI(3),tnsI(3,3),eigI(3,3),tmass
  real(dp),allocatable :: pmat(:,:)
  real(dp),allocatable :: tmpvec(:),  tmpmat(:,:) ! add in tmpmat for any sort of uses
  integer,allocatable :: itmpvec(:), itmpmat(:,:) !
  logical :: cent_read
! TLS stuffs
  real(dp),dimension(:,:),allocatable :: T,L,S
  real(dp),dimension(:,:),allocatable :: Lv,rho
  real(dp),dimension(:),allocatable :: Le
  real(dp),dimension(:), allocatable :: F_sqr ! structure factor squared ie. sum_{i,j}f_i(Q)*f_j(Q)*cos(2pi*Q.r_{ij})
end type 

interface bnm_kirchoff
  module procedure bnm_kirchoff_iso, bnm_kirchoff_pbc
end interface

interface blocks_setup
  module procedure blocks_setup_iso,blocks_setup_pbc
end interface

interface blocks_com
  module procedure blocks_com_iso,blocks_com_pbc
end interface

interface blocks_tnsi
  module procedure blocks_tnsi_iso,blocks_tnsi_pbc
end interface

contains

subroutine bnm_statalloc(blocks)
! DMR 2-19-2010 size printer
use mod_types,       only:sparse
type(block), allocatable, intent(in)    :: blocks(:)
integer :: natoms, natinblock, nixyz, npmat,ntmpvec,ntmpmat,nitmpvec,nitmpmat, &
           nTLS,nLv,nrho,nle,nf_sqr
integer :: i

natoms=0
natinblock=0
nixyz=0
npmat=0
ntmpvec=0
ntmpmat=0
nitmpvec=0
nitmpmat=0
nTLS=0
nLv=0
nrho=0
nle=0
nf_sqr=0

do i = 1, size(blocks)
  natoms = natoms + blocks(i)%natoms
  if(allocated(blocks(i)%atinblock)) natinblock = natinblock + size(blocks(i)%atinblock)
  if(allocated(blocks(i)%ixyz)) nixyz = nixyz + size(blocks(i)%ixyz)
  if(allocated(blocks(i)%pmat)) npmat = npmat + size(blocks(i)%pmat)
  if(allocated(blocks(i)%tmpvec)) ntmpvec = ntmpvec + size(blocks(i)%tmpvec)
  if(allocated(blocks(i)%tmpmat)) ntmpmat = ntmpmat + size(blocks(i)%tmpmat)
  if(allocated(blocks(i)%itmpvec)) nitmpvec = nitmpvec + size(blocks(i)%itmpvec)
  if(allocated(blocks(i)%itmpmat)) nitmpmat = nitmpmat + size(blocks(i)%itmpmat)
  if(allocated(blocks(i)%T)) nTLS = nTLS + size(blocks(i)%T)
  if(allocated(blocks(i)%L)) nTLS = nTLS + size(blocks(i)%L)
  if(allocated(blocks(i)%S)) nTLS = nTLS + size(blocks(i)%S)
  if(allocated(blocks(i)%lv)) nlv = nlv + size(blocks(i)%lv)
  if(allocated(blocks(i)%rho)) nrho = nrho + size(blocks(i)%rho)
  if(allocated(blocks(i)%le)) nle = nle + size(blocks(i)%le)
  if(allocated(blocks(i)%F_sqr)) nF_sqr = nF_sqr + size(blocks(i)%F_sqr)
  
end do

print *, "printing summed dimensions for block data_structure"
print *, "natoms",      natoms
print *, "natinblock",  natinblock
print *, "nixyz",       nixyz
print *, "npmat",       npmat
print *, "ntmpvec",     ntmpvec
print *, "ntmpmat",     ntmpmat
print *, "nitmpvec",    nitmpvec
print *, "nitmpmat",    nitmpmat
print *, "nTLS",        nTLS
print *, "nLv",         nLv
print *, "nrho",        nrho
print *, "nle",         nle
print *, "nf_sqr",      nf_sqr


end subroutine


subroutine bnm_tpbc_blocker(input,asyms,unit_cell,blocks)
! DMR added Oct 1, 2009
use mod_types,     only : inp_par,asym_list_type
type(inp_par),         intent(in out) :: input
type(asym_list_type),  intent(in out) :: asyms, unit_cell
type(block),allocatable, intent(out)  :: blocks(:)
integer :: ialloc,nzeros,nblocks,i,j,natoms
type(sparse) :: bnmhess_coor


if(trim(input%tpbc) .eq. "pbc" .or. trim(input%tpbc) .eq. "bvk") then
  call blocks_setup (input,unit_cell,blocks)
else if(trim(input%tpbc) .eq. "asy") then
  call blocks_setup (input,asyms,blocks)
  natoms = asyms%aunit(1)%natoms
else
  call blocks_setup (input,asyms%aunit(1),blocks)
end if

end subroutine

subroutine blocks_subblocker(input,atoms_in,blocks_in,blocks_out)
! form subblocks (ie. one or few atoms) from bigblocks (all atom) 
! using only subset of atoms original blocks
! this can't be accomplished by forming the smaller collection of atoms
!   and forming the blocks from them because of orthog and centers will be different
!
!
use mod_types,     only: inp_par
use mod_linalg,    only: lapack_eig3
use mod_crysbuild, only: atom_shrinker
type(inp_par),            intent(in)      :: input
type(protein_atom_list_type),     intent(in)      :: atoms_in
type(block), allocatable, intent(in)      :: blocks_in(:)
type(block), allocatable, intent(inout)   :: blocks_out(:)
type(protein_atom_list_type)                      :: subatoms
integer :: i,ii,jj,kk,k,j,ialloc,nat,ier

call atom_shrinker(input,atoms_in,subatoms)
! set up most of the stuff
call blocks_setup (input,subatoms,blocks_out)
print *, 'copy subsection of blocks via subblocker'

! now make blocks_out more like blocks_in, but don't mess up the indices of where atoms
! and coords are, because those are different now
do i = 1, size(blocks_in)
     blocks_out(i)%cent = blocks_in(i)%cent
     blocks_out(i)%com  = blocks_in(i)%com
     blocks_out(i)%rotc = blocks_in(i)%rotc
     blocks_out(i)%momi = blocks_in(i)%momi
     blocks_out(i)%tnsi = blocks_in(i)%tnsi
     blocks_out(i)%eigi = blocks_in(i)%eigi
     blocks_out(i)%tmass = blocks_in(i)%tmass
     blocks_out(i)%cent_read = blocks_in(i)%cent_read
end do

! now copy pmat
do i = 1,size(blocks_in)
  blocks_out(i)%pmat = zero

  ii = 1
  jj = 1
  kk = 1
  do j = 1,size(blocks_in(i)%atinblock)
    if (atoms_in%atom(blocks_in(i)%atinblock(j))%lab .eq. trim(input%atom_anl)) then 
!    set up 3n coordinates
      do k = 0,2
        blocks_out(i)%pmat(:,jj+k) = blocks_in(i)%pmat(:,kk+k) 
      end do
      ii = ii + 1
      jj = jj + 3
    end if
    kk = kk + 3 
  end do
end do

print *, 'exit subblocker'

end subroutine

subroutine bnm_interact_tf(input,atoms1,atoms2,lattvec,ablock1,ablock2,tfint)
use mod_types, only: inp_par
use mod_math,                 only: distance
type(inp_par),          intent(in)  :: input
type(protein_atom_list_type),   intent(in)  :: atoms1,atoms2
real(dp), dimension(3), intent(in)  :: lattvec
type(block),            intent(in)  :: ablock1,ablock2
logical,                intent(out) :: tfint
real(dp) :: dist,rcut,dxyz(3),xyzi(3),xyzj(3)
integer :: i,j

tfint = .false.

do i = 1, ablock1%natoms
  xyzi = atoms1%atom(ablock1%atinblock(i))%x
  do j = 1, ablock2%natoms
    xyzj = atoms2%atom(ablock2%atinblock(j))%x+lattvec
    dxyz = xyzi-xyzj
    dist = distance(dxyz)
    if (dist .le. input%rcut_start) then
      if (dist .gt. zero) tfint = .true.
      return
    end if
  end do
end do

end subroutine

subroutine bnm_kirchoff_iso(input,atoms,blocks,kirchoff)
use mod_types, only: inp_par
use mod_math,                 only: distance
type(inp_par),           intent(in)  :: input
type(protein_atom_list_type),    intent(in)  :: atoms
type(block),allocatable, intent(in)  :: blocks(:)
real(dp), allocatable,  intent(in out) :: kirchoff(:,:)
logical    :: tfint
integer    :: i,j,nblocks
real(dp), dimension(3)               :: lattvec

nblocks = size(blocks)
lattvec = zero 

do i =1,nblocks
  do j = i, nblocks
    call bnm_interact_tf(input,atoms,atoms,lattvec,blocks(i),blocks(j),tfint) 
    if (tfint) then
      kirchoff(i,j) = kirchoff(i,j) - one
      kirchoff(i,i) = kirchoff(i,i) + one
      if(i.ne.j) kirchoff(j,j) = kirchoff(j,j) + one
      kirchoff(j,i) = kirchoff(i,j)
    end if 
  end do
end do

end subroutine

subroutine bnm_kirchoff_pbc(input,asyms,blocks,kirchoff)

use mod_types, only: inp_par,asym_list_type
use mod_math,                 only: distance
type(inp_par),           intent(in)  :: input
type(asym_list_type),    intent(in)  :: asyms
type(block),allocatable, intent(in)  :: blocks(:)
real(dp), allocatable,  intent(in out) :: kirchoff(:,:)
logical    :: tfint
integer     :: i,ii,j,jj,nblocks,ja,jb,jc
real(dp), dimension(3)               :: vec,avec,bvec,cvec

nblocks = size(blocks)

if (asyms%naunits .gt. 1) then
  ii = input%aunit
else
  ii = 1
end if

do i =1,nblocks
  do j = i, nblocks
    do jj =1, asyms%naunits
      do ja=asyms%aneigh(1),asyms%aneigh(2)
      avec = real(ja,kind=dp)*asyms%avec
      do jb=asyms%bneigh(1),asyms%bneigh(2)
      bvec = real(jb,kind=dp)*asyms%bvec
      do jc=asyms%cneigh(1),asyms%cneigh(2)
      cvec = real(jc,kind=dp)*asyms%cvec
      vec = avec + bvec + cvec
      call bnm_interact_tf(input,asyms%aunit(ii),asyms%aunit(jj),vec,blocks(i),blocks(j),tfint)
      if (tfint) then
        kirchoff(i,j) = kirchoff(i,j) - one
        kirchoff(i,i) = kirchoff(i,i) + one
        if (i.ne.j) kirchoff(j,j) = kirchoff(j,j) + one
        kirchoff(j,i) = kirchoff(i,j)
      end if 
      end do
      end do
      end do
    end do
  end do
end do

end subroutine

!subroutine bnm_vcov(input,cell,atoms,unit_cell,kirchoff,hessian,vcov)
!use mod_types,     only : sparse,inp_par,asym_list_type,dos
!use mod_hessian,   only : dvcov
!use mod_inout,     only : full_to_sparse
!use mod_lattdyn,   only : bnm_qdispers
!use crystal_types, only : Crystal_Cell_Type
!type(inp_par),        intent(in out) :: input
!type (Crystal_Cell_Type),intent(in) :: cell
!type(protein_atom_list_type), intent(in out) :: atoms
!type(asym_list_type), intent(in out) :: unit_cell
!type (sparse)       , intent(in)     :: kirchoff,hessian
!real(dp), allocatable , intent(out)  :: vcov(:,:)
!type (sparse) :: bnmhess
!real(dp), allocatable                :: vals(:),vecs(:,:),atvecs(:,:)
!type(block),allocatable              :: blocks(:)
!integer :: ialloc,nzeros
!type(sparse) :: bnmhess_coor
!type(dos) :: states
!
!call blocks_setup (input,atoms,blocks)
!call hess_proj(hessian,blocks,bnmhess)
!
!if (trim(input%tpbc) .ne. "bvk") then
!  call bnm_eigenal_valvec(input,bnmhess_coor,vals,vecs)
!  input%first_mode = 1
!  input%last_mode = size(vals)
!  call vect_proj(blocks,vecs,atvecs)
!  allocate(vcov(size(atvecs(:,1)),size(atvecs(:,1))),stat=ialloc)
!  vcov=zero
!  call dvcov(input,vals,atvecs,vcov,nzeros)
!else
! do it
! call bnm_qdispers(input,cell,blocks,unit_cell,kirchoff,hessian,vcov,states)
!end if
!
!! dmr needs to be tested for mass weighting
!
!end subroutine
!
subroutine block_init(blocks)
type(block), intent(in out) :: blocks

blocks%natoms     = 0
blocks%com        = zero
blocks%rotc       = zero
blocks%rotc       = zero
blocks%tmass      = zero
blocks%momI       = zero
blocks%tnsI       = zero
blocks%eigI       = zero
blocks%cent_read  = .false.

end subroutine

subroutine vect_proj(blocks,bvects,atvects)
! projects 6n vectors of dim 6n to 6n vectors of dim 3N
! this is just quick implementation
! currently, each block has a pmat that is 6x3Nb 
! and a pointer array that gives xyz of the 3Nb
! we may want to make this easier to understand later 
use mod_types,  only : sparse
use mod_inout,  only : sparse_to_full,full_to_sparse
use mod_linalg, only : my_dgemm
type(block),allocatable, intent(in)  :: blocks(:)
real(dp),allocatable,    intent(in)  :: bvects(:,:)
real(dp),allocatable,    intent(out)  :: atvects(:,:)
real(dp),allocatable :: bigproj(:,:)
integer :: row,col,i,j,k,l,ll,ier,nblocks,pi,iblock,fblock
integer :: natoms

nblocks = size(blocks)

natoms = 0
do i = 1, nblocks
  natoms = natoms + blocks(i)%natoms
end do

allocate(atvects(3*natoms,size(bvects,2)),bigproj(3*natoms,6*nblocks),stat=ier)
if(ier /= 0) stop "vect_proj> memory error"
atvects = zero

call gen_bigproj(blocks,bigproj)
!atvects = matmul(bigproj,bvects)
!print *, 'vect_proj dim bigproj', size(bigproj,1),size(bigproj,2)
!print *, 'vect_proj dim bvects' , size(bvects,1),size(bvects,2)
!print *, 'vect_proj dim atvects' , size(atvects,1),size(atvects,2)
call my_dgemm(bigproj,bvects,atvects,'n','n',1.0d0,0.0d0)

end subroutine

subroutine gen_bigproj(blocks,bigproj)
type(block),allocatable, intent(in)  :: blocks(:)
real(dp),allocatable,    intent(out) :: bigproj(:,:)
integer :: i,k,l,ll,pi,iblock,fblock,natoms,nblocks,ier

nblocks = size(blocks)
natoms = 0
do i = 1, nblocks
  natoms = natoms + blocks(i)%natoms
end do

allocate(bigproj(3*natoms,6*nblocks),stat=ier)
if(ier /= 0) stop "vect_proj> memory error"
bigproj = zero

do i = 1, nblocks
  iblock = 6*(i-1)+1
  fblock = iblock + 5
  do k = 1, size(blocks(i)%ixyz)
    ll = 1
    do l = iblock,fblock
      bigproj(blocks(i)%ixyz(k),l) = blocks(i)%pmat(ll,k)
      ll = ll + 1
    end do
  end do
end do

end subroutine

subroutine gen_zbigproj(blocks,zbigproj)
type(block),allocatable, intent(in)  :: blocks(:)
complex(kind=8),allocatable,    intent(out) :: zbigproj(:,:)
integer :: i,k,l,ll,pi,iblock,fblock,natoms,nblocks,ier

nblocks = size(blocks)
natoms = 0
do i = 1, nblocks
  natoms = natoms + blocks(i)%natoms
end do

allocate(zbigproj(3*natoms,6*nblocks),stat=ier)
if(ier /= 0) stop "vect_proj> memory error"
zbigproj = cmplx(zero,zero,kind=8)

do i = 1, nblocks
  iblock = 6*(i-1)+1
  fblock = iblock + 5
  do k = 1, size(blocks(i)%ixyz)
    ll = 1
    do l = iblock,fblock
      zbigproj(blocks(i)%ixyz(k),l) = cmplx(blocks(i)%pmat(ll,k),zero,kind=8)
      ll = ll + 1
    end do
  end do
end do

end subroutine

subroutine hess_proj(shess,blocks,bnmhess)
! projecting the hessian
!   this subrout doesn't compute hessian on the fly
!   leave that for later if needed
!
!  6x6 blocks of the projected hessian constructed as follows
!  pseudo
!   do i=1,nblocks
!     do j=i,nblocks
!       
!       Hb_ij = Pt_i H_ij P_j    ! where Pt_i is the transpose of iblock proj  
!                                ! and H_ij is the block of the hessian corr
!                                ! to proj block.  
!                                ! H_ij dimensions are 3nat_i x 3nat_j
!                                ! P_i  dimensions are 3nat_i x 6
!     
!     end do
!   end do
!
!   since indices are not necessarily continuous, we have debate:
!     1. copy in appropriate H_ij to some tmat with approp dimen
!        and use mkl libs
!     2. write subroutine to carry out multiplication directly
!
!   we'll try 1. first.
!
use mod_types,  only : sparse
use mod_linalg, only : my_dgemm
use mod_inout,  only : sparse_to_full,full_to_sparse
type(sparse)            ,   intent(in)  :: shess
type(block),allocatable,    intent(in)  :: blocks(:)
type(sparse),               intent(out) :: bnmhess
real(dp),allocatable :: fhess(:,:),phess(:,:),rphess(:,:)
real(dp),allocatable :: tmat(:,:)
integer :: i,j,ier,nblocks
integer :: iblock,jblock
integer :: ndb
real(dp) :: time1,time2

nblocks = size(blocks)
ndb = nblocks*6
call sparse_to_full(shess,fhess)

allocate(phess(ndb,ndb),stat=ier)
call cpu_time(time1)

iblock = 1

do i=1,nblocks
  jblock = iblock
  do j=i,nblocks
    call Pti_Hij_Pj(blocks(i),fhess,blocks(j),tmat)   
    phess(iblock:iblock+5,jblock:jblock+5) = tmat
    jblock = jblock+6
  end do
  iblock = iblock+6
end do

call full_to_sparse(phess,bnmhess)
deallocate(phess)
call cpu_time(time2)
print *, "hess_proj> time projecting",time2-time1

end subroutine

subroutine pti_hij_pj(blocki,hess,blockj,mat)
! DMR added December 7,2009 to speed up projection
use mod_linalg, only : my_dgemm
type(block),          intent(in)  :: blocki,blockj
real(dp),allocatable, intent(in)  :: hess(:,:)
real(dp),allocatable :: mat(:,:)
real(dp),allocatable              :: subhess(:,:),rphess(:,:)
integer :: i,j,indim,jndim,ier

indim = size(blocki%ixyz)
jndim = size(blockj%ixyz)
allocate(subhess(indim,jndim),stat=ier)
if(ier /= 0) stop "pti_hij_pj> memory error"

do i = 1, indim
  do j = 1,jndim
    subhess(i,j) = hess(blocki%ixyz(i),blockj%ixyz(j))
  end do
end do 

allocate(rphess(indim,6),stat=ier)
if(ier /= 0) stop "pti_hij_pj> memory error"

call my_dgemm(subhess,blockj%pmat,rphess,'n','t',1.0d0,0.0d0)
deallocate(subhess)
call my_dgemm(blocki%pmat,rphess,mat,'n','n',1.0d0,0.0d0)
deallocate(rphess)

end subroutine

subroutine pti_dij_pj(blocki,dynm,blockj,mat)
! DMR added December 7,2009 to speed up projection
use mod_linalg, only : my_zgemm
type(block),             intent(in)  :: blocki,blockj
complex(dp),allocatable, intent(in)  :: dynm(:,:)
complex(dp),allocatable :: mat(:,:)
complex(dp),allocatable :: subdynm(:,:),rpdynm(:,:),zpmati(:,:),zpmatj(:,:)
integer :: i,j,indim,jndim,ier
complex(dp) :: zone,zzero

zone  = cmplx(one,zero,kind = dp)
zzero = cmplx(zero,zero,kind = dp)

indim = size(blocki%ixyz)
jndim = size(blockj%ixyz)
allocate(subdynm(indim,jndim),stat=ier)
if(ier /= 0) stop "pti_hij_pj> memory error"
allocate(zpmati(indim,6),zpmatj(jndim,6),stat=ier)
if(ier /= 0) stop "pti_hij_pj> memory error"

zpmati = transpose(blocki%pmat) ! fix dimensions
zpmatj = transpose(blockj%pmat)

! copy in the subblock
do i = 1, indim
  do j = 1,jndim
    subdynm(i,j) = dynm(blocki%ixyz(i),blockj%ixyz(j))
!    write (666,*) blocki%ixyz(i),blockj%ixyz(j),real(subdynm(i,j)),dimag(subdynm(i,j))
  end do
end do 

allocate(rpdynm(indim,6),stat=ier)
if(ier /= 0) stop "pti_hij_pj> memory error"

call my_ZGEMM(subdynm,zpmatj,rpdynm,'n','n',zone,zzero)
deallocate(subdynm)
call my_ZGEMM(zpmati,rpdynm,mat,'t','n',zone,zzero)
deallocate(rpdynm)

end subroutine


subroutine hess_proj_old(shess,blocks,bnmhess)
use mod_types,  only : sparse
use mod_linalg, only : my_dgemm
use mod_inout,  only : sparse_to_full,full_to_sparse
type(sparse)            ,   intent(in)  :: shess
type(block),allocatable,    intent(in)  :: blocks(:)
type(sparse),    intent(out)  :: bnmhess
real(dp),allocatable :: fhess(:,:),phess2(:,:),phess1(:,:),phess3(:,:),bigproj(:,:)
integer :: row,col,i,j,ndim,k,l,ll,ier,nblocks,pi,iblock,fblock
integer :: natoms,ndb
real(dp) :: time1,time2,time3

nblocks = size(blocks)
ndim = shess%ndim
ndb = nblocks*6
call sparse_to_full(shess,fhess)

allocate(phess1(ndim,ndb),phess2(ndb,ndb),stat=ier)
if(ier /= 0) stop "hess_proj> memory error"
phess1 = zero
phess2 = zero
print *, "hess_proj> genbigproj"
call cpu_time(time1)
call gen_bigproj(blocks,bigproj)
call cpu_time(time2)
print *, "hess_proj> time generating bigproj,",time2-time1
call my_dgemm(fhess,bigproj,phess1,'n','n',1.0d0,0.0d0)
deallocate(fhess)
call my_dgemm(bigproj,phess1,phess2,'t','n',1.0d0,0.0d0)
deallocate(phess1)
print *, "hess_proj> convert bnmhess to sparse"
call full_to_sparse(phess2,bnmhess)
deallocate(phess2)

call cpu_time(time3)
print *, "hess_proj> time projecting",time3-time1

end subroutine

subroutine dynmat_proj(dynmat,blocks,bnmdyn)
! projecting the dynmat
!
!  6x6 blocks of the projected dynmat constructed as follows
!  pseudo
!   do i=1,nblocks
!     do j=i,nblocks
!       
!       Hb_ij = Pt_i D_ij P_j    ! where Pt_i is the transpose of iblock proj  
!                                ! and D_ij is the block of the dynmat corr
!                                ! to proj block.  
!                                ! H_ij dimensions are 3nat_i x 3nat_j
!                                ! P_i  dimensions are 3nat_i x 6
!     
!     end do
!   end do
!
!   since indices are not necessarily continuous, we have debate:
!     1. copy in appropriate D_ij to some tmat with approp dimen
!        and use mkl libs
!     2. write subroutine to carry out multiplication directly
!
use mod_types,  only : sparse
use mod_linalg, only : my_dgemm
use mod_inout,  only : zsparse_to_full,zfull_to_zsparse,&
                       sparse_upper_to_gen,sparse_deinit
type(sparse)            ,   intent(in)  :: dynmat
type(block),allocatable,    intent(in)  :: blocks(:)
type(sparse),    intent(out)  :: bnmdyn
type(sparse)                  :: tmpdyn
complex(dp),allocatable :: fdyn(:,:),pdynm(:,:)
complex(dp),allocatable :: tmat(:,:)
integer :: i,j,ier,nblocks
integer :: iblock,jblock
integer :: ndb
real(dp) :: time1,time2

nblocks = size(blocks)
ndb = nblocks*6
call zsparse_to_full(dynmat,fdyn)

allocate(pdynm(ndb,ndb),stat=ier)
call cpu_time(time1)

iblock = 1

do i=1,nblocks
  jblock = iblock
  do j=i,nblocks
    call Pti_dij_Pj(blocks(i),fdyn,blocks(j),tmat)   
    pdynm(iblock:iblock+5,jblock:jblock+5) = tmat
    jblock = jblock+6
  end do
  iblock = iblock+6
end do

call zfull_to_zsparse(pdynm,tmpdyn)
deallocate(pdynm)
call  sparse_upper_to_gen(tmpdyn,bnmdyn,"z")
call  sparse_deinit(tmpdyn)
call cpu_time(time2)
print *, "dyn_proj> time projecting",time2-time1

end subroutine

subroutine dynmat_proj_old(dynmat,blocks,bnmdyn)
use mod_types,  only : sparse
use mod_linalg, only : my_zgemm
use mod_inout,  only : zsparse_to_full,zfull_to_zsparse,sparse_upper_to_gen,sparse_deinit
type(sparse)            ,   intent(in out)  :: dynmat
type(block),allocatable,    intent(in)  :: blocks(:)
type(sparse),    intent(out)  :: bnmdyn
INTEGER, PARAMETER :: WP = KIND(1.0D0)
complex(wp),allocatable :: fdyn(:,:),fbnmdyn(:,:),tmp(:,:),zbigproj(:,:)
real(dp),   allocatable :: bigproj(:,:)
integer :: row,col,i,j,ndim,k,l,ll,ier,nblocks,pi,iblock,fblock
integer :: natoms
type(sparse) :: tmpdyn
complex(wp) :: zone,zzero
real(dp) :: time1,time2

zone  = cmplx(one,zero,kind = wp)
zzero = cmplx(zero,zero,kind = wp)

nblocks = size(blocks)
call zsparse_to_full(dynmat,fdyn)
call sparse_deinit(dynmat)

call cpu_time(time1)
call gen_bigproj(blocks,bigproj)

allocate(fbnmdyn(6*nblocks,6*nblocks),stat=ier)
if(ier /= 0) stop "dynmat_proj> memory error"

allocate(tmp(size(fdyn,1),size(bigproj,2)),stat=ier)
if(ier /= 0) stop "dynmat_proj> memory error"

fbnmdyn = cmplx(zero,zero)
tmp     = cmplx(zero,zero)

call dp_zp(bigproj,zbigproj)
call my_ZGEMM(fdyn,zbigproj,tmp,'n','n',zone,zzero)
deallocate (fdyn)
call my_ZGEMM(zbigproj,tmp,fbnmdyn,'t','n',zone,zzero)
!call ZGEMM_MKL95(zbigproj,tmp,fbnmdyn,'t','n',zone,zzero)
deallocate(zbigproj)
deallocate(tmp)
call zfull_to_zsparse(fbnmdyn,tmpdyn)
deallocate(fbnmdyn)
call  sparse_upper_to_gen(tmpdyn,bnmdyn,"z")
call  sparse_deinit(tmpdyn)

end subroutine

subroutine dp_zp(dmat,zmat)
real(dp), allocatable, intent(inout) :: dmat(:,:)
complex(kind=8), allocatable, intent(out) :: zmat(:,:)
integer :: ier,i,j

allocate(zmat(size(dmat,1),size(dmat,2)), stat = ier)
if (ier /= 0) stop "dp_zp> malloc"

do i = 1, size(dmat,1)
  do j = 1, size(dmat,2)
    zmat(i,j) = cmplx(dmat(i,j),kind=8)
  end do
end do

deallocate(dmat)

end subroutine

subroutine vcov_proj(bvcov,blocks,vcov)
use mod_inout
use mod_linalg, only : my_dgemm
real(dp),allocatable,   intent(in)      :: bvcov(:,:)
type(block),allocatable,    intent(in)  :: blocks(:)
real(dp),allocatable,   intent(out)     :: vcov(:,:)
real(dp),allocatable                    :: sub_vcov(:,:)
real(dp),allocatable :: sub_bvcov(:,:)
integer :: natoms, nblocks, iblock,jblock 
integer :: i,j,k,l,ier,cdim 

natoms = 0
nblocks = size(blocks)
do i = 1, nblocks
  natoms = natoms + blocks(i)%natoms
end do
cdim = 3*natoms

allocate(vcov(cdim,cdim),sub_bvcov(6,6),stat=ier)
if(ier /= 0) stop "vcov_proj> memory error"

iblock = 1
do i=1,nblocks

  jblock = iblock
  do j=i,nblocks

    sub_bvcov = bvcov(iblock:iblock+5,jblock:jblock+5)

    call Pi_vcov_Ptj(blocks(i),sub_bvcov,blocks(j),sub_vcov)

!  copy the block into the right spot
    do k =1, size(blocks(i)%ixyz)
      do l =1, size(blocks(j)%ixyz)
        vcov(blocks(i)%ixyz(k),blocks(j)%ixyz(l)) = sub_vcov(k,l)
      end do
    end do

    jblock = jblock+6
  end do
  iblock = iblock+6
end do

! fill up the lower diagonal
do i = 1, cdim
  do j = i+1, cdim
    vcov(j,i) = vcov(i,j)
  end do
end do

end subroutine

subroutine pi_vcov_ptj(blocki,bvcov,blockj,vcov)
use mod_linalg, only : my_dgemm
type(block),          intent(in)   :: blocki,blockj
real(dp),allocatable, intent(in)   :: bvcov(:,:)
real(dp),allocatable, intent(out)  ::  vcov(:,:)
real(dp),allocatable               :: rpvcov(:,:)
real(dp), allocatable              :: pmati(:,:),pmatj(:,:)
integer :: i,j,indim,jndim,ier

indim = size(blocki%ixyz)
jndim = size(blockj%ixyz)
allocate(vcov(indim,jndim),pmati(indim,6),pmatj(jndim,6),stat=ier)
if(ier /= 0) stop "pti_vcov_pj> memory error"

pmati = transpose(blocki%pmat)  ! fix dimensions
pmatj = transpose(blockj%pmat) 

allocate(rpvcov(6,jndim),stat=ier)
if(ier /= 0) stop "pti_hij_pj> memory error"

call my_dgemm(bvcov,pmatj,rpvcov,'n','t',1.0d0,0.0d0)
call my_dgemm(pmati,rpvcov,vcov ,'n','n',1.0d0,0.0d0)
deallocate(rpvcov)

end subroutine

subroutine vcov_proj_old(svcov,blocks,bvcov)
use mod_inout
use mod_linalg, only : my_dgemm
real(dp),allocatable,   intent(in)      :: svcov(:,:)
type(block),allocatable,    intent(in)  :: blocks(:)
real(dp),allocatable,   intent(out)     :: bvcov(:,:)
real(dp),allocatable :: tmp(:,:),bigproj(:,:)
integer :: row,col,i,j,ndim,k,l,ll,ier,nblocks,iblock,fblock
integer :: natoms,cdim,bdim

natoms = 0
nblocks = size(blocks)
do i = 1, nblocks
  natoms = natoms + blocks(i)%natoms
end do
cdim = 3*natoms

allocate(bvcov(cdim,cdim),stat=ier)
allocate(tmp(6*size(blocks),cdim),stat=ier)
if(ier /= 0) stop "vcov_proj> memory error"

call gen_bigproj(blocks,bigproj)
call my_dgemm(svcov,bigproj,tmp,'n','t',1.0d0,0.0d0)
call my_dgemm(bigproj,tmp,bvcov,'n','n',1.0d0,0.0d0)
!bvcov = matmul(bigproj,matmul(svcov,transpose(bigproj)))
!call gemm(bigproj,bvects,atvects,'n','n',1.0d0,0.0d0)

!call matrx2d_write(101,bigproj)
!call matrx2d_write(102,bvcov)
!call matrx2d_write(103,svcov)

end subroutine

subroutine blocks_setup_iso(input,atoms_in,blocks)
use mod_types,  only : inp_par
use mod_linalg, only: lapack_eig3
type(inp_par),            intent(in)    :: input
type(protein_atom_list_type),     intent(inout) :: atoms_in
type(block), allocatable, intent(out)   :: blocks(:)
integer :: i

call set_blocks(input,atoms_in,blocks)
if (input%block_read) then 

  if (size(input%block_cents(:,1)) .ne. size(blocks))  stop "something fishy with block definitions"
  do i = 1, size(input%block_cents(:,1))
    blocks(i)%rotc = input%block_cents(i,:)
  end do
!call rad_block_centers(input,blocks)
end if
call blocks_com(input,atoms_in,blocks)

do i = 1, size(blocks)
  if (blocks(i)%cent_read) then
     blocks(i)%cent = blocks(i)%rotc
  else 
     blocks(i)%cent = blocks(i)%com
  end if
end do

call blocks_tnsi(input,atoms_in,blocks)

do i = 1, size(blocks)
  call lapack_eig3(blocks(i)%tnsi,blocks(i)%momI,blocks(i)%eigi)
end do

  call tr_proj_schmidt(input,atoms_in,blocks,.true.) 
!call tr_proj(input,atoms_in,blocks) 

end subroutine

subroutine blocks_setup_pbc(input,asyms,blocks)
use mod_types,  only : inp_par,asym_list_type
use mod_linalg, only: lapack_eig3
type(inp_par),            intent(in)    :: input
type(Asym_list_Type),     intent(inout) :: asyms
type(block), allocatable, intent(out)   :: blocks(:)
integer :: i

call set_blocks(input,asyms%aunit(1),blocks)
if (input%block_read) then 

  if (size(input%block_cents(:,1)) .ne. size(blocks))  stop "something fishy with block definitions"
  do i = 1, size(input%block_cents(:,1))
    blocks(i)%rotc = input%block_cents(i,:)
  end do
!call rad_block_centers(input,blocks)
end if
call blocks_com(input,asyms,blocks)

do i = 1, size(blocks)
  if (blocks(i)%cent_read) then
     blocks(i)%cent = blocks(i)%rotc
  else 
     blocks(i)%cent = blocks(i)%com
  end if
  !write(655,'(I5,3F10.4)') i,blocks(i)%cent
end do

call blocks_tnsi(input,asyms,blocks)

do i = 1, size(blocks)
  call lapack_eig3(blocks(i)%tnsi,blocks(i)%momI,blocks(i)%eigi)
end do

call tr_proj_schmidt_pbc(input,asyms,blocks,.true.) 
!call tr_proj(input,atoms_in,blocks) 

end subroutine

subroutine tr_proj(input,atoms_in, blocks)
use mod_types,     only : inp_par
use mod_crysbuild, only : atom_mass_constr
use mod_math,      only : dist_sqr
use mod_linalg,    only : root_invrse_mat,schmidt
type(inp_par),            intent(in)    :: input
type(protein_atom_list_type),     intent(in)    :: atoms_in
type(block), allocatable, intent(inout) :: blocks(:)
real(dp),    allocatable :: masses(:)
real(dp) :: shit(6,6),dxyz_com(3),irtmat(3,3)
real(dp) :: x,y,z,xc,yc,zc,momi1,momi2,momi3
integer :: i,j,ii,jj,ier,k,kk,l

call atom_mass_constr (input,atoms_in,masses)

do i = 1, size(blocks)
  call root_invrse_mat(blocks(i)%tnsi,blocks(i)%eigI,irtmat)
  jj = 1
  do j = 1, atoms_in%natoms
    if (atoms_in%atom(j)%iblock .eq. i) then
      dxyz_com  = atoms_in%atom(j)%x - blocks(i)%cent
      xc = dxyz_com(1)
      yc = dxyz_com(2)
      zc = dxyz_com(3)
      blocks(i)%pmat(1,jj)   = dsqrt(masses(j)/blocks(i)%tmass) 
      blocks(i)%pmat(2,jj+1) = dsqrt(masses(j)/blocks(i)%tmass) 
      blocks(i)%pmat(3,jj+2) = dsqrt(masses(j)/blocks(i)%tmass)
      ! compute rotational part as long as nat > 1
      if (blocks(i)%natoms > 1) then
        ! x part
        blocks(i)%pmat(4,jj)   = sqrt(masses(j))*(irtmat(1,2)*zc - irtmat(1,3)*yc)
        blocks(i)%pmat(5,jj)   = sqrt(masses(j))*(irtmat(2,2)*zc - irtmat(2,3)*yc)
        blocks(i)%pmat(6,jj)   = sqrt(masses(j))*(irtmat(3,2)*zc - irtmat(3,3)*yc)
        ! y part
        blocks(i)%pmat(4,jj+1) = sqrt(masses(j))*(irtmat(1,3)*xc - irtmat(1,1)*zc)
        blocks(i)%pmat(5,jj+1) = sqrt(masses(j))*(irtmat(2,3)*xc - irtmat(2,1)*zc)
        blocks(i)%pmat(6,jj+1) = sqrt(masses(j))*(irtmat(3,3)*xc - irtmat(3,1)*zc)
        ! z part
        blocks(i)%pmat(4,jj+2) = sqrt(masses(j))*(irtmat(1,1)*yc - irtmat(1,2)*xc)
        blocks(i)%pmat(5,jj+2) = sqrt(masses(j))*(irtmat(2,1)*yc - irtmat(2,2)*xc)
        blocks(i)%pmat(6,jj+2) = sqrt(masses(j))*(irtmat(3,1)*yc - irtmat(3,2)*xc)
      end if
      jj = jj+3
    end if
  end do
end do
! do i = 1, size(blocks)
!   shit = matmul(blocks(i)%pmat,transpose(blocks(i)%pmat))
!   ier = 665+i
!   do k = 1,6
!     do l = 1,6
!       write(ier,'(2I4,F10.4)') k,l,shit(k,l)
!     end do
!   end do
! end do

end subroutine

subroutine atdxyz_cent_asy (asy_uc,blocks, &
                            atdxyz)
use mod_math,  only: distance
use mod_types, only: asym_list_type
! compute the dx dy dz dist for atom to block center under pbc 
! inout
type(asym_List_Type),intent(in)    :: asy_uc
type(block), allocatable, intent(in) :: blocks(:)
real(dp),allocatable ,intent(out)  :: atdxyz(:,:,:)
real(dp),dimension(3) :: avec,bvec,cvec,vec
integer                           :: i,ii,jj,k,kk,j,h,natoms,nneighs,naunits,ialloc,nblocks
real(sp)                          :: time1, time2
integer :: ja,jb,jc
real(dp)                          :: dist, xyzi(3), xyzj(3),dxyz(3)


nblocks = size(blocks)
nneighs = asy_uc%nneighs
natoms  = asy_uc%aunit(1)%natoms
naunits = asy_uc%naunits

allocate(atdxyz(nblocks,natoms,4),stat=ialloc)
if(ialloc /= 0) stop "atdxyz_cent_asy> malloc"

atdxyz=9999.999d0

do i = 1, size(blocks)
  xyzi = blocks(i)%cent(:)
  do j  = 1, naunits
  do jj = 1, natoms
    !if (asy_uc%aunit(j)%atom(jj)%iblock .eq. i) then
      do ja=asy_uc%aneigh(1),asy_uc%aneigh(2)
      avec = real(ja,kind=dp)*asy_uc%avec
      do jb=asy_uc%bneigh(1),asy_uc%bneigh(2)
      bvec = real(jb,kind=dp)*asy_uc%bvec
      do jc=asy_uc%cneigh(1),asy_uc%cneigh(2)
      cvec = real(jc,kind=dp)*asy_uc%cvec
      vec = avec + bvec + cvec

      xyzj = asy_uc%aunit(j)%atom(jj)%x+vec
      dxyz = xyzj-xyzi
      dist = distance(dxyz)

      if (dist .le. atdxyz(i,jj,4)) then
        atdxyz(i,jj,4)  = dist
        atdxyz(i,jj,1:3) = dxyz 
      end if

      end do
      end do
      end do
    end do
  end do
end do

end subroutine

subroutine tr_proj_schmidt_pbc(input,asyms,blocks,qschmidt)
use mod_types,     only : inp_par,asym_list_type
use mod_crysbuild, only : atom_mass_constr
use mod_math,      only : dist_sqr
use mod_linalg,    only : schmidt
type(inp_par),            intent(in)    :: input
type(Asym_list_Type),     intent(in)    :: asyms
type(block), allocatable, intent(inout) :: blocks(:)
logical, optional,        intent(in)    :: qschmidt      
real(dp),    allocatable :: masses(:),atxyz(:,:,:)
real(dp) :: dx,dy,dz
integer :: i,j,ii,jj,ier,k,kk,l


call atom_mass_constr (input,asyms%aunit(1),masses)
call atdxyz_cent_asy(asyms,blocks,atxyz)

do i = 1, size(blocks)
  jj = 1
  do j = 1, asyms%aunit(1)%natoms
    if (asyms%aunit(1)%atom(j)%iblock .eq. i) then
      dx = atxyz(i,j,1)
      dy = atxyz(i,j,2)
      dz = atxyz(i,j,3)
      ! diagonal mass part
      blocks(i)%pmat(1,jj)   = dsqrt(masses(j)) 
      blocks(i)%pmat(2,jj+1) = dsqrt(masses(j)) 
      blocks(i)%pmat(3,jj+2) = dsqrt(masses(j))
      ! compute rotational part as long as natoms > 1
      if (blocks(i)%natoms .eq. 1) then
        print *, "one atom in block:",i, "setting diag to one"
        blocks(i)%pmat(1,jj)   = one
        blocks(i)%pmat(2,jj+1) = one
        blocks(i)%pmat(3,jj+2) = one
      end if
!It's clear that in calculating the the elements (as done there and in my code) they would be zero, but this could trip up the Schmidt-orthog step.
      if (blocks(i)%natoms .eq. 2) print *, "WARNING>: TWO ATOMS IN BLOCK:",i, "you may want to split this up"
      if (blocks(i)%natoms > 1) then
        ! rot1  
        blocks(i)%pmat(4,jj)     =  zero
        blocks(i)%pmat(4,jj+1)   = -dz*dsqrt(masses(j))
        blocks(i)%pmat(4,jj+2)   =  dy*dsqrt(masses(j)) 
        ! rot2  
        blocks(i)%pmat(5,jj)     =  dz*dsqrt(masses(j))
        blocks(i)%pmat(5,jj+1)   =  zero
        blocks(i)%pmat(5,jj+2)   = -dx*dsqrt(masses(j)) 
        ! rot3  
        blocks(i)%pmat(6,jj)     = -dy*dsqrt(masses(j))
        blocks(i)%pmat(6,jj+1)   =  dx*dsqrt(masses(j))
        blocks(i)%pmat(6,jj+2)   =  zero
      end if
      jj = jj+3
    end if
  end do
  if (present(qschmidt)) then
    if (blocks(i)%natoms > 1) then
      if(qschmidt) call schmidt(blocks(i)%pmat)
    end if
  end if

end do

print *, "projection matrix set"

! do i = 1, size(blocks)
!   ier = 30+i
!   write(ier,'(3F10.5)') blocks(i)%cent
!   shit = matmul(blocks(i)%pmat,transpose(blocks(i)%pmat))
!   do k = 1,6
!     do l = 1,6
!       write(ier,'(2I4,F10.4)') k,l,shit(k,l)
!    end do
!   end do
! end do



end subroutine

subroutine tr_proj_schmidt(input,atoms_in, blocks,qschmidt)
use mod_types,     only : inp_par
use mod_crysbuild, only : atom_mass_constr
use mod_math,      only : dist_sqr
use mod_linalg,    only : schmidt
type(inp_par),            intent(in)    :: input
type(protein_atom_list_type),     intent(in)    :: atoms_in
type(block), allocatable, intent(inout) :: blocks(:)
logical, optional,        intent(in)    :: qschmidt      
real(dp),    allocatable :: masses(:)
real(dp) :: x,y,z,xc,yc,zc,shit(6,6)
integer :: i,j,ii,jj,ier,k,kk,l


call atom_mass_constr (input,atoms_in,masses)

do i = 1, size(blocks)
  xc   = blocks(i)%cent(1)
  yc   = blocks(i)%cent(2)
  zc   = blocks(i)%cent(3)
  jj = 1
  do j = 1, atoms_in%natoms
    if (atoms_in%atom(j)%iblock .eq. i) then
      ! diagonal mass part
      x = atoms_in%atom(j)%x(1)
      y = atoms_in%atom(j)%x(2)
      z = atoms_in%atom(j)%x(3)
      blocks(i)%pmat(1,jj)   = dsqrt(masses(j)) 
      blocks(i)%pmat(2,jj+1) = dsqrt(masses(j)) 
      blocks(i)%pmat(3,jj+2) = dsqrt(masses(j))
      ! compute rotational part as long as natoms > 1
      if (blocks(i)%natoms .eq. 1) then
        print *, "one atom in block:",i, "setting diag to one"
        blocks(i)%pmat(1,jj)   = one
        blocks(i)%pmat(2,jj+1) = one
        blocks(i)%pmat(3,jj+2) = one
      end if
!It's clear that in calculating the the elements (as done there and in my code) they would be zero, but this could trip up the Schmidt-orthog step.
      if (blocks(i)%natoms .eq. 2) print *, "WARNING>: TWO ATOMS IN BLOCK:",i, "you may want to split this up"
      if (blocks(i)%natoms > 1) then
        ! x part
        !blocks(i)%pmat(4,jj)   =  zero
        !blocks(i)%pmat(5,jj)   = -(z-zc)*sqrt(masses(j))
        !blocks(i)%pmat(6,jj)   =  (y-yc)*sqrt(masses(j)) 
        ! y part
        !blocks(i)%pmat(4,jj+1) =  (z-zc)*sqrt(masses(j))
        !blocks(i)%pmat(5,jj+1) =  zero
        !blocks(i)%pmat(6,jj+1) = -(x-xc)*sqrt(masses(j)) 
        ! z part
        !blocks(i)%pmat(4,jj+2) = -(y-yc)*sqrt(masses(j))
        !blocks(i)%pmat(5,jj+2) =  (x-xc)*sqrt(masses(j))
        !blocks(i)%pmat(6,jj+2) =  zero
        ! rot1  
        blocks(i)%pmat(4,jj)     =  zero
        blocks(i)%pmat(4,jj+1)   = -(z-zc)*dsqrt(masses(j))
        blocks(i)%pmat(4,jj+2)   =  (y-yc)*dsqrt(masses(j)) 
        ! rot2  
        blocks(i)%pmat(5,jj)   =  (z-zc)*dsqrt(masses(j))
        blocks(i)%pmat(5,jj+1) =  zero
        blocks(i)%pmat(5,jj+2) = -(x-xc)*dsqrt(masses(j)) 
        ! rot3  
        blocks(i)%pmat(6,jj)   = -(y-yc)*dsqrt(masses(j))
        blocks(i)%pmat(6,jj+1) =  (x-xc)*dsqrt(masses(j))
        blocks(i)%pmat(6,jj+2) =  zero
      end if
      jj = jj+3
    end if
  end do
  if (present(qschmidt)) then
    if (blocks(i)%natoms > 1) then
      if(qschmidt) call schmidt(blocks(i)%pmat)
    end if
  end if

end do

print *, "projection matrix set"

! do i = 1, size(blocks)
!   ier = 30+i
!   write(ier,'(3F10.5)') blocks(i)%cent
!   shit = matmul(blocks(i)%pmat,transpose(blocks(i)%pmat))
!   do k = 1,6
!     do l = 1,6
!       write(ier,'(2I4,F10.4)') k,l,shit(k,l)
!    end do
!   end do
! end do



end subroutine

subroutine read_blocks_atom(input,atoms_in,blockfile)
! this subroutine initiates the blocks
! and assigns block ids to atoms
! currently this information is read in from a file $fileroot.block 
! this will need to be changed, but it's ok for now
use mod_types, only: inp_par
type(inp_par),           intent(in out)  :: input
type (protein_atom_list_type),   intent(in out)  :: atoms_in
character(len=*), optional, intent(in)   :: blockfile
integer :: i,j,blck,at,nat,ier,ires,maxblock,ichain
logical :: pass
character(len=24) :: blck_or_at
character(len=50)  :: blockname 

if (present(blockfile)) then
  blockname = blockfile
else
  blockname = trim(input%fileroot) // ".block"
end if

open(unit=11,file=blockname, &
   status='old',action='read',position='rewind', &
   iostat=ier)

if (ier > 0) then 
  print *, "read_blocks_atom> cannot read block description file"
  print *, "read_blocks_atom> setting all atoms to same block/chain"
  do i = 1, atoms_in%natoms
    atoms_in%atom(i)%iblock = 1
    atoms_in%atom(i)%ichain = 1
    atoms_in%atom(i)%ires   = i
  end do
  return
end if

nat = 0
maxblock = 1 
do
  read(11,*,iostat = ier) blck_or_at,at,ires,blck,ichain
  if (trim(blck_or_at) .eq. "block") cycle
  if (ier > 0) stop "*** BLOCK INPUT ERROR"
  if (ier < 0) exit ! end of file
  nat = nat+1
  if (maxblock .lt. blck) maxblock = blck
  if (nat .gt. atoms_in%natoms) then
   print *, "got up to", nat, "block entries"
   print *, "number of atoms:",atoms_in%natoms
   stop "too many entries in file.block"
  end if
  atoms_in%atom(at)%iblock = blck
  atoms_in%atom(at)%ichain = ichain
  atoms_in%atom(at)%ires   = ires
  atoms_in%atom(at)%resn   = trim(blck_or_at)
end do
close(11)

input%nblocks = maxblock

! let's test the block ids
  ! any zeroes?
  do i = 1, atoms_in%natoms
    if (atoms_in%atom(i)%iblock .eq. 0) then
      if (i .eq. 1) then
        atoms_in%atom(i)%iblock = atoms_in%atom(i+1)%iblock
      else
        atoms_in%atom(i)%iblock = atoms_in%atom(i-1)%iblock
      end if 
      
      print *, 'no block defined for atom:',i
      print *, 'setting block to same as previous one in sequence'
      print *, 'block for atom:',i,' is now:',atoms_in%atom(i)%iblock
    end if
  end do

blocks: do j = 1, maxblock
  pass = .false.
  atoms: do i = 1, atoms_in%natoms
    if (atoms_in%atom(i)%iblock .eq. j) then
      pass = .true.
      exit atoms
    end if
  end do atoms
  if (pass) then
    cycle blocks
  else
    print *, "no block:", j, " something wrong with block def file"
    stop
  end if
end do blocks

if (input%block_read) then
  allocate(input%block_cents(maxblock,3), stat = ier)
  if (ier /= 0) stop "block_cents> malloc"
  input%block_cents = zero

! read in the centers into atoms_in
    call read_block_centers_new(input)
end if

end subroutine

subroutine read_block_centers_new(input)
use mod_types, only: inp_par
type(inp_par), intent(in out)                :: input
integer :: blck,nblocks,ier,i,j,nmultip
character(len=24) :: blck_or_at
real :: center(3)

print *, "reading in the centers"

open(unit=11,file=trim(input%fileroot) // ".block", &
     status='old',action='read',position='rewind', &
     iostat=ier)

if (ier > 0) then
  print *, "read_blocks_atom> cannot read block description file"
  print *, "read_blocks_atom> leave all block reaction centers as is"
  return
end if

nblocks = 0
do
  read(11,*,iostat = ier) blck_or_at,blck,center
  if (ier .eq. 10) cycle ! read error
  if (ier /= 0) exit
  if (trim(blck_or_at) .eq. "block") then
    nblocks = nblocks+1
    if (nblocks .gt. size(input%block_cents(:,1))) then
      print *, "too many block centers defined in file"
      print *, "using the first", size(input%block_cents(:,1)), "blocks"
      nblocks = nblocks - 1
      exit
      return
    end if
    input%block_cents(blck,:) = center
  end if
end do

close(11)

end subroutine

subroutine read_block_centers(input,blocks)
use mod_types, only: inp_par
! this subroutine reads the centers for each block
! it will look for the .block file, if there it will set the reaction centers
! from the file... if not, it will set rotc to center of mass
type(inp_par), intent(in)                :: input
type(block), allocatable, intent(in out) :: blocks(:)
integer :: blck,nblocks,ier,i,j,nmultip
character(len=24) :: blck_or_at
real(dp) :: center(3)
real :: cent(3)
real, allocatable :: orb(:,:)

stop

open(unit=11,file=trim(input%fileroot) // ".block", &
     status='old',action='read',position='rewind', &
     iostat=ier)

if (ier > 0) then 
  print *, "read_blocks_atom> cannot read block description file"
  print *, "read_blocks_atom> leave all block reaction centers as is"
  return
end if

nblocks = 0
do
  read(11,*,iostat = ier) blck_or_at,blck,center
  if (ier .eq. 10) cycle ! read error
  if (ier /= 0) exit 
  if (trim(blck_or_at) .eq. "block") then
    nblocks = nblocks+1
    if (nblocks .gt. size(blocks)) then
      print *, "too many block centers defined in file"
      print *, "using the first", size(blocks), "blocks"
      nblocks = nblocks - 1
      exit
      return
    end if
    blocks(blck)%rotc = center
  end if
end do

! note this will only set cent_read equal to .true. for those blocks found
!  so others will default to the center of mass
do i = 1, nblocks
  blocks(i)%cent_read = .true.
end do

if (nblocks .lt. size(blocks)) then
  print *, "WARNING> fewer centers than blocks defined in file"
  print *, "  those centers not defined will defaul to block center of mass"
end if

end subroutine

subroutine bnm_eigenal_valvec(input,matx,vals,vecs)
use mod_inout,     only: sparse_to_full,sparse_deinit
use mod_arpack,    only: mod_arpack_dsdrv1_sparse
use mod_linalg,    only: lapack_eig,mod_linalg_pseudo
use mod_types,     only: inp_par
type(inp_par), intent(in out) :: input
type(sparse), intent(in)     :: matx
real(dp), allocatable                :: matx_full(:,:)
real(dp), allocatable,intent(out)    :: vals(:),vecs(:,:)
real(dp), allocatable                :: valsthz(:)
integer                              :: ialloc
integer :: nfrq

if (allocated(vals)) deallocate(vals)
if (allocated(vecs)) deallocate(vecs)

allocate(vals(input%nfrqs), &
         vecs(matx%ndim,input%nfrqs),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"
if (trim(input%eigtyp) .eq. "lapack") then
  call sparse_to_full(matx,matx_full)
  call lapack_eig(matx_full, vals, vecs )
else
  call mod_arpack_dsdrv1_sparse(input,input%nfrqs,matx,vals,vecs)
end if

end subroutine

subroutine set_blocks(input,atoms_in,blocks)
! this subroutine initiates the blocks
! and assigns block ids to atoms
use mod_types, only: inp_par
type(inp_par),            intent(in)  :: input
type (protein_atom_list_type),   intent(in out)  :: atoms_in
type(block), allocatable, intent(out) :: blocks(:)
integer :: i,j,k,ii,jj,nblocks,ier,at,nat,blck,max_blck
real(dp) :: center(3), moment_I

allocate(blocks(input%nblocks),stat=ier)

print *, 'number of blocks!',input%nblocks

do i = 1, input%nblocks
  call block_init(blocks(i))
end do 

do i = 1, atoms_in%natoms
  blck = atoms_in%atom(i)%iblock
  blocks(blck)%natoms = blocks(blck)%natoms + 1 
end do

! set up atom ids
do i = 1, size(blocks)
  allocate(blocks(i)%atinblock(blocks(i)%natoms), stat=ier)
  if (ier /= 0) stop "set_blocks> memory err"
  allocate(blocks(i)%ixyz(3*blocks(i)%natoms), stat=ier)
  if (ier /= 0) stop "set_blocks> memory err"
  allocate(blocks(i)%pmat(6,3*blocks(i)%natoms), stat=ier)
  if (ier /= 0) stop "set_blocks> memory err"
  blocks(i)%pmat = zero
  ii = 1
  jj = 1
  do j = 1, atoms_in%natoms
    if (atoms_in%atom(j)%iblock .eq. i) then
      blocks(i)%atinblock(ii) = j
!     set up 3n coordinates
      do k = 0,2
        blocks(i)%ixyz(jj+k) = 3*(j-1)+1+k
      end do
      ii = ii + 1
      jj = jj + 3
    end if
  end do
end do

! 05122009 checked, and works
end subroutine

subroutine blocks_com_iso(input,atoms_in,blocks)
use mod_crysbuild, only : atom_mass_constr
use mod_math,      only : dist_sqr
use mod_types,     only : inp_par
type(inp_par),            intent(in)    :: input
type(protein_atom_list_type),     intent(in)    :: atoms_in
type(block), allocatable, intent(inout) :: blocks(:)
real(dp), allocatable               :: masses(:)
integer   :: i

call atom_mass_constr (input,atoms_in,masses)

! compute center of mass for each atom
do i =1, atoms_in%natoms
  blocks(atoms_in%atom(i)%iblock)%com = blocks(atoms_in%atom(i)%iblock)%com &
                                        + atoms_in%atom(i)%x*masses(i)
  blocks(atoms_in%atom(i)%iblock)%tmass = blocks(atoms_in%atom(i)%iblock)%tmass + masses(i)
end do

do i = 1, size(blocks)
  blocks(i)%com = blocks(i)%com/blocks(i)%tmass
end do

end subroutine



subroutine blocks_com_pbc(input,asyms,blocks)
use mod_crysbuild, only : atom_mass_constr
use mod_math,      only : distance
use mod_types,     only : inp_par,asym_list_type
type(inp_par),            intent(in)    :: input
type(Asym_list_Type),     intent(in)    :: asyms
type(block), allocatable, intent(inout) :: blocks(:)
real(dp), allocatable               :: masses(:),tmpcom(:,:)
real(dp) :: dist,mindist,mindxyz(3),vec(3),avec(3),bvec(3),cvec(3),xyzj(3),dxyz(3)
integer   :: i,j,natoms,naunits,ier
integer   :: ja,jb,jc

natoms  = asyms%aunit(1)%natoms
naunits = asyms%naunits

allocate(tmpcom(size(blocks),3), stat = ier)

call atom_mass_constr (input,asyms%aunit(1),masses)

! we have to find the right block center... maybe not really...  hmm.
! let's compute the center of mass relative to first atom in the block

! compute center of mass for each block which may be split across the unit  
print *, "blocks_com_pbc> correct?"
! set com of each block to the location of the first atom, and save this vector
blockloop: do i = 1, size(blocks)
  do j = 1,natoms 
    if (asyms%aunit(1)%atom(j)%iblock .eq. i) then
     ! tmpcom(i,:)blocks(i)%com = asyms%aunit(1)%atom(1)%x
      tmpcom(i,:) = asyms%aunit(1)%atom(j)%x
     ! write(555,'(I5,3F10.3)') i,j,tmpcom(i,:)
      cycle blockloop
    end if
  end do
end do blockloop

do i =1, natoms

  mindist = 9999.999d0
  do j  = 1, naunits
  do ja=asyms%aneigh(1),asyms%aneigh(2)
  avec = real(ja,kind=dp)*asyms%avec
  do jb=asyms%bneigh(1),asyms%bneigh(2)
  bvec = real(jb,kind=dp)*asyms%bvec
  do jc=asyms%cneigh(1),asyms%cneigh(2)
  cvec = real(jc,kind=dp)*asyms%cvec
  vec = avec + bvec + cvec

  xyzj = asyms%aunit(j)%atom(i)%x+vec
  dxyz = xyzj - tmpcom(asyms%aunit(1)%atom(i)%iblock,:)
  dist = distance(dxyz)

  if (dist .lt. mindist) then
    mindist = dist
    mindxyz = dxyz
  end if

  end do
  end do
  end do
  end do

  blocks(asyms%aunit(1)%atom(i)%iblock)%com = blocks(asyms%aunit(1)%atom(i)%iblock)%com &
                                        + mindxyz*masses(i)
  blocks(asyms%aunit(1)%atom(i)%iblock)%tmass = blocks(asyms%aunit(1)%atom(i)%iblock)%tmass + masses(i)

end do

do i = 1, size(blocks)
  blocks(i)%com = (blocks(i)%com/blocks(i)%tmass) + tmpcom(i,:)
end do

end subroutine

subroutine blocks_tnsi_pbc(input,asyms,blocks)
use mod_crysbuild, only : atom_mass_constr
use mod_math,      only : dist_sqr,distance
use mod_types,     only : inp_par,asym_list_type
type(inp_par),            intent(in)    :: input
type(Asym_list_Type),     intent(in)    :: asyms
type(block), allocatable, intent(inout) :: blocks(:)
real(dp), allocatable               :: masses(:),atxyz(:,:,:)
integer   :: i,j
real(dp)  :: dxyz(3),xxi,yyi,zzi,xyi,xzi,yzi,dist

call atom_mass_constr (input,asyms%aunit(1),masses)

call atdxyz_cent_asy(asyms,blocks,atxyz)

do i=1,asyms%aunit(1)%natoms
  j=asyms%aunit(1)%atom(i)%iblock  
  dxyz = atxyz(j,i,1:3)
  dist = atxyz(j,i,4) 
  if (dist .gt. 50.0d0) then
    print *, 'WARNING> atom ',i, 'is far from block center:', asyms%aunit(1)%atom(i)%iblock
    print *, 'cent:',blocks(asyms%aunit(1)%atom(i)%iblock)%cent
    print *, 'dist:',dist
    print *, 'WARNING> are blocks defined correctly?'
  end if
  xxi   =  masses(i)*(dxyz(2)*dxyz(2) + dxyz(3)*dxyz(3)) 
  yyi   =  masses(i)*(dxyz(1)*dxyz(1) + dxyz(3)*dxyz(3)) 
  zzi   =  masses(i)*(dxyz(2)*dxyz(2) + dxyz(1)*dxyz(1)) 
  xyi   = -masses(i)*(dxyz(1)*dxyz(2)) 
  xzi   = -masses(i)*(dxyz(1)*dxyz(3)) 
  yzi   = -masses(i)*(dxyz(2)*dxyz(3)) 
  blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(1,1) = blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(1,1) + xxi 
  blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(2,2) = blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(2,2) + yyi 
  blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(3,3) = blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(3,3) + zzi 
  blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(1,2) = blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(1,2) + xyi 
  blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(1,3) = blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(1,3) + xzi 
  blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(2,3) = blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(2,3) + yzi 
  blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(2,1) = blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(1,2) 
  blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(3,1) = blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(1,3) 
  blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(3,2) = blocks(asyms%aunit(1)%atom(i)%iblock)%tnsI(2,3)
 
end do

end subroutine

subroutine blocks_tnsi_iso(input,atoms_in,blocks)
use mod_crysbuild, only : atom_mass_constr
use mod_math,      only : dist_sqr,distance
use mod_types,     only : inp_par
type(inp_par),            intent(in)    :: input
type(protein_atom_list_type),     intent(in)    :: atoms_in
type(block), allocatable, intent(inout) :: blocks(:)
real(dp), allocatable               :: masses(:)
integer   :: i
real(dp)  :: dxyz(3),xxi,yyi,zzi,xyi,xzi,yzi,dist

call atom_mass_constr (input,atoms_in,masses)

do i=1,atoms_in%natoms
  dxyz = atoms_in%atom(i)%x(:)-blocks(atoms_in%atom(i)%iblock)%cent
  dist = distance(dxyz)
  if (dist .gt. 50.0d0) then
    print *, 'WARNING> atom ',i, 'is far from block center:', atoms_in%atom(i)%iblock
    print *, 'cent:',blocks(atoms_in%atom(i)%iblock)%cent
    print *, 'dist:',dist
    print *, 'WARNING> are blocks defined correctly?'
  end if
  xxi   =  masses(i)*(dxyz(2)*dxyz(2) + dxyz(3)*dxyz(3)) 
  yyi   =  masses(i)*(dxyz(1)*dxyz(1) + dxyz(3)*dxyz(3)) 
  zzi   =  masses(i)*(dxyz(2)*dxyz(2) + dxyz(1)*dxyz(1)) 
  xyi   = -masses(i)*(dxyz(1)*dxyz(2)) 
  xzi   = -masses(i)*(dxyz(1)*dxyz(3)) 
  yzi   = -masses(i)*(dxyz(2)*dxyz(3)) 
  blocks(atoms_in%atom(i)%iblock)%tnsI(1,1) = blocks(atoms_in%atom(i)%iblock)%tnsI(1,1) + xxi 
  blocks(atoms_in%atom(i)%iblock)%tnsI(2,2) = blocks(atoms_in%atom(i)%iblock)%tnsI(2,2) + yyi 
  blocks(atoms_in%atom(i)%iblock)%tnsI(3,3) = blocks(atoms_in%atom(i)%iblock)%tnsI(3,3) + zzi 
  blocks(atoms_in%atom(i)%iblock)%tnsI(1,2) = blocks(atoms_in%atom(i)%iblock)%tnsI(1,2) + xyi 
  blocks(atoms_in%atom(i)%iblock)%tnsI(1,3) = blocks(atoms_in%atom(i)%iblock)%tnsI(1,3) + xzi 
  blocks(atoms_in%atom(i)%iblock)%tnsI(2,3) = blocks(atoms_in%atom(i)%iblock)%tnsI(2,3) + yzi 
  blocks(atoms_in%atom(i)%iblock)%tnsI(2,1) = blocks(atoms_in%atom(i)%iblock)%tnsI(1,2) 
  blocks(atoms_in%atom(i)%iblock)%tnsI(3,1) = blocks(atoms_in%atom(i)%iblock)%tnsI(1,3) 
  blocks(atoms_in%atom(i)%iblock)%tnsI(3,2) = blocks(atoms_in%atom(i)%iblock)%tnsI(2,3)
 
end do

end subroutine

subroutine atomlist_COM(input,atom_list,iat,fat,com)
! center of mass calc
use mod_types, only: inp_par
use mod_crysbuild, only : atom_mass_constr 
type(inp_par), intent(in) :: input
type(protein_atom_list_type),  intent(in)   :: atom_list
integer, intent(in) :: iat,fat
real(dp),              intent(out)  :: com(3)
real(dp), allocatable               :: masses(:)
real(dp) :: sum_mass
integer   :: i

!print *, 'shit',iat,fat
call atom_mass_constr (input,atom_list,masses)
com = zero
sum_mass = zero

do i =iat,fat
  com = com + atom_list%atom(i)%x*masses(i)
  sum_mass = sum_mass + masses(i)
end do

com = com/sum_mass

end subroutine

subroutine moment_inertia (input,atoms_in,center,iat,fat,moment_I)
use mod_types, only: inp_par
use mod_math, only: dist_sqr 
use mod_crysbuild, only : atom_mass_constr
! DMR 04-27-2009 
!   computes moment of inertia as [sum_{i=iat,fat} m_i (r_i)^2]
!   where r_i is the vector to center   
type(inp_par),           intent(in)  :: input
type (protein_atom_list_type),   intent(in)  :: atoms_in
real (dp) ,              intent(in)  :: center(3)
integer,                 intent(in)  :: iat,fat
real (dp),               intent(out) :: moment_I
integer                              :: i,j, nats_block, ier
real(dp) :: dxyz(3),rsq_i
real(dp),                allocatable :: masses(:)

call atom_mass_constr(input,atoms_in,masses)

moment_I = zero

do i=iat,fat
  dxyz  = atoms_in%atom(i)%x(:)-center(:)
  rsq_i = dist_sqr(dxyz)  
  moment_I = moment_I + masses(i)*rsq_i        
end do

end subroutine

subroutine test_moment_inertia (input,atoms_in,center,moment_I)
use mod_math, only: dist_sqr 
use mod_types, only: inp_par
! DMR 04-27-2009 
!   computes moment of inertia as [sum_{i=iat,fat} m_i (r_i)^2]
!   where r_i is the vector to center
type(inp_par), intent(in) :: input   
type (protein_atom_list_type),   intent(in)  :: atoms_in
real (dp) ,              intent(in)  :: center(3)
real (dp),                intent(out) :: moment_I
integer                              :: iat,fat
logical                              :: tf_mass 
real(dp) :: com(3)

print *,'all atoms' 
tf_mass = .true.
iat = 1
fat = atoms_in%natoms
call atomlist_com(input,atoms_in,iat,fat,com)
call moment_inertia(input,atoms_in,com,iat,fat,moment_I)

print '(A3,1x,3F10.5)', 'com', com
print *, 'moment_I', dsqrt(moment_I)

print *,'atoms 1 to 30'

iat = 1
fat = 8
call atomlist_com(input,atoms_in,iat,fat,com)
call moment_inertia(input,atoms_in,com,iat,fat,moment_I)

print '(A3,1x,3F10.5)', 'com', com
print *, 'moment_I', dsqrt(moment_I)

end subroutine

subroutine isocorr_byblock(input,dr,atcorr,distances,blocks,blockid)
type(inp_par),           intent(in) :: input
real(dp),                intent(in) :: dr
real(dp),allocatable   , intent(in) :: atcorr(:,:),distances(:,:)
type(block),allocatable, intent(inout) :: blocks(:)
character(len=*) ,      intent(in) :: blockid
real(dp),allocatable                :: corr(:),isocorr(:),isocorr_intra(:),isocorr_inter(:)
integer,allocatable                 :: histo(:),histo_intra(:),histo_inter(:)
integer  :: i,j,nbin,ier

open(unit=22,file=trim(input%fileroot)//"-"//trim(input%genm)//"-"// &
                  trim(input%tpbc)//"-"//trim(input%fctyp) //"-"//   &
                  trim(blockid)//"-isocorr.txt", &
      status='unknown',action='write', iostat=ier)
open(unit=23,file=trim(input%fileroot)//"-"//trim(input%genm)//"-"// &
                trim(input%tpbc)//"-"//trim(input%fctyp) //"-"//   &
                trim(blockid)//"-isocorr-intra.txt", &
    status='unknown',action='write', iostat=ier)
open(unit=24,file=trim(input%fileroot)//"-"//trim(input%genm)//"-"// &
                trim(input%tpbc)//"-"//trim(input%fctyp) //"-"//   &
                trim(blockid)//"-isocorr-inter.txt", &
    status='unknown',action='write', iostat=ier)

do i = 1 , size(blocks)
  call byblock_isocorr_calc(dr,atcorr,distances,blocks(i))
end do

nbin = size(blocks(1)%tmpmat,1)
allocate(isocorr(nbin), histo(nbin),&
         isocorr_intra(nbin),histo_intra(nbin),&
         isocorr_inter(nbin),histo_inter(nbin), stat = ier)
isocorr = zero;  histo=0
isocorr_intra = zero;  histo_intra=0
isocorr_inter = zero;  histo_inter=0

do i = 1 , size(blocks)
  do j = 1,nbin
    histo(j)   = histo(j) + blocks(i)%itmpmat(j,1)
    isocorr(j) = isocorr(j) + blocks(i)%tmpmat(j,1)!*real(blocks(i)%itmpmat(j,1),kind=dp)
    histo_intra(j)   = histo_intra(j) + blocks(i)%itmpmat(j,2)
    isocorr_intra(j) = isocorr_intra(j) + blocks(i)%tmpmat(j,2)!*real(blocks(i)%itmpmat(j,2),kind=dp)
    histo_inter(j)   = histo_inter(j) + blocks(i)%itmpmat(j,3)
    isocorr_inter(j) = isocorr_inter(j) + blocks(i)%tmpmat(j,3)!*real(blocks(i)%itmpmat(j,3),kind=dp)
  end do
end do

! average them
do j = 1,nbin
  if(histo(j)>0)        isocorr(j)        = isocorr(j)/real(histo(j),kind=dp)
  if(histo_intra(j)>0)  isocorr_intra(j)  = isocorr_intra(j)/real(histo_intra(j),kind=dp)
  if(histo_inter(j)>0)  isocorr_inter(j)  = isocorr_inter(j)/real(histo_inter(j),kind=dp)
end do

allocate(corr(size(isocorr)), stat=ier)
do i = 1, size(isocorr)
  corr(i) = real(i-1)*dr
end do

do i = 1, size(isocorr)
  if (histo(i) > 0) then
    write(22,'(F10.4,I8,F10.4)') corr(i), histo(i),isocorr(i)
  end if
  if (histo_intra(i) > 0) then
    write(23,'(F10.4,I8,F10.4)') corr(i), histo_intra(i),isocorr_intra(i)
  end if
  if (histo_inter(i) > 0) then
    write(24,'(F10.4,I8,F10.4)') corr(i), histo_inter(i),isocorr_inter(i)
  end if
end do

end subroutine

subroutine byblock_isocorr_calc(dr,atcorr,distances,blockin)
! compute the average correlation (atomic) value as a function of interatomic separation
! basically bin up all things within shell and average their correlation... msd?
! isocorr(:) should start at one (self correlations) and head off to zero
! let's do it blockwise
! we'll use the itmpmat and tmpmat to store the intra, inter, and total block correlations
! so the dim of these matrices will be nbin x 3
real(dp),              intent(in)    :: dr   ! requested resolution... ie. 0.1 \AA\ 
real(dp), allocatable, intent(in)    :: atcorr(:,:),distances(:,:)
type(block),           intent(inout) :: blockin
integer :: maxbin,nbin,ier,i,j,ihist,iat,jat

maxbin = 0
nbin = nint(maxval(distances)/dr) +1

if (allocated(blockin%itmpmat)) deallocate(blockin%itmpmat)
if (allocated(blockin%tmpmat)) deallocate(blockin%tmpmat)
allocate (blockin%tmpmat(nbin,3),blockin%itmpmat(nbin,3), stat=ier)
blockin%tmpmat = zero ; blockin%itmpmat = 0

call block_isocorr_tot(dr,atcorr,distances,blockin)
call block_isocorr_intra(dr,atcorr,distances,blockin)
call block_isocorr_inter(blockin)

! average them up
!do i = 1, maxbin
!  do j = 1,3
!    if(blockin%itmpmat(i,j) > 0)  blockin%tmpmat(i,j) = blockin%tmpmat(i,j)/real(blockin%itmpmat(i,j),kind=dp)
!  end do
!end do

end subroutine

subroutine block_isocorr_inter(blockin)
type(block),           intent(inout) :: blockin
integer :: i,j

if (allocated(blockin%itmpmat) .and. allocated(blockin%tmpmat)) then
  ! do nothing
else
  print *, "you must have values for total and intra correlation to determine the inter block correlations"
  stop
end if

do i = 1, size(blockin%itmpmat,1) 
  blockin%itmpmat(i,3) = blockin%itmpmat(i,1) - blockin%itmpmat(i,2)
  blockin%tmpmat(i,3)  = blockin%tmpmat(i,1)  - blockin%tmpmat(i,2)
end do

end subroutine

subroutine block_isocorr_tot(dr,atcorr,distances,blockin)
! compute the average correlation (atomic) value as a function of interatomic separation
!   bin up within shell and average 
! done blockwise: 
!   split into inter intra and total; using blockin%tmpmat and blockin%itmpmat
real(dp),              intent(in)    :: dr   ! requested resolution... ie. 0.1 \AA\ 
real(dp), allocatable, intent(in)    :: atcorr(:,:),distances(:,:)
type(block),           intent(inout) :: blockin
integer :: maxbin,nbin,ier,i,j,ihist,iat,jat,nwhich

nwhich = 1 !tot is 1, intra is 2, inter is 3
nbin = nint(maxval(distances)/dr) +1

if (allocated(blockin%itmpmat) .and. allocated(blockin%tmpmat)) then
  blockin%tmpmat(:,nwhich) = zero ; blockin%itmpmat(:,nwhich) = 0
else
  allocate (blockin%tmpmat(nbin,3),blockin%itmpmat(nbin,3), stat=ier)
  blockin%tmpmat = zero ; blockin%itmpmat = 0
end if

! self interactions
do i = 1, size(blockin%atinblock,1)
  iat = blockin%atinblock(i)
  ihist = nint(distances(iat,iat)/dr) + 1 ! unnecessary, but keeping it as check
  blockin%itmpmat(ihist,nwhich) = blockin%itmpmat(ihist,nwhich) + 1
  blockin%tmpmat(ihist,nwhich)  = blockin%tmpmat(ihist,nwhich) + atcorr(iat,iat)
end do

! total intra + inter blockin
do i = 1, size(blockin%atinblock,1)
  iat = blockin%atinblock(i)
  do j = iat+1, size(atcorr,1)
    jat = j
    ihist = nint(distances(iat,jat)/dr) + 1
    blockin%itmpmat(ihist,nwhich) = blockin%itmpmat(ihist,nwhich) + 1
    blockin%tmpmat(ihist,nwhich)  = blockin%tmpmat(ihist,nwhich)  + atcorr(iat,jat)
  end do
end do

end subroutine

subroutine block_isocorr_intra(dr,atcorr,distances,blockin)
! compute the average correlation (atomic) value as a function of interatomic separation
!   bin up within shell and average 
! done blockwise: 
!   split into inter intra and total; using blockin%tmpmat and blockin%itmpmat
real(dp),              intent(in)    :: dr   ! requested resolution... ie. 0.1 \AA\ 
real(dp), allocatable, intent(in)    :: atcorr(:,:),distances(:,:)
type(block),           intent(inout) :: blockin
integer :: maxbin,nbin,ier,i,j,ihist,iat,jat,nwhich

nwhich = 2 !intra is 2
nbin = nint(maxval(distances)/dr) +1

if (allocated(blockin%itmpmat) .and. allocated(blockin%tmpmat)) then
  blockin%tmpmat(:,nwhich) = zero ; blockin%itmpmat(:,nwhich) = 0
else
  allocate (blockin%tmpmat(nbin,3),blockin%itmpmat(nbin,3), stat=ier)
  blockin%tmpmat = zero ; blockin%itmpmat = 0
end if

! self interactions, we could just write down the result, but we'll keep it as check
do i = 1, size(blockin%atinblock,1)
  iat = blockin%atinblock(i)
  ihist = nint(distances(iat,iat)/dr) + 1 
  blockin%itmpmat(ihist,nwhich) = blockin%itmpmat(ihist,nwhich) + 1
  blockin%tmpmat(ihist,nwhich)  = blockin%tmpmat(ihist,nwhich) + atcorr(iat,iat)
end do

do i = 1, size(blockin%atinblock,1)
  iat = blockin%atinblock(i)
  do j = i+1, size(blockin%atinblock,1)
    jat = blockin%atinblock(j)
    ihist = nint(distances(iat,jat)/dr) + 1
    blockin%itmpmat(ihist,nwhich) = blockin%itmpmat(ihist,nwhich) + 1
    blockin%tmpmat(ihist,nwhich)  = blockin%tmpmat(ihist,nwhich) + atcorr(iat,jat)
  end do
end do

end subroutine

end module

