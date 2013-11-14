module mod_symop
use cfml_globaldeps,                 only: dp,sp
use mod_constants
use cfml_crystal_metrics,              only: Crystal_Cell_Type
use cfml_crystallographic_symmetry,  only: space_group_type

implicit none

contains

subroutine mat3_rotate(cell,spg,iop,mat)
! rotate matrix:
!     xx xy xz
!     yx yy yz
!     zx zy zz
! in cartesian vec
use mod_linalg,                 only: lapack_eig

type (Crystal_Cell_Type), intent(in)          :: Cell
type (space_group_type),  intent(in)          :: SpG
integer,                  intent(in)          :: iop ! label of the symoperator
real(dp),                 intent(in out)      :: mat(3,3)

mat = matmul(Cell%orth_Cr_cel,matmul(mat,transpose(Cell%orth_Cr_cel)))
mat = matmul(spg%symop(iop)%rot,matmul(mat,transpose(spg%symop(iop)%rot)))
mat = matmul(Cell%Cr_orth_cel,matmul(mat,transpose(Cell%cr_orth_cel))) 

end subroutine

subroutine vecmat3_rotate(cell,spg,iop,matvec)
! rotate matrix:
!     xx xy xz
!        yy yz
!           zz
! in cartesian vec  stored in 6D vector
use mod_linalg,                 only: lapack_eig

type (Crystal_Cell_Type), intent(in)          :: Cell
type (space_group_type),  intent(in)          :: SpG
integer,                  intent(in)          :: iop ! label of the symoperator
real(sp),                 intent(in out)      :: matvec(6)
real(sp) :: mat(3,3)

mat(1,1) = matvec(1)
mat(2,2) = matvec(2)
mat(3,3) = matvec(3)
mat(1,2) = matvec(4)
mat(1,3) = matvec(5)
mat(2,3) = matvec(6)
mat(2,1) = matvec(4)
mat(3,1) = matvec(5)
mat(3,2) = matvec(6)

mat = matmul(Cell%orth_Cr_cel,matmul(mat,transpose(Cell%orth_Cr_cel)))
mat = matmul(spg%symop(iop)%rot,matmul(mat,transpose(spg%symop(iop)%rot)))
mat = matmul(Cell%Cr_orth_cel,matmul(mat,transpose(Cell%cr_orth_cel))) 

matvec(1) = mat(1,1)
matvec(2) = mat(2,2)
matvec(3) = mat(3,3)
matvec(4) = mat(1,2)
matvec(5) = mat(1,3)
matvec(6) = mat(2,3)

end subroutine

subroutine vcov_p1(cell,spg,smatvcov,bigvcov)
type (Crystal_Cell_Type), intent(in)        :: Cell
type (space_group_type),  intent(in)        :: SpG
real(dp), allocatable,    intent(in)        :: smatvcov(:,:)
real(dp), allocatable,    intent(out)       :: bigvcov(:,:)
real(dp), allocatable                       :: trnsvcov(:,:)
integer :: bigdim,ndim,ier,i,ii,jj

ndim = size(smatvcov(1,:))
bigdim = ndim*spg%multip
allocate(bigvcov(bigdim,bigdim), stat=ier)
if (ier /= 0) stop "tls_atcovar_p1> malloc"
bigvcov = 0.0d0

do i = 1, spg%multip
  ii = (i-1)*ndim+1
  call vcov_iop(i,cell,spg,smatvcov,trnsvcov)
  bigvcov(ii:ii+ndim-1,ii:ii+ndim-1) = trnsvcov
end do

end subroutine

subroutine tempfactor_p1(cell,spg,smatdfmat,bigdfmat)
type (Crystal_Cell_Type), intent(in)        :: Cell
type (space_group_type),  intent(in)        :: SpG
real(dp), allocatable,    intent(in)        :: smatdfmat(:)
real(dp), allocatable,    intent(out)       :: bigdfmat(:)
integer :: bigdim,ndim,ier,i,ii,jj

ndim = size(smatdfmat,1)
bigdim = ndim*spg%multip

if (allocated(bigdfmat)) deallocate(bigdfmat)
allocate(bigdfmat(bigdim), stat=ier)
if (ier /= 0) stop "anisotempfactor_p1> malloc"
bigdfmat = 0.0d0

if (bigdim .eq. ndim) then
  bigdfmat = smatdfmat
  return
end if

do i = 1, spg%multip
  ii = (i-1)*ndim+1
  bigdfmat(ii:ii+ndim-1) = smatdfmat
end do

end subroutine

subroutine anisotempfactor_p1(cell,spg,smatdfmat,bigdfmat)
type (Crystal_Cell_Type), intent(in)        :: Cell
type (space_group_type),  intent(in)        :: SpG
real(dp), allocatable,    intent(in)        :: smatdfmat(:,:)
real(dp), allocatable,    intent(out)       :: bigdfmat(:,:)
real(dp), allocatable                       :: trnsdfmat(:,:)
integer :: bigdim,ndim,ier,i,ii,jj

ndim = size(smatdfmat,1)
bigdim = ndim*spg%multip

if (allocated(bigdfmat)) deallocate(bigdfmat)
allocate(bigdfmat(bigdim,3), stat=ier)
if (ier /= 0) stop "anisotempfactor_p1> malloc"
bigdfmat = 0.0d0

if (bigdim .eq. ndim) then
  bigdfmat = smatdfmat
  return
end if

do i = 1, spg%multip
  ii = (i-1)*ndim+1
  call anisotempfactor_iop(i,cell,spg,smatdfmat,trnsdfmat)
  bigdfmat(ii:ii+ndim-1,1:3) = trnsdfmat
end do

end subroutine

subroutine atcor_p1(spg,smatcor,bigatcor,tdim)
! stamps the smatcor along the diagonal
! this is atom-atom correlation so it's isotropic and needs no sym ops
type (space_group_type),  intent(in)        :: SpG
real(dp), allocatable,    intent(in)        :: smatcor(:,:)
real(dp), allocatable,    intent(out)       :: bigatcor(:,:)
integer, optional, intent(in) :: tdim
real(dp), allocatable                       :: trnsvcov(:,:)
integer :: bigdim,ndim,ier,i,ii,jj


ndim = size(smatcor(1,:))

if (present(tdim)) then
  if (tdim .eq. ndim) then
    bigdim = ndim
    allocate(bigatcor(bigdim,bigdim), stat=ier)
    if (ier /= 0) stop "tls_atcor_p1> malloc"
    bigatcor = smatcor
    return
  end if
else
  bigdim = ndim*spg%multip
  allocate(bigatcor(bigdim,bigdim), stat=ier)
  if (ier /= 0) stop "tls_atcor_p1> malloc"
  bigatcor = 0.0d0

  do i = 1, spg%multip
    ii = (i-1)*ndim+1
    print *, 'dimrc',ii,ii+ndim-1
    bigatcor(ii:ii+ndim-1,ii:ii+ndim-1) = smatcor 
  end do
end if

end subroutine

subroutine vcov_iop(iop,cell,spg,smatvcov,trnsvcov)
integer,                  intent(in)        :: iop
type (Crystal_Cell_Type), intent(in)        :: Cell
type (space_group_type),  intent(in)        :: SpG
real(dp), allocatable,    intent(in)        :: smatvcov(:,:)
real(dp), allocatable,    intent(out)       :: trnsvcov(:,:)
real(dp) :: tmat(3,3)
integer :: i,j,ii,jj,ndim,ier,iun

ndim = size(smatvcov(1,:))

allocate(trnsvcov(ndim,ndim),stat =ier)
if (ier /= 0) stop "tls_atcovar_iop> malloc"

trnsvcov = smatvcov

do i = 1,ndim/3
  ii = 3*(i-1)+1
  do j = 1,ndim/3
    jj = 3*(j-1)+1
    tmat = trnsvcov(ii:ii+2,jj:jj+2)
    call mat3_rotate(cell,spg,iop,tmat)
    trnsvcov(ii:ii+2,jj:jj+2) = tmat
  end do
end do

end subroutine

subroutine anisotempfactor_iop(iop,cell,spg,smatvcov,trnsvcov)
integer,                  intent(in)        :: iop
type (Crystal_Cell_Type), intent(in)        :: Cell
type (space_group_type),  intent(in)        :: SpG
real(dp), allocatable,    intent(in)        :: smatvcov(:,:)
real(dp), allocatable,    intent(out)       :: trnsvcov(:,:)
real(dp) :: tmat(3,3)
integer :: i,j,ii,jj,ndim,ier,iun

ndim = size(smatvcov,1)

allocate(trnsvcov(ndim,3),stat =ier)
if (ier /= 0) stop "tls_atcovar_iop> malloc"

trnsvcov = smatvcov

do i = 1,ndim/3
  ii = 3*(i-1)+1
  tmat = smatvcov(ii:ii+2,1:3)
  call mat3_rotate(cell,spg,iop,tmat)
  trnsvcov(ii:ii+2,1:3) = tmat
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


subroutine vec_rotate(cell,spg,iop,vec)
! rotate single vector (x1,y1,z1,.... xn,yn,zn)
! take in cartesian vec
type (Crystal_Cell_Type), intent(in)        :: Cell
type (space_group_type), intent(in)         :: SpG
integer,                 intent(in)         :: iop ! label of the symoperator
real(dp), allocatable,intent(in out)  :: vec(:)
integer :: i,j,k,natoms

  natoms = size(vec)/3      
  call vec_c2f(cell,vec)    
  k=1
  do j = 1,natoms           
    vec(k:k+2) = matmul(spg%symop(iop)%rot,vec(k:k+2))
    k = k + 3            
  end do
  call vec_f2c(cell,vec)

end subroutine

subroutine vecs_rotate(cell,spg,iop,vecs)
! rotate vector of vectors (x1,y1,z1,.... xn,yn,zn)
type (Crystal_Cell_Type), intent(in)        :: Cell
type (space_group_type), intent(in)         :: SpG
integer,                 intent(in)         :: iop ! label of the symoperator
real(dp), allocatable,intent(in out)  :: vecs(:,:)
real(dp), allocatable :: tmpvec(:)
integer :: i,j,k,nrows,ncols,ialloc

nrows=size(vecs(:,1))
ncols=size(vecs(1,:))
allocate(tmpvec(nrows),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

do j = 1,ncols
  tmpvec(:) = vecs(:,j)
  call vec_rotate(cell,spg,iop,tmpvec)
  vecs(:,j) = tmpvec(:)
end do

end subroutine

subroutine vector_combine(vec1,vec2,vec12)
! combines to vectors into a new one with combined dimension
!   pretty sure we'll always be using this for real vectors
real(dp), allocatable, dimension(:), intent(in) :: vec1,vec2
real(dp), allocatable, dimension(:), intent(out) :: vec12
integer :: ialloc, ndim1,ndim2,ndim12

ndim1  = size(vec1)
ndim2  = size(vec2)
ndim12 = ndim1 + ndim2

allocate(vec12(ndim12), stat = ialloc)

vec12(1:ndim1)        = vec1(:)
vec12(ndim1+1:ndim12) = vec2(:)

end subroutine

subroutine vectors_combine(vecs1,vecs2,vecs12)
! takes to sets of vectors and combines them with same indexing
!   pretty sure we'll always be using this for real vectors
real(dp), allocatable, dimension(:,:), intent(in) :: vecs1,vecs2
real(dp), allocatable, dimension(:,:), intent(out) :: vecs12
real(dp), allocatable, dimension(:) :: tmpv1,tmpv2,tmpv12
integer :: ialloc, ndim1,ndim2,nvecs,ndim12
integer :: i

ndim1  = size(vecs1(:,1))
ndim2  = size(vecs2(:,1))
nvecs  = size(vecs1(1,:))
if (nvecs .ne. size(vecs2(1,:))) stop "vectors_combine> error in number of vectors"
ndim12 = ndim1 + ndim2

allocate(vecs12(ndim12,nvecs),tmpv1(ndim1),tmpv2(ndim2), stat = ialloc)

do i=1,nvecs
  tmpv1 = vecs1(:,i) 
  tmpv2 = vecs2(:,i) 
  call vector_combine(tmpv1,tmpv2,tmpv12)       
  vecs12(:,i) = tmpv12(:)
end do

deallocate(tmpv1,tmpv2,tmpv12)

end subroutine

subroutine vectors_grow(addvecs,allvecs)
! adds each of vecs1 to the end of each of vecs12
real(dp), allocatable, dimension(:,:), intent(in)    :: addvecs
real(dp), allocatable, dimension(:,:), intent(inout) :: allvecs
real(dp), allocatable, dimension(:,:)                :: tmpvecs
real(dp), allocatable, dimension(:) :: tmpv1,tmpv2,tmpv12
integer :: ialloc, ndim1,ndim2,nvecs,ndim12
integer :: i

ndim1  = size(allvecs(:,1))
ndim2  = size(addvecs(:,1))
nvecs  = size(addvecs(1,:))
if (nvecs .ne. size(allvecs(1,:))) stop "vectors_combine> error in number of vectors"

ndim12 = ndim1 + ndim2

allocate(tmpvecs(ndim12,nvecs),tmpv1(ndim1),tmpv2(ndim2), stat = ialloc)

! add in the new
do i=1,nvecs
  tmpv1 = allvecs(:,i) 
  tmpv2 = addvecs(:,i) 
  call vector_combine(tmpv1,tmpv2,tmpv12)       
  tmpvecs(:,i) = tmpv12(:)
end do

deallocate(allvecs)
allocate(allvecs(size(tmpvecs(:,1)),size(tmpvecs(1,:))), stat = ialloc)
if (ialloc /= 0) stop "vectors_grow> memory problem"

allvecs = tmpvecs

deallocate(tmpv1,tmpv2,tmpv12,tmpvecs)

end subroutine

subroutine symop_vecbuild(ndim,cell,spg,asym_vecs,p1_vecs)
use cfml_crystallographic_symmetry,only: space_group_type
integer, intent(in)  :: ndim
type (Crystal_Cell_Type),              intent(in)  :: Cell
type (space_group_type),               intent(in)  :: SpG
real(dp), allocatable, dimension(:,:), intent(in)  :: asym_vecs
real(dp), allocatable, dimension(:,:), intent(out) :: p1_vecs
real(dp), allocatable, dimension(:,:)              :: tmpvecs
integer :: i,j, ialloc, iasy,nvec

iasy = size(asym_vecs(:,1))
nvec = size(asym_vecs(1,:))
allocate(tmpvecs(iasy,nvec),p1_vecs(iasy,nvec) , stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"
! copy in the first, no rotation
p1_vecs = asym_vecs

! prep the tmpvecs
tmpvecs = asym_vecs

do j=2,SpG%multip
  if (ndim .eq. 3) call vecs_rotate(cell,spg,j,tmpvecs)
  call vectors_grow(tmpvecs,p1_vecs)
  ! reset the vectors
  tmpvecs = asym_vecs
end do

end subroutine

end module
