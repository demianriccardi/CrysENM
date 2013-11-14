module mod_valvec_store
use cfml_globaldeps,                 only: dp,sp
use mod_constants

implicit none

type :: valvec
integer :: nfrq
real(dp) :: multiplicity  ! multiplicity for a give wave vector... zero,zero,zero = 1, and others =2..
real(dp)       ,allocatable  :: vals(:),rvecs(:,:),ivecs(:,:),qvec(:)
complex(kind=8),allocatable  :: zvecs(:,:)
end type valvec

contains

subroutine zero_store(vals,vecs,branches,qvec)
real(dp), allocatable,        intent(in)    :: vals(:),vecs(:,:)
type (valvec),                intent(out)   :: branches
real(dp), optional,           intent(in)    :: qvec(3)
real(dp) :: testval
integer  :: ialloc,nfrq,ndim
integer  :: i,j

allocate(branches%qvec(3), stat = ialloc)
! test the wavevector
if (present(qvec)) then
  testval = dot_product(qvec,qvec)
  if (testval .ne. zero) then
    print *,"rvec_store> error, calling real storage with nonzero wavevector"
    return
  end if
end if

branches%qvec = zero
branches%multiplicity=one
ndim = size(vecs(:,1))
nfrq = size(vals)
branches%nfrq = nfrq

allocate(branches%vals(nfrq), branches%zvecs(ndim,nfrq), stat = ialloc)
if (ialloc /= 0) stop "rvec_store> memory error"

do i = 1, nfrq
  branches%vals(i) = vals(i)
  do j = 1, ndim
    branches%zvecs(j,i) = cmplx(vecs(j,i),zero,kind=8)
  end do
end do 

end subroutine

subroutine qvec_store(qvec,vals,vecs,branches)
real(dp),                     intent(in)    :: qvec(3)
real(dp), allocatable,        intent(in)    :: vals(:)
complex(kind=8), allocatable, intent(in)    :: vecs(:,:)
type (valvec),                intent(out)   :: branches
real(dp) :: testval
integer  :: ialloc,nfrq,ndim
integer  :: i,j

allocate(branches%qvec(3), stat = ialloc)
branches%qvec = qvec

ndim = size(vecs(:,1))
nfrq = size(vals)
branches%nfrq = nfrq
testval = dot_product(qvec,qvec)
if (testval .ne. zero)  then
  branches%multiplicity=two
else
  branches%multiplicity=one
end if  

allocate(branches%vals(nfrq), branches%zvecs(ndim,nfrq), stat = ialloc)
if (ialloc /= 0) stop "rvec_store> memory error"

do i = 1, nfrq
  branches%vals(i) = vals(i)
  do j = 1, ndim
    branches%zvecs(j,i) = vecs(j,i)
  end do
end do 

end subroutine

end module
