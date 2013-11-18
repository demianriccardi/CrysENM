module mod_linalg
! calls to pseudo and eispack
use cfml_globaldeps,                 only: dp,sp
use mod_math
use mod_constants

implicit none

contains

subroutine lapack_zeig (matin, &
                             zvals,zvecs )
! slower pseudo inverse where identity matrix is converted to pseudo inverse
integer                           :: ndim
complex(kind=8),intent(in), allocatable  :: matin(:,:)
complex(kind=8),intent(out),allocatable  :: zvals(:),zvecs(:,:)
character*1                       :: jobz='V',uplo='U'
integer                           :: info,lda,liwork,lrwork,lwork,n
integer, allocatable              :: iwork(:)
real(dp),     allocatable         :: rwork(:), w(:) 
complex(kind=8),     allocatable  :: A(:,:),   work(:) 
integer                           :: ierr
!lapack

ndim=size(matin(1,:))

if (allocated(zvecs)) deallocate(zvecs)
if (allocated(zvals)) deallocate(zvals)
allocate(zvecs(ndim,ndim),zvals(ndim),stat=ierr)
if (ierr /= 0) STOP "*** not enough memory ***"

n=ndim
lda=n

lwork  = 2*n+n**2
lrwork = 1+5*n+2*n**2
liwork = 3+5*n

allocate(a(ndim,ndim),w(ndim),work(lwork),&
         rwork(lrwork),iwork(liwork),stat=ierr)
if (ierr /= 0) STOP "*** not enough memory ***"

a=matin

!call r4_svd(ndim,ndim,matin,w,matu,u,matv,v,ierr)
call zheevd(jobz,uplo,n,a,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)

zvals=w
zvecs=a

deallocate(a,w,rwork,iwork,work)

!print *, "Time calculating pseudo inverse", time2-time1

end subroutine

subroutine my_dgemm (a, b, c, TransA,TransB,alpha,beta)
real(kind=8),intent(in),  allocatable  :: A(:,:), B(:,:)
real(kind=8),intent(inout), allocatable  :: C(:,:)
character*1 ,intent(in)                :: TransA, TransB
real(kind=8),intent(in)                :: alpha,beta
!
integer                                :: M, N, K, LDA, LDB, LDC
integer                                :: ierr

if (TransA .eq. 'N') then
  M=size(a(:,1))
  K=size(a(1,:))
  LDA = M
else
  M=size(a(1,:))
  K=size(a(:,1))
  LDA = K
endif

if (TransB .eq. 'N') then
  N=size(b(1,:))
  if (K .ne. size(b(:,1))) stop "my_dgemm> dimension mismatch"
  LDB = K
else
  N=size(b(:,1))
  if (K .ne. size(b(1,:))) stop "my_dgemm> dimension mismatch"
  LDB = N
endif

LDC = M

if(allocated(C)) then
  if ( (size(C(:,1)) .ne. M) .or. (size(C(1,:)) .ne. N) ) then 
    stop "my_dgemm> dimension mismatch"
  endif
else
  allocate(C(M,N),stat=ierr)
  C = 0.0d0
endif

call dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )

end subroutine

subroutine my_zgemm (a, b, c, TransA,TransB,alpha,beta)
complex(kind=8),intent(in),    allocatable  :: A(:,:), B(:,:)
complex(kind=8),intent(inout), allocatable  :: C(:,:)
character*1 ,intent(in)                :: TransA, TransB
complex(kind=8),intent(in)                :: alpha,beta
!
integer                                :: M, N, K, LDA, LDB, LDC
integer                                :: ierr

if (TransA .eq. 'N') then
  M=size(a(:,1))
  K=size(a(1,:))
  LDA = M
else
  M=size(a(1,:))
  K=size(a(:,1))
  LDA = K
endif

if (TransB .eq. 'N') then
  N=size(b(1,:))
  if (K .ne. size(b(:,1))) stop "my_dgemm> dimension mismatch"
  LDB = K
else
  N=size(b(:,1))
  if (K .ne. size(b(1,:))) stop "my_dgemm> dimension mismatch"
  LDB = N
endif

LDC = M

if(allocated(C)) then
  if ( (size(C(:,1)) .ne. M) .or. (size(C(1,:)) .ne. N) ) then 
    stop "my_dgemm> dimension mismatch"
  endif
else
  allocate(C(M,N),stat=ierr)
  C = cmplx(0.0d0,0.0d0,kind = wp)
endif

call zgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )

end subroutine

subroutine zheevr_zeig (nfirst,nlast,matin, &
                             vals)
use f95_lapack, only : la_heevr
! slower pseudo inverse where identity matrix is converted to pseudo inverse
integer        ,intent(in)               :: nfirst,nlast
complex(kind=8),intent(in), allocatable  :: matin(:,:)
real(dp),       intent(out),allocatable  :: vals(:)
complex(dp),                    allocatable  :: A(:,:)
real(dp),                       allocatable  :: w(:)
integer                           :: ierr,ndim,nv,il,iu
real(sp) :: time1,time2

il = nfirst
iu = nlast
nv = iu-il+1
call cpu_time(time1)
ndim=size(matin(1,:))
print *, 'using zheevr to compute eigvals from:',il,' to ',iu
print *, 'for matrix dimension:',ndim

if (allocated(vals)) deallocate(vals)
allocate(w(ndim),vals(nv),A(ndim,ndim),stat=ierr)
if (ierr /= 0) STOP "*** not enough memory ***"

A=matin

call la_heevr(a,w,il=il,iu=iu)

vals=w(1:nv)
call cpu_time(time2)
print *, 'time zheevr to compute eigvals from:',il,' to ',iu, ': ',time2-time1


deallocate(a,w)

end subroutine

subroutine lapack_eig (matin, &
                             vals,vecs )
real(dp),intent(in), allocatable           :: matin(:,:)
real(dp),intent(out),allocatable           :: vals(:),vecs(:,:)
character*1                       :: jobz='V',uplo='U' ! jobz='V' for eigvecs
!character*1                       :: jobz='N',uplo='U' ! jobz='N' for no eigvecs
integer                           :: info,lda,liwork,lwork,n
integer, allocatable              :: iwork(:)
real(dp),     allocatable         :: w(:) 
real(dp),     allocatable         :: A(:,:),   work(:) 
integer                           :: ndim,ierr
!lapack

ndim=size(matin(1,:))

if (allocated(vecs)) deallocate(vecs)
if (allocated(vals)) deallocate(vals)
allocate(vecs(ndim,ndim),vals(ndim),stat=ierr)
if (ierr /= 0) STOP "*** not enough memory ***"

n=ndim
lda=n

if (jobz .eq. 'V') then
  lwork  = 1 + 6*n + 2*n**2
  liwork = 3+5*n
else 
  lwork  = 3*n
  liwork = n
end if

allocate(a(ndim,ndim),w(ndim),work(lwork),&
         iwork(liwork),stat=ierr)
if (ierr /= 0) STOP "*** not enough memory ***"

a=matin

call dsyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info)
! print *, "dsyevd info=", info

vals=w
if (jobz .eq. 'V') then
  vecs=a
else 
  vecs=0.0d0
end if

deallocate(a,w,iwork,work)

!print *, "Time calculating pseudo inverse", time2-time1

end subroutine

subroutine lapack_eig3 (matin, &
                             vals,vecs )
real(dp),intent(in)           :: matin(3,3)
real(dp),intent(out)          :: vals(3),vecs(3,3)
character*1                       :: jobz='V',uplo='U' ! jobz='V' for eigvecs
!character*1                       :: jobz='N',uplo='U' ! jobz='N' for no eigvecs
integer                           :: info,lda,liwork,lwork,n
integer, allocatable              :: iwork(:)
real(dp),     allocatable         :: w(:) 
real(dp),     allocatable         :: A(:,:),   work(:) 
integer                           :: ndim,ierr
!lapack

ndim=3

n=ndim
lda=n

if (jobz .eq. 'V') then
  lwork  = 1 + 6*n + 2*n**2
  liwork = 3+5*n
else 
  lwork  = 3*n
  liwork = n
end if

allocate(a(ndim,ndim),w(ndim),work(lwork),&
         iwork(liwork),stat=ierr)
if (ierr /= 0) STOP "*** not enough memory ***"

a=matin

call dsyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info)
! print *, "dsyevd info=", info

vals=w
if (jobz .eq. 'V') then
  vecs=a
else 
  vecs=0.0d0
end if

deallocate(a,w,iwork,work)

!print *, "Time calculating pseudo inverse", time2-time1

end subroutine

subroutine lapack_pseudo (matin, &
                             matout )
real(dp),intent(in), allocatable  :: matin(:,:)
real(dp),intent(out),allocatable  :: matout(:,:)
real(dp),            allocatable  :: tmpmat(:,:),s(:),work(:) 
logical                           :: matu=.true.,matv=.true.
real(dp)                          :: rr=zero
real                              :: time1,time2
integer                           :: ierr,rank,lwork,info,i,ndim
integer, allocatable              :: iwork(:)
character*1                       :: jobz='A'

ndim = size(matin(1,:))

if (allocated(matout)) deallocate(matout)
allocate(matout(ndim,ndim),stat=ierr)
if (ierr /= 0) STOP "*** not enough memory ***"

lwork = ndim*3+2*ndim
allocate(tmpmat(ndim,ndim),s(ndim),work(lwork),stat=ierr)
if (ierr /= 0) STOP "*** not enough memory ***"

tmpmat=matin

! set matout to identity
matout=zero
do i=1,ndim
  matout(i,i)=one
end do

call cpu_time(time1)

call dgelss(ndim, ndim, ndim,tmpmat, ndim, matout,ndim, s,zero,rank, work,lwork,info)

deallocate(tmpmat,s,work)

call cpu_time(time2)

end subroutine

subroutine mod_linalg_pseudo (ndim, matin, &
                             matout,nzeros )

integer, intent(in)               :: ndim
integer, intent(out)              :: nzeros
real(dp),intent(in out), allocatable  :: matin(:,:)
real(dp),intent(out),allocatable  :: matout(:,:)
real(dp),            allocatable  :: s(:),u(:,:),vt(:,:),work(:) !svd stuff
real(dp),            allocatable  :: tmpmat(:,:)
logical                           :: matu=.true.,matv=.true.
real(dp)                          :: rr=zero
real                              :: time1,time2
integer                           :: i,j,ii,jj,ierr,lwork,iworkdim,info
integer, allocatable              :: iwork(:)
!lapack
character*1                       :: jobz='A'

! fill in other elements
do i=1,ndim
  do j=1,ndim
    matin(j,i)= matin(i,j)
  end do
end do


if (allocated(matout)) deallocate(matout)
allocate(matout(ndim,ndim),tmpmat(ndim,ndim),stat=ierr)
if (ierr /= 0) STOP "*** not enough memory ***"
tmpmat=matin

lwork = 8*ndim*ndim+4*ndim
iworkdim=ndim*8
allocate(s(ndim),u(ndim,ndim),vt(ndim,ndim),work(lwork),iwork(iworkdim),stat=ierr)
if (ierr /= 0) STOP "*** not enough memory ***"

call cpu_time(time1)

! DMR: use http://www.netlib.org/lapack/double/dgesdd.f lapack
!      dgesdd return the transpose of v ->vt
!      which is confusingly referred to as v**t in the source code  
call DGESDD( JOBZ, ndim, ndim, tmpmat, ndim, s, u, ndim, vt, ndim, work, &
                        lwork, iwork, INFO )

! computer  pinv(A)=V S^{-1}(nonzero) U^{trans}
! or here:  pinv(A)=U S^{-1} VT
matout = zero
nzeros = 0
do i = 1, ndim
  if (abs(s(i)) .gt. 1.0d-06) then
    matout(i,i) = one/s(i)
  else
    matout(i,i) = zero
    nzeros = nzeros+1
  end if
end do

tmpmat = matmul(u,matout)
matout = matmul(tmpmat,vt)

!tmpmat = matmul(transpose(vt),matout)
!matout = matmul(tmpmat,transpose(u))

deallocate(s,u,vt,work,iwork)

!call mkl_set_num_threads(1)
call cpu_time(time2)

!print *, "Time calculating pseudo inverse with DGESDD", time2-time1

end subroutine

subroutine root_invrse_mat(matin,eigvec,matout)
! uses matrix eigenvectors to compute root inverse
real(dp), dimension(3,3), intent(in)  :: matin,eigvec
real(dp), dimension(3,3), intent(out) :: matout
integer :: i,j
real(dp) :: diag_test

diag_test = 0.0d0
! get diagonal
matout = matmul(matin,eigvec)
matout = matmul(transpose(eigvec),matout)

do i = 1, 3
  do j = i+1,3
   diag_test = diag_test + abs(matout(i,j))
   diag_test = diag_test + abs(matout(j,i))
  end do
end do 

if (diag_test .gt. 1.0d-8) stop "root_invrse_mat> problem with diag mat"

do i = 1,3
  matout(i,i)=one/dsqrt(matout(i,i))
end do
! now, go back
matout = matmul(matout,transpose(eigvec))
matout = matmul(eigvec,matout)

end subroutine

subroutine invrse_mat(matin,matout)
! comput inverse of 3by3 matrix
real(dp), dimension(3,3), intent(in)  :: matin
real(dp), dimension(3,3), intent(out) :: matout
real(dp)             :: evec(3,3),eval(3)
integer :: i,j

call lapack_eig3(matin,eval,evec)
matout = zero
do i=1,3
  matout(i,i)=one/eval(i)
end do
! now, go back
matout = matmul(matout,transpose(evec))
matout = matmul(evec,matout)

end subroutine

subroutine schmidt(c)
! c is m by n matrix.  
! implemented by DMR: 4-24-2009
!     orthogonalization for gram-schmidt
!
! pseudo code for more numerically stable process from wikipedia:
!
! define projection of v on to u:
!   (u dot v) u_vec / (u dot u)          
!
! for j from 1 to m :
!   for i from 1 to j - 1 :   
!     v_j = v_j - projection of v_j onto v_i  ! remove component in direction of v_i
!
!   v_j = v_j/(v_j dot v_j) ! normalize
!
real(dp), allocatable, intent(in out) :: c(:,:)
real(dp) :: anorm,scal
integer  :: i,j,m,n

m = size(c(:,1))
n = size(c(1,:))

! normalize the first vector 
anorm = dot_product(c(1,:),c(1,:))
anorm = 1.d0 / ( dsqrt(anorm) )
c(1,:) = anorm * c(1,:)

do  i=2,m

  do j=1,i-1

    ! compute scalar part of projection
    scal = dot_product(c(i,:),c(j,:))
    ! remove component of c(i,:) that's in direction of c(j,:)
    c(i,:) = c(i,:) - scal*c(j,:)
    ! normalize it
    anorm = dot_product(c(i,:),c(i,:))
    anorm = 1.0d0 / ( dsqrt(anorm) )
    c(i,:) = c(i,:)*anorm  

  enddo

enddo

end subroutine

subroutine test_schmidt()

real(dp),allocatable,dimension(:,:)  :: vf90
integer :: err

allocate(vf90(3,3), stat = err )

vf90(1,:) = (/3.0d0, 0.0d0, 4.0d0/)
vf90(2,:) = (/-6.0d0, -4.0d0, 1.0d0/)
vf90(3,:) = (/5.0d0, 0.0d0, -3.0d0/)

print *, 'input matrix and respective overlaps'
print '(a3,1x,3F12.8)',"v1", vf90(1,:)
print '(a3,1x,3F12.8)',"v2", vf90(2,:)
print '(a3,1x,3F12.8)',"v3", vf90(3,:)
print '(a3,1x,F12.8)','v11', dot_product(vf90(1,:),vf90(1,:))
print '(a3,1x,F12.8)','v12', dot_product(vf90(1,:),vf90(2,:))
print '(a3,1x,F12.8)','v13', dot_product(vf90(1,:),vf90(3,:))
print '(a3,1x,F12.8)','v22', dot_product(vf90(2,:),vf90(2,:))
print '(a3,1x,F12.8)','v23', dot_product(vf90(2,:),vf90(3,:))
print '(a3,1x,F12.8)','v33', dot_product(vf90(3,:),vf90(3,:))

print *, 'using gram-schmidt to orthogonalize'
call schmidt(vf90)

print *, 'all set:'
print '(a3,1x,3F12.8)',"v1", vf90(1,:)
print '(a3,1x,3F12.8)',"v2", vf90(2,:)
print '(a3,1x,3F12.8)',"v3", vf90(3,:)
print '(a3,1x,F12.8)','v11', dot_product(vf90(1,:),vf90(1,:))
print '(a3,1x,F12.8)','v12', dot_product(vf90(1,:),vf90(2,:))
print '(a3,1x,F12.8)','v13', dot_product(vf90(1,:),vf90(3,:))
print '(a3,1x,F12.8)','v22', dot_product(vf90(2,:),vf90(2,:))
print '(a3,1x,F12.8)','v23', dot_product(vf90(2,:),vf90(3,:))
print '(a3,1x,F12.8)','v33', dot_product(vf90(3,:),vf90(3,:))

end subroutine

! keep the following for historical purposes
! DMR: sparse matrix matrix multiplication was stupid and slow.

!subroutine zsparse_matmat ( A,B,C )
!use mod_types, only: sparse
!
!type(sparse)                , intent(IN)  :: A,B
!complex(kind=8),allocatable , intent(out) :: C(:,:)
!complex(kind=8), parameter :: beta =(0.0d0,0.0d0)
!integer :: i,j,ier
!complex(kind=8) :: ab

!allocate(C(A%rdim,B%cdim), stat=ier)
!C = beta

!do i=1,A%nnzero
!  do j =1,B%nnzero 
!    ab = A%zvalues(i)*B%zvalues(j) 
!    C(A%rows(i),B%Columns(j))= C(A%rows(i),B%Columns(j)) + ab
!  end do
!end do

!end subroutine

!subroutine zsparse_matmat_mkl ( transa,A,B,C )
!use mod_types, only: sparse
!
!character(len=1) , intent(in) :: transa 
!type(sparse)                , intent(IN)  :: A,B
!complex(kind=8),allocatable , intent(out) :: C(:,:)
!complex(kind=8), parameter :: beta =(0.0d0,0.0d0),alpha= (1.0d0,0.0d0)
!integer :: i,j,ier
!complex(kind=8) :: ab
!character(len=1) :: matdescra(6)

!allocate(C(A%rdim,B%cdim), stat=ier)
!C = beta

!http://www.intel.com/software/products/mkl/docs/webhelp/bla/bla_IC.html#tbl2-6
!http://www.intel.com/software/products/mkl/docs/webhelp/bla/functn_mkl_dcoomm.html#functn_mkl_dcoomm
!matdescra(1) = 'G'
!matdescra(4) = 'F'
!call mkl_dcoosymv('U', a%ndim, a%values, a%rows, a%columns, a%nnzero, v, w)
!call mkl_zcoomm(transa,a%rdim, b%cdim, a%cdim, alpha, matdescra, &
!                       a%zvalues, a%rows, a%columns, &
!                       a%nnzero, b, ldb, alpha, c, ldc)

!end subroutine


!subroutine sparse_matmat ( A,B,C )
!use mod_types, only: sparse
!
!type(sparse)             , intent(IN)  :: A,B
!real(dp),    allocatable , intent(out) :: C(:,:)
!integer :: i,j,ier
!real(dp) :: ab

!allocate(C(A%rdim,B%cdim), stat=ier)
!C = zero

!do i=1,A%nnzero
!  do j =1,B%nnzero 
!    ab = A%zvalues(i) * B%zvalues(j) 
!    C(A%rows(i),B%Columns(j))= C(A%rows(i),B%Columns(j)) + ab
!  end do
!end do

!end subroutine

end module
