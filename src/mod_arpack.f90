module mod_arpack
! fortran 90 module interface for arpack
! subroutines were made from the corresponding examples
! DMR 01-11-2008
!
! public:
!   subroutine mod_arpack_dsdrv1 (N,NEIG,A,VAL,VEC)
!   subroutine mod_arpack_dsdrv2 (N,NEIG,A,VAL,VEC)
!   subroutine mod_arpack_zndrv1 (N,NEIG,A,VAL,VEC)
!
! private:
!   matrix multiplier real
!     subroutine av_sparse ( n,A, v, w )
!   matrix multiplier complex
!     subroutine zav_sparse ( n,A, v, w )
!     subroutine cmplx_insrt_sort (n,nval,A,vec)
!
! todo:
!       1. is shift-invert mode viable???  it seems way faster, but volatile 
!          it appears to be used by Chao Yang for similar problems:
!          see http://crd.lbl.gov/~chao/nma.html
!          I contacted him about it; he responded initially, but did not respond
!          after more details were given.  Matlab had similar problems via it's eigs
!          command, for the 100x100 circulant matrix.
!         a. if it is viable, bring in zndrv2! for the complex case
!       2. make parallel
!       3. fix up memory handling a bit, if possible
!
!   files
!
!\BeginLib
!
!\Author of arpack and/or example files from which these subroutines were built
!     Richard Lehoucq
!     Danny Sorensen
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!     %-----------------------------%
!     | Define leading dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
!      integer, parameter :: dp = lselected_real_kind(14,80)
!include 'debug.h'
!INTEGER, PARAMETER :: DP = KIND (1.d0) 
use cfml_globaldeps,            only: sp,dp
use mod_types                       
use mod_constants,       only: one,two,zero,pi
implicit none

contains

subroutine mod_arpack_dsdrv1_sparse (input,NEIG,A,VAL,VEC)
! Routines called:
!     dsaupd  ARPACK reverse communication interface routine.
!     dseupd  ARPACK routine that returns Ritz values and Ritz vectors.
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     av_sparse   Matrix vector multiplication routine that computes A*x with sparse storage (coor).
!
      type(inp_par), intent (in)          :: input
      integer, intent(INOUT)              :: NEIG
      type(sparse),intent(IN)             :: A
      real(dp), intent(out), allocatable  :: VAL(:),VEC(:,:)
      integer                             :: n,ldv

!     %--------------%
!     | Local Arrays |
!     %--------------%
      real(dp), pointer             :: pntr1(:),pntr2(:),ax(:)
      real(dp), allocatable, target :: workd(:), v(:,:)
!      real(dp), allocatable               :: v(:,:),workl(:),workd(:),d(:,:),resid(:),ax(:)

      real(dp), allocatable               :: workl(:),d(:,:),resid(:)
      logical,  allocatable               :: lselect(:)
      integer                             :: iparam(11), ipntr(11)
!     %---------------%
!     | Local Scalars |
!     %---------------%
      character ::     bmat*1, which*2
      integer   ::     ido, nev, ncv, lworkl, info, ierr, j, &
                       nconv, maxitr, mode, ishfts
! DMR for allocation
      integer   ::     ialloc,ier,i
      logical   ::     rvec
      real(dp)  ::     tol, sigma
      real      ::     time1,time2
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      real(dp),external  ::     dnrm2  ! function
      external  daxpy                  ! subroutine

!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
      intrinsic ::     abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %----------------------------------------------------%
!     | dimension of the matrix.  A standard eigenvalue    |
!     | problem is solved (BMAT = 'I'). NEV is the number  |
!     | of eigenvalues to be approximated.  The user can   |
!     | modify NEV, NCV, WHICH to solve problems of        |
!     | different sizes, and to get different parts of the |
!     | spectrum.  However, The following conditions must  |
!     | be satisfied:                                      |
!     |                   N <= MAXN,                       | 
!     |                 NEV <= MAXNEV,                     |
!     |             NEV + 1 <= NCV <= MAXNCV               | 
!     %----------------------------------------------------% 
!
      n=A%ndim
      ldv=n
      nev =  neig
      ncv =  neig*2
      lworkl  = ncv*(ncv+8)

      allocate(  v(ldv,ncv), workl(lworkl), &
                  workd(3*n), d(ncv,2), resid(n),lselect(ncv), &
                  ax(n), stat=ialloc)

      if (ialloc /= 0) STOP "*** not enough memory ***"

      call cpu_time(time1)

      open(unit=12,file=trim(input%fileroot)//"-"//trim(input%genm)//"-"//trim(input%tpbc)//"-"// &
                        trim(input%fctyp) //"-"//trim(input%bztyp)//"-arpack.txt", &
      status='unknown',action='write', iostat=ier)

      bmat  = 'I'
      which = 'SM'
!
!   %--------------------------------------------------%
!   | array WORKL is used in DSAUPD as workspace.      |
!   | if TOL<=0, machine precision is used.            | 
!   | IDO is used for rev comm and is init set to 0.   |
!   | INFO=0 random vector in DSAUPD to start the      | 
!   | Arnoldi iteration.                               |
!   %--------------------------------------------------%
!
      tol     = zero !1.0E-08  
      info    = 0
      ido     = 0

!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 1600
      mode   = 1
      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode 

!     M A I N   L O O P (Reverse communication) |
 10   continue
!   Repeatedly call DSAUPD and proceed accord to param IDO  
!   either convergence or maxitr   
!
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, &
                       ncv, v, ldv, iparam, ipntr, workd, workl, &
                       lworkl, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%
!
 
            pntr1 => workd(ipntr(1):ipntr(1)+N)
            pntr2 => workd(ipntr(2):ipntr(2)+N)
            call av_sparse (A,pntr1,pntr2)
            !call av_sparse (n,A,workd(ipntr(1)),workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         end if 
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
!
      else 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%
!           
         rvec = .true.
!
         call dseupd ( rvec, 'All', lselect, d, v, ldv, sigma, & 
             bmat, n, which, nev, tol, resid, ncv, v, ldv,   &
             iparam, ipntr, workd, workl, lworkl, ierr )

!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         if ( ierr .ne. 0) then
!
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%
!
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
!
         else
!
             nconv =  iparam(5)
             do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
                pntr1 => v(1:N,j)          
                call av_sparse(A, pntr1, ax)
                !call av_sparse(n,A, v(1,j), ax)
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
!
 20          continue
!
!            %-------------------------------%
!            | Display computed residuals    |
!            %-------------------------------%
!
             call dmout(12, nconv, 2, d, ncv, -6, &
                 'Ritz values and relative residuals')
         end if
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
      call cpu_time(time2)
      time2=time2-time1
 
      call arp_wrtout (time2,12,"DSDRV1",info,n,nev,ncv,which,nconv,iparam,tol)

      close(12)

!
      end if
!
!     %---------------------------%
!     | Done with program dsdrv1. |
!     %---------------------------%
!
 9000 continue

!DMR allocate the output arrays: VAL and VEC
      if (allocated(val)) deallocate(val)
      if (allocated(vec)) deallocate(vec)
      allocate(val(nconv),vec(n,nconv),stat=ialloc)
      if (ialloc /= 0) STOP "*** not enough memory ***"

      if (nconv .lt. nev) then
        print *, 'hoots mon! fewer converged than requested. see arpack.out'
      end if 

      do i=1,nconv
        val(i)=d(i,1)
        vec(:,i)=v(1:n,i)
      end do

      nev=nconv

end subroutine

subroutine mod_arpack_zndrv1 (input,NEIG,A,VAL,VEC,logical_rvec)
! Routines called:
!     dsaupd  ARPACK reverse communication interface routine.
!     dseupd  ARPACK routine that returns Ritz values and Ritz vectors.
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     av_sparse   Matrix vector multiplication routine that computes A*x with sparse storage (coor).
!
      type(inp_par), intent (in)                  :: input
      integer, intent(INOUT)                      :: NEIG
      type(sparse),intent(IN)                     :: A
      complex(kind=8), intent(out), allocatable   :: VAL(:),VEC(:,:)
      logical,optional,intent(in)                 :: logical_rvec
      integer                                     :: n,ldv

!     %--------------%
!     | Local Arrays |
!     %--------------%
      complex(kind=8), allocatable  :: workl(:),workev(:),d(:),resid(:)
      complex(kind=8), pointer             :: pntr1(:),pntr2(:),ax(:)
      complex(kind=8), allocatable, target :: workd(:), v(:,:)
      real(dp), allocatable         :: rwork(:),rd(:,:)                            
      logical,  allocatable         :: lselect(:)
      integer                       :: iparam(11), ipntr(11)
      complex(kind=8), allocatable  :: arrayin(:),arrayout(:)
!     %---------------%
!     | Local Scalars |
!     %---------------%
      character ::     bmat*1, which*2
      integer   ::     ido, nev, ncv, lworkl, info, ierr, j, &
                       nconv, maxitr, mode, ishfts
! DMR for allocation
      integer         :: ialloc,ier,i
      logical         :: rvec
      complex(kind=8) :: sigma
      real(dp)        :: tol
      real            :: time1,time2
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%

      real(dp),external  :: dznrm2 , dlapy2  ! function
      external   zaxpy                  ! subroutine

!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
      intrinsic ::     abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %----------------------------------------------------%
!     | dimension of the matrix.  A standard eigenvalue    |
!     | problem is solved (BMAT = 'I'). NEV is the number  |
!     | of eigenvalues to be approximated.  The user can   |
!     | modify NEV, NCV, WHICH to solve problems of        |
!     | different sizes, and to get different parts of the |
!     | spectrum.  However, The following conditions must  |
!     | be satisfied:                                      |
!     |                   N <= MAXN,                       | 
!     |                 NEV <= MAXNEV,                     |
!     |             NEV + 1 <= NCV <= MAXNCV               | 
!     %----------------------------------------------------% 
!

      n=A%ndim

      ldv=n
      nev =  neig
      ncv =  neig*2
      lworkl  = 3*ncv*ncv+5*ncv !3*ncv*(ncv+2)

      allocate(  v(ldv,ncv), workl(lworkl),workev(3*ncv), &
                  workd(3*n), d(ncv), resid(n),lselect(ncv), &
                  ax(n),rwork(ncv), rd(ncv,3), stat=ialloc)
      if (ialloc /= 0) STOP "*** not enough memory ***"

! for the matvec      
      allocate( arrayin(A%ndim),arrayout(A%ndim), stat=ialloc)
      if (ialloc /= 0) STOP "*** not enough memory ***"

      call cpu_time(time1)

      open(unit=12,file=trim(input%fileroot)//"-"//trim(input%genm)//"-"//trim(input%tpbc)//"-"// &
                        trim(input%fctyp) //"-"//trim(input%bztyp)//"-arpack.txt", &
      status='unknown',action='write', iostat=ier)

!      print *, 'shit'
      bmat  = 'I'
      which = 'SM'

!
!   %--------------------------------------------------%
!   | array WORKL is used in DSAUPD as workspace.      |
!   | if TOL<=0, machine precision is used.            | 
!   | IDO is used for rev comm and is init set to 0.   |
!   | INFO=0 random vector in DSAUPD to start the      | 
!   | Arnoldi iteration.                               |
!   %--------------------------------------------------%
!
      tol     = zero !1.0E-08  
      info    = 0
      ido     = 0

!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 5600
      mode   = 1
      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode 

!print *, ncv 

!     M A I N   L O O P (Reverse communication) |
 10   continue
!   Repeatedly call ZNAUPD and proceed accord to param IDO  
!   either convergence or maxitr   
!

         call znaupd  ( ido, bmat, n, which, nev, tol, resid, ncv, &
              v, ldv, iparam, ipntr, workd, workl, lworkl,         &
              rwork,info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%
!
            pntr1 => workd(ipntr(1):ipntr(1)+N)
            pntr2 => workd(ipntr(2):ipntr(2)+N)
            call zzav_sparse (A,pntr1,pntr2)
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         end if 
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in ZNAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
!
      else 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%
!           
        if(present(logical_rvec)) then
          rvec = logical_rvec
        else
          rvec = .true.
        end if

!

         call zneupd  (rvec, 'A', lselect, d, v, ldv, sigma, &
              workev, bmat, n, which, nev, tol, resid, ncv, &
              v, ldv, iparam, ipntr, workd, workl, lworkl,  &
              rwork, ierr)

!        %----------------------------------------------%
!        | Eigenvalues are returned in the one          |
!        | dimensional array D.  The corresponding      |
!        | eigenvectors are returned in the first NCONV |
!        | (=IPARAM(5)) columns of the two dimensional  | 
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%         


      if ( ierr .ne. 0) then
!
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%
!
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
!
         else
!
             nconv =  iparam(5)
             do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
                pntr1 => v(1:N,j)          
                call zzav_sparse(A, pntr1, ax)
                !call zzav_sparse(n,A, v(1,j), ax)
                call zaxpy(n, -d(j), v(1,j), 1, ax, 1)
                rd(j,1) = dble (d(j))
                rd(j,2) = dimag (d(j))
                rd(j,3) = dznrm2 (n, ax, 1)
                rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
!
 20          continue
!
!            %-------------------------------%
!            | Display computed residuals    |
!            %-------------------------------%
!
             call dmout(12, nconv, 3, rd, ncv, -6, &
                 'Ritz values (Real, Imag) and relative residuals')
         end if
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
      call cpu_time(time2)
      time2=time2-time1
  
      call arp_wrtout (time2,12,"DSDRV1",info,n,nev,ncv,which,nconv,iparam,tol)
      close(12)

!
      end if
!
!     %---------------------------%
!     | Done with program dsdrv1. |
!     %---------------------------%
!
 9000 continue

!DMR allocate the output arrays: VAL and VEC
      if (allocated(val)) deallocate(val)
      if (allocated(vec)) deallocate(vec)
      allocate(val(nconv),vec(n,nconv),stat=ialloc)
      if (ialloc /= 0) STOP "*** not enough memory ***"

      if (nconv .lt. nev) then
        print *, 'hoots mon! fewer converged than requested. see arpack.out'
      end if 

      do i=1,nconv
        val(i)=d(i)
        vec(:,i)=v(1:n,i)
      end do

      call cmplx_insrt_sort (n,nconv,val,vec)
      nev=nconv

end subroutine

subroutine cmplx_insrt_sort (n,nval,A,vec)
integer,        intent(in)                 :: n,nval
complex(kind=8),intent(in out),allocatable :: A(:),vec(:,:)
complex(kind=8)                            :: tmpval
complex(kind=8),allocatable                :: tmpvec(:)
integer :: i,j,ialloc,ier

if (allocated(tmpvec)) deallocate(tmpvec)
allocate(tmpvec(n),stat=ialloc)
if (ier /= 0) STOP "*** not enough memory ***"

do i = 2,nval
  tmpval = A(i)
  tmpvec(:) = vec(:,i)
  j=i-1
  do while ((j.ge.1) .and. (real(A(j)) .gt. real(tmpval)))
    A(j+1)=A(j)
    vec(:,j+1)=vec(:,j)
    j=j-1
  end do
  A(j+1)=tmpval
  vec(:,j+1)=tmpvec(:)
end do

end subroutine

subroutine restart_vector_write (n,nconv,vecs,vals)
! DMR 02082008.  this restart vector approach did not prove useful
!
! from arpack userguide: http://www.caam.rice.edu/software/ARPACK/UG/
! The parameter info should be set to 0 on the initial call to dsaupd unless 
! the user wants to supply the starting vector that initializes the IRLM. Normally, 
! this default is a reasonable choice. However, if this eigenvalue calculation is one
! of a sequence of closely related problems then convergence may be accelerated 
! if a suitable starting vector is specified. Typical choices in this situation might 
! be to use the final value of the starting vector from the previous eigenvalue 
! calculation (that vector will already be in the first column of V) or to construct 
! a starting vector by taking a linear combination of the computed eigenvectors from 
! the previously converged eigenvalue calculation. If the starting vector is to be 
! supplied, then it should be placed in the array resid and info should be set to 1 
! on entry to dsaupd. On completion, the parameter info may contain the value 0 indicating 
! the iteration was successful or it may contain a nonzero value indicating an error 
! or a warning condition. The meaning of a nonzero value returned in info may 
! be found in the header comments of the subroutine dsaupd.

integer, intent(in)               :: n,nconv
real(dp), allocatable, intent(in) :: vecs(:,:),vals(:)
real(dp)                          :: restart_vec(n)
integer                           :: i,j,ier

open(unit=12,file="eigenvect.restart", status='unknown',action='write', iostat=ier)

restart_vec=zero

do j=1,n
  do i=1,nconv
    if (abs(vals(i)) .gt. 1.0D-07) then
      !linear combo for restart
      restart_vec(j)=restart_vec(j)+vecs(j,i)/vals(i)
    end if
  end do
end do

do i=1,n
  write(12,*) restart_vec(i)
end do 

close(12)

end subroutine

subroutine av_sparse ( A, v, w )
! lapack may have faster ones
!
type(sparse), intent(IN) :: A       
real(dp),pointer,intent(IN)      :: v(:)       
real(dp),pointer,intent(out)     :: w(:)
real(dp),parameter               :: beta = 0.0d0
integer :: i

!call mkl_dcsrsymv('U', n, a%values, a%rows, a%columns, v, w)
!call mkl_dcoosymv('U', a%ndim, a%values, a%rows, a%columns, a%nnzero, v, w)
!call dcoosymv('U', a%ndim, a%values, a%rows, a%columns, a%nnzero, v, w)
w(1:A%ndim)=beta

do i=1,A%nnzero
  w(A%rows(i))=w(A%rows(i))+A%zvalues(i)*v(A%columns(i))
end do

end subroutine



subroutine zzav_sparse ( A, v, w )
!
type(sparse)           , intent(IN) :: A      
complex(kind=8), pointer,intent(IN)     :: v(:)
complex(kind=8), pointer,intent(out)    :: w(:)
complex(kind=8), parameter :: beta =(0.0d0,0.0d0)
integer :: i

w(1:A%ndim)=beta

do i=1,A%nnzero
  w(A%rows(i))=w(A%rows(i))+A%zvalues(i)*v(A%columns(i))
end do

end subroutine


subroutine zav_sparse ( n,A, v, w )
!use mkl95_blas
! lapack may have faster ones
!
integer, intent(IN)      :: N
complex(kind=8), allocatable, intent(IN)     :: A(:,:)
complex(kind=8), pointer,intent(IN)     :: v(:)
complex(kind=8), pointer,intent(out)    :: w(:)
complex(kind=8), parameter :: alpha=(1.0d0,0.0d0)
complex(kind=8), parameter :: beta =(0.0d0,0.0d0)

!call hemv(a,v,w,'u')   ! doesn't work???

!w = matmul(A,v)
call zhemv('U', n, alpha, a, n, v, 1, beta, w, 1)

end subroutine

subroutine arp_wrtout (time,un,atyp,info,n,nev,ncv,which,nconv,iparam,tol)
real(sp)               :: time
integer,   intent(in)  :: un,info,n,nev,ncv,nconv,iparam(11)
real(dp),  intent(in)  :: tol
character, intent(in)  :: which*2,atyp*6
integer :: ier

!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' '
            print *, ' No shifts could be applied during implicit ', &
                      'Arnoldi update, try increasing NCV.'
            print *, ' '
         end if
!
         write(un,*) ' ====== ', trim(atyp), '======='
         write(un,*) ' '
         write(un,*) ' Size of the matrix is ', n
         write(un,*) ' The number of Ritz values requested is ', nev
         write(un,*) ' The number of Arnoldi vectors generated', &
                  ' (NCV) is ', ncv
         write(un,*) ' What portion of the spectrum: ', which
         write(un,*) ' The number of converged Ritz values is ', &
                    nconv
         write(un,*) ' The number of Implicit Arnoldi update', &
                  ' iterations taken is ', iparam(3)
         write(un,*) ' The number of OP*x is ', iparam(9)
         write(un,*) ' The convergence criterion is ', tol
         write(un,*) ' '
         write(un,*) ' Time taken: ', time
         write(un,*) ' '
!
!
end subroutine

end module

