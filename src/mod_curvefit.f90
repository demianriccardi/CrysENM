module mod_curvefit
! DMR: July 16, 2008
! module FOR power fitting of f(x) = a*x^b + c... 
!                             f(x) - c = a*x^b
!                         log(f(x) - c)= log(a)+b*log(x)
!                              
!
!  fits the above analytical expression to a selected region of the data
!
!  adjustment of data:
!    1. curve is normalized by the "theoretical" maximum so it ends at 1.  
!       this constant goes into a
!    2. c set as the first value of f(x): subtracted out 
!
!  fitting of the parameters
!    I. b is searched with a do loop: at each point, least squares is used to find a and the 
!       residual. Best parameters give the lowest residual = sum (func_i - data_i)^2
!   II. linear fit of the log-log function using lsq...  
!

implicit none
type :: datafunc  ! introduced for functions that look like f(x) = a*x^{b}
  ! carries the original data (shifted to zero early on for cases like G(w))
  ! and the fit parameters
  real(8), allocatable,dimension(:) :: coor, func   ! original data adjusted by shift
  real(8), allocatable,dimension(:) :: afunc_min,afunc_log   ! original data adjusted by shift
  integer, dimension(2)             :: pntrs        ! ptrs that can be used for whatever
                                                    ! here they will point to region of orig data 
!
  real(8) :: acoor_offset,coor_offset, func_offset ! pnt at which func rises from zero and shift from zero     
  real(8) :: a,b,a_log,b_log,msd,msd_log ! _log are from the log fit
  integer :: norig,ncoor
  real(8) :: funcbeg,funcend
end type datafunc

real (8), parameter :: zero = 0.0d0, one = 1.0d0

contains

subroutine loglogfit (adfunc)
! use a search of offset and b to minimize residual: a,b,resid all in adfunc 
type(datafunc), intent(in out) :: adfunc
real(8),allocatable,dimension(:) :: logfunc,logcoord
real(8),allocatable,dimension(:) :: analyt_func
real(8) :: a,b,resid
integer :: i,j,ier

if (allocated(adfunc%afunc_log)) deallocate(adfunc%afunc_log)
allocate(logcoord(adfunc%ncoor-1),logfunc(adfunc%ncoor-1),stat=ier) ! minus 1 to avoid log(0)
allocate(analyt_func(adfunc%ncoor),adfunc%afunc_log(adfunc%ncoor),stat=ier) ! minus 1 to avoid log(0)

do i =1,adfunc%ncoor-1
  logcoord(i) = dlog10(adfunc%coor(i+1) - adfunc%acoor_offset) ! use the analytical offset
  logfunc(i)  = dlog10(adfunc%func(i+1))
  !print '(F5.2,4E15.5)',adfunc%acoor_offset, adfunc%coor(i), logcoord(i), adfunc%func(i),logfunc(i)
end do

call lfit_analyt_lsq_offset(logfunc,logcoord,b,a,resid)

adfunc%a_log     = 10.0d0**a
adfunc%b_log     = b

! compute the residual for new analytic function
analyt_func=0.0d0

do i = 1, adfunc%ncoor
  analyt_func(i) = axtotheb(adfunc%coor(i),adfunc%a_log,adfunc%b_log ,adfunc%acoor_offset)
  adfunc%afunc_log(i) = analyt_func(i)
end do                          
                    
adfunc%msd_log = sumsqr(adfunc%func,analyt_func)/real(size(adfunc%func))

deallocate(logcoord,logfunc)

end subroutine

subroutine minsearch (adfunc)
! use a search of offset and b to minimize residual: a,b,resid all in adfunc 
type(datafunc), intent(in out) :: adfunc
real(8),allocatable,dimension(:) :: analyt_func
real(8) :: a,b,resid,minres,offset,bestoffset,besta,bestb,dcor,corend
real(8) :: resol
integer :: i,j,ier

resol = 10.0d0

minres = 99999.99d0

if(adfunc%coor_offset.ne.zero) then

  dcor   = (adfunc%coor(2) - adfunc%coor(1))/resol
  offset = adfunc%coor_offset - dcor*resol/2.0
  corend = adfunc%coor_offset + dcor*resol/2.0 ! step up to halfway between the 1st and 2nd. 
                                             ! don't want to get to close to next point 
else

  dcor   = 9999999.00d0
  offset = 0.0d0 
  corend = offset

end if

if (allocated(adfunc%afunc_min)) deallocate(adfunc%afunc_min)
allocate(analyt_func(adfunc%ncoor),adfunc%afunc_min(adfunc%ncoor), stat = ier)

do while (offset .le. corend)

  b = 1.0d0
  do while (b .le. 16.0d0)


    analyt_func = 0.0d0

    a = 1.0

    do i=1,adfunc%ncoor
      analyt_func(i) = axtotheb(adfunc%coor(i),a,b,offset)
    end do

    call lfit_analyt_lsq(adfunc%func,analyt_func,a,resid)

    if (resid .lt. minres) then
      minres = resid
      besta = a
      bestb = b
      bestoffset = offset
    end if

    b = b + 0.01

  end do

  offset = offset + dcor

end do

adfunc%a              = besta
adfunc%b              = bestb
adfunc%acoor_offset   = bestoffset

do i=1,adfunc%ncoor
  analyt_func(i)      = axtotheb(adfunc%coor(i),adfunc%a,adfunc%b,adfunc%acoor_offset)
  adfunc%afunc_min(i) = analyt_func(i)
end do

adfunc%msd = sumsqr(adfunc%func,analyt_func)/real(size(adfunc%func))


!adfunc%msd       = minres

deallocate(analyt_func)

end subroutine

subroutine data_adjust(coor,func,funcmax,funcbeg,funcend,adfunc)
! 1. if the first value is nonzero (as it is with Gw) substract it off
! 2. find the last zero, this is the offset
! 3. save the part of function needed for fitting
! this is all this should do...  DMR FIX
! make the offset optional... ie. set to zero if not desired
!
! implimented the offset because the curves we're interested in depend on the
! force constant, and since the force constant that we are using is nonphysical, we need
! all for shifts due to diff force constants
!
real(8), allocatable, dimension(:), intent(in) :: coor,func
real(8), intent(in) :: funcmax ,& ! maximum of the function. div by this: 0 to 1 (maybe 1)
                       funcbeg, funcend ! begin and end: choose from 0.0 to 1.0
type(datafunc),intent(out) :: adfunc
logical :: l_offset
real(8), allocatable,dimension(:) :: tmpcor,tmpfnc,scor,sfnc
real(8), parameter :: zero=0.0d0,effzero=1.0D-16 !anything less than this is effectively zero
integer :: ntot,nset,num,ier,i,j,n,nbegin
real(8) :: offval,offset ! tmp vals for offset and shiftval

! set offset for everything
!l_offset = .true.
! set no offset for everything
l_offset = .false.

ntot         = size(coor)
adfunc%norig = ntot

allocate(tmpcor(ntot),tmpfnc(ntot),sfnc(ntot),scor(ntot), stat=ier)
tmpcor = coor 

! if the first value is not zero, set the yoffset to this value; it will be subtracted off

if (func(1) .gt. effzero) then 
  offval = func(1)
else
  offval = zero
end if

! let's make sure zeros are zeros and normalize by funcmax
do i =1, ntot
  tmpfnc(i) = func(i) - offval
  if (abs(tmpfnc(i)) .le. effzero) then
    tmpfnc(i) = zero
  else
    tmpfnc(i) = tmpfnc(i)/(funcmax-offval) ! DMR CHECK THIS
  end if
end do 

if(l_offset) then
! step up to the last zero
  i=1
  do while (tmpfnc(i) .eq. zero)
    offset = tmpcor(i)
    i=i+1
    nbegin = i-1
  end do
else

! set offset to zero and start from 1...  this seems to be the way to do it
nbegin = 1; offset = zero

end if

sfnc = 0.0d0
scor = 0.0d0

! copy in requested section
j=1
do i = nbegin, ntot

  if ((tmpfnc(i)  .le. funcend) .and. (tmpfnc(i) .ge. funcbeg)) then
    if (j .eq. 1) adfunc%pntrs(1)=i
    sfnc(j) = tmpfnc(i)
    scor(j) = tmpcor(i)
    j = j + 1
  end if
end do

nset = j - 1 ; adfunc%pntrs(2) = adfunc%pntrs(1) + nset - 1

num = Adfunc%pntrs(2)-adfunc%pntrs(1)+1 

if (allocated(adfunc%func)) deallocate(adfunc%func)
if (allocated(adfunc%coor)) deallocate(adfunc%coor)
allocate(adfunc%func(num),adfunc%coor(num), stat = ier)

adfunc%func(:)      = sfnc(1:nset)
adfunc%coor(:)      = scor(1:nset)
adfunc%coor_offset  = offset
adfunc%func_offset  = offval
adfunc%funcbeg      = funcbeg
adfunc%funcend      = funcend
adfunc%ncoor        = num

if (offset .ne. zero) then
  print *, 'offset not equal to zero for some reason; stop program'
  stop
end if

! now we have the data we want
deallocate(tmpcor,tmpfnc,sfnc,scor, stat=ier)

end subroutine

subroutine lfit_analyt_lsq_offset(refdat,compdat,m,b,bb)
!a:scales b:offset
real(8), intent(IN),allocatable :: refdat(:),compdat(:)
real(8), allocatable            :: newcompdat(:)
real(8), intent(out)            :: m,b,bb
real(8)                         :: avgref,avgcomp,refsqr,compsqr,ndat,mtop,btop, &
                                bot,dot_rc,dot_cc,sum_c,sum_r,rr,cc,tmpbb
integer :: i,ier

allocate (newcompdat(size(compdat)), stat = ier)

ndat=real(size(refdat))
dot_cc = dot_product(compdat,compdat)
dot_rc = dot_product(refdat,compdat)
sum_c  = sum(compdat)
sum_r  = sum(refdat)
avgref = sum_r/real(size(refdat))

!m
mtop=ndat*dot_rc-sum_r*sum_c
btop=dot_cc*sum_r - dot_rc*sum_c

bot=ndat*dot_cc-sum_c*sum_c
m=mtop/bot

!b
b=btop/bot

do i=1, size(compdat)
  newcompdat(i) =  m*compdat(i)+b
end do

bb = sumsqr(refdat,newcompdat)

deallocate(newcompdat)

return

end subroutine

subroutine lfit_analyt_lsq(refdat,compdat,a,b)
!a:scales b:heuristic offset, set to zero for now
real(8), intent(in out),allocatable  :: refdat(:),compdat(:)
real(8), intent(out)             :: a,b
real(8)                          :: ndat,dot_rc,dot_cc,sum_c,sum_r
real(8)                          :: avgref,avgcomp,tmpb,rr,cc

dot_cc = dot_product(compdat,compdat)
dot_rc = dot_product(refdat,compdat)
avgcomp = sum(compdat)/real(size(compdat))
avgref  = sum(refdat)/real(size(refdat))


a=dot_rc/dot_cc

compdat = compdat * a

b = sumsqr(refdat,compdat)

return

end subroutine

real(8) function sumsqr (refdat,compdat)
real(8), intent(in),allocatable  :: refdat(:),compdat(:)
integer :: i

sumsqr = 0.0d0

do i=1, size(refdat)
sumsqr = sumsqr + (refdat(i)-compdat(i))**2
end do

end function

real(8) function axtotheb (x,a,b,offset)
real(8) :: x,a,b,offset

if (offset .gt. x) then  ! don't want -x**b
  axtotheb = 0.0d0  
else
  axtotheb = a*(x-offset)**b
end if

end function

end
