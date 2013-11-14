module mod_math
!
!
!        

use cfml_globaldeps,                 only: dp,sp
use mod_constants

implicit none

contains

subroutine analyt_lsq(refdat,compdat,a,b)
!a:scales b:heuristic offset, set to zero for now
real(dp), intent(in)          :: refdat(:),compdat(:)
real(dp), intent(out)         :: a,b              
real(dp)                      :: ndat,dot_rc,dot_cc,sum_c,sum_r

dot_cc = dot_product(compdat,compdat)
dot_rc = dot_product(refdat,compdat)

a=dot_rc/dot_cc
b=zero

return

end subroutine


subroutine analyt_lsq_offset(refdat,compdat,a,b)
!a:scales b:offset
real(dp), intent(IN)          :: refdat(:),compdat(:)
real(dp), intent(out)         :: a,b
real(dp)                      :: ndat,atop,btop,bot,dot_rc,dot_cc,sum_c,sum_r

ndat=real(size(refdat))
dot_cc = dot_product(compdat,compdat)
dot_rc = dot_product(refdat,compdat)
sum_c  = sum(compdat)
sum_r  = sum(refdat)

!a
atop=ndat*dot_rc-sum_r*sum_c
btop=dot_cc*sum_r - dot_rc*sum_c

bot=ndat*dot_cc-sum_c*sum_c
a=atop/bot

!b
b=btop/bot

return

end subroutine

real(dp) function dist_sqr(dxyz)
real(dp), intent(in) :: dxyz(3)

dist_sqr = dot_product(dxyz,dxyz)

end function

real(dp) function distance(dxyz)
real(dp) , intent(in) :: dxyz(3)

distance = dsqrt(dot_product(dxyz,dxyz))

end function

subroutine vec_avg_std(vector,avg,std)
real(dp), allocatable, intent(in)  :: vector(:)
real(dp),              intent(out) :: avg, std
integer :: i
real(dp) :: sumsqr

avg=sum(vector)/real(size(vector))

do i = 1,size(vector)
 sumsqr = sumsqr + (vector(i) - avg)**2
end do

std = dsqrt(sumsqr/real(size(vector)-1))

end subroutine


subroutine vec_avg_rmsd(vector,avg,rmsd)
real(dp), allocatable, intent(in)  :: vector(:)
real(dp),              intent(out) :: avg, rmsd
integer :: i
real(dp) :: sumsqr

avg=sum(vector)/real(size(vector))

do i = 1,size(vector)
 sumsqr = sumsqr + (vector(i) - avg)**2
end do

rmsd = dsqrt(sumsqr/real(size(vector)))

end subroutine

subroutine linear_correl(refdat,compdat,correl)
! equation 5 from kundus, phillips Biophys. J. 83(2) 723-732, 2002
!!!!  DMR there's a problem!!
real(dp), intent(in),allocatable :: refdat(:),compdat(:)
real(dp), intent(out) :: correl
integer :: ndat,i
real(dp) :: avg_r,avg_c,top_sum,bot_sum_r,bot_sum_c

! DMR: I found a weird bug.  If the top_sum, etc variables are declared above with
!      real(dp) :: top_sum=zero,bot_sum_r=zero, bot_sum_c=zero
!      when this subroutine or function is called more than once, the values aren't initialized
!      at zero!! it just adds on to what it was from the last time.  Not sure why this is
!      if you are looking at this code, and know why, email me the answer.
!
top_sum=zero
bot_sum_r=zero
bot_sum_c=zero

ndat=size(refdat)

avg_r=sum(refdat)/real(ndat)
avg_c=sum(compdat)/real(ndat)

do i=1,ndat
  top_sum=top_sum + (refdat(i)-avg_r)*(compdat(i)-avg_c)
  bot_sum_r=bot_sum_r + (refdat(i)-avg_r)**2 
  bot_sum_c=bot_sum_c + (compdat(i)-avg_c)**2 
end do

if (abs(bot_sum_r) .lt. 1.0d-8) then
  print *, 'constant reference set?'
  print *, 'linear_correl> sum {(exp-<exp>)^2:}',bot_sum_r
  print *, 'linear_correl> sum {(thr-<thr>)^2:}', bot_sum_c
  print *, 'linear_correl> setting linear-correlation to scale by average'
  correl = avg_r/avg_c
  
else
  correl = top_sum/dsqrt(bot_sum_r*bot_sum_c)
end if

end subroutine

!real(dp) function linear_correl(refdat,compdat)
! equation 5 from kundus, phillips Biophys. J. 83(2) 723-732, 2002

!real(dp), intent(in),allocatable :: refdat(:),compdat(:)
!integer :: ndat,i
!real(dp) :: avg_r,avg_c,top_sum=zero,bot_sum_r=zero,bot_sum_c=zero

!ndat=size(refdat)
!avg_r=sum(refdat)/real(ndat)
!avg_c=sum(compdat)/real(ndat)

!do i=1,ndat
!  top_sum=top_sum + (refdat(i)-avg_r)*(compdat(i)-avg_c)
!  bot_sum_r=bot_sum_r + (refdat(i)-avg_r)**2 
!  bot_sum_c=bot_sum_c + (compdat(i)-avg_c)**2 
!end do

!linear_correl = top_sum/dsqrt(bot_sum_r*bot_sum_c)

!end function 

subroutine realrecip(cell,ang, &
                         rcell,rang)
!GET RECIPROCAL CELL CONSTANTS FROM REAL ONES
!from GNPs lib
real(sp),dimension(3),intent(in)   :: cell,ang
real(dp),dimension(3),intent(out)  :: rcell,rang
real(dp)                            :: a,b,c,alpha,beta,gamma
real(dp)                            :: AS,BS,CS,ALPHAS,BETAS,GAMMAS
real(dp)              :: V,CA,CB,CG,SA,SB,SG

a=cell(1)
b=cell(2)
c=cell(3)
alpha = ang(1)*pi/180.0d0
beta  = ang(2)*pi/180.0d0
gamma = ang(3)*pi/180.0d0

  CA=DCOS(ALPHA)
  CB=DCOS(BETA)
  CG=DCOS(GAMMA)
  SA=DSIN(ALPHA)
  SB=DSIN(BETA)
  SG=DSIN(GAMMA)
  V=A*B*C*DSQRT(one+two*CA*CB*CG-CA**2-CB**2-CG**2)
  AS=B*C*SA/V
  BS=C*A*SB/V
  CS=A*B*SG/V
  ALPHAS=DACOS((CB*CG-CA)/(SB*SG))
  BETAS =DACOS((CG*CA-CB)/(SG*SA))
  GAMMAS=DACOS((CA*CB-CG)/(SA*SB))

rcell(1) = as
rcell(2) = bs
rcell(3) = cs
rang(1)  = alphas*180.0d0/pi
rang(2)  = betas*180.0d0/pi
rang(3)  = gammas*180.0d0/pi

end subroutine


SUBROUTINE INVERT(A,AINVERSE)
!from GNPs lib
  real(DP), intent(in) :: A(3,3)
  real(DP), intent(out) :: AINVERSE(3,3)
  real(DP) :: D

  D=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+ &
            A(2,1)*(A(3,2)*A(1,3)-A(1,2)*A(3,3))+ &
            A(3,1)*(A(1,2)*A(2,3)-A(1,3)*A(2,2))
  IF(D.EQ.0.0) STOP 'MATRIX IS SINGULAR AND CANNOT BE INVERTED'
  AINVERSE(1,1)= (A(2,2)*A(3,3)-A(2,3)*A(3,2))/D
  AINVERSE(1,2)=-(A(1,2)*A(3,3)-A(1,3)*A(3,2))/D
  AINVERSE(1,3)= (A(1,2)*A(2,3)-A(1,3)*A(2,2))/D
  AINVERSE(2,1)=-(A(2,1)*A(3,3)-A(3,1)*A(2,3))/D
  AINVERSE(2,2)= (A(1,1)*A(3,3)-A(1,3)*A(3,1))/D
  AINVERSE(2,3)=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))/D
  AINVERSE(3,1)= (A(2,1)*A(3,2)-A(3,1)*A(2,2))/D
  AINVERSE(3,2)=-(A(1,1)*A(3,2)-A(1,2)*A(3,1))/D
  AINVERSE(3,3)= (A(1,1)*A(2,2)-A(1,2)*A(2,1))/D
  RETURN
end subroutine

subroutine realort(A,B,C,ALPHA,BETA,GAMMA,ORT)
!from GNPs lib
  real(dp),intent(in)  :: A,B,C,ALPHA,BETA,GAMMA
  real(dp),intent(out) :: ORT(3,3)
  real(dp) :: CA,CB,CG,SB,SG

  CA=DCOS(ALPHA)
  CB=DCOS(BETA)
  CG=DCOS(GAMMA)
  SB=DSIN(BETA)
  SG=DSIN(GAMMA)
  ORT(1,1)=A
  ORT(1,2)=B*CG
  ORT(1,3)=C*CB
  ORT(2,1)=ZERO
  ORT(2,2)=B*SG
  ORT(2,3)=C*(CA-CB*CG)/SG
  ORT(3,1)=ZERO
  ORT(3,2)=ZERO
  ORT(3,3)=C*SB*SQRT(one-((CB*CG-CA)/(SB*SG))**2)

end subroutine

real(dp) function determ3(mat)
!determinant of a 3x3 matrix
  real(dp),intent(in)  :: mat(3,3)

  determ3 =  mat(1,1)*(mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)) &
           - mat(1,2)*(mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1)) &
           + mat(1,3)*(mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))

end function

subroutine valvec_inv(vals,vecs,inv,nzeros)
real(dp),        allocatable,intent(in)     :: vals(:),vecs(:,:)
real(dp),        allocatable,intent(in out) :: inv(:,:)
integer,                     intent(out)    :: nzeros
real(dp)        :: tmpval
integer ::   nvals, ndim
integer ::   i,j,k

nvals  = size(vals)
ndim   = size(vecs(:,1))
nzeros = 0

do i=1,nvals
  if (abs(vals(i)).gt.1.0D-08) then
    do j=1,ndim
      do k=1,ndim
        tmpval    = vecs(k,i)*vecs(j,i)
        inv(k,j) = inv(k,j) + tmpval/vals(i)
      end do
    end do
  else
    nzeros = nzeros + 1
  end if
end do

end subroutine


!GET RESOLUTION OF A REFLECTION
!(REQUIRES MATRIX FROM RECIPORT)
!RETURNS VALUE DSTAR IN RECIPROCAL ANGSTROMS
!
!FROM PAGE 19 IN PRINCE
!
real(DP) function dstar(v,iort)
!from GNPs lib
  real(DP),intent(in) :: v(3),iort(3,3)
  real(DP)            :: a(3),tmp

  a=matmul(v,iort)
  tmp=dot_product(a,a)
  dstar=dsqrt(tmp)

end function

! for the dstar squared
real(dp) function dstar_sqr(v,iort)
real(dp),intent(in) :: v(3),iort(3,3)
real(dp)            :: a(3)

  a=matmul(v,iort)
  dstar_sqr=dot_product(a,a)

end function

!
!
!      ________________________________________________________
!     |                                                        |
!     | GIVEN THE SINGULAR VALUE DECOMPOSITION QDP TRANSPOSE OF|
!     | A GENERAL MATRIX, COMPUTE ITS REGULARIZED PSEUDOINVERSE|
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         LA    --LEADING (ROW) DIMENSION OF ARRAY A     |
!     |                                                        |
!     |         Q     --FIRST FACTOR IN THE SVD DECOMPOSITION  |
!     |                                                        |
!     |         LQ    --LEADING (ROW) DIMENSION OF ARRAY Q     |
!     |                                                        |
!     |         MQ    --ROW DIMENSION OF MATRIX STORED IN Q    |
!     |                                                        |
!     |         D     --SINGULAR VALUES                        |
!     |                                                        |
!     |         P     --LAST FACTOR IN THE SVD DECOMPOSITION   |
!     |                                                        |
!     |         LP    --LEADING (ROW) DIMENSION OF ARRAY P     |
!     |                                                        |
!     |         MP    --ROW DIMENSION OF MATRIX STORED IN P    |
!     |                                                        |
!     |         R     --REGULARIZATION PARAMETER               |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |                                                        |
!     |         A     --THE REGULARIZED PSEUDOINVERSE          |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS,MIN0,dsqrt                    |
!     |________________________________________________________|
!
        SUBROUTINE PSEUDO(A,LA,Q,LQ,MQ,D,P,LP,MP,R)
        INTEGER I,J,K,LA,LP,LQ,MP,MQ,N
        REAL(dp) A(LA,1),D(1),Q(LQ,1),P(LP,1),R,S,T,Y,Z
        N = MIN0(MQ,MP)
        Z = DABS(R)
        Y = dsqrt(Z)
        DO 10 J = 1,MQ
             DO 10 I = 1,MP
  10              A(I,J) = ZERO
        DO 50 K = 1,N
             IF ( D(K) .LE. 1.0D-10 ) GOTO 50
             IF ( D(K) .GT. Y ) GOTO 20
             S = D(K)/(D(K)**2+Z)
             GOTO 30
  20         S = 1./(D(K)+Z/D(K))
  30         DO 40 J = 1,MQ
                  T = S*Q(J,K)
                  DO 40 I = 1,MP
  40                   A(I,J) = A(I,J) + T*P(I,K)
  50    CONTINUE
        RETURN
end subroutine

! approximate integration using trapezoid method.  
! index starts at zero!
subroutine trapezoidal_int (coord,func,intfunc)
real(dp), allocatable,intent(in out)  :: coord(:), func(:)
real(dp), allocatable,intent(out)     :: intfunc(:)
integer :: i,n,ier

n = size(coord)

! make sure everything is in the right order
call dual_vectsort(coord,func)

allocate (intfunc(n), stat=ier)

do i = 1, n

  intfunc(i) = intfunc(i-1) + (coord(i)-coord(i-1))*(func(i)+func(i-1))/2.0d0

end do

end subroutine

subroutine dual_vectsort(coor,func)
real(dp), allocatable,intent(in out)  :: coor(:), func(:)
integer :: nval, i,j,ier
real(dp) :: tmpvalc,tmpvalf

nval = size(coor)-1

do i = 1,nval
  tmpvalc = coor(i)
  tmpvalf = func(i)
  j=i-1
  do while ((j.ge.0) .and. (coor(j) .gt. tmpvalc) )
    coor(j+1) = coor(j)
    func(j+1) = func(j)
    j=j-1
  end do
  coor(j+1) = tmpvalc
  func(j+1) = tmpvalf
end do

end subroutine


end module 

