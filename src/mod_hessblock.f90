module mod_hessblock
!
! DMR 05-13-2008
!
!  provides subroutines for constructing 3x3 hessian blocks for atom-atom interactions
!
!
!CrysFML modules
use cfml_globaldeps,      only: dp,sp
!external modules
use mod_constants, only: one,two,zero,pi,joule2cal

implicit none

! public subroutines: 
public :: sprnghessblock,reach_hessblock,reach_hessblock_temperat
! private subroutines 

contains
subroutine reach_hessblock_temperat(input,i,j,chaini,chainj,dist,dxyz,hess)
use mod_types,  only: inp_par
! quick and dirty version of REACH taken from 2009 table... tetramer
! input: distance and dxyz vector between two point
! outpt: 3x3 hessian for the spring interaction
!        with various forceconstants requested in input
! i and j are the residues
type(inp_par)         :: input
integer,  intent(in)  :: i,j,chaini,chainj
real(dp), intent(in)  :: dist,dxyz(3)
real(dp), intent(out) :: hess(3,3)
real(dp)              :: fcnst,dstsqr, xixj,xiyj,xizj,yiyj,yizj,zizj
integer :: i_minus_j
logical :: chains
real(dp) :: b,nonbond_fcnst

fcnst = zero
i_minus_j = abs(i-j)

if (chaini .eq. chainj) then
  chains = .true.
else
  chains = .false.
end if

!nonbond fcnst
b = 0.0011*input%temperature + 0.41d0
nonbond_fcnst = expo_fcnst(1040.0d0*joule2cal,b,dist)

    select case (chains) ! intra or inter
      case(.true.)

        if (input%temperature .gt. 170.0d0) then
          select case (i_minus_j)
            case(0) ! self
              fcnst = zero
            case(1) ! k12
              fcnst = (-0.955d0*input%temperature +1020.0d0) *joule2cal
            case(2) ! k13
              fcnst = (-0.155d0*input%temperature+55.3d0)*joule2cal
            case(3) ! k14
              fcnst = (-0.0541d0*input%temperature + 43.7d0)*joule2cal
            case default
              fcnst = nonbond_fcnst
          end select  
        else
          select case (i_minus_j)
            case(0) ! self
              fcnst = zero
            case(1) ! k12
              fcnst = 847.0d0 *joule2cal
            case(2) ! k13
              fcnst = 29.0d0*joule2cal
            case(3) ! k14
              fcnst = (-0.0541d0*input%temperature + 43.7d0)*joule2cal
            case default
              fcnst = nonbond_fcnst
          end select  
        end if
      case default
        fcnst = nonbond_fcnst
    end select

call hessian_block(dist,dxyz,fcnst,hess)

end subroutine

subroutine reach_hessblock(input,i,j,chaini,chainj,dist,dxyz,hess)
use mod_types,  only: inp_par
! quick and dirty version of REACH taken from 2009 table... tetramer
! input: distance and dxyz vector between two point
! outpt: 3x3 hessian for the spring interaction
!        with various forceconstants requested in input
! i and j are the residues
type(inp_par)         :: input
integer,  intent(in)  :: i,j,chaini,chainj
real(dp), intent(in)  :: dist,dxyz(3)
real(dp), intent(out) :: hess(3,3)
real(dp)              :: fcnst,dstsqr, xixj,xiyj,xizj,yiyj,yizj,zizj
integer :: i_minus_j
logical :: chains
real(dp) :: b, nonbond_fcnst

fcnst = zero
i_minus_j = abs(i-j)

if (chaini .eq. chainj) then
  chains = .true.
else
  chains = .false.
end if

!nonbond fcnst
  select case (trim(input%tpbc))
    case("iso") ! iso refers to solvated myoglobin from biophys j 2009 97, 1158
      nonbond_fcnst = expo_fcnst(1980.0d0*joule2cal,0.749d0,dist)
    case default
      if (chains) then
        nonbond_fcnst = expo_fcnst(6000.0d0*joule2cal,0.896d0,dist)
      else
        nonbond_fcnst = expo_fcnst(2650.0d0*joule2cal,0.750d0,dist)
      end if
  end select

select case (chains) ! intra or inter
case(.true.)
  select case (trim(input%tpbc))
    case("iso") ! iso refers to solvated myoglobin from biophys j 2009 97, 1158
      select case (i_minus_j)
        case(0) ! self
          fcnst = zero
        case(1) ! k12
          fcnst = 735.0d0 *joule2cal
        case(2) ! k13
          fcnst = 5.6d0   *joule2cal
        case(3) ! k14
          fcnst = 31.9d0  *joule2cal
        case default
          fcnst = nonbond_fcnst
      end select  
    case default ! refers to p222 myoglobin from biophys j 2009 97, 1158 
      select case (i_minus_j)
        case(1)
          fcnst = 723.0d0 *joule2cal
        case(2)
          fcnst = 13.8d0  *joule2cal
        case(3)
          fcnst = 29.6d0  *joule2cal
        case default
          fcnst = nonbond_fcnst
      end select
  end select  
case default
  fcnst = nonbond_fcnst
end select

call hessian_block(dist,dxyz,fcnst,hess)

end subroutine

real(dp) function expo_fcnst(a,b,r)
real(dp) , intent(in) :: a,b,r

expo_fcnst = a*dexp(-b*r)

end function

subroutine hessian_block(dist,dxyz,fcnst,hess)
real(dp), intent(in)  :: fcnst,dist,dxyz(3)
real(dp), intent(out) :: hess(3,3)
real(dp)              :: dstsqr, xixj,xiyj,xizj,yiyj,yizj,zizj

dstsqr = dist*dist

xixj=-(fcnst/dstsqr)*(dxyz(1))**2
yiyj=-(fcnst/dstsqr)*(dxyz(2))**2
zizj=-(fcnst/dstsqr)*(dxyz(3))**2
xiyj=-(fcnst/dstsqr)*(dxyz(1))*(dxyz(2))
xizj=-(fcnst/dstsqr)*(dxyz(1))*(dxyz(3))
yizj=-(fcnst/dstsqr)*(dxyz(2))*(dxyz(3))

hess(1,1) = xixj
hess(1,2) = xiyj
hess(1,3) = xizj

hess(2,1) = xiyj
hess(2,2) = yiyj
hess(2,3) = yizj

hess(3,1) = xizj
hess(3,2) = yizj
hess(3,3) = zizj

end subroutine


subroutine sprnghessblock (input, dist, dxyz, hess)
use mod_types,  only: inp_par
! input: input contains requests
! input: distance and dxyz vector between two point
! outpt: 3x3 hessian for the spring interaction
!        with various forceconstants requested in input
type(inp_par)         :: input
real(dp), intent(in)  :: dist,dxyz(3)
real(dp), intent(out) :: hess(3,3)
real(dp)              :: fcnst,dstsqr, xixj,xiyj,xizj,yiyj,yizj,zizj


if (trim(input%fctyp) .eq. "all1") then
  if (dist .le. 4.0d0) then
    fcnst = input%fcnst*1.0d01 !205.d0
  else
    fcnst = input%fcnst*2.0D03/(dist**6)  ! 1968
  end if
else if (trim(input%fctyp) .eq. "all2") then
  if (dist .le. 3.0d0) then
    fcnst = input%fcnst*1050.0d0
  else if (dist .gt. 3.0d0 .and. dist .lt. 9.0d0 ) then
    fcnst = input%fcnst*37.8d02/(dist**2)
  else
    fcnst = input%fcnst*305.9D03/(dist**6)  ! went with hinsen rather then 
  end if
else if (trim(input%fctyp) .eq. "all3") then
  if (dist .le. 3.0d0) then
    fcnst = input%fcnst*250.0d0
  else
    fcnst = input%fcnst*1139.0d0/(dist**6)  ! went with hinsen rather then 
  end if
else if (trim(input%fctyp) .eq. "all4") then
  fcnst = input%fcnst*300.0D00/(dist**3)  ! went with hinsen rather then 

else if (trim(input%fctyp) .eq. "all5") then
  if (dist .le. 3.0d0) then
    fcnst = input%fcnst*37.2d0 !205.d0
  else
    fcnst = input%fcnst*357.4D0/(dist**6)  ! 1968
  end if
else if (trim(input%fctyp) .eq. "all6") then
  if (dist .le. 3.0d0) then
    fcnst = input%fcnst*372.0d0 !205.d0
  else
    fcnst = input%fcnst*357.4D0/(dist**6)  ! 1968
  end if

else if (trim(input%fctyp) .eq. "hinsen") then
  if (dist .le. 4.0d0) then
    fcnst = input%fcnst*(205.5d0*dist - 571.2d0) !hinsen amber fit, calpha
  else 
    fcnst = input%fcnst*(305.9D03/(dist**6)) !hinsen amber fit, calpha
  end if
else ! if (trim(input%fctyp) .eq. "heavy")
  if (input%fcpower .ne. 0) then
  !if (input%fcpower .ne. zero) then
    fcnst = input%fcnst/(dist**input%fcpower)
  else
    fcnst = input%fcnst
  end if
end if


call hessian_block(dist,dxyz,fcnst,hess)

end subroutine

end module

