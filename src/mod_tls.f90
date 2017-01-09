module mod_tls
! DMR 05-19-2009
!
! contains subroutines for reading in TLS tensors and constructing the projection matrix 
! uses the block datatype      
use mod_types,                        only: sparse, protein_atom_list_type
use mod_constants
use mod_bnm
use cfml_crystal_metrics,            only: Crystal_Cell_Type
use cfml_crystallographic_symmetry,  only: space_group_type


implicit none

! found in mod_bnm.f90
!type block
!  integer :: natoms
!  integer, allocatable :: atinblock(:) ! this is contain "pointers" of each atom in the block
!                                   ! useful if blocks contain atoms that are not in order
!  integer, allocatable :: ixyz(:)  ! give the indices corresponding to ats, ie. the location
!                                   ! of block
!  real(dp) :: com(3),rotc(3),cent(3),momI(3),tnsI(3,3),eigI(3,3),tmass
!  real(dp),allocatable :: pmat(:,:)
!  logical :: cent_read
! TLS stuffs
!  real(dp),dimension(:,:),allocatable :: T,L,S
!  real(dp),dimension(:,:),allocatable :: Lv,rho
!  real(dp),dimension(:),allocatable :: Le
!end type

contains

subroutine tls_vcov(input,asym_un,vcov)
use mod_types,  only: inp_par
type(inp_par),            intent(in)     :: input
type(protein_atom_list_type),     intent(in out) :: asym_un 
real(dp),    allocatable, intent(out)    :: vcov(:,:)
type(block), allocatable                 :: blocks(:)

  call tls_block_setup(input,asym_un,blocks)
  call tls_atcovar(blocks,vcov)

end subroutine

subroutine tls_block_setup(input,atoms,blocks)
use mod_types,  only: inp_par
use mod_linalg, only: invrse_mat,lapack_eig3
type(inp_par),           intent(in)    :: input
type(protein_atom_list_type),     intent(inout)    :: atoms
type(block), allocatable,intent(out)   :: blocks(:)
integer :: blck,i,j
real(dp), allocatable :: atvcov(:,:)
real(dp) :: dxyz(3),uij(6),tf

call tls_read(input,atoms,blocks)

call tls_transform(blocks)
do i = 1, size(blocks)
  blocks(i)%cent = blocks(i)%rotc
end do
call tr_proj_schmidt(input,atoms,blocks)
!call tls_atcovar(blocks,atvcov)
!do i = 1, atoms%natoms
!  dxyz = atoms%atom(i)%x-blocks(atoms%atom(i)%iblock)%rotc
!  call calc_uij(dxyz,blocks(atoms%atom(i)%iblock),UIJ)
!  tf = zero
!  do j = 1,3
!    tf = tf+uij(j)
!  end do
!  tf = tf/3.0d0
!  write (55,'(I4,7E12.4)') i,tf,uij 
!end do

call tlsblock_writeout(input,atoms,blocks)

end subroutine

subroutine tls_atcovar(blocks,atvcov)
type(block), allocatable, intent(in)  :: blocks(:)
real(dp),    allocatable, intent(out) :: atvcov(:,:)
real(dp),    allocatable :: tmpmat(:,:)
integer :: i,j,k,ier,natoms

natoms = 0
do i = 1, size(blocks)
  natoms = natoms+blocks(i)%natoms
end do

allocate(atvcov(3*natoms,3*natoms),stat=ier)
if (ier /= 0) stop "tls_atvcovar> malloc"

atvcov = zero

do i = 1, size(blocks)
  call tls_atcovar_block(blocks(i),tmpmat)
  do j = 1, size(blocks(i)%ixyz)
    do k = 1, size(blocks(i)%ixyz) 
      atvcov(blocks(i)%ixyz(j),blocks(i)%ixyz(k)) = tmpmat(j,k)
    end do
  end do
end do

end subroutine

subroutine tls_atcovar_block(ablock,atvcov,iun)
type(block),              intent(in)  :: ablock
real(dp),    allocatable, intent(out) :: atvcov(:,:)
integer,     intent(in), optional     :: iun
real(dp),    allocatable :: tlsmat(:,:),tmpmat(:,:)
integer :: i,ii,jj,j,k,l,ier,iun2
real(dp) :: tf,shit(6,6)


if (allocated(atvcov)) deallocate(atvcov)
allocate(atvcov(3*ablock%natoms,3*ablock%natoms), & 
         tmpmat(6,3*ablock%natoms),tlsmat(6,6), stat=ier)
if(ier /= 0) stop "tls_atcovar> malloc"

atvcov = zero
tmpmat = zero

tlsmat(1:3,1:3) = ablock%T 
tlsmat(4:6,4:6) = ablock%L!*rad_per_deg*rad_per_deg 
tlsmat(1:3,4:6) = ablock%S!*rad_per_deg
tlsmat(4:6,1:3) = transpose(ablock%S)!*rad_per_deg

do i = 1, 6
  do j = 1, 3*ablock%natoms
    do k = 1,6
      tmpmat(i,j) = tmpmat(i,j)+tlsmat(i,k)*ablock%pmat(k,j)
    end do
  end do
end do

do i = 1, 3*ablock%natoms
  do j = 1, 3*ablock%natoms
    do k = 1,6
      atvcov(i,j) = atvcov(i,j)+ablock%pmat(k,i)*tmpmat(k,j)
    end do
  end do
end do

atvcov = atvcov*u2b

!atvcov = matmul(transpose(ablock%pmat),matmul(tlsmat,ablock%pmat))

if (present(iun)) then
!do i = 1, 3*ablock%natoms
 ! do j = 1, 3*ablock%natoms
    !write(iun,'(2I4,F16.4)'), i,j,atvcov(i,j)
 ! end do
!end do
iun2 = iun+100
ii = 1
do i = 1, ablock%natoms
  tf = (atvcov(ii,ii)+atvcov(ii+1,ii+1)+atvcov(ii+2,ii+2))/3.0
  write(iun2,'(I4,7E12.4)'), ablock%atinblock(i),tf, &
                      atvcov(ii,ii),atvcov(ii+1,ii+1),atvcov(ii+2,ii+2),&
                      atvcov(ii,ii+1),atvcov(ii,ii+2),atvcov(ii+1,ii+2)
  ii = ii+3
end do

end if

end subroutine 

subroutine calc_uij(dxyz,ablock,UIJ)
real(dp),    intent(in)  :: dxyz(3)  ! atom to center of reaction
type(block), intent(in)  :: ablock
real(dp),    intent(out) :: uij(6) ! u11,u22,u33,u12,u13,u23
real(dp) :: T(3,3),L(3,3),S(3,3)  
real(dp) :: x,y,z,xx,yy,zz,xy,xz,yz

T = ablock%T
L = ablock%L
S = ablock%S

x = dxyz(1)
y = dxyz(2)
z = dxyz(3)
xy = x*y
xz = x*z
yz = y*z
xx = x*x
yy = y*y
zz = z*z

uij(1) =  L(2,2)*zz + L(3,3)*yy - two*L(2,3)*yz + two*S(2,1)*z - two*S(3,1)*y + T(1,1)
uij(2) =  L(1,1)*zz + L(3,3)*xx - two*L(3,1)*xz - two*S(1,2)*z + two*S(3,2)*x + T(2,2)
uij(3) =  L(1,1)*yy + L(2,2)*xx - two*L(1,2)*xy - two*S(2,3)*x + two*S(1,3)*y + T(3,3)

uij(4) = -L(3,3)*xy + L(2,3)*xz + L(1,3)*yz - L(1,2)*zz  &
         - S(1,1)*z + S(2,2)*z + S(3,1)*x - S(3,2)*y + T(1,2)
uij(5) = -L(2,2)*xz + L(2,3)*xy - L(1,3)*yy + L(1,2)*yz  &
         + S(1,1)*y - S(3,3)*y + S(2,3)*z - S(2,1)*x + T(1,3) 
uij(6) = -L(1,1)*yz - L(2,3)*xx + L(1,3)*xy + L(1,2)*xz  &
         - S(2,2)*x + S(3,3)*x + S(1,2)*y - S(1,3)*z + T(2,3)

uij = uij*u2b

end subroutine

subroutine tls_read(input,atoms,blocks)
use mod_types, only: inp_par
use mod_linalg, only: lapack_eig
use mod_math,   only: determ3
type(inp_par),           intent(in)    :: input
type(protein_atom_list_type),     intent(inout)    :: atoms
type(block), allocatable,intent(out)   :: blocks(:)
character(len=3) :: chartls
character(len=20) :: charints,crap
integer :: j,i,ii,jj,k,nblocks,ier,iblock,inti,intj,natoms
integer :: good_val(3),ngood
real(dp) :: smat(3,3),tmprotc(3),bigP(3,3)
real(dp) :: detval,imat(3,3)
real(dp) :: l1,l2,l3,l12,l13,l23

print *, "tls_read> reading tls information, the block information will be overwritten"
print *, "tls_read> with that found in", trim(input%fileroot) //".tlsin"

open(unit=11,file=trim(input%fileroot) // ".tlsin", &
     status='old',action='read',position='rewind', &
     iostat=ier)

if (ier > 0) then
  print *, "tls_read> cannot tls file"
  print *, "tls_read> stopping program"
  stop
end if

! read: REFMAC $nblocks
!       $\n
  read(11,*,iostat = ier) crap,nblocks

allocate(blocks(nblocks),stat=ier)
if( ier /= 0) stop "tls_read> malloc"

do i = 1, nblocks
  call block_init(blocks(i))
  blocks(i)%cent_read = .true.
  allocate(blocks(i)%T(3,3),blocks(i)%L(3,3),blocks(i)%S(3,3),stat = ier)
  allocate(blocks(i)%Lv(3,3),blocks(i)%rho(3,3),blocks(i)%Le(3),stat = ier)
  if( ier /= 0) stop "tls_read> malloc"
  blocks(i)%T   =zero
  blocks(i)%L   =zero
  blocks(i)%S   =zero
  blocks(i)%Lv  =zero
  blocks(i)%rho =zero
end do
iblock=1
do
  read(11,*,iostat = ier) chartls, charints
  call charints_to_ints(charints,inti,intj)
  ! overwrite block info for atoms, note that since there should be 2 lists of atoms
  ! p1 and nonp1 sym, this information is still out there
  call res_to_atom_block(atoms,inti,intj,iblock,natoms)
  blocks(iblock)%natoms=natoms
  print *, 'tls_read> block:',iblock,'; residues range:', inti, " to ", intj
  if (ier < 0) stop "tls_read> nblocks and tls entries disagree" 
  read(11,*,iostat = ier) crap, charints
  read(11,*,iostat = ier) crap, blocks(iblock)%rotc
  read(11,*,iostat = ier) crap, blocks(iblock)%T(1,1),blocks(iblock)%T(2,2), &
                                blocks(iblock)%T(3,3),blocks(iblock)%T(1,2), &
                                blocks(iblock)%T(1,3),blocks(iblock)%T(2,3)  
  read(11,*,iostat = ier) crap, blocks(iblock)%L(1,1),blocks(iblock)%L(2,2), &
                                blocks(iblock)%L(3,3),blocks(iblock)%L(1,2), &
                                blocks(iblock)%L(1,3),blocks(iblock)%L(2,3)  
! refmac ordering <s22-s11> <s11-s33> s12 s13 s23 s21 s31 s32
! http://www.ccp4.ac.uk/html/refmac5/files/tls.html
!   s11+s22+s33=0 used to determine entries
  read(11,*,iostat = ier) crap, blocks(iblock)%S(2,2),blocks(iblock)%S(1,1), &
                                blocks(iblock)%S(1,2),blocks(iblock)%S(1,3), &
                                blocks(iblock)%S(2,3),blocks(iblock)%S(2,1), &  
                                blocks(iblock)%S(3,1),blocks(iblock)%S(3,2)  
  call matfill(blocks(iblock)%T)
  call matfill(blocks(iblock)%L)
  call sdiag(blocks(iblock)%S(2,2),blocks(iblock)%S(1,1),&
             blocks(iblock)%S(1,1),blocks(iblock)%S(2,2),blocks(iblock)%S(3,3))

  ! convert degrees to radians
  blocks(iblock)%S = blocks(iblock)%S*Rad_Per_Deg
  blocks(iblock)%L = blocks(iblock)%L*Rad_Per_Deg*Rad_Per_Deg

  call lapack_eig(blocks(iblock)%L,blocks(iblock)%Le,blocks(iblock)%Lv)

! Lv right handed?
  detval = determ3(blocks(iblock)%Lv)
  if (detval .lt. zero) then
    imat = zero
    imat(1,1) = -one
    imat(2,2) =  one
    imat(3,3) =  one
    blocks(iblock)%Lv = matmul(blocks(iblock)%Lv,imat)
  end if

 
  print '(A6,1x,3F10.4)', 'ORIGIN', blocks(iblock)%rotc
  print '(A6,1x,3F10.4)', 'leig', blocks(iblock)%Le
  do j = 1,3
    print '(A3,3F12.6)','T',blocks(iblock)%T(j,:)
  end do
  do j = 1,3
    print '(A3,3F12.6)','L',blocks(iblock)%L(j,:)
  end do
  do j = 1,3
    print '(A3,3F12.6)','S',blocks(iblock)%S(j,:)
  end do
  do j = 1,3
   print '(A3,3F12.6)','Lv',blocks(iblock)%Lv(j,:)
  end do
  
  iblock = iblock+1
  if (iblock > nblocks ) exit 
end do

close(11)

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
  do j = 1, atoms%natoms
    if (atoms%atom(j)%iblock .eq. i) then
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

end subroutine

subroutine tls_transform(blocks)

use mod_linalg, only: lapack_eig
use mod_math,   only: determ3
type(block), allocatable,intent(inout)   :: blocks(:)
real(dp) :: smat(3,3),tmprotc(3),bigP(3,3)
real(dp) :: detval,imat(3,3)
real(dp) :: l1,l2,l3,l12,l13,l23
integer :: iblock,i,j,k

print *, 'starting with transform'


do iblock = 1, size(blocks)

  do i = 1, 3
    if (blocks(iblock)%Le(i) .le. 1.0d-08) then
      print *, 'block', iblock, 'is funky; eig:',blocks(iblock)%Le(i)
      print *, 'expect fishiness'
    end if
  end do

! Lv right handed?
  detval = determ3(blocks(iblock)%Lv)
  if (detval .lt. zero) then
    imat = zero
    imat(1,1) = -one
    imat(2,2) =  one
    imat(3,3) =  one
    blocks(iblock)%Lv = matmul(blocks(iblock)%Lv,imat)
  end if

  smat = matmul(transpose(blocks(iblock)%Lv),matmul(blocks(iblock)%S,blocks(iblock)%Lv))

  l1 = blocks(iblock)%Le(1)
  l2 = blocks(iblock)%Le(2)
  l3 = blocks(iblock)%Le(3)
  l12 = l1 + l2
  l13 = l1 + l3
  l23 = l2 + l3

  tmprotc(1) = (smat(2,3)-smat(3,2))/(blocks(iblock)%Le(2)+blocks(iblock)%Le(3))
  tmprotc(2) = (smat(3,1)-smat(1,3))/(blocks(iblock)%Le(3)+blocks(iblock)%Le(1))
  tmprotc(3) = (smat(1,2)-smat(2,1))/(blocks(iblock)%Le(1)+blocks(iblock)%Le(2))
  tmprotc = matmul(blocks(iblock)%Lv,tmprotc) 

  bigP = zero
  bigP(1,2) =  tmprotc(3) 
  bigP(2,1) = -tmprotc(3) 
  bigP(1,3) = -tmprotc(2) 
  bigP(3,1) =  tmprotc(2) 
  bigP(2,3) =  tmprotc(1) 
  bigP(3,2) = -tmprotc(1)  
 
  blocks(iblock)%T = blocks(iblock)%T + matmul(bigP,blocks(iblock)%S) + &
                     matmul(transpose(blocks(iblock)%S),transpose(bigP)) + & 
                     matmul(transpose(bigP),matmul(blocks(iblock)%L,bigP))
  blocks(iblock)%S = blocks(iblock)%S + matmul(blocks(iblock)%L,transpose(bigP))

  smat = matmul(transpose(blocks(iblock)%Lv),matmul(blocks(iblock)%S,blocks(iblock)%Lv))

  blocks(iblock)%rho = zero 

  blocks(iblock)%rho(:,1) = (/ zero, -smat(1,3)/l1, smat(1,2)/l1/)
  blocks(iblock)%rho(:,2) = (/ smat(2,3)/l2, zero, -smat(2,1)/l2/)

  blocks(iblock)%rho(:,3) = (/-smat(3,2)/l3, smat(3,1)/l3, zero/)

  blocks(iblock)%rho =  matmul(blocks(iblock)%Lv,blocks(iblock)%rho) 
  blocks(iblock)%rotc = blocks(iblock)%rotc+tmprotc
 
  print '(A6,1x,3F10.4)', 'ORIGIN', blocks(iblock)%rotc
  print '(A6,1x,3F10.4)', 'cordsp', tmprotc
  print '(A11,1x,3F10.4)', 'leig(L-deg)', deg_per_rad*deg_per_rad*blocks(iblock)%Le
  do j = 1,3
    print '(A13,3F12.6)','T(Angs)',blocks(iblock)%T(j,:)
  end do
  do j = 1,3
    print '(A13,3F12.6)','L(deg^2)',blocks(iblock)%L(j,:)*deg_per_rad*deg_per_rad
  end do
  do j = 1,3
    print '(A13,3F12.6)','S(angs*deg)',blocks(iblock)%S(j,:)*deg_per_rad
 end do
  do j = 1,3
   print '(A3,3F12.6)','Lv',blocks(iblock)%Lv(j,:)
  end do
  do j = 1,3
   print '(A3,1x,3F12.6)','rho',blocks(iblock)%rho(j,:)
  end do
  
end do

print *, 'done with transform'

end subroutine
  
subroutine res_to_atom_block(atoms,inti,intj,iblock,natoms)
! take in residue range and set the blocks for atoms within that residue
! tls lists by residue
type(protein_atom_list_type), intent(inout) :: atoms
integer, intent(in)  :: inti,intj,iblock
integer, intent(out) :: natoms
integer :: i

natoms = 0

do i = 1, atoms%natoms
  if ((atoms%atom(i)%ires .ge. inti ) .and. (atoms%atom(i)%ires .le. intj)) then
    atoms%atom(i)%iblock=iblock
    natoms=natoms+1
  end if
end do

end subroutine

subroutine charints_to_ints(charints,int1,int2)
use cfml_string_utilities,         only: getnum
! takes range of numbers written as a character:
!      1-34  and returns the numbers
character(len=20) , intent(in) :: charints
integer           , intent(out) :: int1,int2
real(sp)  :: crap(10)
integer   :: icrap(10),icrap2
integer :: i,j

i = len_trim(charints)
j = INDEX(charints, "-")
call getnum(charints(1:j-1),crap,icrap,icrap2)
int1 = icrap(1)
call getnum(charints(j+1:i),crap,icrap,icrap2)
int2 = icrap(1)

end subroutine

subroutine matfill(mat)
! takes symmetric mat with upper diag entries and spits out full matrix
real(dp),allocatable, intent(inout) :: mat(:,:)
integer :: i,j

do i = 1, size(mat(:,1))
  do j = i+1, size(mat(1,:))
    mat(j,i) = mat(i,j)
  end do
end do

end subroutine

subroutine sdiag (diff21,diff13,v11,v22,v33)
! using v11+v22+v33=0 and input vals for v22-v11 and v11-v33
real(dp), intent(in)  :: diff21,diff13
real(dp), intent(out) :: v11,v22,v33

v22 = (2.0d0*diff21+diff13)/3.0d0
v11 = v22-diff21
v33 = v11-diff13

end subroutine

subroutine rmatrixz(vec,R)
! converted from mmlib
real(dp), dimension(3), intent(in)   :: vec
real(dp), dimension(3,3), intent(out) :: R
real(dp), dimension(3,3)             :: Rxz,Rxz2z
real(dp) :: u,v,w,d
real(dp), dimension(3) :: nvec
integer :: i

nvec = vec/dsqrt(dot_product(vec,vec))

u = nvec(1)
v = nvec(2)
w = nvec(3)

d = dsqrt(u*u+v*v)

if (d /= zero) then
  Rxz(1,:) = (/  u/d,  v/d, zero/)
  Rxz(2,:) = (/ -v/d,  u/d, zero/)
  Rxz(3,:) = (/ zero, zero, one /)
else
  Rxz = zero
  do i = 1,3
    Rxz(i,i) = one
  end do
end if

Rxz2z(1,:) = (/    w,  zero,   -d/)
Rxz2z(2,:) = (/ zero,   one, zero/)
Rxz2z(3,:) = (/    d,  zero,    w/)

R = matmul(Rxz2z,Rxz)

end subroutine

subroutine cross_product(veca,vecb,vcross)
real(dp), dimension(3), intent(in)   :: veca,vecb
real(dp), dimension(3), intent(out)  :: vcross

vcross(1) = veca(2)*vecb(3)-veca(3)*vecb(2)
vcross(2) = veca(3)*vecb(1)-veca(1)*vecb(3)
vcross(3) = veca(1)*vecb(2)-veca(2)*vecb(1)

end subroutine

subroutine alloc_cross_product(veca,vecb,vcross)
real(dp), dimension(:),allocatable, intent(in)   :: veca,vecb
real(dp), dimension(3), intent(out)  :: vcross

vcross(1) = veca(2)*vecb(3)-veca(3)*vecb(2)
vcross(2) = veca(3)*vecb(1)-veca(1)*vecb(3)
vcross(3) = veca(1)*vecb(2)-veca(2)*vecb(1)

end subroutine

subroutine tlsblock_writeout(input,atoms,blocks)
use mod_types,  only: inp_par
type(inp_par) ,           intent(in) :: input
type(protein_atom_list_type),     intent(in) :: atoms
type(block), allocatable, intent(in) :: blocks(:)
integer :: i,ier

open(unit=11,file=trim(input%fileroot) // ".tlsblock", &
     status='unknown',action='write',iostat=ier)

do i = 1, size(blocks)
  write(11,'(A5,1x,I4,3F10.5)') 'block',i,blocks(i)%rotc
end do

do i = 1, atoms%natoms
  write(11,'(A4,1x,3I5)') trim(atoms%atom(i)%resn), i, atoms%atom(i)%ires,atoms%atom(i)%iblock
end do

close(11)

end subroutine

end module



