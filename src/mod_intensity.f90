module mod_intensity
! DMR note to self... I think there may be additional games to play with atom-atom correlation functions... slightly
!       diff from var-cov mats at they have 1's down the diagonal...  holy shit, I can take the experimental bfactors
!       and multiply them by the atom-atom correlations to get something interesting!!!!!!!!!!  YES!
!
! DMR 11-10-2007
! 
! modified for use with crysFML defined types
! compared with old code and cctbx and it seems to perform well
! checked for 1c1g, 2cba, 2tma it various hkls
!
! takes in protein_atom_list_type in fractional coordinates
!       
!use cfml_Atom_typedef,                only: Atom_List_Type 
use cfml_Scattering_Chemical_Tables, only: xray_form_type,Xray_Form,set_xray_form
use cfml_crystal_metrics,              only: Crystal_Cell_Type
use cfml_String_Utilities,           only: l_case
use cfml_Reflections_Utilities
use mod_constants
use mod_types,                  only: inp_par,protein_atom_list_type

implicit none
integer                           :: hstart,hstop,kstart,kstop,lstart,lstop
integer, parameter                :: natom_entry=214
type(Xray_Form_Type),allocatable  :: smllatlib(:) ! global subset of atomlib5(:) 
real(dp), allocatable             :: attscat(:,:,:,:) ! global lookup table for atomic scatt factors; hkl indices  
integer, allocatable              :: ityp(:)
!DMR: 01-29-2009  added reflection datatypes
!  this way, a one d vector of 3d reflections can be constructed
!  and the intensities or sfactor can be stored for each reflection
type :: reflx
real(dp) :: vec(3)
end type reflx

type :: reflections
integer :: nbragg,noffbragg
type(reflx), allocatable :: qBragg_list(:)
type(reflx), allocatable :: hBragg_list(:)
type(reflx), allocatable :: qvec_list(:) ! list of qvecs, these are used relative to bragg spots 
end type reflections

type :: scatter
integer :: nscatt ! number of reflections
real(dp), allocatable :: scatt(:)
end type scatter

type :: intensity_store
 ! store however many sets of intensities, ie. as a function of qvec
  type(scatter),    allocatable :: qintns(:) ! dim1: reflection, ! dim2 qvec_label
  real(dp),         allocatable :: qvecs(:,:)  ! dim1: qvec_label, ! dim2 3d qvec
  integer :: nqvecs
end type intensity_store

type :: atscatter
! collection of atomic scatters.  Tabulated beforehand to speed things up
!   ntypes is the number of uniqe atom 
integer  :: ntypes ! number of unique atom types, eg. C,N,O,H,S => 5
integer  :: natoms ! number of atoms 
real(dp) :: qvec(3)
type(scatter),    allocatable :: atype(:)
integer ,         allocatable :: itype(:) ! array indexed with atoms that gives the numeric atom type, eg. C=1,N=2 ..
character(len=6), allocatable :: SfacSymb(:) ! store the unique symbols 
end type atscatter

public :: hstart,hstop,kstart,kstop,lstart,lstop,intensity_calc,atscatter

contains

subroutine intensity_calc(input,cell,atoms,vcov,int_thrx,int_thr,int_exp)
use mod_crysbuild,              only: cart_to_fract, fract_to_cart
type(inp_par),           intent(in)   :: input
type (Crystal_Cell_Type),intent(in)   :: Cell
type(protein_atom_list_type)    ,intent(in out)   :: atoms
real(dp),allocatable,    intent(in)   :: vcov(:,:)
real(dp),allocatable,    intent(out)  :: int_thrx(:,:,:),int_thr(:,:,:),int_exp(:,:,:)
integer                     :: reflcts(6)!hstart,hstop,kstart,kstop,lstart,lstop
integer :: i,j,h,k,l !,hstart,hstop,kstart,kstop,lstart,lstop

call cart_to_fract(cell,atoms)
call mod_intensity_setup(input%hihfkikflilf,atoms,cell)
if (trim(input%genm) .eq. "enm") then
  call direct_aniso(input,cell,atoms,vcov,int_thrx,int_thr,int_exp)
else if (trim(input%genm) .eq. "gnm") then
  call direct_iso(input,cell,atoms,vcov,int_thrx,int_thr,int_exp)
end if
call fract_to_cart(cell,atoms)

end subroutine

subroutine mod_intensity_setup(reflcts,atoms,cell) 

integer                 ,intent(in)   :: reflcts(6)
type(protein_atom_list_type)    ,intent(in)   :: atoms
type (Crystal_Cell_Type),intent(in)   :: Cell
real(sp)                              :: iort(3,3)
integer :: i,j,h,k,l !,hstart,hstop,kstart,kstop,lstart,lstop

hstart=reflcts(1)
hstop =reflcts(2)
kstart=reflcts(3)
kstop =reflcts(4)
lstart=reflcts(5)
lstop =reflcts(6)

!set the fract_cart transformation matrix
iort=Cell%orth_cr_cel

call uniqat(atoms,ityp)
! set up lookup table for ssq4 and atomic scattering factors
call atscattsm(iort)


end subroutine

subroutine direct_iso(input,cell,atoms,vcov_in,int_thrx,int_thr,int_exp)
! computes the intensity for the globally defined recip latt indices
! takes in vcov as Angstroms ^2.  
use mod_crysbuild,              only: cart_to_fract, fract_to_cart
type (inp_par),intent(in)             :: input
type (Crystal_Cell_Type),intent(in)   :: Cell
type(protein_atom_list_type),    intent(inout)   :: atoms !fractional cartesian
real(DP), intent(in), allocatable     :: vcov_in(:,:)
real(DP), intent(out), allocatable    :: int_thrx(:,:,:),int_thr(:,:,:),int_exp(:,:,:)
real(DP), allocatable                 :: coor(:,:),occ(:),&
                                         bfact(:), vcov(:,:)
real(DP)                              :: twpih,twpik,twpil,cdmult
real(DP)                              :: dcosterm1
integer                               :: i,j,hh,kk,ll,iat,jat
real(DP)                              :: xi,yi,zi,xj,yj,zj
real(DP)                              :: iatsct,jatsct,iatsct_thr,jatsct_thr,atsct_thrx
real(DP)                              :: h,k,l,ssq4,ssq2
integer                               :: ialloc,natom,ier, eghtpisqr
real(DP)                              :: q(3),qc(3)
real                                  :: time1,time2

natom=atoms%natoms

allocate(vcov(size(vcov_in(1,:)),size(vcov_in(:,1))), stat = ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"
vcov= (8.0d0*pi*pi)*vcov_in

allocate(coor(atoms%natoms,3), &
         occ(atoms%natoms),bfact(atoms%natoms), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"


open(unit=22,file=trim(input%fileroot)//"-intns-"//trim(input%fctyp)//"-"// &
  trim(input%tpbc)//"-"// trim(input%genm)// ".txt", &
     status='unknown',action='write', &
     iostat=ier)

do i=1,natom
  coor(i,:) = atoms%atom(i)%x(:)
  bfact(i)  = atoms%atom(i)%Biso
  occ(i)    = atoms%atom(i)%occ
!  write(99,'(A3,1x,3F12.5,2F6.2)') atoms%atom(i)%ChemSymb,atoms%atom(i)%x(:),occ(i),bfact(i)
end do

if (allocated(int_thr))  deallocate(int_thr)
if (allocated(int_thrx)) deallocate(int_thrx)
if (allocated(int_exp))  deallocate(int_exp)
allocate(int_thr(hstart:hstop,kstart:kstop,lstart:lstop),& 
         int_thrx(hstart:hstop,kstart:kstop,lstart:lstop),&
         int_exp(hstart:hstop,kstart:kstop,lstart:lstop),&
         stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"


int_thr  = zero
int_thrx = zero
int_exp  = zero

call cpu_time(time1)

do ll=lstart,lstop 
  l=real(ll)
  twpil=two*pi*l
  do kk=kstart,kstop 
    k=real(kk)
    twpik=two*pi*k
   do hh=hstart,hstop 
      h=real(hh)
      twpih=two*pi*h
      ssq4=attscat(0,hh,kk,ll)
      ssq2=two*attscat(0,hh,kk,ll)
      q=(/h,k,l/)
      call recart(q,Cell%orth_cr_cel,qc)
! self and cross terms
      do i=1,natom
        ! atom i atomic,occupancy,dwf
        iatsct    =attscat(ityp(i),hh,kk,ll)*occ(i)*dexp(-bfact(i)*ssq4)
        iatsct_thr=attscat(ityp(i),hh,kk,ll)*occ(i)*dexp(-vcov(i,i)*ssq4)
        !iatsct_thr=attscat(ityp(i),hh,kk,ll)*dexp(-vcov(i,i)*ssq4)
        xi=coor(i,1)
        yi=coor(i,2)
        zi=coor(i,3)
        do j=i,natom
          jatsct     = attscat(ityp(j),hh,kk,ll)*occ(j)*dexp(-bfact(j)*ssq4)
          jatsct_thr = attscat(ityp(j),hh,kk,ll)*occ(j)*dexp(-vcov(j,j)*ssq4)
          !print *, 'shit', bfact(j), vcov(j,j)*eghtpisqr 
        !jatsct_thr = attscat(ityp(j),hh,kk,ll)*dexp(-vcov(j,j)*ssq4)
          atsct_thrx = dexp(vcov(i,j)*ssq2)
          if (i.eq.j) then
            cdmult = one
          else
            cdmult=two
          end if
          xj=coor(j,1)
          yj=coor(j,2)
          zj=coor(j,3)
         
!              dexp(-0.5*( twpih2*(vcovar(xi,xi)+vcovar(xj,xj)) + &
!                          twpik2*(vcovar(yi,yi)+vcovar(yj,yj)) + &
!                          2.0*twpihk*(vcovar(xi,yi) + vcovar(xj,yj) ) ) + &
!                    xmult * ( twpih2*vcovar(xi,xj) + twpihk*vcovar(xi,yj) + &
!                         twpihk*vcovar(yi,xj) +twpik2*vcovar(yi,yj) ) &
!                   )
 
          dcosterm1 = cdmult * &
                      dcos(twpih * (xi-xj) + twpik * (yi - yj)+twpil * (zi-zj))
!          dcosterm1=cdmult*iatsct*jatsct*dexp(-ssq4*(bfact(j)+bfact(i))*(cdmult-one)) * &
!                  dcos(twpih*(xi-xj)+twpik*(yi - yj)+twpil*(zi-zj))

          int_thr(hh,kk,ll)  = int_thr(hh,kk,ll) + dcosterm1*iatsct_thr*jatsct_thr 
          int_thrx(hh,kk,ll) = int_thrx(hh,kk,ll) + dcosterm1*iatsct_thr*jatsct_thr*atsct_thrx
          int_exp(hh,kk,ll)  = int_exp(hh,kk,ll) + dcosterm1*iatsct*jatsct
      

        end do
      end do
          write (22,'(3I4,3F10.4,4F16.3)') hh,kk,ll,qc,int_exp(hh,kk,ll), &
            int_thr(hh,kk,ll),int_thrx(hh,kk,ll), int_thrx(hh,kk,ll)-int_thr(hh,kk,ll)
    end do
  end do
end do
call cpu_time(time2)
print *, 'time in direct sum', time2-time1, 'seconds'

end subroutine 

subroutine direct_aniso(input,cell,atoms,vcov_in,int_thrx,int_thr,int_exp)
! computes the intensity for the globally defined recip latt indices
! DMR: I think that can be sped up.  The attscat lib may be slowing it down
!      due to the dimensionality.
type(inp_par),           intent(in)   :: input
type (Crystal_Cell_Type),intent(in)   :: Cell
type(protein_atom_list_type),    intent(in)   :: atoms
real(DP), intent(in), allocatable     :: vcov_in(:,:)
real(DP), intent(out), allocatable    :: int_thrx(:,:,:),int_thr(:,:,:),int_exp(:,:,:)
real(DP), allocatable                 :: coor(:,:),occ(:),bfact(:),bfact_thr(:)
real(dp), allocatable                 :: vcov(:,:)
real(DP)                              :: twpih,twpik,twpil,cdmult
real(DP), allocatable                 :: int_thrbf(:,:,:)
real(DP)                              :: dcosterm1
integer                               :: i,j,hh,kk,ll,iat,jat
real(DP)                              :: ai,bi,ci,aj,bj,cj
real(DP)                              :: iatsct,jatsct,iatsct_thr,jatsct_thr,atsct_thrx
real(DP)                              :: jatsct_thrbf,ithr_bf,ijthr_bf
real(DP)                              :: h,k,l,ssq4,ssq2,iexp_exp,ijexp_exp,iexp_thr,ijexp_thr
real(dp)                              :: q(3),qc(3),qc_dot_qc,ddot
integer                               :: ialloc,ier,natom,icor,jcor,incr
real                                  :: time1,time2
external :: ddot
!DMR this will be 1 or 3 depending on iso or aniso
incr=3

allocate(vcov(size(vcov_in(1,:)),size(vcov_in(:,1))), stat = ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"
vcov= (8.0d0*pi*pi)*vcov_in
allocate (coor(atoms%natoms,3),occ(atoms%natoms),bfact(atoms%natoms),bfact_thr(atoms%natoms), stat = ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"


if(trim(input%tpbc).eq."bvk") then
  open(unit=22,file=trim(input%fileroot)//"-intns-"//trim(input%fctyp)//"-"// &
    trim(input%tpbc)//"-"// trim(input%genm)//"-"//trim(input%bztyp)// ".txt", &
     status='unknown',action='write', &
     iostat=ier)
else
  open(unit=22,file=trim(input%fileroot)//"-intns-"//trim(input%fctyp)//"-"// &
    trim(input%tpbc)//"-"// trim(input%genm)// ".txt", &
     status='unknown',action='write', &
     iostat=ier)
end if


natom=atoms%natoms
j = 1
do i=1,natom
  coor(i,:) = atoms%atom(i)%x(:)
  bfact(i)  = atoms%atom(i)%Biso
  bfact_thr(i)  = (vcov(j,j)+vcov(j+1,j+1)+vcov(j+2,j+2))/3.0d0
  occ(i)    = atoms%atom(i)%occ
!  write(99,'(A3,1x,3F12.5,2F6.2)') atoms%atom(i)%ChemSymb,atoms%atom(i)%x(:),occ(i),bfact(i)
  j = j+3
end do

if (allocated(int_thr))  deallocate(int_thr)
if (allocated(int_thrx)) deallocate(int_thrx)
if (allocated(int_exp))  deallocate(int_exp)
allocate(int_thr(hstart:hstop,kstart:kstop,lstart:lstop),& 
         int_thrx(hstart:hstop,kstart:kstop,lstart:lstop),&
         int_exp(hstart:hstop,kstart:kstop,lstart:lstop),&
         int_thrbf(hstart:hstop,kstart:kstop,lstart:lstop),&
         stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"


int_thr  = zero
int_thrbf  = zero
int_thrx = zero
int_exp  = zero

call cpu_time(time1)

do ll=lstart,lstop 
  l=real(ll,kind=dp)
  twpil=two*pi*l
  do kk=kstart,kstop 
    k=real(kk,kind=dp)
    twpik=two*pi*k
   do hh=hstart,hstop 
      h=real(hh,kind=dp)
      twpih=two*pi*h
      !ssq4=attscat(0,hh,kk,ll)
      !ssq2=two*attscat(0,hh,kk,ll)
      q=(/h,k,l/)
      call recart(q,Cell%orth_cr_cel,qc)
      qc_dot_qc = (qc(1)*qc(1)+qc(2)*qc(2)+qc(3)*qc(3) )/4.0d0
      !qc_dot_qc = dot_product(qc,qc)
! self and cross terms
      icor=1
      do i=1,natom
        ! atom i atomic,occupancy,dwf
!        iexp_exp=-bfact(i)*dot_product(qc,qc)/4.0d0
        iexp_exp=-bfact(i)*qc_dot_qc
        ithr_bf = -bfact_thr(i)*qc_dot_qc
        !iexp_thr=-vcov(i,i)*dot_product(qc,qc)/4.0d0

        iexp_thr=-(qc(1)*dot_product(qc,vcov(icor:icor+2,icor)) + &
                   qc(2)*dot_product(qc,vcov(icor:icor+2,icor+1)) + &
                   qc(3)*dot_product(qc,vcov(icor:icor+2,icor+2)))/4.0d0 

    !    iexp_thr=-(qc(1)*( qc(1)*vcov(icor  ,icor) +     &
    !                       qc(2)*vcov(icor+1,icor) +     &
    !                       qc(3)*vcov(icor+2,icor) ) +   &
    !               qc(2)*( qc(1)*vcov(icor  ,icor+1) +   &
    !                       qc(2)*vcov(icor+1,icor+1) +   &
    !                       qc(3)*vcov(icor+2,icor+1) ) + &
    !               qc(3)*( qc(1)*vcov(icor  ,icor+2) +   &
    !                       qc(2)*vcov(icor+1,icor+2) +   &
    !                       qc(3)*vcov(icor+2,icor+2) ) )/4.0d0
        
        iatsct    =attscat(ityp(i),hh,kk,ll)*occ(i)
        ai=coor(i,1)
        bi=coor(i,2)
        ci=coor(i,3)
        jcor=icor
        do j=i,natom
         ! ijexp_exp = iexp_exp-bfact(j)*dot_product(qc,qc)/4.0d0
          ijexp_exp = iexp_exp-bfact(j)*qc_dot_qc
          ijthr_bf = ithr_bf-bfact_thr(j)*qc_dot_qc
          ijexp_thr=iexp_thr-(qc(1)*dot_product(qc,vcov(jcor:jcor+2,jcor)) + &
                              qc(2)*dot_product(qc,vcov(jcor:jcor+2,jcor+1)) + &
                              qc(3)*dot_product(qc,vcov(jcor:jcor+2,jcor+2)))/4.0d0 
        
      !   ijexp_thr=iexp_thr - (qc(1)*( qc(1)*vcov(jcor  ,jcor) +     &
      !                     qc(2)*vcov(jcor+1,jcor) +     &
      !                     qc(3)*vcov(jcor+2,jcor) ) +   &
      !             qc(2)*( qc(1)*vcov(jcor  ,jcor+1) +   &
      !                     qc(2)*vcov(jcor+1,jcor+1) +   &
      !                     qc(3)*vcov(jcor+2,jcor+1) ) + &
      !             qc(3)*( qc(1)*vcov(jcor  ,jcor+2) +   &
      !                     qc(2)*vcov(jcor+1,jcor+2) +   &
      !                     qc(3)*vcov(jcor+2,jcor+2) ) )/4.0d0

          jatsct     = attscat(ityp(j),hh,kk,ll)*occ(j)*dexp(ijexp_exp)
          jatsct_thrbf     = attscat(ityp(j),hh,kk,ll)*occ(j)*dexp(ijthr_bf)
          jatsct_thr = attscat(ityp(j),hh,kk,ll)*occ(j)*dexp(ijexp_thr)

          atsct_thrx = (qc(1)*dot_product(qc,vcov(icor:icor+2,jcor)) + &
                        qc(2)*dot_product(qc,vcov(icor:icor+2,jcor+1)) + &
                        qc(3)*dot_product(qc,vcov(icor:icor+2,jcor+2)))/2.0d0
          atsct_thrx = dexp(atsct_thrx)

          if (i.eq.j) then
            cdmult = one
          else
            cdmult=two
          end if
          aj=coor(j,1)
          bj=coor(j,2)
          cj=coor(j,3)
 
          dcosterm1=cdmult* &
                  dcos(twpih*(ai-aj)+twpik*(bi - bj)+twpil*(ci-cj))

          int_thr(hh,kk,ll)  = int_thr(hh,kk,ll)  + dcosterm1*iatsct*jatsct_thr 
          int_thrx(hh,kk,ll) = int_thrx(hh,kk,ll) + dcosterm1*iatsct*jatsct_thr*atsct_thrx
          int_exp(hh,kk,ll)  = int_exp(hh,kk,ll)  + dcosterm1*iatsct*jatsct
          int_thrbf(hh,kk,ll)  = int_thrbf(hh,kk,ll)  + dcosterm1*iatsct*jatsct_thrbf

          jcor=jcor+incr
        end do
        icor=icor+incr
      end do
      write (22,'(3I4,3F10.4,5F16.3)') hh,kk,ll,qc,int_exp(hh,kk,ll),int_thrbf(hh,kk,ll), &
            int_thr(hh,kk,ll),int_thrx(hh,kk,ll), int_thrx(hh,kk,ll)-int_thr(hh,kk,ll)

    end do
  end do
end do
call cpu_time(time2)
print *, 'time in direct sum', time2-time1, 'seconds'

close (22)

end subroutine 

subroutine zdirect_aniso(input,cell,atoms,zvcov,int_thrx,int_thr,int_exp)
! computes the intensity for the globally defined recip latt indices
! ! old subroutine kept for future reference
! uses complex vcov matrix..  vcov is not complex if the BZ is sampled correctly
type(inp_par),           intent(in)   :: input
type (Crystal_Cell_Type),intent(in)   :: Cell
type(protein_atom_list_type),    intent(in)   :: atoms
complex(kind=8), intent(in), allocatable  :: zvcov(:,:)
real(DP), intent(out), allocatable    :: int_thrx(:,:,:),int_thr(:,:,:),int_exp(:,:,:)
real(DP)                              :: coor(atoms%natoms,3),occ(atoms%natoms),bfact(atoms%natoms)
real(DP)                              :: twpih,twpik,twpil,cdmult
real(DP)                              :: dcosterm1,dcostermx,dcostermexp
real(DP)                              :: dsinterm1
integer                               :: i,j,hh,kk,ll,iat,jat
real(DP)                              :: ai,bi,ci,aj,bj,cj
real(DP)                              :: iatsct,jatsct,iatsct_thr,jatsct_thr,jatsct_thrx,atsct_thrx
real(DP)                              :: cos_thr,sin_thr
real(DP)                              :: cos_thrx,sin_thrx
real(DP)                              :: h,k,l,ssq4,ssq2,iexp_exp,ijexp_exp,iexp_thr,ijexp_thr
real(DP)                              :: real_ithr,real_ijthr,cmpx_ithr,cmpx_ijthr
real(DP)                              :: real_thrx,cmpx_thrx
real(dp)                              :: q(3),qc(3),qc_dot_qc,ddot
integer                               :: ialloc,ier,natom,icor,jcor,incr
real                                  :: time1,time2
external :: ddot
!DMR this will be 1 or 3 depending on iso or aniso
incr=3

print *, 'old subroutine zdirect_aniso called, stopping program'
stop

open(unit=22,file=trim(input%fileroot)//"-zintns-"//trim(input%fctyp)//"-"// &
  trim(input%tpbc)//"-"// trim(input%genm)// ".txt", &
     status='unknown',action='write', &
     iostat=ier)

natom=atoms%natoms
do i=1,natom
  coor(i,:) = atoms%atom(i)%x(:)
  bfact(i)  = atoms%atom(i)%Biso
  occ(i)    = atoms%atom(i)%occ
!  write(99,'(A3,1x,3F12.5,2F6.2)') atoms%atom(i)%ChemSymb,atoms%atom(i)%x(:),occ(i),bfact(i)
end do

if (allocated(int_thr))  deallocate(int_thr)
if (allocated(int_thrx)) deallocate(int_thrx)
if (allocated(int_exp))  deallocate(int_exp)
allocate(int_thr(hstart:hstop,kstart:kstop,lstart:lstop),& 
         int_thrx(hstart:hstop,kstart:kstop,lstart:lstop),&
         int_exp(hstart:hstop,kstart:kstop,lstart:lstop),&
         stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"


int_thr  = zero
int_thrx = zero
int_exp  = zero

call cpu_time(time1)

do ll=lstart,lstop 
  l=real(ll,kind=dp)
  twpil=two*pi*l
  do kk=kstart,kstop 
    k=real(kk,kind=dp)
    twpik=two*pi*k
   do hh=hstart,hstop 
      h=real(hh,kind=dp)
      twpih=two*pi*h
      ssq4=attscat(0,hh,kk,ll)
      ssq2=two*attscat(0,hh,kk,ll)
      q=(/h,k,l/)
      call recart(q,Cell%orth_cr_cel,qc)
      qc_dot_qc = (qc(1)*qc(1)+qc(2)*qc(2)+qc(3)*qc(3) )/4.0d0
      !qc_dot_qc = dot_product(qc,qc)
! self and cross terms
      icor=1
      do i=1,natom
! DEBYE-WALLER Terms
        iexp_exp=-bfact(i)*qc_dot_qc       ! experimental
        ! next theoretical: real and imaginary components
        real_ithr=-(qc(1)*dot_product(qc,dble(zvcov(icor:icor+2,icor))) + &
                   qc(2)*dot_product(qc,dble(zvcov(icor:icor+2,icor+1) )) + &
                   qc(3)*dot_product(qc,dble(zvcov(icor:icor+2,icor+2))))/4.0d0 
    ! complex part of the diagonal yields zero
    !    cmpx_ithr=-(qc(1)*dot_product(qc,dimag(zvcov(icor:icor+2,icor))) + &
    !                qc(2)*dot_product(qc,dimag(zvcov(icor:icor+2,icor+1))) + &
    !                qc(3)*dot_product(qc,dimag(zvcov(icor:icor+2,icor+2))))/4.0d0 

       ! atomic scattering factor and occupancy for atom i 
        iatsct    =attscat(ityp(i),hh,kk,ll)*occ(i)
        ai=coor(i,1)
        bi=coor(i,2)
        ci=coor(i,3)
        jcor=icor

        do j=i,natom
! DEBYE-WALLER Terms
          ijexp_exp = iexp_exp - bfact(j)*qc_dot_qc
        ! next theoretical: real and imaginary components
          real_ijthr=real_ithr - (qc(1)*dot_product(qc,dble(zvcov(jcor:jcor+2,jcor))) + &
                                  qc(2)*dot_product(qc,dble(zvcov(jcor:jcor+2,jcor+1))) + &
                                  qc(3)*dot_product(qc,dble(zvcov(jcor:jcor+2,jcor+2))))/4.0d0 
    ! complex part of the diagonal yields zero
          !cmpx_ijthr=cmpx_ithr - (qc(1)*dot_product(qc,dimag(zvcov(jcor:jcor+2,jcor))) + &
          !                        qc(2)*dot_product(qc,dimag(zvcov(jcor:jcor+2,jcor+1))) + &
          !                        qc(3)*dot_product(qc,dimag(zvcov(jcor:jcor+2,jcor+2))))/4.0d0 
        
        ! next theoretical: real and imaginary cross terms
          real_thrx = (qc(1)*dot_product(qc,dble(zvcov(icor:icor+2,jcor))) + &
                        qc(2)*dot_product(qc,dble(zvcov(icor:icor+2,jcor+1))) + &
                        qc(3)*dot_product(qc,dble(zvcov(icor:icor+2,jcor+2))))/2.0d0
    ! complex off diag part may yield nonzero vals 
          cmpx_thrx = (qc(1)*dot_product(qc,dimag(zvcov(icor:icor+2,jcor))) + &
                        qc(2)*dot_product(qc,dimag(zvcov(icor:icor+2,jcor+1))) + &
                        qc(3)*dot_product(qc,dimag(zvcov(icor:icor+2,jcor+2))))/2.0d0
       ! atomic scattering factor and occupancy for atom i 
          jatsct     = attscat(ityp(j),hh,kk,ll)*occ(j)*dexp(ijexp_exp)
          jatsct_thr = attscat(ityp(j),hh,kk,ll)*occ(j)*dexp(real_ijthr)
          jatsct_thrx = attscat(ityp(j),hh,kk,ll)*occ(j)*dexp(real_ijthr+real_thrx)

          cos_thrx = dcos(cmpx_ijthr+cmpx_thrx)
          sin_thrx = dsin(cmpx_ijthr+cmpx_thrx)

          if (i.eq.j) then
            cdmult = one
          else
            cdmult=two
          end if
          aj=coor(j,1)
          bj=coor(j,2)
          cj=coor(j,3)
 
          dcosterm1=cdmult* &
                  dcos(twpih*(ai-aj)+twpik*(bi - bj)+twpil*(ci-cj) )
          
          dsinterm1=cdmult* &
                  dsin(twpih*(ai-aj)+twpik*(bi - bj)+twpil*(ci-cj))

          int_thr(hh,kk,ll)  = int_thr(hh,kk,ll)  + iatsct*jatsct_thr*dcosterm1
          int_thrx(hh,kk,ll)  = int_thrx(hh,kk,ll)  + iatsct*jatsct_thrx*(dcosterm1*cos_thrx - &
                                                                       dsinterm1*sin_thrx )
          int_exp(hh,kk,ll)  = int_exp(hh,kk,ll)  + dcosterm1*iatsct*jatsct
          jcor=jcor+incr
        end do
        icor=icor+incr
      end do
      write (22,'(3I4,3F10.4,4F16.3)') hh,kk,ll,qc,int_exp(hh,kk,ll), &
            int_thr(hh,kk,ll),int_thrx(hh,kk,ll),int_thrx(hh,kk,ll)-int_thr(hh,kk,ll)

    end do
  end do
end do
call cpu_time(time2)
print *, 'time in direct sum', time2-time1, 'seconds'

close (22)

end subroutine 

subroutine uniqat(atoms,ityp)
! uses global smllatlib for now
! 11-10-2007: works
type(protein_atom_list_type),intent(in)       :: atoms
integer, allocatable,    intent(out)  :: ityp(:)
Type(Xray_Form_Type),    allocatable  :: tmpatlib(:)
character(5)                          :: tmpat1,tmpat2
integer                               :: ialloc, i,ii, j,smnum=0,natom 
real                                  :: time1,time2

natom = atoms%natoms

! load up the atomic scattering factors
call set_xray_form()
if (allocated(tmpatlib))  deallocate(tmpatlib)
if (allocated(ityp))  deallocate(ityp)
allocate(tmpatlib(natom_entry),ityp(natom),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

js: do j=1,natom
  tmpat1=trim(l_case(atoms%atom(j)%SfacSymb))
  do i=1,smnum
    tmpat2=trim(l_case(tmpatlib(i)%symb))
    if (tmpat1.eq.tmpat2) then
      ityp(j)=i
      cycle js
    end if
  end do 

  smnum=smnum+1
  ityp(j)=smnum

  do ii=1,natom_entry
    tmpat2=trim(l_case(xray_form(ii)%symb))
    if (tmpat1 .eq. tmpat2) then
      tmpatlib(smnum)=xray_form(ii)
    end if
  end do
end do js

if (allocated(smllatlib)) deallocate(smllatlib)
allocate(smllatlib(smnum),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"
 

do i=1,smnum
  smllatlib(i)=tmpatlib(i)
end do

deallocate(tmpatlib)
deallocate(xray_form)

end subroutine

subroutine atscattsm(iort)
! load up the atscatt lookup table
real(sp),intent(in)           :: iort(3,3)
integer                       :: nat,ialloc
real(sp)                      :: q(3),a(4),b(4),c,ssq4
integer                       :: i,hh,kk,ll
real(sp)                      :: h,k,l
real                          :: time1,time2

call cpu_time(time1)

nat=size(smllatlib)

if (allocated(attscat))  deallocate(attscat)
allocate(attscat(0:nat,hstart:hstop,kstart:kstop,lstart:lstop),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

! lookup entry for the dstar_sqr
do ll=lstart,lstop
  l=real(ll)
  do kk=kstart,kstop
    k=real(kk)
    do hh=hstart,hstop
      h=real(hh)
      q(1)=h ; q(2)=k; q(3)=l
      attscat(0,hh,kk,ll)=dstar_sqr(q,iort)/four
    end do
  end do
end do

do i=1,nat
  a=smllatlib(i)%a
  b=smllatlib(i)%b
  c=smllatlib(i)%c
  do ll=lstart,lstop 
    l=real(ll)
    do kk=kstart,kstop
      k=real(kk)
      do hh=hstart,hstop
        h=real(hh)
        q(1)=h ; q(2)=k; q(3)=l
        ssq4=attscat(0,hh,kk,ll)
        attscat(i,hh,kk,ll)=fat4(a,b,c,ssq4)
      end do
    end do
  end do
end do
call cpu_time(time2)
print *, 'time to set up look up table', time2-time1, 'seconds'

end subroutine 

real(sp) function fat4(a,b,c,ssq4)
real(sp), intent(in) :: a(4),b(4),c,ssq4
real(sp) :: qsqr

fat4   = a(1)*exp(-b(1)*ssq4)+ &
         a(2)*exp(-b(2)*ssq4)+ &
         a(3)*exp(-b(3)*ssq4)+ &
         a(4)*exp(-b(4)*ssq4)+ &
         c

end function


real(DP) function fat5(a,b,c,ssq4)
real(DP), intent(in) :: a(5),b(5),c,ssq4
real(DP) :: qsqr

fat5   = a(1)*exp(-b(1)*ssq4)+ &
         a(2)*exp(-b(2)*ssq4)+ &
         a(3)*exp(-b(3)*ssq4)+ &
         a(4)*exp(-b(4)*ssq4)+ &
         a(5)*exp(-b(5)*ssq4)+ &
         c

end function

!GET RESOLUTION OF A REFLECTION
!(REQUIRES MATRIX FROM RECIPORT)
!RETURNS VALUE DSTAR IN RECIPROCAL ANGSTROMS
!
!FROM PAGE 19 IN PRINCE
!
real(SP) function dstar(v,iort)
  real(SP),intent(in) :: v(3),iort(3,3)
  real(SP)            :: a(3),tmp

  a=matmul(v,iort)
  tmp=dot_product(a,a)
  dstar=sqrt(tmp)

end function

! for the dstar squared
real(SP) function dstar_sqr(v,iort)
  real(SP),intent(in) :: v(3),iort(3,3)
  real(SP)            :: a(3)

  a=matmul(v,iort)
  dstar_sqr=dot_product(a,a)

end function

subroutine recart(v,iort,a)
  real(DP),intent(in)  :: v(3)
  real(SP),intent(in)  :: iort(3,3)
  real(DP),intent(out) :: a(3)

  a=matmul(v,iort)

end subroutine


! below this is new as of 012009
subroutine reflection_setup (input,cell,reflects)
type(inp_par),     intent(in)  :: input
type(reflections), intent(out) :: reflects
type (Crystal_Cell_Type), intent(in)        :: Cell
integer :: h,k,l,hi,hf,ki,kf,li,lf,nrflx,nh,nk,nl
integer :: ih,ik,il
integer :: ialloc,ndh,ndk,ndl
real(dp) :: vec(3),dh,dk,dl,rh,rk,rl


hi=input%hihfkikflilf(1)
hf=input%hihfkikflilf(2)
dh = input%dh
ndh = npoints(dh)
ki=input%hihfkikflilf(3)
kf=input%hihfkikflilf(4)
dk = input%dk
ndk = npoints(dk)
li=input%hihfkikflilf(5)
lf=input%hihfkikflilf(6)
dl = input%dl
ndl = npoints(dl)

! determine number of reflections
h = total_points(hi,hf,dh)
k = total_points(ki,kf,dk)
l = total_points(li,lf,dl)
!h = total_integers(hi,hf)
!k = total_integers(ki,kf)
!l = total_integers(li,lf)

nrflx = h*k*l

print *, "number of reciprocal space points:",nrflx
print *, "hdivs:",ndh,"total h:", h
print *, "kdivs:",ndk,"total k:", k
print *, "ldivs:",ndl,"total l:", l

reflects%nbragg = nrflx
allocate(reflects%hBragg_list(reflects%nbragg),reflects%qBragg_list(reflects%nbragg),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

nrflx = 0
do l = li,lf
do il = 0,ndl-1
  rl = real(l,kind=8) + dl*real(il,kind=8)
  if(rl > real(lf,kind=8)) cycle

  do k = ki,kf
  do ik = 0,ndk-1
    rk = real(k,kind=8) + dk*real(ik,kind=8)
    if(rk > real(kf,kind=8)) cycle

    do h = hi,hf
    do ih = 0,ndh-1
      rh = real(h,kind=8) + dh*real(ih,kind=8)
      if(rh > real(hf,kind=8)) cycle
      nrflx = nrflx +1
      vec = (/rh,rk,rl/) 
      !vec = (/real(h,kind=8),real(k,kind=8),real(l,kind=8)/)
      reflects%hBragg_list(nrflx)%vec = two * pi * vec
      reflects%qBragg_list(nrflx)%vec = two * pi * matmul(vec,Cell%orth_cr_cel)
!      write(777,'(I5,3F10.4)') nrflx, reflects%hBragg_list(nrflx)%vec/(two*pi)
    end do
    end do

  end do
  end do

end do
end do

end subroutine

subroutine reflection_setup_resolution (input,res_angs,cell,reflects)
type(inp_par),            intent(in)   :: input
real(dp),                 intent(in)   :: res_angs
type (Crystal_Cell_Type), intent(in)   :: Cell
type(reflections),        intent(out)  :: reflects
integer :: ialloc,ncnt,i
real(sp) :: vec(3),rh,rk,rl,rki,rli,testres
real(dp) :: maxh,maxk,maxl,res
real(dp), allocatable :: tmp(:,:)

res = one/res_angs

maxh = anint(maxvalres(res,cell,"h"))
maxk = anint(maxvalres(res,cell,"k"))
maxl = anint(maxvalres(res,cell,"l"))
print *, 'resolution', res, 'in angs',res_angs
print *, 'maxh', maxh/input%dh
print *, 'maxk', maxk/input%dk
print *, 'maxl', maxl/input%dl

rh = -maxh
rk = -maxk
rl = -maxl

! we only need half of the intensities we'll go from zero - max of biggest
!       and then -val to val for the other 2
if (maxh .gt. maxk .and. maxh .gt. maxl) rh = zero
if (maxk .gt. maxh .and. maxk .gt. maxl) rk = zero
if (maxl .gt. maxh .and. maxl .gt. maxk) rl = zero

rki = rk
rli = rl

reflects%nbragg = nint((maxh-rh +one)*(maxk-rk+one)*(maxl-rl+one))
allocate(tmp(reflects%nbragg,3),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

!print *, "max count",reflects%nbragg,input%dh,input%dk,input%dl
ncnt = 1
do while (rh .le. maxh)
  rk = rki
  do while (rk .le. maxk)
    rl = rli
    do while (rl .le. maxl)
      vec = (/rh,rk,rl/) 
      testres = dstar(vec,Cell%orth_cr_cel) 
      if (testres .le. res) then
        tmp(ncnt,:) = vec
        ncnt = ncnt+1
      end if
      rl = rl + input%dl
    end do
    rk = rk + input%dk
  end do
  rh = rh + input%dh
end do

reflects%nbragg = ncnt
allocate(reflects%hBragg_list(reflects%nbragg),reflects%qBragg_list(reflects%nbragg),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

do i = 1,ncnt
   vec = tmp(i,:) 
   reflects%hBragg_list(i)%vec = two * pi * vec
   reflects%qBragg_list(i)%vec = two * pi * matmul(vec,Cell%orth_cr_cel)
end do

end subroutine

subroutine unique_reflect_crysfml (res_angs,spg,cell,reflects)
use cfml_crystallographic_symmetry,  only: space_group_type
real(dp),                          intent(in)   :: res_angs
type (Space_Group_Type) ,          intent(in)   :: spg
type (Crystal_Cell_Type),          intent(in)   :: Cell
type(reflections),        intent(out)  ::          reflects
type (Reflect_Type)                             :: reflex(1000000)
integer                                         :: num_ref,ialloc,i
real(sp) :: res
real(dp) :: vec(3)

! this is how we do it! range of sin(th)/lam  = 1/2 \AA-1
res = one/(2.0d0*res_angs)
call Hkl_Gen(cell,Spg,.true.,0.001_sp,res,Num_Ref,Reflex)

reflects%nbragg = num_ref
allocate(reflects%hBragg_list(reflects%nbragg), &
         reflects%qBragg_list(reflects%nbragg),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

do i = 1,num_ref
   vec = real(reflex(i)%h, kind =dp) 
   reflects%hBragg_list(i)%vec = two * pi * vec
   reflects%qBragg_list(i)%vec = two * pi * matmul(vec,Cell%orth_cr_cel)
end do

end subroutine

real(dp) function maxvalres(res,cell,hkl)
real(dp),                 intent(in)   :: res
type (Crystal_Cell_Type), intent(in)   :: Cell
character(len=1)        , intent(in)   :: hkl
real(dp)                :: vec(3),qvec(3),tmp
real(dp) :: vrsres

vec = zero
if(hkl .eq. "h") vec(1) = one
if(hkl .eq. "k") vec(2) = one
if(hkl .eq. "l") vec(3) = one

qvec = matmul(vec,Cell%orth_cr_cel)
! this is inv \AA of a bragg peak
tmp  = dsqrt(dot_product(qvec,qvec))
maxvalres = res/tmp

end function

subroutine uniqat_new(atoms,smllatlib_new,atomic_scatt)
! DMR 01-21-2009
type(protein_atom_list_type),intent(in)       :: atoms
type(Xray_Form_Type),allocatable, intent(out)  :: smllatlib_new(:) ! global subset of atomlib5(:) 
type(atscatter),intent(out)           :: atomic_scatt 
Type(Xray_Form_Type),    allocatable  :: tmpatlib(:)
character(5)                          :: tmpat1,tmpat2
integer                               :: ialloc, i,ii, j,smnum=0,natom 
real                                  :: time1,time2
integer, allocatable  :: ityp(:)

natom = atoms%natoms
atomic_scatt%natoms = natom
! load up the atomic scattering factors
call set_xray_form()
if (allocated(tmpatlib))  deallocate(tmpatlib)
if (allocated(ityp))  deallocate(ityp)
allocate(tmpatlib(natom_entry),atomic_scatt%itype(natom),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

js: do j=1,natom
  tmpat1=trim(l_case(atoms%atom(j)%SfacSymb))
  do i=1,smnum
    tmpat2=trim(l_case(tmpatlib(i)%symb))
    if (tmpat1.eq.tmpat2) then
      atomic_scatt%itype(j)=i
      cycle js
    end if
  end do 

  smnum=smnum+1
  atomic_scatt%itype(j)=smnum

  do ii=1,natom_entry
    tmpat2=trim(l_case(xray_form(ii)%symb))
    if (tmpat1 .eq. tmpat2) then
      tmpatlib(smnum)=xray_form(ii)
    end if
  end do
end do js

if (allocated(smllatlib_new)) deallocate(smllatlib_new)
allocate(smllatlib_new(smnum),atomic_scatt%SfacSymb(smnum), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

do i = 1, atoms%natoms
  atomic_scatt%SfacSymb(atomic_scatt%itype(i)) = trim(atoms%atom(j)%SfacSymb)
end do

do i=1,smnum
  smllatlib_new(i)=tmpatlib(i)
end do

deallocate(tmpatlib)
deallocate(xray_form)

end subroutine

subroutine atscatt_oew(reflects,smllatlib_new,atomic_scatt)
type(Xray_Form_Type),allocatable, intent(in) :: smllatlib_new(:) 
type(reflections),                intent(in) :: reflects
type(atscatter),                  intent(inout) :: atomic_scatt 
! load up the atscatt lookup table
integer                       :: nat,ialloc
real(sp)                      :: q(3),a(4),b(4),c,ssq4
integer                       :: i,hh,kk,ll
real(sp)                      :: h,k,l
real                          :: time1,time2

call cpu_time(time1)

nat=size(smllatlib_new)
atomic_scatt%ntypes = nat

if (allocated(atomic_scatt%atype))  deallocate(atomic_scatt%atype)
allocate(atomic_scatt%atype(nat),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"
do i =1, nat
  allocate(atomic_scatt%atype(i)%scatt(reflects%nbragg),stat=ialloc)
  if (ialloc /= 0) STOP "*** not enough memory ***"
end do

do i=1,nat
  a=smllatlib_new(i)%a
  b=smllatlib_new(i)%b
  c=smllatlib_new(i)%c
  do kk = 1, reflects%nbragg
    q = reflects%qBragg_list(kk)%vec
    ssq4=dot_product(q,q)/(16.0d0*pi**2) !4.0d0
    atomic_scatt%atype(i)%scatt(kk) = fat4(a,b,c,ssq4)
  end do
end do
call cpu_time(time2)
print *, 'time to set up look up table', time2-time1, 'seconds'

end subroutine 



subroutine atscatt_new(qvec,reflects,smllatlib_new,atomic_scatt)
real(dp), intent(in) :: qvec(3)
type(Xray_Form_Type),allocatable, intent(in) :: smllatlib_new(:) 
type(reflections),                intent(in) :: reflects
type(atscatter),                  intent(inout) :: atomic_scatt 
! load up the atscatt lookup table
integer                       :: nat,ialloc
real(sp)                      :: q(3),a(4),b(4),c,ssq4
integer                       :: i,hh,kk,ll
real(sp)                      :: h,k,l
real                          :: time1,time2
real(dp) :: qmag


call cpu_time(time1)

nat=size(smllatlib_new)
atomic_scatt%ntypes = nat

if (allocated(atomic_scatt%atype))  deallocate(atomic_scatt%atype)
allocate(atomic_scatt%atype(nat),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"
do i =1, nat
  allocate(atomic_scatt%atype(i)%scatt(reflects%nbragg),stat=ialloc)
  if (ialloc /= 0) STOP "*** not enough memory ***"
end do

qmag = dot_product(qvec,qvec)

if (qmag /= zero)  print *, "adding qvec: ",qvec," to all reciprocal points"

do i=1,nat
  a=smllatlib_new(i)%a
  b=smllatlib_new(i)%b
  c=smllatlib_new(i)%c
  do kk = 1, reflects%nbragg
    q = reflects%qBragg_list(kk)%vec+qvec
    ssq4=dot_product(q,q)/(16.0d0*pi**2) 
    atomic_scatt%atype(i)%scatt(kk) = fat4(a,b,c,ssq4)
  end do
end do
call cpu_time(time2)
print *, 'time to set up look up table', time2-time1, 'seconds'

end subroutine 

integer function total_integers(inti,intf)
! determine number of integers between to integers including those integers
integer, intent(in) :: inti, intf 
integer :: ntmp,first,last

if (inti .le. intf) then
  first = inti
  last  = intf
else
  first = intf
  last  = inti
end if 

total_integers = last - first + 1

end function

integer function total_points(inti,intf,dh)
! determine number of points between to numbers using some step size
! eg. -1.0 -0.75 -0.5 -0.25, 0.0 0.25 0.5 0.75 1.0
! what about :
integer,  intent(in) :: inti, intf 
real(dp), intent(in) :: dh
real(dp) :: hh
integer :: ntmp,first,last,h,nval

nval = npoints(dh)

if (inti .le. intf) then
  first = inti
  last  = intf
else
  first = intf
  last  = inti
end if 

total_points = nval*(last - first) + 1

end function

integer function npoints(dh)
real(dp),intent(in) :: dh
real(dp) :: hh
! count up number of points separated by dh between integers ! when dh isn't sum divisible thing, it's ugly, but I think
! it's ok 
hh = one 
npoints = 1
! now add the rest as long as dh > zero
if (dh > zero .and. dh < one ) then
  do while(hh+dh + 1.0d-08 < two ) ! add error to be sure we're not rounded low
    npoints = npoints+1
    hh = hh + dh
  end do
end if

end function

subroutine reflections_atscatt_setup(input,atoms,cell,reflects,atomic_scatt)
type(inp_par), intent(in)  :: input
type(protein_atom_list_type)    ,intent(in)   :: atoms
type (Crystal_Cell_Type),intent(in)   :: Cell
type(reflections),       intent(out)  :: reflects
type(atscatter),         intent(out)  :: atomic_scatt 
type(Xray_Form_Type),allocatable      :: smllatlib_new(:) 
integer :: i,j,h,k,l !,hstart,hstop,kstart,kstop,lstart,lstop
real (dp) :: val,qvec(3)

qvec = zero

call reflection_setup (input,cell,reflects)
call uniqat_new(atoms,smllatlib_new,atomic_scatt)
call atscatt_new(qvec,reflects,smllatlib_new,atomic_scatt)
!call atscatt_oew(reflects,smllatlib_new,atomic_scatt)
deallocate(smllatlib_new)
!do i = 1, reflects%nbragg
!  val = dsqrt(dot_product(reflects%qBragg_list(i)%vec,reflects%qBragg_list(i)%vec)) 
!  print '(2F12.6)', val, atomic_scatt%atype(1)%Bragg(i)
!end do
end subroutine

subroutine reflections_atscatt_setres(input,res,spg,atoms,cell,reflects,atomic_scatt)
use cfml_crystallographic_symmetry,  only: space_group_type
type(inp_par),            intent(in)   :: input
real(dp),                 intent(in)   :: res
type (Space_Group_Type) , intent(in)   :: spg
type(protein_atom_list_type)    , intent(in)   :: atoms
type (Crystal_Cell_Type), intent(in)   :: Cell
type(reflections),        intent(out)  :: reflects
type(atscatter),          intent(out)  :: atomic_scatt
type(Xray_Form_Type),     allocatable  :: smllatlib_new(:)
integer :: i,j,h,k,l !,hstart,hstop,kstart,kstop,lstart,lstop
real (dp) :: val,qvec(3)

qvec = zero

call unique_reflect_crysfml (res,spg,cell,reflects)
!call reflection_setup_resolution (input,res,cell,reflects)
call uniqat_new(atoms,smllatlib_new,atomic_scatt)
call atscatt_new(qvec,reflects,smllatlib_new,atomic_scatt)
!call atscatt_oew(reflects,smllatlib_new,atomic_scatt)
deallocate(smllatlib_new)
!do i = 1, reflects%nbragg
!  val = dsqrt(dot_product(reflects%qBragg_list(i)%vec,reflects%qBragg_list(i)%vec)) 
!  print '(2F12.6)', val, atomic_scatt%atype(1)%Bragg(i)
!end do
end subroutine

subroutine reflections_atscatt_setbig(input,res,spg,atoms,cell,reflects,atomic_scatt)
use cfml_crystallographic_symmetry,  only: space_group_type
type(inp_par),            intent(in)   :: input
real(dp),                 intent(in)   :: res
type (Space_Group_Type) , intent(in)   :: spg
type(protein_atom_list_type)    , intent(in)   :: atoms
type (Crystal_Cell_Type), intent(in)   :: Cell
type(reflections),        intent(out)  :: reflects
type(atscatter),          intent(out)  :: atomic_scatt
type(Xray_Form_Type),     allocatable  :: smllatlib_new(:)
integer :: i,j,h,k,l !,hstart,hstop,kstart,kstop,lstart,lstop
real (dp) :: val,qvec(3)

qvec = zero

!call unique_reflect_crysfml (res,spg,cell,reflects)
call reflection_setup_resolution (input,res,cell,reflects)
call uniqat_new(atoms,smllatlib_new,atomic_scatt)
call atscatt_new(qvec,reflects,smllatlib_new,atomic_scatt)
!call atscatt_oew(reflects,smllatlib_new,atomic_scatt)
deallocate(smllatlib_new)
!do i = 1, reflects%nbragg
!  val = dsqrt(dot_product(reflects%qBragg_list(i)%vec,reflects%qBragg_list(i)%vec)) 
!  print '(2F12.6)', val, atomic_scatt%atype(1)%Bragg(i)
!end do
end subroutine



subroutine atscatt_setup(qvec,atoms,reflects,atomic_scatt)
real(dp)                ,intent(in)   :: qvec(3)
type(protein_atom_list_type)    ,intent(in)   :: atoms
type(reflections),       intent(in)   :: reflects
type(atscatter),         intent(out)  :: atomic_scatt 
type(Xray_Form_Type),allocatable      :: smllatlib_new(:) 
integer :: i,j,h,k,l !,hstart,hstop,kstart,kstop,lstart,lstop
real (dp) :: val

call uniqat_new(atoms,smllatlib_new,atomic_scatt)
call atscatt_new(qvec,reflects,smllatlib_new,atomic_scatt)
deallocate(smllatlib_new)

!do i = 1, reflects%nbragg
!  val = dsqrt(dot_product(reflects%qBragg_list(i)%vec,reflects%qBragg_list(i)%vec)) 
!  print '(2F12.6)', val, atomic_scatt%atype(1)%Bragg(i)
!end do
end subroutine

subroutine exp_structure_factor(input,cell,atoms,reflects,atomic_scatt,struct_fact)
! DMR 012109 tested the following against cctbx and 
!    our old intensity calc... all good
!    the phases work too, off by some rotes but seemed ok via cctbx
! uses bfactor which 8pi^2*flucts in angs
! |F_hkl| = dsqrt(costerm*costerm + sinterm*sinterm)
! theta   = atan(sinterm/costerm)*180.d0/pi
type(inp_par),              intent(in)      :: input
type (Crystal_Cell_Type),   intent(in)      :: Cell
type(protein_atom_list_type),       intent(inout)   :: atoms
type(reflections),          intent(in)      :: reflects
type(atscatter),            intent(in)      :: atomic_scatt
complex(kind=8),allocatable,intent(out)     :: struct_fact(:)
real(dp) :: eightpisqr ,qc_dot_qc,tfactj,costerm,sinterm,qdotr,twopi,vec(3)
real(dp) :: time1,time2
integer  :: i,j,k,l,ialloc

eightpisqr = 8.0d0 * pi *pi
twopi = 2.0d0*pi
call cpu_time(time1)
allocate(struct_fact(reflects%nbragg), stat=ialloc)

do i = 1, reflects%nbragg
  !vec(:) = real(reflects%qbragg_list(i)%vec(:),kind=8)
  !qc_dot_qc = dot_product(vec,vec)/(2.0d0*eightpisqr)
  qc_dot_qc = dot_product(reflects%qbragg_list(i)%vec,reflects%qbragg_list(i)%vec)/(2.0d0*eightpisqr)
  costerm = zero
  sinterm = zero
 do j = 1, atoms%natoms
    tfactj = dexp(-real(atoms%atom(j)%biso,kind=8)*qc_dot_qc)*real(atoms%atom(j)%occ,kind=8) * &
                       atomic_scatt%atype(atomic_scatt%itype(j))%scatt(i)
    !qdotr = twopi*dot_product(reflects%qbragg_list(i)%vec, real(atoms%atom(j)%x,kind=8))
    qdotr = dot_product(reflects%qbragg_list(i)%vec, real(atoms%atom(j)%x,kind=8))
    costerm = costerm + dcos(qdotr)*tfactj
    sinterm = sinterm + dsin(qdotr)*tfactj
 end do
 struct_fact(i) = cmplx(costerm,sinterm,kind=8)
 !write (33, '(6F16.6)'), reflects%hbragg_list(i)%vec/twopi,dsqrt(costerm*costerm+sinterm*sinterm)
end do

call cpu_time(time2)
print *, 'Structure_factor-directsum', time2-time1, 'seconds'

end subroutine 

subroutine intensity_zero_aniso(atoms,reflects,aniso,atomic_scatt, &
                           intensity)
! DMR 012109 
!    our old intensity calc... all good
!    the phases work too, off by some rotes but seemed ok via cctbx
use mod_vcov_store
type(protein_atom_list_type),       intent(in)      :: atoms
type(reflections),          intent(in)      :: reflects
real(dp)        , allocatable,  intent(in)  :: aniso(:,:)
type(atscatter),            intent(in)      :: atomic_scatt
real(dp),allocatable,       intent(out)     :: intensity(:)
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,qdotr,qc(3)
real(dp) :: time1,time2
integer  :: i,j,ialloc,jcor


call cpu_time(time1)
allocate(intensity(reflects%nbragg), stat=ialloc)

do i = 1, reflects%nbragg
  qc = reflects%qbragg_list(i)%vec
  costerm = zero
  sinterm = zero
  jcor = 1
  do j = 1, atoms%natoms
    ! DMR, checked some timing with this, and there doesn't seem to be much gain by enumerating
    !      but using matmuls does slow things down
    jexpon  =-(qc(1)*dot_product(qc,aniso(jcor:jcor+2,1)) + &
               qc(2)*dot_product(qc,aniso(jcor:jcor+2,2)) + &
               qc(3)*dot_product(qc,aniso(jcor:jcor+2,3)))/2.0d0

! some thinking required on the occupancy stuff
    tfactj  = dexp(jexpon)*real(atoms%atom(j)%occ,kind=8)*atomic_scatt%atype(atomic_scatt%itype(j))%scatt(i)
    qdotr   = dot_product(reflects%qbragg_list(i)%vec, real(atoms%atom(j)%x,kind=8))
    costerm = costerm + dcos(qdotr)*tfactj
    sinterm = sinterm + dsin(qdotr)*tfactj
    jcor = jcor + 3
  end do
  intensity(i) = costerm*costerm+sinterm*sinterm
end do

call cpu_time(time2)

print *, 'intensity_zero_aniso>', time2-time1, 'seconds'

end subroutine 

subroutine blocks_Fsqr(atoms,reflects,atomic_scatt,blocks)
! DMR added Oct 1, 2009
! what does it do asshole!
use mod_bnm,  only: block
type(protein_atom_list_type),       intent(in)      :: atoms
type(reflections),          intent(in)      :: reflects
type(atscatter),            intent(in)      :: atomic_scatt
type(block),    allocatable,intent(inout)   :: blocks(:)
real(dp) :: qc(3),qdotr,atscti,atsctj,tmpvec(3),cdmult
real(dp) :: time1,time2,costerm,sinterm
integer  :: i,j,ier,ialloc,ati,atj,ii,jj

do i = 1, size(blocks,1)
  if (allocated(blocks(i)%F_sqr)) deallocate(blocks(i)%F_sqr)
  allocate(blocks(i)%F_sqr(reflects%nbragg), stat=ialloc)
  blocks(i)%F_sqr=zero
end do

do i = 1, reflects%nbragg
  qc = reflects%qbragg_list(i)%vec

  do j = 1, size(blocks,1)
    costerm = zero
    sinterm = zero
    do ii=1,size(blocks(j)%atinblock,1)
      ati = blocks(j)%atinblock(ii)
      atscti = real(atoms%atom(ati)%occ,kind=8) * atomic_scatt%atype(atomic_scatt%itype(ati))%scatt(i)
      qdotr = dot_product(reflects%qbragg_list(i)%vec, real(atoms%atom(ati)%x,kind=8))
      costerm = costerm + dcos(qdotr)*atscti
      sinterm = sinterm + dsin(qdotr)*atscti
    end do
    blocks(j)%F_sqr(i) = costerm*costerm + sinterm*sinterm 
  end do

end do    

end subroutine

subroutine atself_diffscat_aniso(input,atoms,reflects,aniso,atomic_scatt)
! DMR added Oct 1, 2009
use mod_vcov_store
type(inp_par), intent(in)                   :: input
type(protein_atom_list_type),       intent(in)      :: atoms
type(reflections),          intent(in)      :: reflects
real(dp)        , allocatable,  intent(in)  :: aniso(:,:)
type(atscatter),            intent(in)      :: atomic_scatt
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,qdotr,qc(3)
real(dp) :: time1,time2
integer  :: i,j,ier,ialloc,jcor


open(unit=23,file=trim(input%fileroot)//"-diffscat-"//trim(input%tpbc)//"-"// &
  trim(input%genm)// "-" //trim(input%fctyp)//"-atself.txt",& 
     status='unknown',action='write', &
     iostat=ier)

call cpu_time(time1)
!allocate(intensity(reflects%nbragg), stat=ialloc)

do i = 1, reflects%nbragg
  qc = reflects%qbragg_list(i)%vec
  costerm = zero
  jcor = 1
  do j = 1, atoms%natoms
    ! DMR, checked some timing with this, and there doesn't seem to be much gain by enumerating
    !      but using matmuls does slow things down
    jexpon  =-(qc(1)*dot_product(qc,aniso(jcor:jcor+2,1)) + &
               qc(2)*dot_product(qc,aniso(jcor:jcor+2,2)) + &
               qc(3)*dot_product(qc,aniso(jcor:jcor+2,3)))/2.0d0

! some thinking required on the occupancy stuff
    tfactj  = (one-dexp(2.0d0*jexpon))*real(atoms%atom(j)%occ,kind=8)*(atomic_scatt%atype(atomic_scatt%itype(j))%scatt(i)**2)
    costerm = costerm + tfactj
    jcor = jcor + 3
  end do
  write(23, '(6F10.5,F16.5)'), &
                            reflects%hbragg_list(i)%vec/(two*pi), &
                            reflects%qbragg_list(i)%vec/(two*pi), &
                            costerm
  !intensity(i) = costerm
end do

call cpu_time(time2)

print *, 'intensity_atself_aniso>', time2-time1, 'seconds'

end subroutine 

subroutine blockself_diffscat_iso(input,blocks,reflects,iso)
! DMR added Oct 1, 2009
use mod_bnm, only: block
type(inp_par), intent(in)                   :: input
type(block)  ,allocatable,  intent(in)      :: blocks(:)
type(reflections),          intent(in)      :: reflects
real(dp)        , allocatable,  intent(in)  :: iso(:)
real(dp) :: tfactj,costerm,qc_dot_qc
real(dp) :: time1,time2
integer  :: i,j,ier,ialloc,jcor


open(unit=23,file=trim(input%fileroot)//"-diffscat-"//trim(input%tpbc)//"-"// &
  trim(input%genm)// "-" //trim(input%fctyp)//"-isoblockself.txt",& 
     status='unknown',action='write', &
     iostat=ier)

call cpu_time(time1)
!allocate(intensity(reflects%nbragg), stat=ialloc)

do i = 1, reflects%nbragg
  qc_dot_qc = dot_product(reflects%qbragg_list(i)%vec,reflects%qbragg_list(i)%vec)/2.0d0
  costerm = zero
  jcor = 1
  do j = 1, size(blocks,1)
    tfactj  = (one-dexp(-qc_dot_qc*iso(j)))*blocks(j)%F_sqr(i)
    costerm = costerm + tfactj
  end do

  write(23, '(6F10.5,F16.5)'), &
                            reflects%hbragg_list(i)%vec/(two*pi), &
                            reflects%qbragg_list(i)%vec/(two*pi), &
                            costerm
  !intensity(i) = costerm
end do

call cpu_time(time2)

print *, 'intensity_atself_aniso>', time2-time1, 'seconds'

end subroutine 

subroutine blockself_diffscat_aniso(input,blocks,reflects,aniso)
! DMR added Oct 1, 2009
use mod_bnm, only: block
type(inp_par), intent(in)                   :: input
type(block)  ,allocatable,  intent(in)      :: blocks(:)
type(reflections),          intent(in)      :: reflects
real(dp)        , allocatable,  intent(in)  :: aniso(:,:)
real(dp) :: tfactj,jexpon,costerm,sinterm,qdotr,qc(3)
real(dp) :: time1,time2
integer  :: i,j,ier,ialloc,jcor


open(unit=23,file=trim(input%fileroot)//"-diffscat-"//trim(input%tpbc)//"-"// &
  trim(input%genm)// "-" //trim(input%fctyp)//"-blockself.txt",& 
     status='unknown',action='write', &
     iostat=ier)

call cpu_time(time1)
!allocate(intensity(reflects%nbragg), stat=ialloc)

do i = 1, reflects%nbragg
  qc = reflects%qbragg_list(i)%vec
  costerm = zero
  jcor = 1
  do j = 1, size(blocks,1)
    ! DMR, checked some timing with this, and there doesn't seem to be much gain by enumerating
    !      but using matmuls does slow things down
    jexpon  =-(qc(1)*dot_product(qc,aniso(jcor:jcor+2,1)) + &
               qc(2)*dot_product(qc,aniso(jcor:jcor+2,2)) + &
               qc(3)*dot_product(qc,aniso(jcor:jcor+2,3)))/2.0d0

! some thinking required on the occupancy stuff
    tfactj  = (one-dexp(jexpon))*blocks(j)%F_sqr(i)
    costerm = costerm + tfactj
    jcor = jcor + 3
  end do
  write(23, '(6F10.5,F16.5)'), &
                            reflects%hbragg_list(i)%vec/(two*pi), &
                            reflects%qbragg_list(i)%vec/(two*pi), &
                            costerm
  !intensity(i) = costerm
end do

call cpu_time(time2)

print *, 'intensity_atself_aniso>', time2-time1, 'seconds'

end subroutine 

subroutine intensity_vcov(atoms,reflects,vcov_cross,atomic_scatt, &
                           intensity,anisoin)
use mod_vcov_store
type(protein_atom_list_type),       intent(in)           :: atoms
type(reflections),          intent(in)           :: reflects
real(dp), allocatable,       intent(in)          :: vcov_cross(:,:)
type(atscatter),            intent(in)           :: atomic_scatt
real(dp),allocatable,       intent(out)          :: intensity(:)
real(dp), allocatable,      intent(in), optional :: anisoin(:,:)
real(dp), allocatable :: aniso(:,:)
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,qdotr,qc(3)
real(dp) :: time1,time2
real(dp) :: dwfi,dwfj,dwfij,varcovar,cdmult,atscti,tmpvec(3)
integer  :: i,ii,jj,ialloc,icor,jcor,ndim


call cpu_time(time1)

ndim = atoms%natoms*3
i = size(vcov_cross(:,1))
if (ndim .ne. i)  then
  print *, "intensity_vcov>",ndim, i
  stop "intensity_vcov> something wrong with dimensions"
end if
allocate(aniso(ndim,3), stat=ialloc)
if (ialloc /= 0) stop "intensity_vcov> memory problem"

! construct self blocks
if (present(anisoin)) then
  aniso = anisoin
else

  icor = 1
  do i=1, atoms%natoms
    aniso(icor:icor+2,:) = vcov_cross(icor:icor+2,icor:icor+2)
    icor = icor + 3
  end do
end if

allocate(intensity(reflects%nbragg), stat=ialloc)
if (ialloc /= 0) stop "intensity_vcov> memory problem"

do i = 1, reflects%nbragg
  qc = reflects%qbragg_list(i)%vec
  costerm = zero
  icor = 1
  do ii = 1, atoms%natoms
    dwfi =-(qc(1)*dot_product(qc,aniso(icor:icor+2,1)) + &
            qc(2)*dot_product(qc,aniso(icor:icor+2,2)) + &
            qc(3)*dot_product(qc,aniso(icor:icor+2,3)))/2.0d0
    atscti = real(atoms%atom(ii)%occ,kind=8)*atomic_scatt%atype(atomic_scatt%itype(ii))%scatt(i)
    jcor = icor

    do jj = ii, atoms%natoms
      
      dwfj  =-(qc(1)*dot_product(qc,aniso(jcor:jcor+2,1)) + &
               qc(2)*dot_product(qc,aniso(jcor:jcor+2,2)) + &
               qc(3)*dot_product(qc,aniso(jcor:jcor+2,3)))/2.0d0
      dwfij = (qc(1)*dot_product(qc,vcov_cross(icor:icor+2,jcor)) + &
               qc(2)*dot_product(qc,vcov_cross(icor:icor+2,jcor+1)) + &
               qc(3)*dot_product(qc,vcov_cross(icor:icor+2,jcor+2)))      
      varcovar = dexp(dwfi+dwfj+dwfij)*atscti*real(atoms%atom(jj)%occ,kind=8)* &
                 atomic_scatt%atype(atomic_scatt%itype(jj))%scatt(i)

      if (ii.eq.jj) then
        cdmult = one
      else
        cdmult=two
      end if

      tmpvec  = real(atoms%atom(ii)%x,kind=8)-real(atoms%atom(jj)%x,kind=8)
      qdotr   = dot_product(reflects%qbragg_list(i)%vec, tmpvec)
      costerm = costerm + cdmult*dcos(qdotr)*varcovar
      jcor = jcor +3
    end do
    icor = icor +3
  end do
  intensity(i) = costerm
end do

call cpu_time(time2)

print *, 'intensity_vcov>', time2-time1, 'seconds'

end subroutine 

subroutine diffscat_vcov(input,atoms,reflects,vcov_cross,atomic_scatt, &
                           anisoin)
use mod_vcov_store
type(inp_par),              intent(in)           :: input
type(protein_atom_list_type),       intent(in)           :: atoms
type(reflections),          intent(in)           :: reflects
real(dp), allocatable,       intent(in)          :: vcov_cross(:,:)
type(atscatter),            intent(in)           :: atomic_scatt
real(dp), allocatable,      intent(in), optional :: anisoin(:,:)
real(dp), allocatable :: aniso(:,:)
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,qdotr,qc(3)
real(dp) :: time1,time2
real(dp) :: dwfi,dwfj,dwfij,varcovar,cdmult,atscti,tmpvec(3)
integer  :: i,ii,jj,ialloc,icor,jcor,ndim,ier


call cpu_time(time1)

open(unit=23,file=trim(input%fileroot)//"-diffscat-"//trim(input%tpbc)//"-"// &
  trim(input%genm)// "-" //trim(input%fctyp)//"-vcov.txt",& 
     status='unknown',action='write', &
     iostat=ier)

ndim = atoms%natoms*3
i = size(vcov_cross(:,1))
if (ndim .ne. i)  then
  print *, "intensity_vcov>",ndim, i
  stop "intensity_vcov> something wrong with dimensions"
end if
allocate(aniso(ndim,3), stat=ialloc)
if (ialloc /= 0) stop "intensity_vcov> memory problem"

! construct self blocks
if (present(anisoin)) then
  aniso = anisoin
else

  icor = 1
  do i=1, atoms%natoms
    aniso(icor:icor+2,:) = vcov_cross(icor:icor+2,icor:icor+2)
    icor = icor + 3
  end do
end if

do i = 1, reflects%nbragg
  qc = reflects%qbragg_list(i)%vec
  costerm = zero
  icor = 1
  do ii = 1, atoms%natoms
    dwfi =-(qc(1)*dot_product(qc,aniso(icor:icor+2,1)) + &
            qc(2)*dot_product(qc,aniso(icor:icor+2,2)) + &
            qc(3)*dot_product(qc,aniso(icor:icor+2,3)))/2.0d0
    atscti = real(atoms%atom(ii)%occ,kind=8)*atomic_scatt%atype(atomic_scatt%itype(ii))%scatt(i)
    jcor = icor

    do jj = ii, atoms%natoms
      
      dwfj  =-(qc(1)*dot_product(qc,aniso(jcor:jcor+2,1)) + &
               qc(2)*dot_product(qc,aniso(jcor:jcor+2,2)) + &
               qc(3)*dot_product(qc,aniso(jcor:jcor+2,3)))/2.0d0
      dwfij = (qc(1)*dot_product(qc,vcov_cross(icor:icor+2,jcor)) + &
               qc(2)*dot_product(qc,vcov_cross(icor:icor+2,jcor+1)) + &
               qc(3)*dot_product(qc,vcov_cross(icor:icor+2,jcor+2)))      
      !varcovar = dexp(dwfi+dwfj+dwfij)*atscti*real(atoms%atom(jj)%occ,kind=8)* &
      !           atomic_scatt%atype(atomic_scatt%itype(jj))%scatt(i)
      varcovar = atscti*real(atoms%atom(jj)%occ,kind=8)* &
                 atomic_scatt%atype(atomic_scatt%itype(jj))%scatt(i) * &
                 (dexp(dwfi+dwfj+dwfij)-dexp(dwfi+dwfj))

      if (ii.eq.jj) then
        cdmult = one
      else
        cdmult=two
      end if

      tmpvec  = real(atoms%atom(ii)%x,kind=8)-real(atoms%atom(jj)%x,kind=8)
      qdotr   = dot_product(reflects%qbragg_list(i)%vec, tmpvec)
      costerm = costerm + cdmult*dcos(qdotr)*varcovar
      jcor = jcor +3
    end do
    icor = icor +3
  end do
  !intensity(i) = costerm
  write(23, '(6F10.5,F16.5)'), &
                            reflects%hbragg_list(i)%vec/(two*pi), &
                            reflects%qbragg_list(i)%vec/(two*pi), &
                            costerm
end do

call cpu_time(time2)

print *, 'intensity_vcov>', time2-time1, 'seconds'

end subroutine 

subroutine block_diffscat_vcov(input,blocks,atoms,reflects,vcov_cross,atomic_scatt, &
                           anisoin)
use mod_vcov_store
use mod_bnm, only: block
type(inp_par),              intent(in)           :: input
type(block),allocatable,    intent(in)           :: blocks(:)
type(protein_atom_list_type),       intent(in)           :: atoms
type(reflections),          intent(in)           :: reflects
real(dp), allocatable,       intent(in)          :: vcov_cross(:,:)
type(atscatter),            intent(in)           :: atomic_scatt
real(dp), allocatable,      intent(in), optional :: anisoin(:,:)
real(dp), allocatable :: aniso(:,:)
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,qdotr,qc(3)
real(dp) :: time1,time2
real(dp) :: dwfi,dwfj,dwfij,varcovar,cdmult,atscti,tmpvec(3)
integer  :: i,ii,jj,ati,atj,iblk,ialloc,icor,jcor,ndim,ier


call cpu_time(time1)

open(unit=23,file=trim(input%fileroot)//"-diffscat-"//trim(input%tpbc)//"-"// &
  trim(input%genm)// "-" //trim(input%fctyp)//"-blockvcov.txt",& 
     status='unknown',action='write', &
     iostat=ier)

ndim = atoms%natoms*3
i = size(vcov_cross(:,1))
if (ndim .ne. i)  then
  print *, "intensity_vcov>",ndim, i
  stop "intensity_vcov> something wrong with dimensions"
end if
allocate(aniso(ndim,3), stat=ialloc)
if (ialloc /= 0) stop "intensity_vcov> memory problem"

! construct self blocks
if (present(anisoin)) then
  aniso = anisoin
else

  icor = 1
  do i=1, atoms%natoms
    aniso(icor:icor+2,:) = vcov_cross(icor:icor+2,icor:icor+2)
    icor = icor + 3
  end do
end if

do i = 1, reflects%nbragg
  qc = reflects%qbragg_list(i)%vec
  
! loop over blocks
  costerm = zero
  do iblk = 1, size(blocks)

  do ii = 1, size(blocks(iblk)%atinblock,1) 
    ati = blocks(iblk)%atinblock(ii)
    icor = 3*(ati-1)+1
    dwfi =-(qc(1)*dot_product(qc,aniso(icor:icor+2,1)) + &
            qc(2)*dot_product(qc,aniso(icor:icor+2,2)) + &
            qc(3)*dot_product(qc,aniso(icor:icor+2,3)))/2.0d0

    atscti =  real(atoms%atom(ati)%occ,kind=8) * &
                   atomic_scatt%atype(atomic_scatt%itype(ati))%scatt(i)

    do jj = ii, size(blocks(iblk)%atinblock,1)
      atj = blocks(iblk)%atinblock(jj)
      jcor = 3*(atj-1)+1

      dwfj  =-(qc(1)*dot_product(qc,aniso(jcor:jcor+2,1)) + &
               qc(2)*dot_product(qc,aniso(jcor:jcor+2,2)) + &
               qc(3)*dot_product(qc,aniso(jcor:jcor+2,3)))/2.0d0

      dwfij = (qc(1)*dot_product(qc,vcov_cross(icor:icor+2,jcor)) + &
               qc(2)*dot_product(qc,vcov_cross(icor:icor+2,jcor+1)) + &
               qc(3)*dot_product(qc,vcov_cross(icor:icor+2,jcor+2)))      
      !varcovar = dexp(dwfi+dwfj+dwfij)*atscti*real(atoms%atom(jj)%occ,kind=8)* &
      !           atomic_scatt%atype(atomic_scatt%itype(jj))%scatt(i)
      varcovar = atscti*real(atoms%atom(atj)%occ,kind=8)* &
                 atomic_scatt%atype(atomic_scatt%itype(atj))%scatt(i) * &
                 (dexp(dwfi+dwfj+dwfij)-dexp(dwfi+dwfj))

      if (ati.eq.atj) then
        cdmult = one
      else
        cdmult=two
      end if

      tmpvec  = real(atoms%atom(ati)%x,kind=8)-real(atoms%atom(atj)%x,kind=8)
      qdotr   = dot_product(reflects%qbragg_list(i)%vec, tmpvec)
      costerm = costerm + cdmult*dcos(qdotr)*varcovar
    end do
  end do
  end do
  !intensity(i) = costerm
  write(23, '(6F10.5,F16.5)'), &
                            reflects%hbragg_list(i)%vec/(two*pi), &
                            reflects%qbragg_list(i)%vec/(two*pi), &
                            costerm
end do

call cpu_time(time2)

print *, 'intensity_vcov>', time2-time1, 'seconds'

end subroutine 

subroutine block_diffscat_vcov_iso(input,blocks,atoms,reflects,isovcov,atomic_scatt, &
                           isoin,blockn)
use mod_vcov_store
use mod_bnm, only: block
type(inp_par),              intent(in)           :: input
type(block),allocatable,    intent(in)           :: blocks(:)
type(protein_atom_list_type),       intent(in)           :: atoms
type(reflections),          intent(in)           :: reflects
real(dp), allocatable,       intent(in)          :: isovcov(:,:)
type(atscatter),            intent(in)           :: atomic_scatt
real(dp), allocatable,      intent(in), optional :: isoin(:)
character(len=25)    ,      intent(in), optional :: blockn
real(dp), allocatable :: iso(:)
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,qdotr,qc(3)
real(dp) :: time1,time2
real(dp) :: dwfi,dwfj,dwfij,varcovar,cdmult,atscti,tmpvec(3)
integer  :: i,j,ii,jj,ati,atj,iblk,ialloc,icor,jcor,ndim,ier
character(len=25)     :: nameit

if (present(blockn)) then
  nameit = blockn
else
  nameit = "block"
end if

call cpu_time(time1)

open(unit=23,file=trim(input%fileroot)//"-diffscat-"//trim(input%tpbc)//"-"// &
  trim(input%genm)// "-" //trim(input%fctyp)//"-"//trim(nameit)//"-isov.txt",& 
     status='unknown',action='write', &
     iostat=ier)

ndim = atoms%natoms
i = size(isovcov(:,1))
if (ndim .ne. i) stop "intensity_vcov> something wrong with dimensions"
allocate(iso(ndim), stat=ialloc)
if (ialloc /= 0) stop "intensity_vcov> memory problem"

! construct self blocks
if (present(isoin)) then
  
  iso = isoin

else

  do i=1, atoms%natoms
    iso(i) = isovcov(i,i)
  end do

end if

do i = 1, reflects%nbragg
  qc_dot_qc = dot_product(reflects%qbragg_list(i)%vec,reflects%qbragg_list(i)%vec)/2.0d0
  
! loop over blocks
  costerm = zero
  do iblk = 1, size(blocks)
  do ii = 1, size(blocks(iblk)%atinblock,1) 
    ati = blocks(iblk)%atinblock(ii)
    dwfi =-qc_dot_qc*iso(ati) 

    atscti =  real(atoms%atom(ati)%occ,kind=8) * &
                   atomic_scatt%atype(atomic_scatt%itype(ati))%scatt(i)

    do jj = ii, size(blocks(iblk)%atinblock,1)
      atj = blocks(iblk)%atinblock(jj)

      dwfj = -qc_dot_qc*iso(atj) 
      dwfij= qc_dot_qc*isovcov(ati,atj)*2.0d0

      varcovar = atscti*real(atoms%atom(atj)%occ,kind=8)* &
                 atomic_scatt%atype(atomic_scatt%itype(atj))%scatt(i) * &
                 (dexp(dwfi+dwfj+dwfij)-dexp(dwfi+dwfj))
                 !dexp(dwfi+dwfj)*(dexp(dwfij)-one) 

      if (ati.eq.atj) then
        cdmult = one
      else
        cdmult=two
      end if

      tmpvec  = real(atoms%atom(ati)%x,kind=8)-real(atoms%atom(atj)%x,kind=8)
      qdotr   = dot_product(reflects%qbragg_list(i)%vec, tmpvec)
      costerm = costerm + cdmult*dcos(qdotr)*varcovar
    end do
  end do
  end do
  !intensity(i) = costerm
  write(23, '(6F10.5,F16.5)'), &
                            reflects%hbragg_list(i)%vec/(two*pi), &
                            reflects%qbragg_list(i)%vec/(two*pi), &
                            costerm
end do

call cpu_time(time2)

print *, 'intensity_vcov>', time2-time1, 'seconds'

end subroutine 

subroutine intensity_vcov_iso(atoms,reflects,isovcov,atomic_scatt, &
                           intensity,isoin)
use mod_vcov_store
type(protein_atom_list_type),       intent(in)           :: atoms
type(reflections),          intent(in)           :: reflects
real(dp), allocatable,       intent(in)          :: isovcov(:,:)
type(atscatter),            intent(in)           :: atomic_scatt
real(dp),allocatable,       intent(out)          :: intensity(:)
real(dp), allocatable,      intent(in), optional :: isoin(:)
real(dp), allocatable :: iso(:)
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,qdotr,qc(3)
real(dp) :: time1,time2
real(dp) :: dwfi,dwfj,dwfij,varcovar,cdmult,atscti,tmpvec(3)
integer  :: i,ii,jj,ialloc,icor,jcor,ndim


call cpu_time(time1)

ndim = atoms%natoms
i = size(isovcov(:,1))
if (ndim .ne. i) stop "intensity_vcov> something wrong with dimensions"
allocate(iso(ndim), stat=ialloc)
if (ialloc /= 0) stop "intensity_vcov> memory problem"

! construct self blocks
if (present(isoin)) then
  iso = isoin
else

  do i=1, atoms%natoms
    iso(i) = isovcov(i,i)
  end do

end if

allocate(intensity(reflects%nbragg), stat=ialloc)
if (ialloc /= 0) stop "intensity_vcov> memory problem"

do i = 1, reflects%nbragg
  qc_dot_qc = dot_product(reflects%qbragg_list(i)%vec,reflects%qbragg_list(i)%vec)/2.0d0
  costerm = zero
  icor = 1
  do ii = 1, atoms%natoms
    dwfi   = -qc_dot_qc*iso(ii)
    atscti = real(atoms%atom(ii)%occ,kind=8)*atomic_scatt%atype(atomic_scatt%itype(ii))%scatt(i)
    jcor   = icor

    do jj = ii, atoms%natoms
      
      dwfj = -qc_dot_qc*iso(jj)
      dwfij= qc_dot_qc*isovcov(ii,jj)*2.0d0
      varcovar = dexp(dwfi+dwfj+dwfij)*atscti*real(atoms%atom(jj)%occ,kind=8)* &
                 atomic_scatt%atype(atomic_scatt%itype(jj))%scatt(i)

      if (ii.eq.jj) then
        cdmult = one
      else
        cdmult=two
      end if

      tmpvec  = real(atoms%atom(ii)%x,kind=8)-real(atoms%atom(jj)%x,kind=8)
      qdotr   = dot_product(reflects%qbragg_list(i)%vec, tmpvec)
      costerm = costerm + cdmult*dcos(qdotr)*varcovar
      jcor = jcor +3
    end do
    icor = icor +3
  end do
  intensity(i) = costerm
end do

call cpu_time(time2)

print *, 'intensity_vcov_iso>', time2-time1, 'seconds'

end subroutine 

subroutine diffscat_vcov_iso(input,atoms,reflects,isovcov,atomic_scatt, &
                           isoin)
use mod_vcov_store
type(inp_par), intent(in) :: input
type(protein_atom_list_type),       intent(in)           :: atoms
type(reflections),          intent(in)           :: reflects
real(dp), allocatable,       intent(in)          :: isovcov(:,:)
type(atscatter),            intent(in)           :: atomic_scatt
!real(dp),allocatable,       intent(out)          :: intensity(:)
real(dp), allocatable,      intent(in), optional :: isoin(:)
real(dp), allocatable :: iso(:)
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,qdotr,qc(3)
real(dp) :: time1,time2
real(dp) :: dwfi,dwfj,dwfij,varcovar,cdmult,atscti,tmpvec(3)
integer  :: i,j,ii,jj,ialloc,ndim,ier


open(unit=23,file=trim(input%fileroot)//"-diffscat-"//trim(input%tpbc)//"-"// &
  trim(input%genm)// "-" //trim(input%fctyp)//"-isovcov.txt",& 
     status='unknown',action='write', &
     iostat=ier)

call cpu_time(time1)

ndim = atoms%natoms
i = size(isovcov(:,1))
if (ndim .ne. i) stop "intensity_vcov> something wrong with dimensions"
allocate(iso(ndim), stat=ialloc)
if (ialloc /= 0) stop "intensity_vcov> memory problem"

! construct self blocks
if (present(isoin)) then
  iso = isoin
else

  do i=1, atoms%natoms
    iso(i) = isovcov(i,i)
  end do

end if

!allocate(intensity(reflects%nbragg), stat=ialloc)
!if (ialloc /= 0) stop "intensity_vcov> memory problem"

do i = 1, reflects%nbragg
  qc_dot_qc = dot_product(reflects%qbragg_list(i)%vec,reflects%qbragg_list(i)%vec)/2.0d0
  costerm = zero
  do ii = 1, atoms%natoms
    dwfi   = -qc_dot_qc*iso(ii)
    atscti = real(atoms%atom(ii)%occ,kind=8)*atomic_scatt%atype(atomic_scatt%itype(ii))%scatt(i)

    do jj = ii, atoms%natoms
      dwfj = -qc_dot_qc*iso(jj)
      dwfij= qc_dot_qc*isovcov(ii,jj)*2.0d0
      varcovar = atscti*real(atoms%atom(jj)%occ,kind=8)* &
                 atomic_scatt%atype(atomic_scatt%itype(jj))%scatt(i)*(dexp(dwfi+dwfj+dwfij)-dexp(dwfi+dwfj)) !morestable numwise, but slower

     ! varcovar = atscti*real(atoms%atom(jj)%occ,kind=8)* &
    !             atomic_scatt%atype(atomic_scatt%itype(jj))%scatt(i)*dexp(dwfi+dwfj)*(dexp(dwfij)-one) !morestable numwise, but slower
      if (ii.eq.jj) then
        cdmult = one
      else
        cdmult=two
      end if

      tmpvec  = real(atoms%atom(ii)%x,kind=8)-real(atoms%atom(jj)%x,kind=8)
      qdotr   = dot_product(reflects%qbragg_list(i)%vec, tmpvec)
      costerm = costerm + cdmult*dcos(qdotr)*varcovar
    end do
  end do
!  intensity(i) = costerm
  write(23, '(6F10.5,F16.5)'), &
                            reflects%hbragg_list(i)%vec/(two*pi), &
                            reflects%qbragg_list(i)%vec/(two*pi), &
                            costerm
end do

call cpu_time(time2)

print *, 'intensity_vcov_iso>', time2-time1, 'seconds'

end subroutine 

subroutine atself_diffscat_iso(input,atoms,reflects,iso,atomic_scatt)
! DMR added Oct 1, 2009
! takes atomic fluctuations (isoin) in angstrom
use mod_vcov_store
type(inp_par), intent(in) :: input
type(protein_atom_list_type),       intent(in)           :: atoms
type(reflections),          intent(in)           :: reflects
real(dp), allocatable,       intent(in)          :: iso(:)
type(atscatter),            intent(in)           :: atomic_scatt
!real(dp),allocatable,       intent(out)          :: intensity(:)
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,qdotr,qc(3)
real(dp) :: time1,time2
real(dp) :: one_minus_dwfi,varcovar,cdmult,atscti,tmpvec(3)
integer  :: i,ii,jj,ialloc,icor,jcor,ndim,ier


open(unit=23,file=trim(input%fileroot)//"-diffscat-"//trim(input%tpbc)//"-"// &
  trim(input%genm)// "-" //trim(input%fctyp)//"-isoatself.txt",& 
     status='unknown',action='write', &
     iostat=ier)

call cpu_time(time1)

! not keeping intensity in memory for now, just writing out
!allocate(intensity(reflects%nbragg), stat=ialloc)
!if (ialloc /= 0) stop "intensity_vcov> memory problem"

do i = 1, reflects%nbragg
  qc_dot_qc = dot_product(reflects%qbragg_list(i)%vec,reflects%qbragg_list(i)%vec)/2.0d0
  costerm = zero
  do ii = 1, atoms%natoms
    one_minus_dwfi   = (one - dexp(-qc_dot_qc*iso(ii)*2.0d0))
    costerm = costerm + one_minus_dwfi*real(atoms%atom(ii)%occ,kind=8)*(atomic_scatt%atype(atomic_scatt%itype(ii))%scatt(i)**2)
  end do
!  intensity(i) = costerm
  write(23, '(6F10.5,F16.5)'), &
                            reflects%hbragg_list(i)%vec/(two*pi), &
                            reflects%qbragg_list(i)%vec/(two*pi), &
                            costerm
end do

call cpu_time(time2)

print *, 'intensity_atself>', time2-time1, 'seconds'

end subroutine 

subroutine intensity_vcov_test(atoms,reflects,vcov_self,vcov_blocks,atomic_scatt, &
                           intensity)
use mod_vcov_store
type(protein_atom_list_type),       intent(in)      :: atoms
type(reflections),          intent(in)      :: reflects
real(dp), allocatable,       intent(in)      :: vcov_self(:,:),vcov_blocks(:,:)
type(atscatter),            intent(in)      :: atomic_scatt
real(dp),allocatable,       intent(out)     :: intensity(:)
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,qdotr,qc(3)
real(dp) :: time1,time2
real(dp) :: dwfi,dwfj,dwfij,varcovar,cdmult,atscti,tmpvec(3)
integer  :: i,ii,jj,ialloc,icor,jcor


call cpu_time(time1)
allocate(intensity(reflects%nbragg), stat=ialloc)

do i = 1, reflects%nbragg
  qc = reflects%qbragg_list(i)%vec
  costerm = zero
  icor = 1
  do ii = 1, atoms%natoms
    dwfi =-(qc(1)*dot_product(qc,vcov_self(icor:icor+2,1)) + &
            qc(2)*dot_product(qc,vcov_self(icor:icor+2,2)) + &
            qc(3)*dot_product(qc,vcov_self(icor:icor+2,3)))/2.0d0
    atscti = real(atoms%atom(ii)%occ,kind=8)*atomic_scatt%atype(atomic_scatt%itype(ii))%scatt(i)
    jcor = icor

    do jj = ii, atoms%natoms
      
      dwfj =-(qc(1)*dot_product(qc,vcov_self(jcor:jcor+2,1)) + &
            qc(2)*dot_product(qc,vcov_self(jcor:jcor+2,2)) + &
            qc(3)*dot_product(qc,vcov_self(jcor:jcor+2,3)))/2.0d0
      dwfij =(qc(1)*dot_product(qc,vcov_blocks(icor:icor+2,jcor)) + &
            qc(2)*dot_product(qc,vcov_blocks(icor:icor+2,jcor+1)) + &
            qc(3)*dot_product(qc,vcov_blocks(icor:icor+2,jcor+2)))      
      varcovar = dexp(dwfi+dwfj+dwfij)*atscti*real(atoms%atom(jj)%occ,kind=8)* &
                 atomic_scatt%atype(atomic_scatt%itype(jj))%scatt(i)

      if (ii.eq.jj) then
        cdmult = one
      else
        cdmult=two
      end if

      tmpvec  = real(atoms%atom(ii)%x,kind=8)-real(atoms%atom(jj)%x,kind=8)
      qdotr   = dot_product(reflects%qbragg_list(i)%vec, tmpvec)
      costerm = costerm + cdmult*dcos(qdotr)*varcovar
      jcor = jcor +3

    end do
    icor = icor +3
  end do
  intensity(i) = costerm
end do

call cpu_time(time2)

print *, 'intensity_vcov>', time2-time1, 'seconds'

end subroutine 

subroutine intensity_zero_iso(atoms,reflects,flucts,atomic_scatt, &
                           intensity)
! DMR 012109 tested the following against cctbx and 
!    our old intensity calc... all good
!    the phases work too, off by some rotes but seemed ok via cctbx
! |F_hkl| = dsqrt(costerm*costerm + sinterm*sinterm)
! theta   = atan(sinterm/costerm)*180.d0/pi
type(protein_atom_list_type),       intent(in)      :: atoms
type(reflections),          intent(in)      :: reflects
real(dp)   ,  allocatable,  intent(in)      :: flucts(:)
type(atscatter),            intent(in)      :: atomic_scatt
real(dp),allocatable,       intent(out)     :: intensity(:)
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,qdotr,qc(3)
real(dp) :: time1,time2
integer  :: i,j,k,l,ialloc

call cpu_time(time1)
allocate(intensity(reflects%nbragg), stat=ialloc)
! now let's calculate some intensities
do i = 1, reflects%nbragg
  qc_dot_qc = dot_product(reflects%qbragg_list(i)%vec,reflects%qbragg_list(i)%vec)/2.0d0
  costerm = zero
  sinterm = zero
 do j = 1, atoms%natoms
    tfactj  = dexp(-flucts(j)*qc_dot_qc)*real(atoms%atom(j)%occ,kind=8) &
                  * atomic_scatt%atype(atomic_scatt%itype(j))%scatt(i)
    qdotr   =  dot_product(reflects%qbragg_list(i)%vec, real(atoms%atom(j)%x,kind=8))
    costerm = costerm + dcos(qdotr)*tfactj
    sinterm = sinterm + dsin(qdotr)*tfactj
 end do
 intensity(i) = costerm*costerm+sinterm*sinterm
end do

call cpu_time(time2)
print *, 'intensity_zero_iso>', time2-time1, 'seconds'

end subroutine 

subroutine intensity_one(input, atoms,reflects,branches,aniso, &
                           intensity)
! need to bring in the eigvecs/vals and associated qvec... datatype: bybranch_vals_vecs 
! construct the atomic_scatt with wavevector from bybranch 
! 
use cfml_crystal_metrics,        only: Crystal_Cell_Type
use mod_valvec_store
use mod_constants, only: Rkcal,one,two,pi
type(inp_par),              intent(in)      :: input
type(protein_atom_list_type),       intent(in)      :: atoms
type(reflections),          intent(in)      :: reflects
type(valvec),               intent(in)      :: branches
real(dp)        , allocatable,  intent(in)  :: aniso(:,:)
real(dp),allocatable,       intent(out)     :: intensity(:)
type(atscatter)                             :: atomic_scatt
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,hdotr,Qdotu,qc(3)
real(dp) :: val,bigQ(3), qvec(3),scatvec(3), rQdotu, iQdotu,tmp_im,tmp_re 
real(dp) :: time1,time2,kt
integer  :: i,j,k,ialloc,jcor

qvec  = branches%qvec
kT    = input%temperature*Rkcal

call cpu_time(time1)
call atscatt_setup(qvec,atoms,reflects,atomic_scatt)

allocate(intensity(reflects%nbragg), stat=ialloc)
intensity = zero

do i = 1, reflects%nbragg
  qc    = reflects%qbragg_list(i)%vec
  bigQ  = qc-qvec
  do k = input%first_mode, input%last_mode
    if (abs(branches%vals(k)) .gt. 1.0d-08) then
      val = (kt/branches%vals(k))! kt/lambda

      jcor = 1
      costerm = zero
      sinterm = zero
    do j = 1, atoms%natoms
      ! mass term
      jexpon  =-(qc(1)*dot_product(qc,aniso(jcor:jcor+2,1)) + &
                 qc(2)*dot_product(qc,aniso(jcor:jcor+2,2)) + &
                 qc(3)*dot_product(qc,aniso(jcor:jcor+2,3)))/2.0d0

! some thinking required on the occupancy stuff
      tfactj = dexp(jexpon)*real(atoms%atom(j)%occ,kind=8)*&
               atomic_scatt%atype(atomic_scatt%itype(j))%scatt(i)
      hdotr =  dot_product(reflects%qbragg_list(i)%vec, real(atoms%atom(j)%x,kind=8))
      rQdotu = dot_product(bigQ,dble(branches%zvecs(jcor:jcor+2,k))) 
      iQdotu = dot_product(bigQ,dimag(branches%zvecs(jcor:jcor+2,k))) 
      tmp_re = tfactj*(rQdotu*dcos(hdotr) - iQdotu*dsin(hdotr))
      tmp_im = tfactj*(rQdotu*dsin(hdotr) + iQdotu*dcos(hdotr))
      costerm = costerm + tmp_re 
      sinterm = sinterm + tmp_im
      jcor = jcor + 3
    end do
    intensity(i) = intensity(i) + val*(costerm*costerm+sinterm*sinterm)
    end if
  end do
end do

call cpu_time(time2)

print *, 'intensity_one>', time2-time1, 'seconds'

end subroutine 

subroutine intensity_one_test(input, atoms,reflects,branches,aniso, &
                           intensity)
! need to bring in the eigvecs/vals and associated qvec... datatype: bybranch_vals_vecs 
! construct the atomic_scatt with wavevector from bybranch 
! 
use cfml_crystal_metrics,        only: Crystal_Cell_Type
use mod_valvec_store
use mod_constants, only: Rkcal,one,two,pi
type(inp_par),              intent(in)      :: input
type(protein_atom_list_type),       intent(in)      :: atoms
type(reflections),          intent(in)      :: reflects
type(valvec),               intent(in)      :: branches
real(dp)        , allocatable,  intent(in)  :: aniso(:,:)
real(dp),allocatable,       intent(out)     :: intensity(:)
type(atscatter)                             :: atomic_scatt
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,hdotr,Qdotu,qc(3)
real(dp) :: val,bigQ(3), qvec(3),scatvec(3), rQdotu, iQdotu,tmp_im,tmp_re 
real(dp) :: time1,time2,kt
integer  :: i,j,k,ialloc,jcor

qvec  = branches%qvec
kT    = input%temperature*Rkcal

call cpu_time(time1)
call atscatt_setup(qvec,atoms,reflects,atomic_scatt)

allocate(intensity(reflects%nbragg), stat=ialloc)
intensity = zero

do i = 1, reflects%nbragg
  qc    = reflects%qbragg_list(i)%vec
  bigQ  = qc+qvec
  do k = input%first_mode, input%last_mode
    if (abs(branches%vals(k)) .gt. 1.0d-08) then
      val = (kt/branches%vals(k))! kt/lambda

      jcor = 1
      costerm = zero
      sinterm = zero
    do j = 1, atoms%natoms
      ! mass term
      jexpon  =-(bigQ(1)*dot_product(BigQ,aniso(jcor:jcor+2,1)) + &
                 bigQ(2)*dot_product(BigQ,aniso(jcor:jcor+2,2)) + &
                 bigQ(3)*dot_product(BigQ,aniso(jcor:jcor+2,3)))/2.0d0
! some thinking required on the occupancy stuff
      tfactj = dexp(jexpon)*real(atoms%atom(j)%occ,kind=8)*&
               atomic_scatt%atype(atomic_scatt%itype(j))%scatt(i)
      hdotr =  dot_product(qc, real(atoms%atom(j)%x,kind=8))
      rQdotu = dot_product(bigQ,dble(branches%zvecs(jcor:jcor+2,k))) 
      iQdotu = dot_product(bigQ,dimag(branches%zvecs(jcor:jcor+2,k))) 
      tmp_re = tfactj*(rQdotu*dcos(hdotr) + iQdotu*dsin(hdotr))
      tmp_im = tfactj*(iQdotu*dcos(hdotr) - rQdotu*dsin(hdotr))
      costerm = costerm + tmp_re 
      sinterm = sinterm + tmp_im
      jcor = jcor + 3
    end do
    intensity(i) = intensity(i) + val*(costerm*costerm+sinterm*sinterm)
    end if
  end do
end do

call cpu_time(time2)

print *, 'intensity_one>', time2-time1, 'seconds'

end subroutine 

subroutine intensity_one_meinhold(input, atoms,reflects,branches,aniso, &
                           intensity)
! need to bring in the eigvecs/vals and associated qvec... datatype: bybranch_vals_vecs 
! construct the atomic_scatt with wavevector from bybranch 
! 
use cfml_crystal_metrics,        only: Crystal_Cell_Type
use mod_valvec_store
use mod_constants, only: Rkcal,one,two,pi
type(inp_par),              intent(in)      :: input
type(protein_atom_list_type),       intent(in)      :: atoms
type(reflections),          intent(in)      :: reflects
type(valvec),               intent(in)      :: branches
real(dp)        , allocatable,  intent(in)  :: aniso(:,:)
real(dp),allocatable,       intent(out)     :: intensity(:)
type(atscatter)                             :: atomic_scatt
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,hdotr,Qdotu,qc(3)
real(dp) :: val,bigQ(3), qvec(3),scatvec(3), rQdotu, iQdotu,tmp_im,tmp_re 
real(dp) :: time1,time2,kt
integer  :: i,j,k,ialloc,jcor

qvec  = branches%qvec
kT    = input%temperature*Rkcal

call cpu_time(time1)
call atscatt_setup(qvec,atoms,reflects,atomic_scatt)

allocate(intensity(reflects%nbragg), stat=ialloc)
intensity = zero

do i = 1, reflects%nbragg
  qc    = reflects%qbragg_list(i)%vec
  bigQ  = qc-qvec
  do k = input%first_mode, input%last_mode
    if (abs(branches%vals(k)) .gt. 1.0d-08) then
      val = (kt/branches%vals(k))! kt/lambda

      jcor = 1
      costerm = zero
      sinterm = zero
    do j = 1, atoms%natoms
      ! mass term
      jexpon  =-(qc(1)*dot_product(qc,aniso(jcor:jcor+2,1)) + &
                 qc(2)*dot_product(qc,aniso(jcor:jcor+2,2)) + &
                 qc(3)*dot_product(qc,aniso(jcor:jcor+2,3)))/2.0d0

! some thinking required on the occupancy stuff
      tfactj = dexp(jexpon)*real(atoms%atom(j)%occ,kind=8)*&
               atomic_scatt%atype(atomic_scatt%itype(j))%scatt(i)
      !hdotr =  dot_product(reflects%qbragg_list(i)%vec, real(atoms%atom(j)%x,kind=8))
      hdotr =  dot_product(bigQ,real(atoms%atom(j)%x,kind=8))
      rQdotu = dot_product(bigQ,dble(branches%zvecs(jcor:jcor+2,k))) 
      iQdotu = dot_product(bigQ,dimag(branches%zvecs(jcor:jcor+2,k))) 
      tmp_re = tfactj*(rQdotu*dcos(hdotr) + iQdotu*dsin(hdotr))
      tmp_im = tfactj*(iQdotu*dcos(hdotr) - rQdotu*dsin(hdotr))
      costerm = costerm + tmp_re 
      sinterm = sinterm + tmp_im
      jcor = jcor + 3
    end do
    intensity(i) = intensity(i) + val*(costerm*costerm+sinterm*sinterm)
    end if
  end do
end do

call cpu_time(time2)

print *, 'intensity_one>', time2-time1, 'seconds'

end subroutine 

subroutine intensity_one_testold(input, atoms,reflects,branches,aniso, &
                           intensity)
! need to bring in the eigvecs/vals and associated qvec... datatype: bybranch_vals_vecs 
! construct the atomic_scatt with wavevector from bybranch 
! 
use cfml_crystal_metrics,        only: Crystal_Cell_Type
use mod_valvec_store
use mod_constants, only: Rkcal,one,two,pi
type(inp_par),              intent(in)      :: input
type(protein_atom_list_type),       intent(in)      :: atoms
type(reflections),          intent(in)      :: reflects
type(valvec),               intent(in)      :: branches
real(dp)        , allocatable,  intent(in)  :: aniso(:,:)
real(dp),allocatable,       intent(out)     :: intensity(:)
type(atscatter)                             :: atomic_scatt
real(dp) :: qc_dot_qc,tfactj,jexpon,costerm,sinterm,hdotr,Qdotu,qc(3)
real(dp) :: val,bigQ(3), qvec(3),rscatvec(3),iscatvec(3), rQdotu, iQdotu,tmp_im,tmp_re 
real(dp) :: time1,time2,kt,tval,itval
integer  :: i,j,k,ialloc,jcor

qvec  = branches%qvec
kT    = input%temperature*Rkcal

call cpu_time(time1)
call atscatt_setup(qvec,atoms,reflects,atomic_scatt)

allocate(intensity(reflects%nbragg), stat=ialloc)
intensity = zero

do i = 1, reflects%nbragg
  qc    = two*pi*reflects%qbragg_list(i)%vec
  bigQ  = qc-two*pi*qvec
  do k = input%first_mode, input%last_mode
    if (abs(branches%vals(k)) .gt. 1.0d-08) then
      val = (kt/branches%vals(k))! kt/lambda

      costerm = zero
      sinterm = zero
      rscatvec = zero
      iscatvec = zero
      tval  = zero
      itval = zero
      jcor = 1
    do j = 1, atoms%natoms
      ! mass term
      !jexpon  =-(qc(1)*dot_product(qc,aniso(jcor:jcor+2,1)) + &
      !           qc(2)*dot_product(qc,aniso(jcor:jcor+2,2)) + &
      !           qc(3)*dot_product(qc,aniso(jcor:jcor+2,3)))/2.0d0

! some thinking required on the occupancy stuff
      tfactj = dexp(jexpon)*real(atoms%atom(j)%occ,kind=8)*&
                atomic_scatt%atype(atomic_scatt%itype(j))%scatt(i)
      hdotr =  dot_product(qc, real(atoms%atom(j)%x,kind=8))

      rQdotu = dot_product(bigQ,dble(branches%zvecs(jcor:jcor+2,k))) 
      iQdotu = dot_product(bigQ,dimag(branches%zvecs(jcor:jcor+2,k))) 
      tmp_re = tfactj*(rQdotu*dcos(hdotr) - iQdotu*dsin(hdotr))
      tmp_im = tfactj*(rQdotu*dsin(hdotr) + iQdotu*dcos(hdotr))
      costerm = costerm + tmp_re 
      sinterm = sinterm + tmp_im
      jcor = jcor + 3
    end do
    !print '(A4,2I4,2F10.4)', 'shit',i,k,val,rscatvec
    !intensity(i) = intensity(i) + val*((dot_product(rscatvec,rscatvec))**2 + (dot_product(iscatvec,iscatvec))**2)
    intensity(i) = intensity(i) + val*(costerm*costerm+sinterm*sinterm)
    end if
  end do
    print '(A8,1x,I4,4F16.5)', 'shitshit',i,val,intensity(i), costerm,sinterm
end do

call cpu_time(time2)

print *, 'intensity_one>', time2-time1, 'seconds'

end subroutine 

real(dp) function rval(intens1,intens2)
real(dp), allocatable,dimension(:), intent(in) :: intens1, intens2

rval = sum(abs(intens1-intens2))/sum(intens1) 

end function

end module 

