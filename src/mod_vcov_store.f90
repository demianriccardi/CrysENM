module mod_vcov_store
! contain the subroutine to make bfactor comparisons 
! this code needs to be cleaned up and commented
use cfml_globaldeps,                 only: dp,sp
use mod_constants
use mod_types,                        only: protein_atom_list_type

implicit none

! DMR 01302009 vcov_store is what is used now.  the byatom storage is not fast enough
type :: vcov_store
!DMR ok... new vcov_storage needed.  Damn-atom based storage was too slow when the 2D matrix of 3x3 blocks was used
  integer :: natoms,dzeros,zzeros
  !                       dimension = nat ; 3nat,3;    nat,nat;       3nat,3nat
  real(dp), allocatable        ::    iso(:),anisotropy(:),aniso(:,:), isovcov(:,:) , vcov(:,:)
  complex(kind=8), allocatable ::    zaniso(:,:), zvcov(:,:)
  real(dp) :: kt 
end type

!type :: vcov_entry
! store a block of the vcov matrix as well as the trace/3 of that block
!  real(dp), pointer            :: iso                            
!  real(dp), allocatable        :: vcov(:,:)     ! anisotropic(6), vcov(3,3)
!  complex(kind=8),allocatable  :: zvcov(:,:)    ! zvcov(3,3)
!end type

public ::  vcov_writer,svcov_massadjust,vcov_store,svcov_x_scalar,dvec_store_svcov,zvec_store_svcov, &
          zvcov_store_svcov,vcov_store_svcov,svcov_iso,svcov_aniso,svcov_isovcov,valvecs_svcov

contains

subroutine vcov_writer(svcov)
type(vcov_store)    , intent(in) :: svcov
integer                          :: i,j

do i =1, 3*svcov%natoms
  do j = i, 3*svcov%natoms
    write(666,*) i, j, svcov%vcov(i,j)
  end do
end do

end subroutine

subroutine p1vcov_calculator(input,cell,spg,asym_un,asyms,unit_cell,svcov)
! so ugly!!!
use cfml_crystal_metrics,            only: Crystal_Cell_Type
use cfml_crystallographic_symmetry,only: space_group_type
use mod_types,                only : inp_par, asym_list_type,dos,sparse
use mod_hessian,              only : kirchess_new,kirchess_run_new,kirchess_valvec,dvcov
use mod_correl_dist,            only: corrtab_dexp
use mod_tls
use mod_bnm
use mod_crysbuild, only: atom_shrinker,asym_shrinker
use mod_tempfact, only:vcov_Bscale
use mod_symop,                only : vcov_p1,atcor_p1
use mod_lattdyn,              only : bnm_vcov,qdispers
use mod_inout,                only : vcov_atshrink
type(inp_par) ,       intent(inout)  :: input
type(crystal_cell_type), intent(in)  :: cell
type(space_group_type),  intent(in)  :: spg
type(protein_atom_list_type), intent(in out)  :: asym_un
type(asym_list_type), intent(in out)  :: asyms,unit_cell
type(vcov_store)    , intent(out) :: svcov
type(protein_atom_list_type)                  :: tmpats
type(asym_list_type)                  :: subun_cell
real(dp), allocatable :: vcov(:,:),tmp_vcov(:,:)
type(sparse) :: kirchoff,hessian
type(dos) :: states
integer   :: ialloc,i,j,nviab
real(dp) :: correl,fcnst

! TLS
if (trim(input%genm) .eq. "tls") then
  call tls_vcov(input,asym_un,tmp_vcov)
  !do i =1, asym_un%natoms
  !  j = 3*(i-1)+1
  !  write(777,*) i,asym_un%atom(i)%biso,(tmp_vcov(j,j)+tmp_vcov(j+1,j+1)+tmp_vcov(j+2,j+2))/3.0d0
  !end do
  tmp_vcov = tmp_vcov*b2u
! ISOCORR
else if (trim(input%genm) .eq. "isocorr") then
  !do something that is somewhere
  if  (trim(input%tpbc) .eq. "asy") then
    print *, "sending in the asyms"
    call corrtab_dexp(input,asyms,tmp_vcov)
  else if ( (trim(input%tpbc) .eq. "pbc") .or. &
            (trim(input%tpbc) .eq. "bvk") ) then 
    call corrtab_dexp(input,unit_cell,vcov)
    print *, "sending in the asyms"
  else
    print *, "sending in the atoms"
    call corrtab_dexp(input,asym_un,tmp_vcov)
  end if
! GNM,ENM,etc
else
  !      print *, 'tshit',trim(input%genm), "one of these print things!"
  call kirchess_new(input,asyms,unit_cell,kirchoff,hessian)
  !      print *, 'tshit',trim(input%genm)
  if ((trim(input%tpbc) .eq. "bvk") .or. (trim(input%tpbc) .eq. "pbc")) then
    if (trim(input%genm) .eq. "bnm") then
      print *, 'going into bnm_vcov'
      !call bnm_vcov(input,cell,unit_cell%aunit(1),unit_cell,kirchoff,hessian,vcov)
      call bnm_vcov(input,cell,asyms,unit_cell,kirchoff,hessian,vcov)

    else 
      if (trim(input%tpbc) .eq. "bvk") then
        call qdispers(input,cell,unit_cell,kirchoff,hessian,tmp_vcov,states)
      else
        print *, 'shitshit'
        call kirchess_run_new(input,unit_cell,kirchoff,hessian,&
                          tmp_vcov,states)
        print *, 'shitshit'
      end if
!     get only the CA atoms
      call vcov_atshrink(input,unit_cell%aunit(1),tmp_vcov,vcov)
      deallocate(tmp_vcov)
    end if
  else
    if (trim(input%genm) .eq. "bnm") then
      !call bnm_vcov(input,cell,asyms%aunit(1),unit_cell,kirchoff,hessian,tmp_vcov)
      call bnm_vcov(input,cell,asyms,unit_cell,kirchoff,hessian,tmp_vcov)
    else
      call kirchess_run_new(input,unit_cell,kirchoff,hessian,&
                          vcov,states)
      call vcov_atshrink(input,asyms%aunit(1),vcov,tmp_vcov)
      deallocate(vcov)
    end if
  end if
end if

  !do i =1, asym_un%natoms
  !  j = 3*(i-1)+1
  !  if (allocated(vcov))     write(777,*) i,asym_un%atom(i)%biso,u2b*(vcov(j,j)+vcov(j+1,j+1)+vcov(j+2,j+2))/3.0d0
  !  if (allocated(tmp_vcov)) write(778,*) i,asym_un%atom(i)%biso,u2b*(tmp_vcov(j,j)+tmp_vcov(j+1,j+1)+tmp_vcov(j+2,j+2))/3.0d0
  !end do

if (allocated(tmp_vcov)) then
  if ((trim(input%genm) .ne. "gnm") .and. (trim(input%genm) .ne. "isocorr")) then
    print *, 'calling vcov_p1'
    call vcov_p1(cell,spg,tmp_vcov,vcov)
  else
    call atcor_p1(spg,tmp_vcov,vcov)
  end if
  deallocate(tmp_vcov)
end if

  !do i =1, asym_un%natoms
  !  j = 3*(i-1)+1
  !  write(779,*) i,asym_un%atom(i)%biso,u2b*(vcov(j,j)+vcov(j+1,j+1)+vcov(j+2,j+2))/3.0d0
  !end do

call atom_shrinker(input,asym_un,tmpats)
call asym_shrinker(input,unit_cell,subun_cell)
input%natoms = tmpats%natoms*asyms%naunits
call vcov_Bscale(input,tmpats,vcov,correl,fcnst,nviab)
input%nviab = nviab

print *, "all clear for svcov!"

!  vcov = vcov*u2b

print *, "all clear for svcov!"

if ((trim(input%genm) .ne. "gnm") .and. (trim(input%genm) .ne. "isocorr")) then
  call vcov_store_svcov(input,vcov,svcov)
  call svcov_aniso(svcov)
  call svcov_isovcov(svcov)
  if (trim(input%isovcov_scale) .eq. "bfact") call atcorr_mat(subun_cell%aunit(1),svcov%isovcov)
  call svcov_iso(svcov)
else
  call atcor_store_svcov(input,vcov,svcov)
  if (trim(input%isovcov_scale) .eq. "bfact") call atcorr_mat(subun_cell%aunit(1),svcov%isovcov)
  call svcov_iso(svcov)
end if

!do i = 1, size(svcov%isovcov,1)
!  do j = i, size(svcov%isovcov,1)
!    if(svcov%isovcov(i,j) .ne. zero) write(655,*) i,j,svcov%isovcov(i,j)
!  end do
!end do

end subroutine

subroutine p1vcov_bigcalc(input,cell,spg,asym_un,asyms,unit_cell,svcov,blocks)
! so ugly!!!
use cfml_crystal_metrics,            only: Crystal_Cell_Type
use cfml_crystallographic_symmetry,only: space_group_type
use mod_types,                only : inp_par, asym_list_type,dos,sparse
use mod_hessian,              only : kirchess_new,kirchess_run_new,kirchess_valvec,dvcov
use mod_correl_dist,            only: corrtab_dexp
use mod_tls
use mod_bnm
use mod_tempfact, only:vcov_Bscale
use mod_symop,                only : vcov_p1,atcor_p1
use mod_lattdyn,              only : bnm_bigvcov,qdispers
type(inp_par) ,       intent(inout)  :: input
type(crystal_cell_type), intent(in)  :: cell
type(space_group_type),  intent(in)  :: spg
type(protein_atom_list_type), intent(in out)  :: asym_un
type(asym_list_type), intent(in out)  :: asyms,unit_cell
type(vcov_store)    , intent(out) :: svcov
type(block),allocatable,optional, intent(out) :: blocks(:)
real(dp), allocatable :: vcov(:,:),tmp_vcov(:,:)
type(sparse) :: kirchoff,hessian
type(dos) :: states
integer   :: ialloc,i,j,nviab
real(dp) :: correl,fcnst

! TLS
if (trim(input%genm) .eq. "tls") then
  call tls_vcov(input,asym_un,tmp_vcov)
  !do i =1, asym_un%natoms
  !  j = 3*(i-1)+1
  !  write(777,*) i,asym_un%atom(i)%biso,(tmp_vcov(j,j)+tmp_vcov(j+1,j+1)+tmp_vcov(j+2,j+2))/3.0d0
  !end do
  tmp_vcov = tmp_vcov*b2u
! ISOCORR
else if (trim(input%genm) .eq. "isocorr") then
  !do something that is somewhere
  if  (trim(input%tpbc) .eq. "asy") then
    print *, "sending in the asyms"
    call corrtab_dexp(input,asyms,tmp_vcov)
  else if ( (trim(input%tpbc) .eq. "pbc") .or. &
            (trim(input%tpbc) .eq. "bvk") ) then 
    call corrtab_dexp(input,unit_cell,vcov)
    print *, "sending in the asyms"
  else
    print *, "sending in the atoms"
    call corrtab_dexp(input,asym_un,tmp_vcov)
  end if
! GNM,ENM,etc
else
  !      print *, 'tshit',trim(input%genm), "one of these print things!"
  call kirchess_new(input,asyms,unit_cell,kirchoff,hessian)
  !      print *, 'tshit',trim(input%genm)
  if ((trim(input%tpbc) .eq. "bvk") .or. (trim(input%tpbc) .eq. "pbc")) then
    if (trim(input%genm) .eq. "bnm") then
      print *, 'going into bnm_vcov'
     ! call bnm_bigvcov(input,cell,asymsunit_cell%aunit(1),unit_cell,kirchoff,hessian,vcov)
      call bnm_bigvcov(input,cell,asyms,unit_cell,kirchoff,hessian,vcov,states,blocks)

    else 
      if (trim(input%tpbc) .eq. "bvk") then
        call qdispers(input,cell,unit_cell,kirchoff,hessian,tmp_vcov,states)
      else
        call kirchess_run_new(input,unit_cell,kirchoff,hessian,&
                          tmp_vcov,states)
      end if
    end if
  else
    if (trim(input%genm) .eq. "bnm") then
      !call bnm_bigvcov(input,cell,asyms%aunit(1),unit_cell,kirchoff,hessian,tmp_vcov)
      call bnm_bigvcov(input,cell,asyms,unit_cell,kirchoff,hessian,tmp_vcov,states,blocks)
    else
      call kirchess_run_new(input,unit_cell,kirchoff,hessian,&
                          tmp_vcov,states)
    end if
  end if
end if

print *, 'shitgoingwell'

if (allocated(tmp_vcov)) then
  if ((trim(input%genm) .ne. "gnm") .and. (trim(input%genm) .ne. "isocorr")) then
    print *, 'calling vcov_p1'
    call vcov_p1(cell,spg,tmp_vcov,vcov)
  else
    call atcor_p1(spg,tmp_vcov,vcov)
  end if
  deallocate(tmp_vcov)
end if

input%natoms = unit_cell%aunit(1)%natoms
print *, 'shitgoingwell',input%natoms
if (input%fcnst .eq. one) then
print *, 'shitshitgwell',input%natoms
  call vcov_Bscale(input,asyms%aunit(1),vcov,correl,fcnst,nviab)
print *, 'shitshitgwell',input%natoms
end if

input%nviab = nviab

print *, "all clear for svcov!",(vcov(1,1)+vcov(2,2)+vcov(3,3)/3.0d0)

!vcov = vcov*u2b

print *, "all clear for svcov!"

if ((trim(input%genm) .ne. "gnm") .and. (trim(input%genm) .ne. "isocorr")) then
  call vcov_store_svcov(input,vcov,svcov)
  call svcov_aniso(svcov)
  call svcov_isovcov(svcov)
  if (trim(input%isovcov_scale) .eq. "bfact") call atcorr_mat(unit_cell%aunit(1),svcov%isovcov)
  call svcov_iso(svcov)
else
  call atcor_store_svcov(input,vcov,svcov)
  if (trim(input%isovcov_scale) .eq. "bfact") call atcorr_mat(unit_cell%aunit(1),svcov%isovcov)
  call svcov_iso(svcov)
end if

if(present(blocks)) then

  print *, "p1vcov_bigcalc> setting up block info from  fileroot.block file for p1 unit cell"   
  if (allocated(blocks)) deallocate(blocks) 
  input%nblocks = unit_cell%nblocks
  call blocks_setup (input,unit_cell%aunit(1),blocks)
  print *, "p1vcov_bigcalc> blocks generated"   

end if

!do i = 1, size(svcov%isovcov,1)
!  do j = i, size(svcov%isovcov,1)
!    if(svcov%isovcov(i,j) .ne. zero) write(655,*) i,j,svcov%isovcov(i,j)
!  end do
!end do

end subroutine

subroutine vcov_calc(input,cell,asym_un,asyms,unit_cell,svcov,states,blocks)
use cfml_crystal_metrics,            only: Crystal_Cell_Type
use mod_types,                only : inp_par, asym_list_type,dos,sparse
use mod_hessian,              only : kirchess_new,kirchess_run_new,kirchess_valvec,dvcov
use mod_correl_dist,            only: corrtab_dexp
use mod_tls
use mod_bnm
use mod_lattdyn,              only : bnm_bigvcov,qdispers
type(inp_par) ,       intent(in out)   :: input
type(crystal_cell_type) ,    intent(in) :: cell
type(protein_atom_list_type), intent(in out) :: asym_un
type(asym_list_type), intent(in out) :: asyms,unit_cell
type(vcov_store)    , intent(out)    :: svcov
type(dos),            intent(out)    :: states
type(block),allocatable,optional, intent(out) :: blocks(:)
real(dp),  allocatable :: vcov(:,:)
type(sparse) :: kirchoff,hessian
integer   :: ialloc,i,j,nviab,icor
real(dp) :: correl,fcnst

! TLS
if (trim(input%genm) .eq. "tls") then
  input%natoms = asym_un%natoms
  print *, 'number of atoms in input', input%natoms
  call tls_vcov(input,asym_un,vcov)
  vcov = vcov*b2u
print *, 'number of atoms in input', input%natoms
! ISOCORR
else if (trim(input%genm) .eq. "isocorr") then
  !do something that is somewhere
  if  (trim(input%tpbc) .eq. "asy") then
    print *, "sending in the asyms"
    call corrtab_dexp(input,asyms,vcov)
  else if ( (trim(input%tpbc) .eq. "pbc") .or. &
            (trim(input%tpbc) .eq. "bvk") ) then 
    call corrtab_dexp(input,unit_cell,vcov)
    print *, "sending in the asyms"
  else
    print *, "sending in the atoms"
    call corrtab_dexp(input,asym_un,vcov)
  end if
! GNM,ENM,etc
else
  !      print *, 'tshit',trim(input%genm), "one of these print things!"
  call kirchess_new(input,asyms,unit_cell,kirchoff,hessian)
  !      print *, 'tshit',trim(input%genm)
  if ((trim(input%tpbc) .eq. "bvk") .or. (trim(input%tpbc) .eq. "pbc")) then
    if (trim(input%genm) .eq. "bnm") then
      print *, 'going into bnm_vcov'
      call bnm_bigvcov(input,cell,asyms,unit_cell,kirchoff,hessian,vcov,states)
      !call bnm_bigvcov(input,cell,unit_cell%aunit(1),unit_cell,kirchoff,hessian,vcov)

    else 
      if (trim(input%tpbc) .eq. "bvk") then
        call qdispers(input,cell,unit_cell,kirchoff,hessian,vcov,states)
      else
        call kirchess_run_new(input,unit_cell,kirchoff,hessian,&
                          vcov,states)
      end if
    end if
  else
    if (trim(input%genm) .eq. "bnm") then
      !call bnm_bigvcov(input,cell,asyms%aunit(1),unit_cell,kirchoff,hessian,vcov)
      call bnm_bigvcov(input,cell,asyms,unit_cell,kirchoff,hessian,vcov,states)
    else
      call kirchess_run_new(input,unit_cell,kirchoff,hessian,&
                          vcov,states)
    end if
  end if
end if

if ((trim(input%genm) .ne. "gnm") .and. (trim(input%genm) .ne. "isocorr")) then
  call vcov_store_svcov(input,vcov,svcov)
  call svcov_aniso(svcov)
  call svcov_isovcov(svcov)
  call svcov_iso(svcov)
  call svcov_anisotropy_calc(svcov)
else
  call atcor_store_svcov(input,vcov,svcov)
  call svcov_iso(svcov)
end if

if(present(blocks)) then

  print *, "vcov_calc> setting up block info from  fileroot.block file"
  if (allocated(blocks)) deallocate(blocks)
  !input%nblocks = unit_cell%nblocks       ! for p1vcov
  if(trim(input%tpbc) .eq. "pbc" .or. trim(input%tpbc) .eq. "bvk") then
    call blocks_setup (input,unit_cell%aunit(1),blocks)
  else
    call blocks_setup (input,asyms%aunit(1),blocks)
  end if
  print *, "vcov_calc> blocks generated"

end if


end subroutine

subroutine svcov_x_scalar(scalar,svcov)
!takes in complex eigenvectors
real(dp),                     intent(in) :: scalar
type(vcov_store),            intent(in out) :: svcov

if (allocated(svcov%iso))     svcov%iso     =scalar * svcov%iso
if (allocated(svcov%aniso))   svcov%aniso   =scalar * svcov%aniso
if (allocated(svcov%isovcov)) svcov%isovcov =scalar * svcov%isovcov
if (allocated(svcov%vcov))    svcov%vcov    =scalar * svcov%vcov
if (allocated(svcov%zvcov))   svcov%zvcov   =scalar * svcov%zvcov
if (allocated(svcov%zaniso))  svcov%zaniso  =scalar * svcov%zaniso

end subroutine

subroutine svcov_massadjust(input,atoms,svcov)
use mod_crysbuild, only: invrt_rtmass_vector 
use mod_types,     only: inp_par 
type(inp_par),               intent(in)     :: input
type(protein_atom_list_type),        intent(in)     :: atoms
type(vcov_store),            intent(in out) :: svcov
real(dp), allocatable                       :: thr_rtm(:), one_rtm(:)
integer :: i,j,icor,jcor,ii,jj

if (trim(input%weigh) .eq. "mass" ) then
  call invrt_rtmass_vector(atoms, one_rtm)

icor = 1
do i = 1, svcov%natoms
  if (allocated(svcov%iso))    svcov%iso(i)      = one_rtm(i) * svcov%iso(i)      * one_rtm(i)
  do ii = 0,2
      if (allocated(svcov%aniso))  svcov%aniso(icor+ii,:)  = one_rtm(i) * svcov%aniso(icor+ii,:)  * one_rtm(i)
      if (allocated(svcov%zaniso)) svcov%zaniso(icor+ii,:) = one_rtm(i) * svcov%zaniso(icor+ii,:) * one_rtm(i)
  end do
  icor = icor + 3
end do

icor = 1
do i = 1, svcov%natoms
  jcor=1
  do j = 1, svcov%natoms
    if (allocated(svcov%isovcov))   svcov%isovcov(i,j)  = one_rtm(i) * svcov%isovcov(i,j) * one_rtm(j)
    do ii = 0,2
      do jj = 0,2
        if (allocated(svcov%vcov))    svcov%vcov(icor+ii,jcor+jj)  = one_rtm(i) &
                                             * svcov%vcov(icor+ii,jcor+jj) * one_rtm(j)
        if (allocated(svcov%zvcov))   svcov%zvcov(icor+ii,jcor+jj) = one_rtm(i) &
                                             * svcov%zvcov(icor+ii,jcor+jj) * one_rtm(j)
      end do
    end do
    if (allocated(svcov%isovcov)) svcov%isovcov(i,j)  = one_rtm(i) * svcov%isovcov(i,j) * one_rtm(j)
    jcor = jcor + 3
  end do
  icor = icor + 3
end do

else
  return
end if

end subroutine

subroutine dvec_store_svcov(input,vals,vecs,svcov)
use mod_types,     only: inp_par
use mod_constants, only: Rkcal
! take in the eigenvalues: units rad^2 kcal mol^-1 amu^-1 angs^2
! 
type(inp_par) :: input
real(dp),        allocatable,intent(in)     :: vals(:),vecs(:,:)
type(vcov_store),            intent(in out) :: svcov
real(dp)        :: tmpval
real(dp)        :: kT
integer ::   nvals, ndim
integer ::   i,j,k
logical :: tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zvcov,tf_zaniso
logical, parameter :: T=.true., F=.false.

tf_iso     = F 
tf_aniso   = F
tf_isovcov = F
tf_vcov    = T
tf_zvcov   = F
tf_zaniso  = F

kT = input%temperature*Rkcal
svcov%kt = kt
call vcov_store_init(input%natoms,tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zaniso,tf_zvcov, svcov)

nvals  = size(vals)
ndim   = size(vecs(:,1))
k = ndim/3
if (input%natoms .ne. k) print *, "dvec_store_svcov> warning, something weird with the number of atoms", input%natoms, k
svcov%dzeros = 0
do i=input%first_mode,input%last_mode
  if (abs(vals(i)).gt.1.0D-08) then
    do j=1,ndim
      do k=1,ndim
        tmpval    = vecs(k,i)*vecs(j,i)
        svcov%vcov(k,j) = svcov%vcov(k,j) + tmpval*kt/vals(i)
      end do
    end do
  else
    svcov%dzeros = svcov%dzeros + 1
  end if
end do

end subroutine

subroutine valvecs_svcov(input,valvecs,svcov)
use mod_types,     only: inp_par
use mod_constants, only: Rkcal
use mod_valvec_store
type(inp_par) :: input
type(valvec),    allocatable,intent(in)     :: valvecs(:)
type(vcov_store),            intent(out) :: svcov
real (dp) :: kt,tmpval,summult
logical, parameter :: T=.true., F=.false.
logical :: tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zvcov,tf_zaniso
integer :: i,j,k,l,m,n,ndim

tf_iso     = F 
tf_aniso   = F
tf_isovcov = F
tf_vcov    = T
tf_zvcov   = F
tf_zaniso  = F

kT = input%temperature*Rkcal
svcov%kt = kt

call vcov_store_init(input%natoms,tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zaniso,tf_zvcov, svcov)

ndim   = size(valvecs(1)%zvecs(:,1))
svcov%zzeros = 0

do l = 1, size(valvecs)
  summult = summult + valvecs(l)%multiplicity
  do i=input%first_mode,input%last_mode
    if (abs(valvecs(l)%vals(i)).gt.1.0D-08) then
      do j=1,ndim
        do k=1,ndim
          tmpval    = dble(valvecs(l)%zvecs(k,i)*conjg(valvecs(l)%zvecs(j,i)))
          svcov%vcov(k,j) = svcov%vcov(k,j) + valvecs(l)%multiplicity*tmpval*kt/valvecs(l)%vals(i)
        end do
      end do
    else
      svcov%zzeros = svcov%zzeros + 1
    end if
  end do
end do

! note the use of multiplicity above.  If it's a nonzero wavevector, we multiply by two.  if zero, mult by 1
print *, 'total number of wavevectors:', summult
svcov%vcov = svcov%vcov/summult 

call svcov_aniso(svcov)

end subroutine

subroutine zvec_store_svcov(input,vals,vecs,svcov)
use mod_types,     only: inp_par
use mod_constants, only: Rkcal
! take in the eigenvalues: units rad^2 kcal mol^-1 amu^-1 angs^2
! 
type(inp_par) :: input
real(dp),        allocatable,intent(in)     :: vals(:)
complex(kind=8), allocatable,intent(in)     :: vecs(:,:)
type(vcov_store),            intent(in out) :: svcov
complex(kind=8) :: tmpval
real(dp)        :: kT
integer ::   nvals, ndim
integer ::   i,j,k
logical :: tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zvcov,tf_zaniso
logical, parameter :: T=.true., F=.false.

tf_iso     = F 
tf_aniso   = F
tf_isovcov = F
tf_vcov    = F
tf_zvcov   = T
tf_zaniso  = F

kT = input%temperature*Rkcal
svcov%kt = kt
call vcov_store_init(input%natoms,tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zaniso,tf_zvcov, svcov)

nvals  = size(vals)
ndim   = size(vecs(:,1))
k = ndim/3
if (input%natoms .ne. k) print *, "dvec_store_svcov> warning, something weird with the number of atoms",input%natoms, k

svcov%zzeros = 0
do i=input%first_mode,input%last_mode
  if (abs(vals(i)).gt.1.0D-08) then
    do j=1,ndim
      do k=1,ndim
        tmpval    = vecs(k,i)*conjg(vecs(j,i))
        svcov%zvcov(k,j) = svcov%zvcov(k,j) + tmpval*kt/vals(i)
      end do
    end do
  else
    svcov%zzeros = svcov%zzeros + 1
  end if
end do

end subroutine

subroutine zvcov_store_svcov(input,vcov,svcov)
use mod_types,     only: inp_par
use mod_constants, only: Rkcal
type(inp_par) :: input
real(dp),        allocatable,intent(in)     :: vcov(:,:)
type(vcov_store),            intent(in out) :: svcov
logical :: tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zvcov,tf_zaniso
logical, parameter :: T=.true., F=.false.
real(dp) :: kt

tf_iso     = F 
tf_aniso   = F
tf_isovcov = F
tf_vcov    = T
tf_zvcov   = F
tf_zaniso  = F
call vcov_store_init(input%natoms,tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zaniso,tf_zvcov, svcov)

kT = input%temperature*Rkcal
svcov%kt = kt
svcov%zvcov = vcov

end subroutine

subroutine z2dvcov_store_svcov(svcov)
type(vcov_store),            intent(in out) :: svcov
logical :: tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zvcov,tf_zaniso
logical, parameter :: T=.true., F=.false.
real(dp) :: kt

tf_iso     = F 
tf_aniso   = F
tf_isovcov = F
tf_vcov    = T
tf_zvcov   = F
tf_zaniso  = F
call vcov_store_init(svcov%natoms,tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zaniso,tf_zvcov, svcov)

svcov%vcov = dble(svcov%zvcov)

!deallocate(svcov%zvcov)

end subroutine

subroutine vcov_store_svcov(input,vcov,svcov)
use mod_types,     only: inp_par
use mod_constants, only: Rkcal
type(inp_par) :: input
real(dp),        allocatable,intent(in)     :: vcov(:,:)
type(vcov_store),            intent(in out) :: svcov
logical :: tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zvcov,tf_zaniso
logical, parameter :: T=.true., F=.false.
real(dp) :: kt

tf_iso     = F 
tf_aniso   = F
tf_isovcov = F
tf_vcov    = T
tf_zvcov   = F
tf_zaniso  = F

call vcov_store_init(input%natoms,tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zaniso,tf_zvcov, svcov)
kT = input%temperature*Rkcal
svcov%kt = kt
svcov%vcov = vcov
svcov%natoms = input%natoms

end subroutine

subroutine atcor_store_svcov(input,atcor,svcov)
use mod_types,     only: inp_par
use mod_constants, only: Rkcal
type(inp_par) :: input
real(dp),        allocatable,intent(in)     :: atcor(:,:)
type(vcov_store),            intent(in out) :: svcov
logical :: tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zvcov,tf_zaniso
logical, parameter :: T=.true., F=.false.
real(dp) :: kt

tf_iso     = F 
tf_aniso   = F
tf_isovcov = T
tf_vcov    = F
tf_zvcov   = F
tf_zaniso  = F

call vcov_store_init(input%natoms,tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zaniso,tf_zvcov, svcov)
kT = input%temperature*Rkcal
svcov%kt = kt
svcov%isovcov = atcor
svcov%natoms = input%natoms

end subroutine

subroutine svcov_iso(svcov)
type(vcov_store), intent(in out) :: svcov
integer :: i,icor
logical :: tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zvcov,tf_zaniso
logical, parameter :: T=.true., F=.false.

tf_iso     = T
tf_aniso   = F
tf_isovcov = F
tf_vcov    = F
tf_zvcov   = F
tf_zaniso  = F

call vcov_store_init(svcov%natoms,tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zaniso,tf_zvcov, svcov)
if (allocated(svcov%isovcov)) then
  print *, "taking isotropic temp factors from isovcov"
  do i = 1,svcov%natoms
    svcov%iso(i)=svcov%isovcov(i,i)
  end do
return
end if

icor = 1
if (allocated(svcov%vcov)) then
  do i=1, svcov%natoms
    svcov%iso(i) = (svcov%vcov(icor,icor)+svcov%vcov(icor+1,icor+1)+svcov%vcov(icor+2,icor+2))/3.0d0
    icor = icor + 3
  end do
  return
else if (allocated(svcov%zvcov)) then
  do i=1, svcov%natoms
    svcov%iso(i) = dble((svcov%zvcov(icor,icor)+svcov%zvcov(icor+1,icor+1)+svcov%zvcov(icor+2,icor+2))/3.0d0)
    icor = icor + 3
  end do
  return
end if

end subroutine

subroutine atcorr_mat(atoms,matrx)
! take natom x natom correlation matrx and 
type(protein_atom_list_type), intent(in) :: atoms
real(dp), allocatable, intent(inout) :: matrx(:,:)
real(dp), allocatable :: tmpmat(:,:)
integer :: i,j
real(dp) :: episqr

print *, 'atcorr_mat has been called. the isovcov will be normalized by the experimental bfactors'
print *, 'make sure this is what you want'
allocate(tmpmat(atoms%natoms,atoms%natoms), stat=i)
tmpmat = 0.0d0

  episqr = 8.0d0*pi*pi
  do i = 1,atoms%natoms
    do j = 1,atoms%natoms
! make sure we're in angstroms!
      tmpmat(i,j) = real(sqrt(atoms%atom(i)%biso*atoms%atom(j)%biso)/episqr,kind=dp)* &
                           matrx(i,j)/dsqrt(matrx(i,i)*matrx(j,j))
    end do
  end do

matrx = tmpmat

end subroutine

subroutine atcorr_mat_isomult(iso,matrx)
! take natom x natom correlation matrx and 
real(dp), allocatable,intent(in) :: iso(:)
real(dp), allocatable, intent(inout) :: matrx(:,:)
real(dp), allocatable :: tmpmat(:,:)
integer :: i,j
real(dp) :: episqr

print *, 'atcorr_mat_isomat>  Renormalizing by Bfactors read in'
allocate(tmpmat(size(iso),size(iso)), stat=i)
tmpmat = 0.0d0


  do i = 1,size(matrx,1)
    do j = 1,size(matrx,2)
! make sure we're in angstroms!
      tmpmat(i,j) = sqrt(iso(i)*iso(j))*matrx(i,j)/dsqrt(matrx(i,i)*matrx(j,j))
    end do
  end do

matrx = tmpmat

end subroutine



subroutine svcov_anisotropy_calc(svcov)
use mod_tempfact, only: aniso_eigratio3
type(vcov_store),intent(in out) :: svcov
real(dp) :: atmat(3,3),aniso
integer :: i,icor

allocate(svcov%anisotropy(svcov%natoms),stat=icor)
if(icor /= 0) stop "svcov_anisotropy_calc> malloc"

icor = 1
do i=1, svcov%natoms
  atmat = svcov%aniso(icor:icor+2,:) 
  call aniso_eigratio3(atmat,aniso)
  svcov%anisotropy(i)=aniso
  icor = icor + 3
end do

end subroutine

subroutine svcov_aniso(svcov)
type(vcov_store), intent(in out) :: svcov
integer :: i,icor
logical :: tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zvcov,tf_zaniso
logical, parameter :: T=.true., F=.false.

tf_iso     = F
tf_isovcov = F
tf_vcov    = F
tf_zvcov   = F
tf_aniso   = T
tf_zaniso  = F
call vcov_store_init(svcov%natoms,tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zaniso,tf_zvcov, svcov)

if (allocated(svcov%vcov)) then

  icor = 1
  do i=1, svcov%natoms
    svcov%aniso(icor:icor+2,:) = svcov%vcov(icor:icor+2,icor:icor+2)
    icor = icor + 3
  end do
end if
!if (allocated(svcov%zvcov)) then
!  icor = 1
!  do i=1, svcov%natoms
!    svcov%aniso(icor:icor+2,:) = dble(svcov%zvcov(icor:icor+2,icor:icor+2))
!    icor = icor + 3
!  end do
!end if

end subroutine

subroutine svcov_isovcov(svcov)
type(vcov_store), intent(in out) :: svcov
integer :: i,icor,j,jcor
logical :: tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zvcov,tf_zaniso
logical, parameter :: T=.true., F=.false.

tf_iso     = F
tf_aniso   = F
tf_isovcov = T
tf_vcov    = F
tf_zaniso  = F
tf_zvcov   = F
call vcov_store_init(svcov%natoms,tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zaniso,tf_zvcov, svcov)

if (allocated(svcov%vcov)) then

  icor = 1
  do i = 1, svcov%natoms
    jcor = 1
    do j = 1, svcov%natoms
      svcov%isovcov(i,j) = ( svcov%vcov(icor,jcor)     + &
                             svcov%vcov(icor+1,jcor+1) + &
                             svcov%vcov(icor+2,jcor+2))/3.0d0
      jcor = jcor + 3
    end do
    icor = icor + 3
  end do
end if
if (allocated(svcov%zvcov)) then
  icor = 1
  do i = 1, svcov%natoms
    jcor = 1
    do j = 1, svcov%natoms
      svcov%isovcov(i,j) = dble((svcov%zvcov(icor,jcor)     + &
                                 svcov%zvcov(icor+1,jcor+1) + &
                                 svcov%zvcov(icor+2,jcor+2))/3.0d0)
      jcor = jcor + 3
    end do
    icor = icor + 3
  end do
end if

end subroutine


subroutine vcov_store_init(natoms,tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zaniso,tf_zvcov, svcov)
! initiator
integer,          intent(in)     :: natoms
logical,          intent(in)     :: tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zvcov,tf_zaniso
type(vcov_store), intent(in out) :: svcov
integer :: ialloc

svcov%natoms = natoms

if(tf_iso) then 
  if (allocated(svcov%iso)) deallocate(svcov%iso)
  allocate(svcov%iso(natoms), stat=ialloc)
  if (ialloc /= 0) stop "vcov_store_init> memory problem"
  svcov%iso = zero
end if
if(tf_aniso) then
  if (allocated(svcov%aniso)) deallocate(svcov%aniso)
  allocate(svcov%aniso(3*natoms,3), stat=ialloc)
  if (ialloc /= 0) stop "vcov_store_init> memory problem"
  svcov%aniso = zero
end if
if(tf_isovcov) then
  if (allocated(svcov%isovcov)) deallocate(svcov%isovcov)
  allocate(svcov%isovcov(natoms,natoms), stat=ialloc)
  if (ialloc /= 0) stop "vcov_store_init> memory problem"
  svcov%isovcov = zero
end if
if(tf_vcov) then
  if (allocated(svcov%vcov)) deallocate(svcov%vcov)
  allocate(svcov%vcov(3*natoms,3*natoms), stat=ialloc)
  if (ialloc /= 0) stop "vcov_store_init> memory problem"
  svcov%vcov = zero
end if
if(tf_zvcov) then
  if (allocated(svcov%zvcov)) deallocate(svcov%zvcov)
  allocate(svcov%zvcov(3*natoms,3*natoms), stat=ialloc)
  if (ialloc /= 0) stop "vcov_store_init> memory problem"
  svcov%zvcov = cmplx(zero,zero,kind=8)
end if
if(tf_zaniso) then
  if (allocated(svcov%zaniso)) deallocate(svcov%zaniso)
  allocate(svcov%zaniso(3*natoms,3), stat=ialloc)
  if (ialloc /= 0) stop "vcov_store_init> memory problem"
  svcov%aniso = cmplx(zero,zero,kind=8)
end if

end subroutine

subroutine vcov_store_deinit(natoms,tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zvcov,tf_zaniso, svcov)
! deinitiator
integer,intent(in) :: natoms
logical, intent(in) :: tf_iso,tf_aniso,tf_isovcov,tf_vcov,tf_zvcov,tf_zaniso
type(vcov_store), intent(out) :: svcov
integer :: ialloc

if(tf_iso) then 
  if (allocated(svcov%iso)) deallocate(svcov%iso)
end if
if(tf_aniso) then
  if (allocated(svcov%aniso)) deallocate(svcov%aniso)
end if
if(tf_isovcov) then
  if (allocated(svcov%isovcov)) deallocate(svcov%isovcov)
end if
if(tf_vcov) then
  if (allocated(svcov%vcov)) deallocate(svcov%vcov)
end if
if(tf_zvcov) then
  if (allocated(svcov%zvcov)) deallocate(svcov%zvcov)
end if
if(tf_zaniso) then
  if (allocated(svcov%zaniso)) deallocate(svcov%zaniso)
end if

end subroutine

end module
