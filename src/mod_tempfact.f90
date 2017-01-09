module mod_tempfact
! contain the subroutine to make bfactor comparisons 
! this code needs to be cleaned up and commented
use cfml_globaldeps,                 only: dp,sp
use mod_constants
use mod_types, only : inp_par, protein_atom_list_type

implicit none
type :: fluct_aniso
  real(dp) :: uxyz(6)
end type
type :: fluct_vcov
  real(dp) :: vcov(3,3)
end type
type :: fluct_zvcov
  complex(kind=8) :: zvcov(3,3)
end type

type :: fluctuations
  integer :: natoms
  real(dp),         allocatable :: iso(:)
  type(fluct_aniso),allocatable :: aniso(:)
  type(fluct_vcov), allocatable :: atii(:),atij(:,:) 
  type(fluct_zvcov), allocatable :: zatii(:),zatij(:,:) 
end type

contains

!DMR Jan 26, 2009
!  added:
!     type        fluctuations and subtypes
!     subroutine  fluct_real_init
!     subroutine  fluct_cmplx_init
!     subroutine  fluct_x_scalar
!     subroutine  z_nm_flucts
subroutine fluct_real_init(input,flucts)
type(inp_par), intent(in) :: input
type(fluctuations), intent(out) :: flucts
integer :: i,j,k,ialloc
flucts%natoms = input%natoms
if(allocated(flucts%iso)) deallocate(flucts%iso)
if(allocated(flucts%aniso)) deallocate(flucts%aniso)
if(allocated(flucts%atii)) deallocate(flucts%atii)
if(allocated(flucts%atij)) deallocate(flucts%atij)
if(allocated(flucts%zatii)) deallocate(flucts%zatii)
if(allocated(flucts%zatij)) deallocate(flucts%zatij)
allocate(flucts%iso(flucts%natoms),flucts%aniso(flucts%natoms),&
         flucts%atii(flucts%natoms),flucts%zatii(flucts%natoms),&
         flucts%atij(flucts%natoms,flucts%natoms),& 
         flucts%zatij(flucts%natoms,flucts%natoms), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

do i =1, flucts%natoms
  flucts%iso(i)        = zero  
  flucts%aniso(i)%uxyz = zero  
  flucts%atii(i)%vcov  = zero  

end do

do i =1, flucts%natoms
  do j = 1, flucts%natoms
    flucts%atij(i,j)%vcov  = zero  
  end do
end do

end subroutine

subroutine fluct_x_scalar(scalar,flucts)
real(dp), intent(in) :: scalar
type(fluctuations), intent(inout) :: flucts
integer :: i,j,natoms

do i =1, flucts%natoms
  flucts%iso(i)   = flucts%iso(i)   * scalar  
  if(allocated(flucts%aniso)) flucts%aniso(i)%uxyz  = flucts%aniso(i)%uxyz  * scalar  
  if(allocated(flucts%atii))  flucts%atii(i)%vcov   = flucts%atii(i)%vcov   * scalar
  if(allocated(flucts%zatii)) flucts%zatii(i)%zvcov = flucts%zatii(i)%zvcov * scalar
end do

do i =1, flucts%natoms
  do j = 1, flucts%natoms
    if(allocated(flucts%atij))  flucts%atij(i,j)%vcov   = flucts%atij(i,j)%vcov  * scalar  
    if(allocated(flucts%zatij)) flucts%zatij(i,j)%zvcov = flucts%zatij(i,j)%zvcov * scalar  
  end do
end do

end subroutine

subroutine fluct_cmplx_init(flucts)
! it will be rarer to need the complex components so we'll separate this
use mod_types, only : inp_par
type(fluctuations), intent(inout) :: flucts
integer :: i,j,k,ialloc

if(allocated(flucts%zatii)) deallocate(flucts%zatii)
if(allocated(flucts%zatij)) deallocate(flucts%zatij)
allocate(flucts%zatii(flucts%natoms),&
         flucts%zatij(flucts%natoms,flucts%natoms), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

do i =1, flucts%natoms
  flucts%zatii(i)%zvcov  = cmplx(zero,zero,kind=8) 
end do

do i =1, flucts%natoms
  do j = 1, flucts%natoms
    flucts%zatij(i,j)%zvcov  = cmplx(zero,zero,kind=8) 
  end do
end do

end subroutine

subroutine z_nm_flucts(input,vals,vecs,flucts)
! take complex vectors and store the corresponding fluctuations
use mod_types, only : inp_par
use mod_constants, only: Rkcal
type(inp_par), intent(in) :: input
real(dp),        allocatable, intent(in) :: vals(:)
complex(kind=8), allocatable, intent(in) :: vecs(:,:)
type(fluctuations), intent(in out) :: flucts
real(dp) :: val,kt
complex(kind=8) :: mat(3,3)
integer :: i,j,k,ii,jj,m,n,l

kT = input%temperature*Rkcal

mat = cmplx(zero,zero,kind = 8)
! self terms
do i = 1, size(vals)
  if (abs(vals(i)) .gt. 1.0d-08) then
    val = kT/vals(i)
    j = 1
    do k =1, flucts%natoms
      do ii=0,2
        do jj=0,2
         mat(ii+1,jj+1) = vecs(j+jj,i)*CONJG(vecs(j+ii,i))
        end do
      end do

      if (allocated(flucts%zatii)) flucts%zatii(k)%zvcov  = flucts%zatii(k)%zvcov + val*mat
      if (allocated(flucts%atii))  flucts%atii(k)%vcov    = flucts%atii(k)%vcov   + val*dble(mat)
      flucts%iso(k) = flucts%iso(k) + val*(dble(mat(1,1))+dble(mat(2,2))+dble(mat(3,3)))/3.0d0
      j = j+3
    end do
  end if
end do

if ((allocated(flucts%atii)) .and. (allocated(flucts%aniso))) then

  do k =1, flucts%natoms
    flucts%aniso(k)%uxyz(1) = flucts%atii(k)%vcov(1,1) 
    flucts%aniso(k)%uxyz(2) = flucts%atii(k)%vcov(2,2) 
    flucts%aniso(k)%uxyz(3) = flucts%atii(k)%vcov(3,3) 
    flucts%aniso(k)%uxyz(4) = flucts%atii(k)%vcov(1,2) 
    flucts%aniso(k)%uxyz(5) = flucts%atii(k)%vcov(1,3) 
    flucts%aniso(k)%uxyz(6) = flucts%atii(k)%vcov(2,3) 
  end do

else if (allocated(flucts%aniso)) then
  print *, 'z_nm_flucts> Warning: vcov not allocated, aniso set to zero'
end if

!cross terms

if (allocated(flucts%atij)) then
  mat = cmplx(zero,zero,kind=8)
  do i = 1, size(vals)
    if (abs(vals(i)) .gt. 1.0d-08) then
      val = kT/vals(i)
      m = 1
      do k = 1, flucts%natoms
        n = 1
        do l =1, flucts%natoms
          do ii=0,2
            do jj=0,2
              mat(ii+1,jj+1) = vecs(m+jj,i)*CONJG(vecs(n+ii,i))
            end do
          end do
          if (allocated(flucts%zatij)) flucts%zatij(k,l)%zvcov  = flucts%zatij(k,l)%zvcov + val*mat
          if (allocated(flucts%atij))  flucts%atij(k,l)%vcov    = flucts%atij(k,l)%vcov   + val*dble(mat)
          n = n+ 3
        end do
        m = m + 3
      end do
    end if
  end do
end if

end subroutine

subroutine vcov_Bscale(input,mol_biso,vcov,correl,fcnst,nviab)
use mod_types,                only: inp_par,isoaniso
use mod_math,                 only: analyt_lsq,analyt_lsq_offset,linear_correl
! DMR: 02-20-2008
! take in only the asym_unit, since the vcov is periodic with asym units
! use bfactors for atoms with 1.0 occupancy as did dkon
type (inp_par)   ,        intent(in out)    :: input
type(protein_atom_list_type) ,intent(in)  :: mol_biso
real(dp), allocatable,intent(in out)  :: vcov(:,:)
real(dp)             ,intent(out) :: correl,fcnst
integer              ,intent(out) :: nviab
real(dp), allocatable             :: bfexp(:),bfthr(:)
real(dp), allocatable             :: tmpbfexp(:),tmpbfthr(:)
integer,  allocatable             :: resid(:)
integer                           :: ndim,i,j,k,ialloc,ier,nzeros
real(dp)                          :: b, mean_bfact,orig_fcnst
type(isoaniso)                    :: compars

orig_fcnst = input%fcnst
call compisoaniso(input,mol_biso,vcov,compars)
call analyt_lsq(compars%exp_flucu,compars%thr_flucu,fcnst,b)
!call analyt_lsq_offset(compars%exp_flucu,compars%thr_flucu,fcnst,b)

if (input%fcnst .eq. 1.0d0) then
  !print *, 'scaling the vcov with the fcnst from lsq w/o offset:', fcnst,b
  vcov = vcov*fcnst+b ! get vcov on same scale as bfactor
  compars%thr_flucu=compars%thr_flucu*fcnst+b
  fcnst  = one/fcnst
  input%fcnst = fcnst
else
  fcnst  = one/fcnst
end if

if (input%print_level .ge. 1) then

  write(6,'(A8,1x,F7.3,1x,A8,1x,F7.3,1x,A20,1x,E14.7,1x,A14,1x,I5)') &
           'isocor =',compars%isocorr,':: aniscor =',compars%anisocorr,':: fcnst_fit =', &
            fcnst

      open(unit=22,file=trim(input%fileroot)//"-"//trim(input%genm)//"-"//trim(input%tpbc)//"-"// &
                        trim(input%fctyp) //"-"//trim(input%bztyp)//"-isoaniso.txt", &
      status='unknown',action='write', iostat=ier)

do i=1,compars%nviab
  if (input%fcpower > 0 .and. orig_fcnst .ne. 1.0d0) then
    write (22,'(I5,4F12.5)') compars%resid(i),compars%exp_flucu(i),compars%thr_flucu(i)/fcnst +b ,&
                           compars%exp_aniso(i),compars%thr_aniso(i)  ! had to add cuz looks like pfanm needs high fcnst
                                                                      ! to run due to small nums
  else
    write (22,'(I5,4F12.5)') compars%resid(i),compars%exp_flucu(i),compars%thr_flucu(i),&
                           compars%exp_aniso(i),compars%thr_aniso(i)
  end if

end do

close(22)

end if

correl = compars%isocorr
nviab  = compars%nviab

end subroutine

subroutine bfactor_scale(input,mol_biso,thr_iso,&
                            bfexp,bfthr,resid,correl,fcnst,b,nviab,offsetin)
use mod_math,                 only: analyt_lsq,analyt_lsq_offset,linear_correl
use mod_types,                only: inp_par
! DMR: 11-17-2009
! take in asym_unit
! use bfactors for atoms with 1.0 occupancy 
type(inp_par),        intent(inout)  :: input
type(protein_atom_list_type) ,intent(in)  :: mol_biso
real(dp), allocatable,intent(in)  :: thr_iso(:)
real(dp)             ,intent(out) :: correl,fcnst,b
integer              ,intent(out) :: nviab
logical, optional,    intent(in)  :: offsetin ! if a fit with offset is desired
real(dp), allocatable,intent(out) :: bfexp(:),bfthr(:)
integer,  allocatable,intent(out) :: resid(:)
logical :: offset

if (present(offsetin)) then
  offset = offsetin
else
  offset = .false.
end if

! first let's get all occ one's into arrays
call getbf_occ1(input,mol_biso,thr_iso,bfexp,bfthr,resid,nviab)
input%nviab = nviab

if (offset) then
  call analyt_lsq_offset(bfexp,bfthr,fcnst,b)
else 
  call analyt_lsq(bfexp,bfthr,fcnst,b)
end if

bfthr = bfthr*fcnst+b
call linear_correl(bfexp,bfthr,correl)
if (correl > 1.0) bfthr = bfthr*correl  ! dumb fix

end subroutine

subroutine aniso_analysis(input,mol_biso,thr_aniso, &
            resid,expanisotropy,thranisotropy,dotab,ccmod,suij, &
            avg_expaniso,avg_thraniso,avg_dotab,avg_ccmod,avg_suij,&
            std_expaniso,std_thraniso,std_dotab,std_ccmod,std_suij,navg)
use mod_math,                 only: analyt_lsq,analyt_lsq_offset,linear_correl
use mod_types,                only: inp_par
! DMR: 11-17-2009
! take in asym_unit
! use bfactors for atoms with 1.0 occupancy 
type(inp_par),        intent(inout)  :: input
type(protein_atom_list_type) ,intent(in)     :: mol_biso
real(dp), allocatable,intent(in)     :: thr_aniso(:,:)
integer,  allocatable ,intent(out)   :: resid(:)
real(dp), allocatable ,intent(out)   :: expanisotropy(:),thranisotropy(:),&
                                     dotab(:),ccmod(:),suij(:)
real(dp),              intent(out)   :: avg_expaniso,avg_thraniso,avg_dotab,&
                                        avg_ccmod,avg_suij,std_expaniso, &
                                        std_thraniso,std_dotab,std_ccmod,std_suij
integer, intent(out) :: navg
real(dp), allocatable                :: aniso_thr(:,:),aniso_exp(:,:)
real(dp), allocatable                :: sub_expan(:),sub_thran(:),&
                                        sub_dotab(:),sub_ccmod(:),sub_suij(:)
integer                              :: i,nviab,icor
integer, allocatable :: sub_resid(:)

! first let's get all occ one's into arrays
call getaniso_occ1(input,mol_biso,thr_aniso,aniso_exp,aniso_thr,resid,nviab)
input%nviab = nviab

call aniso_vectors(aniso_exp,aniso_thr,&
                         expanisotropy,thranisotropy,dotab,ccmod,suij)

call mask_int(zero,0.5d0,expanisotropy,resid,sub_resid)
call mask_it(zero,0.5d0,expanisotropy,expanisotropy,sub_expan)
call mask_it(zero,0.5d0,expanisotropy,thranisotropy,sub_thran)
call mask_it(zero,0.5d0,expanisotropy,dotab,sub_dotab)
call mask_it(zero,0.5d0,expanisotropy,ccmod,sub_ccmod)
call mask_it(zero,0.5d0,expanisotropy,suij,sub_suij)

do i = 1, size(expanisotropy)
  !print '(A6,I7,5F10.4)',"aniso>",sub_resid(i),sub_expan(i),sub_thran(i),sub_dotab(i),sub_ccmod(i),sub_suij(i)
  print '(A6,I7,5F10.4)',"aniso>",resid(i),expanisotropy(i),thranisotropy(i),dotab(i),ccmod(i),suij(i)
end do

call aniso_analysis_avgstd(expanisotropy,thranisotropy,dotab,ccmod,suij, &
                             avg_expaniso,avg_thraniso,avg_dotab,avg_ccmod,avg_suij, &
                             std_expaniso,std_thraniso,std_dotab,std_ccmod,std_suij,navg)

end subroutine

subroutine getbf_occ1(input,mol_biso,thr_iso,bfexp,bfthr,resid,nviab)
use mod_types,                only: inp_par
! first let's get all occ one's into arrays
type(inp_par),        intent(in)  :: input
type(protein_atom_list_type) ,intent(in)  :: mol_biso
real(dp), allocatable,intent(in)  :: thr_iso(:)
real(dp), allocatable,intent(out) :: bfexp(:),bfthr(:)
integer,  allocatable,intent(out) :: resid(:)
integer,              intent(out) :: nviab
real(dp), allocatable             :: tmpbfexp(:),tmpbfthr(:)
integer,  allocatable             :: tmpresid(:)
integer                           :: i,j,ier

allocate(tmpbfexp(mol_biso%natoms),tmpbfthr(mol_biso%natoms), &
         tmpresid(mol_biso%natoms),stat =ier)

j = 0
do i = 1, mol_biso%natoms
  if (mol_biso%atom(i)%occ .eq. 1.0) then
    if (trim(input%atom_anl) .eq. "all") then
      j = j+1
      tmpbfexp(j) = mol_biso%atom(i)%biso * b2u ! b2u is in constants
      tmpbfthr(j) = thr_iso(i)
      tmpresid(j) = mol_biso%atom(i)%ires
    else
      if (trim(input%atom_anl) .eq. trim(mol_biso%atom(i)%lab)) then
        j = j+1
        tmpbfexp(j) = mol_biso%atom(i)%biso * b2u
        tmpbfthr(j) = thr_iso(i)
        tmpresid(j) = mol_biso%atom(i)%ires

      end if
    end if
  end if
end do


nviab = j

allocate(bfexp(nviab),bfthr(nviab),resid(nviab),stat =ier)

bfexp = tmpbfexp(1:nviab) ; deallocate(tmpbfexp)
bfthr = tmpbfthr(1:nviab) ; deallocate(tmpbfthr)
resid = tmpresid(1:nviab) ; deallocate(tmpresid)

end subroutine

subroutine aniso_analysis_avgstd(expaniso,thraniso,dotab,ccmod,suij, &
                             avg_expaniso,avg_thraniso,avg_dotab,avg_ccmod,avg_suij, &
                             std_expaniso,std_thraniso,std_dotab,std_ccmod,std_suij,navg)
use mod_types,                only: inp_par
use mod_math,                 only: vec_avg_std
! DMR: 11-17-2009
! return vectors of anisotropys, dots, ccmod, and suij for viable atoms (occ=1)
real(dp), allocatable,intent(in)  :: expaniso(:),thraniso(:),dotab(:),ccmod(:),suij(:)
real(dp),             intent(out) :: avg_expaniso,avg_thraniso,avg_dotab,&
                                     avg_ccmod,avg_suij
real(dp),             intent(out) :: std_expaniso,std_thraniso,std_dotab,&
                                     std_ccmod,std_suij
integer, intent(out) :: navg
real(dp), allocatable             :: sub_dotab(:),sub_ccmod(:),sub_suij(:)
integer :: i

print *, '<dots>,<ccmod>,<suij> analyzed for atoms with experiment anisotropy <= 0.5'
print *, 'all are retained for printing or whatever'
call mask_it(zero,0.5d0,expaniso,dotab,sub_dotab)
call mask_it(zero,0.5d0,expaniso,ccmod,sub_ccmod)
call mask_it(zero,0.5d0,expaniso,suij,sub_suij)
navg = size(sub_dotab)
call vec_avg_std(expaniso,avg_expaniso,std_expaniso)
call vec_avg_std(thraniso,avg_thraniso,std_thraniso)
call vec_avg_std(sub_dotab,avg_dotab,std_dotab)
call vec_avg_std(sub_ccmod,avg_ccmod,std_ccmod)
call vec_avg_std(sub_suij,avg_suij,std_suij)

end subroutine


subroutine mask_it(val_low,val_high,testvec,modvec,subvec)
! take subselection of modvec using val and testvec
real(dp),              intent(in)  :: val_low,val_high
real(dp), allocatable, intent(in)  :: testvec(:),modvec(:)
real(dp), allocatable, intent(out) :: subvec(:)
real(dp), allocatable              :: tmpvec(:)
integer :: ncnt, i,ier

if (size(testvec) .ne. size(modvec)) stop "mask_it> sizes not same"

allocate(tmpvec(size(testvec)),stat=ier)
if(ier/=0) stop "mask_it> malloc"

ncnt=0
do i = 1, size(testvec)
  if (testvec(i) .gt. val_low .and. testvec(i) .le. val_high) then
    ncnt = ncnt+1
    tmpvec(ncnt) = modvec(i)
  end if  
end do

allocate(subvec(ncnt), stat=ier)
if(ier/=0) stop "mask_it> malloc"

subvec = tmpvec(1:ncnt) ; deallocate(tmpvec)

end subroutine

subroutine mask_int(val_low,val_high,testvec,modvec,subvec)
! take subselection of modvec using val and testvec
real(dp),              intent(in)  :: val_low,val_high
real(dp), allocatable, intent(in)  :: testvec(:)
integer, allocatable, intent(in)  :: modvec(:)
integer, allocatable, intent(out) :: subvec(:)
integer, allocatable              :: tmpvec(:)
integer :: ncnt, i,ier

if (size(testvec) .ne. size(modvec)) stop "mask_it> sizes not same"

allocate(tmpvec(size(testvec)),stat=ier)
if(ier/=0) stop "mask_it> malloc"

ncnt=0
do i = 1, size(testvec)
  if (testvec(i) .gt. val_low .and. testvec(i) .le. val_high) then
    ncnt = ncnt+1
    tmpvec(ncnt) = modvec(i)
  end if
end do

allocate(subvec(ncnt), stat=ier)
if(ier/=0) stop "mask_it> malloc"

subvec = tmpvec(1:ncnt) ; deallocate(tmpvec)

end subroutine


subroutine getaniso_occ1(input,mol_biso,thr_aniso,aniso_exp,aniso_thr,resid,nviab)
use mod_types,                only: inp_par
! first let's get all occ one's into arrays
type(inp_par),        intent(in)  :: input
type(protein_atom_list_type) ,intent(in)  :: mol_biso
real(dp), allocatable,intent(in)  :: thr_aniso(:,:)
real(dp), allocatable,intent(out) :: aniso_exp(:,:),aniso_thr(:,:)
integer,  allocatable,intent(out) :: resid(:)
integer,              intent(out) :: nviab
real(dp), allocatable             :: tmpexp(:,:),tmpthr(:,:)
real(dp) :: threebyexp(3,3)
integer,  allocatable             :: tmpresid(:)
integer                           :: i,j,icor,jcor,ier

allocate(tmpexp(3*mol_biso%natoms,3),tmpthr(3*mol_biso%natoms,3), &
         tmpresid(mol_biso%natoms),stat =ier)

j = 0
icor=1
jcor=1
do i = 1, mol_biso%natoms
  if (mol_biso%atom(i)%occ .eq. 1.0) then
    if (trim(input%atom_anl) .eq. "all") then
      j = j+1
      threebyexp(1,1) = mol_biso%atom(i)%u(1)
      threebyexp(2,2) = mol_biso%atom(i)%u(2)
      threebyexp(3,3) = mol_biso%atom(i)%u(3)
      threebyexp(1,2) = mol_biso%atom(i)%u(4)
      threebyexp(1,3) = mol_biso%atom(i)%u(5)
      threebyexp(2,3) = mol_biso%atom(i)%u(6)
      threebyexp(2,1) = mol_biso%atom(i)%u(4)
      threebyexp(3,1) = mol_biso%atom(i)%u(5)
      threebyexp(3,2) = mol_biso%atom(i)%u(6)
      tmpexp(jcor:jcor+2,:)   = threebyexp/10000.0d0  ! pdb aniso tempfactors are in (\AA^2)*10000. 
      tmpthr(jcor:jcor+2,:)   = thr_aniso(icor:icor+2,:)
      tmpresid(j) = mol_biso%atom(i)%ires
      jcor = jcor+3
    else
      if (trim(input%atom_anl) .eq. trim(mol_biso%atom(i)%lab)) then
        j = j+1
        threebyexp(1,1) = mol_biso%atom(i)%u(1)
        threebyexp(2,2) = mol_biso%atom(i)%u(2)
        threebyexp(3,3) = mol_biso%atom(i)%u(3)
        threebyexp(1,2) = mol_biso%atom(i)%u(4)
        threebyexp(1,3) = mol_biso%atom(i)%u(5)
        threebyexp(2,3) = mol_biso%atom(i)%u(6)
        threebyexp(2,1) = mol_biso%atom(i)%u(4)
        threebyexp(3,1) = mol_biso%atom(i)%u(5)
        threebyexp(3,2) = mol_biso%atom(i)%u(6)
        tmpexp(jcor:jcor+2,:)   = threebyexp/10000.0d0  ! pdb aniso tempfactors are in (\AA^2)*10000. 
        tmpthr(jcor:jcor+2,:)   = thr_aniso(icor:icor+2,:)
        tmpresid(j) = mol_biso%atom(i)%ires
        jcor = jcor+3
      end if
    end if
  end if
  icor=icor+3
end do


nviab = j

allocate(aniso_exp(3*nviab,3),aniso_thr(3*nviab,3),resid(nviab),stat =ier)

aniso_exp = tmpexp(1:3*nviab,:) ; deallocate(tmpexp)
aniso_thr = tmpthr(1:3*nviab,:) ; deallocate(tmpthr)
resid     = tmpresid(1:3*nviab) ; deallocate(tmpresid)

end subroutine

subroutine get_expaniso(input,mol_biso,aniso_exp)
use mod_types,                only: inp_par
! first let's get all occ one's into arrays
type(inp_par),        intent(in)  :: input
type(protein_atom_list_type) ,intent(in)  :: mol_biso
real(dp), allocatable,intent(out) :: aniso_exp(:,:)
real(dp), allocatable             :: tmpexp(:,:)
real(dp) :: threebyexp(3,3)
integer                           :: i,j,icor,jcor,ier

allocate(aniso_exp(3*mol_biso%natoms,3),stat =ier)

jcor=1
do i = 1, mol_biso%natoms
  if (trim(input%atom_anl) .eq. "all") then
    threebyexp(1,1) = mol_biso%atom(i)%u(1)
    threebyexp(2,2) = mol_biso%atom(i)%u(2)
    threebyexp(3,3) = mol_biso%atom(i)%u(3)
    threebyexp(1,2) = mol_biso%atom(i)%u(4)
    threebyexp(1,3) = mol_biso%atom(i)%u(5)
    threebyexp(2,3) = mol_biso%atom(i)%u(6)
    threebyexp(2,1) = mol_biso%atom(i)%u(4)
    threebyexp(3,1) = mol_biso%atom(i)%u(5)
    threebyexp(3,2) = mol_biso%atom(i)%u(6)
    aniso_exp(jcor:jcor+2,:)   = threebyexp/10000.0d0  ! pdb aniso tempfactors are in (\AA^2)*10000. 
    jcor = jcor+3
  else
    if (trim(input%atom_anl) .eq. trim(mol_biso%atom(i)%lab)) then
      threebyexp(1,1) = mol_biso%atom(i)%u(1)
      threebyexp(2,2) = mol_biso%atom(i)%u(2)
      threebyexp(3,3) = mol_biso%atom(i)%u(3)
      threebyexp(1,2) = mol_biso%atom(i)%u(4)
      threebyexp(1,3) = mol_biso%atom(i)%u(5)
      threebyexp(2,3) = mol_biso%atom(i)%u(6)
      threebyexp(2,1) = mol_biso%atom(i)%u(4)
      threebyexp(3,1) = mol_biso%atom(i)%u(5)
      threebyexp(3,2) = mol_biso%atom(i)%u(6)
      aniso_exp(jcor:jcor+2,:)   = threebyexp/10000.0d0  ! pdb aniso tempfactors are in (\AA^2)*10000. 
      jcor = jcor+3
    end if
  end if
end do


end subroutine

subroutine get_expiso(input,mol_biso,iso_exp)
use mod_types,                only: inp_par
! first let's get all occ one's into arrays
type(inp_par),        intent(in)  :: input
type(protein_atom_list_type) ,intent(in)  :: mol_biso
real(dp), allocatable,intent(out) :: iso_exp(:)
integer                           :: i,j,icor,jcor,ier

allocate(iso_exp(mol_biso%natoms),stat =ier)

jcor=1
do i = 1, mol_biso%natoms
  if (trim(input%atom_anl) .eq. "all") then
    iso_exp(i) = mol_biso%atom(i)%biso*b2u
  else
    if (trim(input%atom_anl) .eq. trim(mol_biso%atom(i)%lab)) then
      iso_exp(i) = mol_biso%atom(i)%biso*b2u
      jcor = jcor+1
    end if
  end if
end do

end subroutine

subroutine aniso_vectors(expaniso,thraniso,&
                         expanisotropy,thranisotropy,dotab,ccmod,suij)
real(dp), allocatable,intent(in)  :: expaniso(:,:),thraniso(:,:)
real(dp), allocatable,intent(out) :: expanisotropy(:),thranisotropy(:),&
                                     dotab(:),ccmod(:),suij(:)
integer  :: i,icor, ier,nviab
real(dp),allocatable :: mata(:,:),matb(:,:)
real(dp) :: sh1,sh2,sh3,sh4,sh5

nviab = size(expaniso,1)/3

allocate(mata(3,3),matb(3,3))

allocate(expanisotropy(nviab),thranisotropy(nviab),&
                                     dotab(nviab),ccmod(nviab),suij(nviab), stat=ier)
if (ier/=0) stop "aniso_vectors>malloc"

icor = 1

do i = 1,nviab
  mata = expaniso(icor:icor+2,:)  
  matb = thraniso(icor:icor+2,:)  
  call thrxthr_aniso(mata,matb, &
                       expanisotropy(i),thranisotropy(i),&
                       dotab(i),ccmod(i),suij(i))
  icor = icor+3
end do

end subroutine

subroutine thrxthr_aniso(thrxexp,thrxthr, &
                              expanisotropy,thranisotropy,dotab,ccmod,suij)
real(dp), allocatable, intent(in)  :: thrxexp(:,:),thrxthr(:,:)
real(dp), intent(out) :: expanisotropy,thranisotropy,dotab,ccmod,suij
            
call aniso_eigratio(thrxexp,expanisotropy)
call aniso_eigratio(thrxthr,thranisotropy)
call aniso_vecdot (thrxexp,thrxthr,dotab)
call aniso_ccmod (thrxexp,thrxthr,ccmod,suij)

end subroutine

subroutine compisoaniso(input,mol_biso,vcov,comps)
use mod_math,                 only: linear_correl
use mod_types,                only: inp_par,isoaniso
use mod_inout,                only: threed,oned
! DMR
! take in var-covar matrix and pump out some comparisons
! DMR 12-04-2008: carry out calcs for just CA atoms... for when all atom is used
type(inp_par), intent(in)             :: input
type(protein_atom_list_type) ,intent(in)      :: mol_biso
real(dp), allocatable,intent(in)      :: vcov(:,:)
type(isoaniso),       intent(out)     :: comps
real(dp), dimension(:),allocatable    :: tmpbfexp,tmpbfthr,tmpaniso,tmpanexp
real(dp), dimension(:,:),allocatable  :: threebythr,threebyexp
integer,  allocatable                 :: resid(:)
integer :: ndim, nviab, i,j,k,ialloc
real(dp)                              :: correl
real(dp)                              :: dotab,sumdotab,avgdotab
real(dp)                              :: ccmod,sumccmod,avgccmod
real(dp)                              :: suij,sumsuij,avgsuij
real(dp)                              :: anisovalexp,anisovalthr
real(dp)                              :: sumanisoexp,avganisoexp,sumanisothr,avganisothr
integer :: naniso,ndotcc

! number of anisos used
naniso   = 0
ndotcc     = 0
sumdotab = zero
sumccmod = zero
sumanisoexp = zero
sumanisothr = zero

if (trim(input%atom_anl) .eq. "all") then
ndim=mol_biso%natoms
! DMR 12-04-08 
! find out how many CA atoms there are
else
ndim = 0
  do i=1,mol_biso%natoms
    if(trim(mol_biso%atom(i)%lab) .eq. trim(input%atom_anl)) then
      ndim = ndim + 1
    end if
  end do
end if

if (allocated(tmpbfexp)) deallocate(tmpbfexp)
if (allocated(tmpbfthr)) deallocate(tmpbfthr)
if (allocated(tmpaniso)) deallocate(tmpaniso)
if (allocated(tmpanexp)) deallocate(tmpanexp)
if (allocated(resid)) deallocate(resid)
allocate(tmpbfexp(ndim),tmpbfthr(ndim),resid(ndim),tmpaniso(ndim),tmpanexp(ndim),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"
allocate(threebythr(3,3),threebyexp(3,3),stat=ialloc)

if (input%print_level .ge. 1) then
  !if ((input%aniso) .and. (trim(input%genm) .eq. "enm")) then
   if((input%aniso) .and. (trim(input%genm) .ne. "gnm" ) .and. (trim(input%genm) .ne. "isocorr")) then
    print '(A6,4x,A4,4x,A5,4x,A5,5x,A5,6x,A5,4x,A6,4x,A5)',"aniso ", "rcut", "resid",  "exper", &
            "theor", "dotpr",  "ccmod",  "suij"
  end if
end if

  j=0
  k=1
  do i=1,mol_biso%natoms
    ! use only occ=1.0 bfactors
!  CARRY OUT ANALYSIS FOR CA 
!  if(trim(mol_biso%atom(i)%lab) .eq. "CA") then
    if (mol_biso%atom(i)%occ .eq. one) then
      j=j+1
      resid(j)=mol_biso%atom(i)%ires
      if (threed(input)) then
        tmpbfthr(j) = (vcov(k,k)+vcov(k+1,k+1)+vcov(k+2,k+2))/3.0d0
        if(input%aniso) then
          threebythr=vcov(k:k+2,k:k+2)
          call aniso_eigratio(threebythr,anisovalthr)
          tmpaniso(j) = anisovalthr
          if (sum(mol_biso%atom(i)%u(1:3)) .gt. zero) then
            threebyexp(1,1) = mol_biso%atom(i)%u(1)
            threebyexp(2,2) = mol_biso%atom(i)%u(2)
            threebyexp(3,3) = mol_biso%atom(i)%u(3)
            threebyexp(1,2) = mol_biso%atom(i)%u(4)
            threebyexp(1,3) = mol_biso%atom(i)%u(5)
            threebyexp(2,3) = mol_biso%atom(i)%u(6)
            threebyexp(2,1) = mol_biso%atom(i)%u(4)
            threebyexp(3,1) = mol_biso%atom(i)%u(5)
            threebyexp(3,2) = mol_biso%atom(i)%u(6)
            call aniso_eigratio(threebyexp,anisovalexp)
            tmpanexp(j) = anisovalexp
          if (anisovalexp .gt. 0.0) then
            !compare two ellipsoids, normalized by trace: dot of the major axis
!1eb6        if ((i .eq. 33) .or. (i .eq. 38) .or. (i .eq. 47) .or. (i .eq. 50))  then
!3ltz           if (i .eq. 106)  then
!1lni
!           if (i .eq. 298)  then
              naniso = naniso + 1
              sumanisoexp = sumanisoexp+anisovalexp
              sumanisothr = sumanisothr+anisovalthr
              avganisoexp = sumanisoexp/real(naniso)
              avganisothr = sumanisothr/real(naniso)
              comps%avganisoexp = avganisoexp
              comps%avganisothr = avganisothr
          ! compute running average of dotproduct for major axis 
           if (anisovalexp .lt. 0.5) then
              ndotcc = ndotcc + 1
              call aniso_vecdot (threebyexp,threebythr,dotab)
              sumdotab       = sumdotab + dotab
              avgdotab       = sumdotab/real(ndotcc)
              comps%avgdotab = avgdotab
          ! compute ~overlap between two ellipsoids
              call aniso_ccmod (threebyexp,threebythr,ccmod,suij)
              sumccmod       = sumccmod + ccmod
              avgccmod       = sumccmod/real(ndotcc)
              comps%avgccmod = avgccmod
              sumsuij       = sumsuij + suij
              avgsuij       = sumsuij/real(ndotcc)
              comps%avgsuij = avgsuij
              if(input%print_level .ge. 1) then 
                print '(A6,1x,F8.2,1x,I6,1x,5F10.4)', &
                      "aniso>", input%rcut_start, i, &
                       anisovalexp,anisovalthr,dotab,ccmod,suij
              end if
!residue if above end if
             end if
            end if

          else
            tmpanexp(j)    = zero
          end if
        else
          tmpaniso(j)      = zero
          comps%avgdotab   = zero
          comps%avgccmod   = zero
          comps%anisocorr  = zero
          comps%avganisoexp = zero
          comps%avganisothr = zero
        end if
      else if((trim(input%genm) .eq. "gnm" ) .or. (trim(input%genm) .eq. "isocorr")) then
        tmpbfthr(j) = vcov(i,i)
        tmpaniso(j) = zero
      end if
      !tmpbfexp(j)= 3.0d0*mol_biso%atom(i)%biso/(8.0d0*pi**2)
      tmpbfexp(j)= mol_biso%atom(i)%biso/(8.0d0*pi**2)
! DMR carry out for CA
end if
    k=k+3
!    end if
  end do
nviab=j

!fill up the comps data structure
if (allocated(comps%exp_flucu)) deallocate(comps%exp_flucu)
if (allocated(comps%thr_flucu)) deallocate(comps%thr_flucu)
if (allocated(comps%exp_aniso)) deallocate(comps%exp_aniso)
if (allocated(comps%thr_aniso)) deallocate(comps%thr_aniso)
allocate( comps%exp_flucu(nviab),comps%thr_flucu(nviab), &
          comps%thr_aniso(nviab),comps%exp_aniso(nviab),comps%resid(nviab), stat=ialloc)

comps%nviab    = nviab
comps%exp_flucu=tmpbfexp(1:nviab)
comps%thr_flucu=tmpbfthr(1:nviab)
comps%resid=resid(1:nviab)
call linear_correl(comps%exp_flucu,comps%thr_flucu,correl)
comps%isocorr = correl

if((input%aniso) .and. (trim(input%genm) .ne. "gnm" ) .and. (trim(input%genm) .ne. "isocorr")) then
  comps%naniso   = naniso  
  comps%ndotcc   = ndotcc
  comps%thr_aniso=tmpaniso(1:nviab)
  comps%exp_aniso=tmpanexp(1:nviab)
  call linear_correl(comps%exp_aniso,comps%thr_aniso,correl)
  comps%anisocorr = correl
  print '(A13,2x,A4,2x,A6,3x,A3,4x,A3,3x,A8,1x,A5,2x,A5,2x,A4)', &
  "scan> <aniso>",  "rcut", "#aniso", "exp", "thr", "#dot-ccs", "dotpr", "ccmod", "suij"
  print '(A13,1x,F5.2,I7,1x,2F7.3,1x,I7,1x,3F7.3)',"SCAN> <ANISO>", input%rcut_start,comps%naniso,comps%avganisoexp, &
             comps%avganisothr,comps%ndotcc,comps%avgdotab,comps%avgccmod,comps%avgsuij
else
  comps%thr_aniso=zero
  comps%exp_aniso=zero
end if

deallocate(tmpbfexp,tmpbfthr,tmpaniso,resid)

end subroutine

subroutine aniso_eigratio (matrx,anisoval)
use mod_linalg,               only: lapack_eig
real(dp), intent(in),allocatable  :: matrx(:,:)
real(dp), intent(out) :: anisoval
real(dp), allocatable :: vals(:), vecs(:,:)

call lapack_eig(matrx,vals,vecs)

anisoval=vals(1)/vals(3)

end subroutine

subroutine aniso_eigratio3 (matrx,anisoval)
use mod_linalg,               only: lapack_eig3
real(dp), intent(in)  :: matrx(3,3)
real(dp), intent(out) :: anisoval
real(dp) :: vals(3), vecs(3,3)

call lapack_eig3(matrx,vals,vecs)

anisoval=vals(1)/vals(3)

end subroutine

subroutine aniso_vecdot (matrxa,matrxb,dotab)
use mod_linalg,               only: lapack_eig
real(dp), intent(in),allocatable  :: matrxa(:,:),matrxb(:,:)
real(dp), intent(out) :: dotab
real(dp), allocatable :: norm_matrxa(:,:),norm_matrxb(:,:)
real(dp), allocatable :: valsa(:), vecsa(:,:)
real(dp), allocatable :: valsb(:), vecsb(:,:)
integer :: ialloc,maxa,maxb

allocate( norm_matrxa(3,3), norm_matrxb(3,3), stat=ialloc)

norm_matrxa=matrxa/(matrxa(1,1)+matrxa(2,2)+matrxa(3,3))
norm_matrxb=matrxb/(matrxb(1,1)+matrxb(2,2)+matrxb(3,3))

call lapack_eig(norm_matrxa,valsa,vecsa)
call lapack_eig(norm_matrxb,valsb,vecsb)

maxa = maxloc(valsa, dim=1)
maxb = maxloc(valsb, dim=1)

!dotab=valsa(maxa)*valsb(maxb)*abs(dot_product(vecsa(:,maxa),vecsb(:,maxb)))
dotab=abs(dot_product(vecsa(:,maxa),vecsb(:,maxb)))

end subroutine

subroutine aniso_ccmod (matU,matV,ccmod,suij)
! see Acta Crys D55,1997-2004 (1999) and Dkon's structure paper for the
! ellipsoidal comparison
use mod_linalg,             only: lapack_eig,lapack_pseudo
use mod_math,               only: determ3,valvec_inv
real(dp), intent(in),allocatable  :: matU(:,:),matV(:,:)
real(dp), intent(out) :: ccmod,suij
real(dp), allocatable :: norm_matU(:,:),norm_matV(:,:)
real(dp), allocatable :: matUeq(:,:),matVeq(:,:)
real(dp), allocatable :: valsU(:), vecsU(:,:),tmpvalsU(:)
real(dp), allocatable :: valsV(:), vecsV(:,:),tmpvalsV(:)
real(dp) :: cc,ccUiso,ccViso, ccstr,Ueq,Veq
integer :: ialloc,nzeros,i

allocate( norm_matU(3,3), norm_matV(3,3),matUeq(3,3),matVeq(3,3), stat=ialloc)

matUeq = zero
matVeq = zero

Ueq=(matU(1,1)+matU(2,2)+matU(3,3))/3.0d0
Veq=(matV(1,1)+matV(2,2)+matV(3,3))/3.0d0

do i=1,3
  matUeq(i,i)=Ueq
  matVeq(i,i)=Veq
end do
! comparison of U with itself if it were isotropic 
! used in merrit's comp of suij = ccuvnorm/(ccuiso*ccviso)
! denominator
call aniso_cc(matU,matUeq,ccUiso)
call aniso_cc(matV,matVeq,ccViso)
! numerator
norm_matV=(Ueq/Veq)*matV
call aniso_cc(matU,norm_matV,cc)

suij = cc/(ccUiso*ccViso)

! now for dkon's comparison
! first divide both by trace
norm_matU=matU/(3.0d0*Ueq)
norm_matV=matV/(3.0d0*Veq)

!print '(A5,1x,3F10.5)', 'nmatU', matU(:,1)
!print '(A5,1x,3F10.5)', 'nmatU', matU(:,2)
!print '(A5,1x,3F10.5)', 'nmatU', matU(:,3)
!print '(A5,1x,3F10.5)', 'nmatV', matV(:,1)
!print '(A5,1x,3F10.5)', 'nmatV', matV(:,2)
!print '(A5,1x,3F10.5)', 'nmatV', matV(:,3)

call aniso_cc(norm_matU,norm_matV,cc)
! compute ccstr where the major and minor are switched for one of the
! eigenvalues
call aniso_ccstr(norm_matU,norm_matV,ccstr)

ccmod = (cc - ccstr)/(one-ccstr)

end subroutine

subroutine aniso_cc (matU,matV,cc)
use mod_linalg,             only: lapack_eig,lapack_pseudo
use mod_math,               only: determ3,valvec_inv
real(dp), intent(in),allocatable  :: matU(:,:),matV(:,:)
real(dp), intent(out) :: cc
real(dp), allocatable :: invU(:,:),invV(:,:)
real(dp), allocatable :: valsU(:), vecsU(:,:)
real(dp), allocatable :: valsV(:), vecsV(:,:)
integer :: ialloc,nzeros

allocate( invU(3,3),invV(3,3), stat=ialloc)

invU=zero
invV=zero

call lapack_eig(matU,valsU,vecsU)
call lapack_eig(matV,valsV,vecsV)
call valvec_inv(valsU,vecsU,invU,nzeros)
!call valvec_inv(valsU,vecsV,invU,nzeros)
call valvec_inv(valsV,vecsV,invV,nzeros)

cc    = dsqrt(dsqrt(determ3(invU)*determ3(invV))) &
      / dsqrt((one/8.0d0)*determ3(invU+invV))

end subroutine

subroutine aniso_ccstr (matU,matV,ccstr)
use mod_linalg,             only: lapack_eig,lapack_pseudo
use mod_math,               only: determ3,valvec_inv
real(dp), intent(in),allocatable  :: matU(:,:),matV(:,:)
real(dp), intent(out) :: ccstr
real(dp), allocatable :: invU(:,:),invV(:,:), ordered(:,:)
real(dp), allocatable :: valsU(:), vecsU(:,:),tmpvalsU(:)
real(dp), allocatable :: valsV(:), vecsV(:,:),tmpvalsV(:)
integer :: i,ialloc,nzeros,minv,maxv,medv,minu,maxu,medu

allocate( invU(3,3),invV(3,3),tmpvalsV(3),tmpvalsU(3),ordered(3,3), stat=ialloc)

invU=zero
invV=zero

call lapack_eig(matU,valsU,vecsU)
call lapack_eig(matV,valsV,vecsV)
maxv = maxloc(valsv, dim=1)
minv = minloc(valsv, dim=1)
medv = 6 - maxv - minv

maxu = maxloc(valsu, dim=1)
minu = minloc(valsu, dim=1)
medu = 6 - maxu - minu

!return
! switch major with minor for V
! making sure the order is correct
tmpvalsV(1)  = valsV(maxv)
tmpvalsV(2)  = valsV(medv)
tmpvalsV(3)  = valsV(minv)
tmpvalsU(1)  = valsU(minu)
tmpvalsU(2)  = valsU(medu)
tmpvalsU(3)  = valsU(maxu)

ordered(:,1) = vecsV(:,minv)
ordered(:,2) = vecsV(:,medv)
ordered(:,3) = vecsV(:,maxv)

call valvec_inv(tmpvalsU,ordered,invU,nzeros)
call valvec_inv(tmpvalsV,ordered,invV,nzeros)

!print '(A5,1x,3F10.5)', 'invU', invU(:,1)
!print '(A5,1x,3F10.5)', 'invU', invU(:,2)
!print '(A5,1x,3F10.5)', 'invU', invU(:,3)
!print '(A5,1x,3F10.5)', 'invV', invV(:,1)
!print '(A5,1x,3F10.5)', 'invV', invV(:,2)
!print '(A5,1x,3F10.5)', 'invV', invV(:,3)

ccstr    = dsqrt(dsqrt(determ3(invU)*determ3(invV))) &
      / dsqrt((one/8.0d0)*determ3(invU+invV))

end subroutine

end module
