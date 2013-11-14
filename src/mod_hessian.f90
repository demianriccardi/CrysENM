module mod_hessian
!
! updated: DMR 05-11-2008
!
! super trimmed down version of the hessian module.  uses the kirchoff matrix to build the hessian
!
! subroutines to construct hessian from kirchoff
!
!CrysFML modules
use cfml_globaldeps,                   only: dp,sp
use cfml_Atom_typedef,                only: Atom_List_Type
use cfml_crystal_metrics,              only: Crystal_Cell_Type
!external modules
use mod_types,                only: sparse,uc_list_type,asym_list_type,Neigh_List_Type,inp_par
!external modules
use mod_constants,            only: one,two,zero,pi
!use mod_inout
use mod_inout,                only: full_to_sparse,sparse_to_full,matrx2d_init

implicit none

! public subroutines: 
public :: kirchess_run_new,kirchess 
! private subroutines 
private :: wrt_hess

interface fullhess_massweigh
  module procedure fullhess_massweigh_d,fullhess_massweigh_z
end interface

contains

subroutine kirchess_new(input,asyms,unit_cell, &
                            kirchoff_coor,hessian_coor)
use mod_inout,  only: full_to_sparse,sparse_to_full,matrx2d_init
use mod_kirchoff,  only : kirchoff_run
use mod_crysbuild, only : mass_weightit
use mod_hessblock, only : sprnghessblock
use mod_math,      only : distance
type(inp_par)                 ,intent(in out)   :: input
type(asym_List_Type)              ,intent(in)   :: asyms
type(asym_list_Type)              ,intent(in)   :: unit_cell
type(sparse)                  ,intent(out)  :: kirchoff_coor,hessian_coor
type(sparse)                                :: mhessian_coor
real(dp), allocatable                       :: hessian(:,:)
! internal
real(dp)                          :: rcut,dist,hessblock(3,3),xyzi(3),xyzj(3),dxyz(3)
integer                           :: i,ii,jj,ll,iii,jjj,kkk,lll,mm,k,kk,h,natom,ndim,ialloc
real(sp)                          :: time1,time2
integer :: nneighs,naunits,natoms,chaini,chainj,iaunit,jaunit,ichain,jchain
real(dp), dimension(3) :: avec,bvec,cvec,vec
real(dp) :: dottest
! aunit and chain are used in concert to determine whether the atom is indeed on the same chain

chaini = 1
chainj = 1
call cpu_time(time1)

nneighs = asyms%nneighs
naunits = asyms%naunits
natoms  = asyms%aunit(1)%natoms


call kirchoff_run(input,asyms,unit_cell,kirchoff_coor)

ndim=3*kirchoff_coor%ndim
rcut = input%rcut_start
if (allocated(hessian) ) deallocate(hessian)
allocate(hessian(ndim,ndim), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

hessian=zero

do i = 1, kirchoff_coor%nnzero

  ii   = kirchoff_coor%rows(i)
  jj   = kirchoff_coor%columns(i)

  if(input%multint) then
    
    if ((trim(input%tpbc) .eq. "pbc" ) .or. (trim(input%tpbc) .eq. "bvk")) then
      iaunit = unit_cell%aunit(1)%atom(ii)%iaunit       
      jaunit = unit_cell%aunit(1)%atom(jj)%iaunit       
      ichain = unit_cell%aunit(1)%atom(ii)%ichain       
      jchain = unit_cell%aunit(1)%atom(jj)%ichain       

      xyzi = unit_cell%aunit(1)%atom(ii)%x

      do iii=unit_cell%aneigh(1),unit_cell%aneigh(2)
      avec = real(iii,kind=dp)*unit_cell%avec
      do jjj=unit_cell%bneigh(1),unit_cell%bneigh(2)
      bvec = real(jjj,kind=dp)*unit_cell%bvec
      do kkk=unit_cell%cneigh(1),unit_cell%cneigh(2)
        cvec = real(kkk,kind=dp)*unit_cell%cvec
        vec = avec + bvec + cvec
! this is for REACH intra vs inter
        dottest = dot_product(vec,vec)
        if (abs(dottest) .lt. 1.0d-06) then
          if(iaunit .eq. jaunit .and. ichain .eq. jchain) then
            chaini = 1; chainj = 1
          else
            chaini = 1; chainj = 2
          end if
        else
              chaini = 1; chainj = 2
        end if

        xyzj = unit_cell%aunit(1)%atom(jj)%x + vec
        dxyz = xyzj-xyzi
        dist  = distance(dxyz)

        if (dist .le. rcut) then
          if (dist .gt. zero) then
            call hess_common(input, dist,ii,jj,dxyz,   &
                                 hessian,chaini,chainj)
          end if
        end if

      end do
      end do
      end do

    else if (trim(input%tpbc) .eq. "asy" ) then

! note that this uses input%aunit.  this will use a specific asymmetric unit for the calculation
! this is the dumb way of generating the var-covar matrix for a molecule with diff rots   
      xyzi = asyms%aunit(input%aunit)%atom(ii)%x
      iaunit = asyms%aunit(input%aunit)%atom(ii)%iaunit       
      ichain = asyms%aunit(input%aunit)%atom(ii)%ichain       

      do h = 1, asyms%naunits
        jchain = asyms%aunit(h)%atom(jj)%ichain       
        jaunit = asyms%aunit(h)%atom(jj)%iaunit       
      do iii=asyms%aneigh(1),asyms%aneigh(2)
      avec = real(iii,kind=dp)*asyms%avec
      do jjj=asyms%bneigh(1),asyms%bneigh(2)
      bvec = real(jjj,kind=dp)*asyms%bvec
      do kkk=asyms%cneigh(1),asyms%cneigh(2)
        cvec = real(kkk,kind=dp)*asyms%cvec
        vec = avec + bvec + cvec
! this is for REACH intra vs inter
        dottest = dot_product(vec,vec)
        if (abs(dottest) .lt. 1.0d-06) then
          if(iaunit .eq. jaunit .and. ichain .eq. jchain) then
            chaini = 1; chainj = 1
          else
            chaini = 1; chainj = 2
          end if
        else
              chaini = 1; chainj = 2
        end if

        xyzj = asyms%aunit(h)%atom(jj)%x + vec
        dxyz = xyzj-xyzi
        dist  = distance(dxyz)

        if (dist .le. rcut) then
          if (dist .gt. zero) then
            call hess_common(input, dist,ii,jj,dxyz,   &
                                 hessian,chaini,chainj)
          end if
        end if

      end do
      end do
      end do
      end do
    else
      stop 'using PBC for non PBC case, please use: asy or pbc'
    end if

  else
    xyzi = asyms%aunit(input%aunit)%atom(ii)%x
    xyzj = asyms%aunit(input%aunit)%atom(jj)%x
    dxyz = xyzj-xyzi
    dist  = distance(dxyz)
          
    if (dist .le. rcut) then
       if (dist .gt. zero) then
          call hess_common(input, dist,ii,jj,dxyz,   &
                                 hessian,1,1)
       end if
    end if

  end if

end do

call full_to_sparse(hessian,hessian_coor)
deallocate(hessian)

call sparsehess_massweigh(input,unit_cell%aunit(1),hessian_coor,"d")

input%natoms=kirchoff_coor%ndim
input%enm_sparsity = hessian_coor%sparsity
! compute the number of modes!
call nfrq_comp(input)

if (trim(input%genm) .eq."enm") input%sparsity=input%enm_sparsity

call wrt_hess(input,hessian_coor)

call cpu_time(time2)
print *, "time constructing hessian", time2-time1

end subroutine

subroutine sparsehess_massweigh(input,atoms,matin,zd)
! atom_mass_constr will return masses of atoms, which may be one, mass_atom, or mass_residue
! depends on input so this can always be run
use mod_types,     only : inp_par,sparse
use mod_crysbuild, only : atom_mass_constr
type(inp_par)       , intent(in)     :: input
type(atom_list_type), intent(in)     :: atoms
type(sparse)        , intent(inout)  :: matin
character(len=1)    , intent(in)     :: zd
real(dp), allocatable :: atmasses(:),xyzmasses(:) 
real(dp) :: weight
integer :: i

call atom_mass_constr(input,atoms,atmasses)
call mass_at2xyz(atmasses,xyzmasses)

xyzmasses = one/dsqrt(xyzmasses)

do i=1,matin%nnzero
  weight = xyzmasses(matin%rows(i))*xyzmasses(matin%columns(i))
  if (zd .eq. "d") then
    matin%values(i)  = weight*matin%values(i)
  else if (zd .eq. "z") then
    matin%zvalues(i) = weight*matin%zvalues(i)
  end if
end do

end subroutine

subroutine fullhess_massweigh_d(input,atoms,matin)
! atom_mass_constr will return masses of atoms, which may be one, mass_atom, or mass_residue
! depends on input so this can always be run
use mod_types,     only : inp_par,sparse
use mod_crysbuild, only : atom_mass_constr
type(inp_par)       , intent(in)     :: input
type(atom_list_type), intent(in)     :: atoms
real(dp),allocatable, intent(inout)  :: matin(:,:)
real(dp), allocatable :: atmasses(:),xyzmasses(:)
real(dp) :: weight
integer :: i,j

call atom_mass_constr(input,atoms,atmasses)
call mass_at2xyz(atmasses,xyzmasses)

xyzmasses = one/dsqrt(xyzmasses)

do i=1,size(matin,1)
  do j=i,size(matin,2)
    weight = xyzmasses(i)*xyzmasses(j)
    matin(i,j) = weight*matin(i,j)
    matin(j,i) = matin(i,j)
  end do
end do

end subroutine

subroutine fullhess_massweigh_z(input,atoms,matin)
! atom_mass_constr will return masses of atoms, which may be one, mass_atom, or mass_residue
! depends on input so this can always be run
use mod_types,     only : inp_par,sparse
use mod_crysbuild, only : atom_mass_constr
type(inp_par)       , intent(in)     :: input
type(atom_list_type), intent(in)     :: atoms
complex(dp),allocatable, intent(inout)  :: matin(:,:)
real(dp), allocatable :: atmasses(:),xyzmasses(:)
real(dp) :: weight
integer :: i,j

call atom_mass_constr(input,atoms,atmasses)
call mass_at2xyz(atmasses,xyzmasses)

xyzmasses = one/dsqrt(xyzmasses)

do i=1,size(matin,1)
  do j=i,size(matin,2)
    weight = xyzmasses(i)*xyzmasses(j)
    matin(i,j) = weight*matin(i,j)
    matin(j,i) = matin(i,j)
  end do
end do

end subroutine

subroutine mass_at2xyz(atmasses,xyzmasses)
real(dp), allocatable, intent(in) :: atmasses(:)
real(dp), allocatable, intent(out) :: xyzmasses(:)
integer :: i,j,nat,nxyz,iat,ier

nat = size(atmasses,1)
nxyz = 3*nat
allocate(xyzmasses(nxyz), stat=ier)
if (ier /= 0) stop "mass_at2xyz> malloc"

iat = 1
do i = 1, nat
  do j = 1, 3
    xyzmasses(iat) = atmasses(i)
    iat = iat+1
  end do
end do

end subroutine

subroutine hess_weigh_byres (atoms,matin,matout,zd)
use mod_types,     only: sparse
use mod_inout,     only: sparse_refinit
use mod_crysbuild, only: residue_masses
type(atom_list_type), intent(in)     :: atoms
type(sparse),intent(in)           :: matin
type(sparse),intent(out)          :: matout
real(dp),allocatable                 :: atmasses(:),masses(:)
real(dp) :: invrtmassi,invrtmassj,weight
character(len=1)                     :: zd
integer                              :: i,j,nnzero,ialloc,ndim

call sparse_refinit(matin,matout,zd)
call residue_masses (atoms,atmasses)
call mass_at2xyz(atmasses,masses)
deallocate(atmasses)

do i=1,matin%nnzero
  invrtmassi=one/dsqrt(masses(matin%rows(i)))
  invrtmassj=one/dsqrt(masses(matin%columns(i)))
  weight = invrtmassi*invrtmassj
if (zd .eq. "d") then
    matout%values(i)  = weight*matin%values(i)
else if (zd .eq. "z") then
    matout%zvalues(i) = weight*matin%zvalues(i)
end if

end do

end subroutine

subroutine hess_invrt_byres(atoms,invrtmasses)
use mod_crysbuild, only: residue_masses
type(atom_list_type), intent(in)     :: atoms
real(dp),allocatable, intent(out)    :: invrtmasses(:)
real(dp),allocatable                 :: atmasses(:)

call residue_masses (atoms,atmasses)
call mass_at2xyz(atmasses,invrtmasses)

invrtmasses = one/dsqrt(invrtmasses)

end subroutine

subroutine kirchess(input,molecule,asy_im,crystal,&
                          kirchoff_coor,hessian_coor)
use mod_inout,  only: full_to_sparse,sparse_to_full,matrx2d_init
use mod_kirchoff,  only : kirchoff_run
use mod_crysbuild, only : mass_weightit
use mod_hessblock, only : sprnghessblock
use mod_math,      only : distance
type(inp_par)                 ,intent(in out)   :: input
type(atom_List_Type)          ,intent(in)   :: molecule
type(uc_List_Type)            ,intent(in)   :: asy_im
type(neigh_List_Type)         ,intent(in)   :: crystal
type(sparse)                  ,intent(out)  :: kirchoff_coor,hessian_coor
type(sparse)                                :: mhessian_coor
real(dp), allocatable                       :: hessian(:,:)
! internal
real(dp)                          :: rcut,dist,hessblock(3,3),xyzi(3),xyzj(3),dxyz(3)
integer                           :: i,ii,jj,ll,mm,k,kk,h,natom,ndim,ialloc
real(sp)                          :: time1,time2
integer :: nneighs,naunits,natoms
type(crystal_cell_type) :: cell

call cpu_time(time1)

nneighs = asy_im%nneighs
naunits = asy_im%neigh(0)%naunits
natoms  = asy_im%neigh(0)%aunit(1)%natoms
kirchoff_coor%ndim = 0

print *, 'outdated'
stop
return
! build kirchoff matrix
!call kirchoff_run(input,molecule,asy_im,crystal,kirchoff_coor)
!call kirchoff_run(input,asy_im%neigh(0),kirchoff_coor)

ndim=3*kirchoff_coor%ndim
rcut = input%rcut_start
if (allocated(hessian) ) deallocate(hessian)
allocate(hessian(ndim,ndim), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

hessian=zero

do i = 1, kirchoff_coor%nnzero

  ii   = kirchoff_coor%rows(i)
  jj   = kirchoff_coor%columns(i)

  if(input%multint) then
    
    if ((trim(input%tpbc) .eq. "pbc" ) .or. (trim(input%tpbc) .eq. "bvk")) then

      xyzi = crystal%neigh(0)%atom(ii)%x

      do k=0, crystal%nneighs

    !  print *, 'multint'
        xyzj = crystal%neigh(k)%atom(jj)%x
        dxyz = xyzj-xyzi
        dist  = distance(dxyz)

        if (dist .le. rcut) then
          if (dist .gt. zero) then
            call hess_common(input, dist,ii,jj,dxyz,   &
                                 hessian)
          end if
        end if

      end do

    else if (trim(input%tpbc) .eq. "asy" ) then

! note that this uses input%aunit.  this will use a specific asymmetric unit for the calculation
! this is the dumb way of generating the var-covar matrix for a molecule with diff rots   
      xyzi = asy_im%neigh(0)%aunit(input%aunit)%atom(ii)%x

      do k=0, nneighs
        do h=1, naunits

          !xyzj = crystal%neigh(k)%atom(jj+(h-1)*kirchoff_coor%ndim)%x
          xyzj = asy_im%neigh(k)%aunit(h)%atom(jj)%x
          dxyz = xyzj-xyzi
          dist  = distance(dxyz)

          if (dist .le. rcut) then
            if (dist .gt. zero) then
              call hess_common(input, dist,ii,jj,dxyz,   &
                                 hessian)
            end if
          end if

        end do
      end do

    else
      stop 'using PBC for non PBC case, please use: asy or pbc'
    end if

  else
    xyzi = molecule%atom(ii)%x
    xyzj = molecule%atom(jj)%x
    dxyz = xyzj-xyzi
    dist  = distance(dxyz)
          
    if (dist .le. rcut) then
       if (dist .gt. zero) then
          call hess_common(input, dist,ii,jj,dxyz,   &
                                 hessian)
       end if
    end if

  end if

end do

call full_to_sparse(hessian,hessian_coor)
deallocate(hessian)

if(trim(input%weigh) .eq. 'mass') then
    call mass_weightit("enm",crystal%neigh(0),hessian_coor,mhessian_coor,"d") 
    hessian_coor%values  = mhessian_coor%values
end if

input%natoms=kirchoff_coor%ndim
input%enm_sparsity = hessian_coor%sparsity
! compute the number of modes!
call nfrq_comp(input)

if (trim(input%genm) .eq."enm") input%sparsity=input%enm_sparsity

call wrt_hess(input,hessian_coor)

call cpu_time(time2)

end subroutine

subroutine kirchess_run_new(input,unit_cell,kirchoff_coor,hessian_coor,&
                          vcov,states)
use mod_types,     only: sparse,dos,dos_init
type(inp_par)             ,intent(in out)   :: input
type(asym_list_type)             , intent(in)      :: unit_cell
type(sparse)                  ,intent(in)  :: kirchoff_coor,hessian_coor
real(dp), allocatable         ,intent(out)  :: vcov(:,:)
type(dos)                     ,intent(out)  :: states
integer  :: ialloc,i,j
real(dp), allocatable :: rtmasses(:)
! internal

if (trim(input%tpbc) .ne. "bvk") then

  if (trim(input%genm) .eq. "gnm") call vcov_an(input,kirchoff_coor,vcov,states)
  if (trim(input%genm) .eq. "enm") then 
    call vcov_an(input,hessian_coor,vcov,states)
    call fullhess_massweigh(input,unit_cell%aunit(1),vcov)
  end if
end if

end subroutine

subroutine kirchess_valvec(input,kirchoff_coor,hessian_coor,&
                                 val,vec)
use mod_types,     only: uc_List_Type,neigh_List_Type,sparse
use mod_kirchoff,  only: kirchoff_run
type(inp_par)             ,intent(in out)   :: input
type(sparse)                  ,intent(in)   :: kirchoff_coor,hessian_coor
real(dp), allocatable         ,intent(out)  :: val(:),vec(:,:)
integer  :: ialloc
! internal


  if (trim(input%genm) .eq. "gnm")  call eigenal(input,kirchoff_coor,val,vec)
  if (trim(input%genm) .eq. "enm")  call eigenal(input,hessian_coor,val,vec)

end subroutine

subroutine vcov_an(input,matx,vcov,states)
use mod_types,     only: inp_par,sparse,dos,dos_init,dos_histgen,dos_append
use mod_inout,     only: sparse_to_full
use mod_arpack,    only: mod_arpack_dsdrv1_sparse
use mod_linalg,    only: lapack_eig,mod_linalg_pseudo
type(inp_par),        intent(in out) :: input
type(sparse),         intent(in)     :: matx
real(dp), allocatable                :: matx_full(:,:)
real(dp), allocatable,intent(out)    :: vcov(:,:)
type(dos),             intent(out)   :: states
real(dp), allocatable                :: vals(:),valsthz(:),vecs(:,:)
integer                              :: i,nzeros,ialloc 

if (allocated(vcov)) deallocate(vcov)
if (allocated(vals)) deallocate(vals)
if (allocated(vecs)) deallocate(vecs)
allocate(vals(input%nfrqs), &
         vecs(matx%ndim,input%nfrqs),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"
! ok...  soon be done with density of states!

    if (trim(input%eigtyp) .eq. "lapack") then
      call sparse_to_full(matx,matx_full)
      call lapack_eig(matx_full, vals, vecs )
    else
      call mod_arpack_dsdrv1_sparse(input,input%nfrqs,matx,vals,vecs)
    end if

     if (trim(input%eigtyp) .eq. "lapack") then
      call dos_calculator(75,vals,states)
    else
      call dos_calculator(input%nfrqs/2,vals,states)
    end if

call dvcov(input,vals,vecs,vcov,nzeros)

if ((nzeros .ne. 1) .and. (nzeros .ne. input%nnzeros )) then
  print '(A57)', &
         'WARNING> number of zero eigenvalues greater than expected'
  print '(A16,1x,I4,1x,A9,1x,I2,1x,A5,1x,F10.3)', &
         'WARNING> #zeros:',nzeros, 'expected:', input%nnzeros, 'rcut:',input%rcut_start
end if

deallocate(vals,vecs)

end subroutine

subroutine dos_calculator(nbin,vals,states)
use mod_types,     only: dos,dos_init,dos_histgen,dos_append
integer              , intent(in)    :: nbin               
real(dp), allocatable, intent(in)    :: vals(:)
type(dos),             intent(out)   :: states
real(dp), allocatable                :: valsthz(:)

call dos_init(nbin,states)
call eig_to_THz(vals,valsthz)
call dos_append(valsthz,states)
call dos_histgen(states)

end subroutine

subroutine eigenal(input,matx,vals,vecs)
use mod_types,     only: inp_par,sparse
use mod_inout,     only: sparse_to_full
use mod_arpack,    only: mod_arpack_dsdrv1_sparse
use mod_linalg,    only: lapack_eig,mod_linalg_pseudo
type(inp_par),        intent(in out) :: input
type(sparse),         intent(in)     :: matx
real(dp), allocatable                :: matx_full(:,:)
real(dp), allocatable,intent(out)    :: vals(:),vecs(:,:)
real(dp), allocatable                :: valsthz(:)
integer                              :: nzeros,ialloc 

if (allocated(vals)) deallocate(vals)
if (allocated(vecs)) deallocate(vecs)
allocate(vals(input%nfrqs), &
         vecs(matx%ndim,input%nfrqs),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

if (trim(input%eigtyp) .eq. "lapack") then
  call sparse_to_full(matx,matx_full)
  call lapack_eig(matx_full, vals, vecs )
else
  call mod_arpack_dsdrv1_sparse(input,input%nfrqs,matx,vals,vecs)
end if

end subroutine

subroutine eig_to_THz(eig,frq)
! convert eigenvalues (square of circ freq) to freqs in THz 
use mod_constants, only : CONVFRQ, CLIGHT
real(dp),        allocatable, intent(in)   :: eig(:)
real(dp),        allocatable, intent(out)  :: frq(:)
integer :: i,ialloc

if (allocated(frq)) deallocate (frq)
allocate(frq(size(eig)),stat=ialloc)

do i=1,size(eig)
    if(real(eig(i)) .lt. 1.0D-08) then
      frq(i)=zero
    else
      frq(i)=dsqrt(eig(i))*CONVFRQ*CLIGHT*1.0D02/1.0D12
      ! CLIGHT is in m/s 1.0D2 converts CLIGHT to cm/s
      ! this leaves hertz...  /1.0D12 converts to terahertz 
    end if
end do

end subroutine

subroutine nfrq_comp(input)
type(inp_par) , intent(in out) :: input
integer                        :: ndim, natoms 
real(sp)                       :: val

!DMR adding option to compute specific number of modes (up through last_mode)
!    if last mode corresponds to a mode within the first 10% we'll switch to arpack
!    and calculate through last_mode
!    do the same as before if:
!    1. first or last < 0
!    2. last  < first 
!    3. first = last = 0

print *, "eigtyp requested:", trim(input%eigtyp)

natoms = input%natoms
call expected_nzeros(input)

if (trim(input%genm) .eq."enm") then 
  ndim = 3*natoms
else if (trim(input%genm) .eq."gnm") then
  ndim          = natoms
else if (trim(input%genm) .eq."bnm") then
  ndim = 6*input%nblocks
end if

if (input%first_mode .le. 0 .or. input%first_mode .gt. input%last_mode) input%first_mode = 1
if (input%last_mode .le. 0 .or. input%last_mode .lt. input%first_mode .or. input%last_mode .gt. ndim) input%last_mode = ndim

if (trim(input%eigtyp) .eq. "lapack_lapack") then
  print *,"using lapack and modes:",input%first_mode,"to ",input%last_mode
  input%eigtyp = "lapack" ; return
end if


if( (input%first_mode .le. 0 .or. input%last_mode .le. 0) .or. &
     (input%first_mode .gt. input%last_mode) .or.  &
     (trim(input%eigtyp) .eq. "arpack_arpack") ) then
  
  if (trim(input%eigtyp) .eq. "arpack_arpack") input%eigtyp= "arpack"
  if (trim(input%eigtyp) .eq. "arpack") print *, 'using percentage to calculate nfrqs for arpack:', input%perc

  if (trim(input%eigtyp) .eq. "lapack") then

    input%nfrqs=ndim
    input%last_mode = ndim

  print *, 'lapack number of freqs=',input%nfrqs, 'first mode'
  print *, 'first mode: ', input%first_mode
  print *, 'last mode: ', input%last_mode
    return
  end if

!take percentage of nonzero eigenvals
  input%nfrqs=nint(real(ndim-input%nnzeros)*input%perc)
  input%nfrqs=input%nfrqs+input%nnzeros ! shift percentage above the zero eigenvals 
  input%nfrqs=min(input%maxfrqs,input%nfrqs)
!set last mode
  input%last_mode = input%nfrqs

else

! if le 10% go with arpack, else use lapack and all the frqs
  val = real(input%last_mode)/real(ndim)
!switch to arpack if val is < 0.1
  if ((val .le. 0.10) .and. (trim(input%eigtyp) .ne. "lapack_lapack"))  then
    input%eigtyp = "arpack"
    if(input%last_mode .le. input%nnzeros) then
      input%nfrqs = input%nnzeros
    else
      input%nfrqs = input%last_mode
    end if
    print '(A33,1x,F10.2)', "nfrq_comp> fraction of the modes:", val
    print '(A27,1x,A9,I4)', "nfrq_comp> eigtyp and nfrq:", input%eigtyp,input%nfrqs
  else
    input%eigtyp = "lapack"
    input%nfrqs = ndim
    print '(A33,1x,F10.2)', "nfrq_comp> fraction of the modes:", val
    print '(A27,1x,A9,I4)', "nfrq_comp> eigtyp and nfrq:", input%eigtyp,input%nfrqs
  end if 

 
end if

print *, 'number of freqs=',input%nfrqs

end subroutine

subroutine expected_nzeros(input)
type(inp_par) , intent(in out) :: input

if (trim(input%genm) .eq."enm") then 

  if ((trim(input%tpbc) .eq. "pbc") .or. &
      (trim(input%tpbc) .eq. "bvk") .or. &
      (trim(input%tpbc) .eq. "asy")) then 
    input%nnzeros = 3
  else
    input%nnzeros = 6
  end if


else if (trim(input%genm) .eq."gnm") then

  input%nnzeros = 1

end if


end subroutine
 
subroutine dvcov(input,vals,vecs,vcov,nzeros)
use mod_constants, only: Rkcal
use mkl95_blas, only : gemm
! take in the eigenvalues: units rad^2 kcal mol^-1 amu^-1 angs^2
! 
type(inp_par) , intent(in) :: input
real(dp),        allocatable,intent(in)     :: vals(:),vecs(:,:)
real(dp),        allocatable,intent(out) :: vcov(:,:)
integer,                     intent(out)    :: nzeros
real(dp)                  :: temp
real(dp)        :: tmpval
real(dp)        :: kT
integer ::   nvals, ndim, nowvals
integer ::   i,j,k,nfirst,nlast,ier
real(dp),allocatable :: mat1(:,:),mat2(:,:)


temp = input%temperature
kT = temp*Rkcal
nvals  = size(vals)
ndim   = size(vecs(:,1))
nzeros = 0

do i = 1, nvals
  if (abs(vals(i)).lt.1.0D-07) then
    nzeros = nzeros+1
  end if
end do

nfirst = input%first_mode
if (nzeros >= nfirst) nfirst = nzeros+1
nlast = input%last_mode

nowvals = nlast - nfirst + 1

! somewhere else and added to here, so to avoid problems, we'll use matv temp
allocate(mat2(ndim,nowvals),vcov(ndim,ndim),stat = ier)
!allocate(mat1(ndim,nowvals),mat2(ndim,nowvals),vcov(ndim,ndim),stat = ier)
if (ier /= 0 ) stop "dvcov> malloc"

!mat1 = vecs(:,nfirst:nlast) 
mat2 = vecs(:,nfirst:nlast) 

! first step: mult mat2 by kt/eig
j = 1
do i = nfirst,nlast
  mat2(:,j) = mat2(:,j)*dsqrt(kt/vals(i))
  j = j+1
end do

print *, 'computing the variance-covariance matrix'
print *, 'from first mode: ', nfirst, 'to last mode: ', nlast  
call gemm(mat2,mat2,vcov,'n','t',1.0d0,0.0d0)
print *, 'done computing variance-covariance matrix'

!if(trim(input%weigh) .eq. 'mass' .or. trim(input%weigh) .eq. 'byatom')  then
!    call vcov_mass_unweigh(unit_cell%aunit(1),vcov)
!elseif(trim(input%weigh) .eq. 'byres') then
!    print *, "multiplying vcov by masses of residues"
!    call hess_weigh_byres(unit_cell%aunit(1),vcov,"d")
!end if

end subroutine

subroutine dvcov_old(input,vals,vecs,vcov,nzeros)
use mod_constants, only: Rkcal
! take in the eigenvalues: units rad^2 kcal mol^-1 amu^-1 angs^2
! 
type(inp_par) , intent(in) :: input
real(dp),        allocatable,intent(in)     :: vals(:),vecs(:,:)
real(dp),        allocatable,intent(in out) :: vcov(:,:)
integer,                     intent(out)    :: nzeros
real(dp)                  :: temp
real(dp)        :: tmpval
real(dp)        :: kT
integer ::   nvals, ndim
integer ::   i,j,k

temp = input%temperature
kT = temp*Rkcal
nvals  = size(vals)
ndim   = size(vecs(:,1))
nzeros = 0

print *, 'computing the variance-covariance matrix'
print *, 'from first mode: ', input%first_mode, 'to last mode: ', input%last_mode  

! this is really sloooowwww
! we should rewrite this as:
!    1. determine nonzero eig vecs
!    2. matmat multiplications using blas

do i=input%first_mode,input%last_mode  
  if (abs(vals(i)).gt.1.0D-08) then
    do j=1,ndim
      do k=1,ndim
        tmpval    = vecs(k,i)*vecs(j,i)
        vcov(k,j) = vcov(k,j) + tmpval*kt/vals(i)
      end do
    end do
  else
    nzeros = nzeros + 1
  end if
end do

print *, 'done computing the variance-covariance matrix'

end subroutine

subroutine hess_common(input,dist,ii,jj,dxyz, &
                          hessian,chaini,chainj)
use mod_hessblock, only: sprnghessblock,reach_hessblock_temperat,reach_hessblock
type(inp_par),        intent(in)            :: input
integer ,             intent(in)            :: ii,jj
integer , optional,   intent(in)            :: chaini,chainj
real(dp),             intent(in)            :: dist,dxyz(3)
real(dp),allocatable, intent(in out)        :: hessian(:,:)
integer                                     :: i,j,ll,mm
real(dp)                                    :: hessblock(3,3)


i     = (ii-1)*3+1
j     = (jj-1)*3+1

select case (trim(input%fctyp))
  case("reach")
    if (present(chaini) .and. present(chainj)) then
      call reach_hessblock(input,ii,jj,chaini,chainj,dist,dxyz,&   ! this call is not perfect
                         hessblock)                              ! for the case of missing resids
    else                                                         ! ii,jj should be ires,jres
      print *, "hess_common called REACH without chain info"
    end if
  case("treach")
    if (present(chaini) .and. present(chainj)) then
      call reach_hessblock_temperat(input,ii,jj,chaini,chainj,dist,dxyz,&   ! this call is not perfect
                         hessblock)                              ! for the case of missing resids
    else                                                         ! ii,jj should be ires,jres
      print *, "hess_common called REACH without chain info"
    end if
  case default
    call sprnghessblock(input,dist,dxyz, &
                         hessblock)
end select

! hess dyn
!above
do ll=0,2
  do mm=0,2
    hessian(i+ll,j+mm)  = hessian(i+ll,j+mm) + hessblock(1+ll,1+mm) 
  end do
end do

!diagonal
do ll=0,2
  do mm=0,2
    hessian(i+ll,i+mm)                = hessian(i+ll,i+mm) - hessblock(1+ll,1+mm)
    if (ii.ne.jj)  hessian(j+ll,j+mm) = hessian(j+ll,j+mm) - hessblock(1+ll,1+mm)
  end do
end do

end subroutine



subroutine wrt_hess(input,hessian)
type(inp_par), intent(in) :: input
type(sparse),  intent(in) :: hessian
real(dp)             :: summ,trace
integer :: ii,jj,ier,degfred

summ=0.0d0

if (input%print_level .ge. 3) then

open(unit=13,file=trim(input%fileroot)//"-hessian.txt", &
     status='unknown',action='write', &
     iostat=ier)

write(13,'(I6,F10.2,F10.4)') hessian%ndim,input%rcut_start,input%fcnst

trace = zero

do ii = 1, hessian%nnzero
    
 if (hessian%rows(ii).eq.hessian%columns(ii)) trace=trace+hessian%values(ii)

end do

write(13,*) "Trace:",trace

do ii = 1, hessian%nnzero
    write(13,'(2I5,F15.5)') hessian%rows(ii),hessian%columns(ii),hessian%values(ii)
    

!    if (hessian%rows(ii).ne.hessian%columns(ii))  write(13,'(2I5,F15.5)')  &
!                                          hessian%columns(ii),hessian%rows(ii),hessian%values(ii)
    summ = summ+hessian%values(ii) 
    if (hessian%rows(ii).ne.hessian%columns(ii)) summ=summ+hessian%values(ii)
end do

print *, 'sum of all nonzero hessian els', summ

close(13)
end if

end subroutine

subroutine hess_crys(input,crystal,kirchoff_coor, &
                                          hessian_coor)
! Compute the crystal hessian: takes in kirchoff(sparse) output hessian(sparse)
use mod_inout,  only: full_to_sparse,sparse_to_full,matrx2d_init
use mod_crysbuild, only : mass_weightit
use mod_hessblock, only : sprnghessblock
use mod_math,      only : distance
type(inp_par)                 ,intent(in out) :: input
type(neigh_List_Type)         ,intent(in)     :: crystal
type(sparse)                  ,intent(in)     :: kirchoff_coor
type(sparse)                  ,intent(out)    :: hessian_coor
type(sparse)                                  :: mhessian_coor
real(dp), allocatable                         :: hessian(:,:)
! internal
real(dp)                          :: rcut,dist,hessblock(3,3),xyzi(3),xyzj(3),dxyz(3)
integer                           :: i,ii,jj,k,ndim,ialloc,tmp_print

ndim=3*kirchoff_coor%ndim
rcut = input%rcut_start
if (allocated(hessian) ) deallocate(hessian)
allocate(hessian(ndim,ndim), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

hessian=zero

do i = 1, kirchoff_coor%nnzero

  ii   = kirchoff_coor%rows(i)
  jj   = kirchoff_coor%columns(i)

  xyzi = crystal%neigh(0)%atom(ii)%x

  do k=0, crystal%nneighs

    !print *, 'multint'
    xyzj = crystal%neigh(k)%atom(jj)%x
    dxyz = xyzj-xyzi
    dist  = distance(dxyz)

    if (dist .le. rcut) then
      if (dist .gt. zero) then
        call hess_common(input, dist,ii,jj,dxyz,   &
                            hessian)
      end if
    end if

  end do
end do

call full_to_sparse(hessian,hessian_coor)
deallocate(hessian)

if(trim(input%weigh) .eq. 'mass') then
    call mass_weightit("enm",crystal%neigh(0),hessian_coor,mhessian_coor,"d")
    hessian_coor%values  = mhessian_coor%values
end if

input%enm_sparsity = hessian_coor%sparsity

if (trim(input%genm) .eq."enm") input%sparsity=input%enm_sparsity

tmp_print = input%print_level
input%print_level = 5
call wrt_hess(input,hessian_coor)
input%print_level = tmp_print

end subroutine

end module

