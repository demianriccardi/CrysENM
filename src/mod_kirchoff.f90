module mod_kirchoff
!
! last update: DMR 04-07-2008
! subroutines to construct kirchoff matrix
! 
!   input: 
!       molecule config : atom_list_type, neigh_list_type, uc_list_type
!       pbc type        : def, asy, pbc
!       fcnst,rcut
!   
!   output:
!       kirchoff matrix (sparse storage), optional( vcov, g(w) ) if vibrational analysis requested
!
!  The kirchoff matrix is the basic connectivity matrix; there are 3 types of "boundary conditions"
!     1. default: single molecule with no crystal environment
!     2. asy:     periodic with respect to assymetric unit
!     3. pbc:     periodic with respect to unit cell
!
!
!       1. denoted by NAME_un
!       2. denoted by NAME_asy
!       3. denoted by NAME_crys; periodic boundary conditions where primary unit cell is surrounded by 
!          nearest neighbors as constructed with crysFML.  The matrix dimensions correspond to the
!          unit cell; the elements of the matrix for two atoms are added for both intra- and inter- unit
!          cell interactions
!
!       ie. kirchoff for a short one dimensional chain of 4 atoms:
!             1 -1  0  0 
!            -1  2 -1  0
!             0 -1  2 -1
!             0  0 -1  1
!       becomes, with PBC:
!             2 -1  0 -1 
!            -1  2 -1  0
!             0 -1  2 -1
!            -1  0 -1  2
!
!  arpack is used to compute limited number of eigenvals and vecs
!  ie.    mod_arpack_dsdrv1(degfred,znfrq,dynmat,zvals,zvecs)
!
!CrysFML modules
use cfml_globaldeps,                   only: dp,sp
use cfml_Atom_typedef,                only: Atom_List_Type
use cfml_crystal_metrics,              only: Crystal_Cell_Type
!external modules
use mod_types,                only: inp_par,sparse,asym_list_type,uc_list_type,Neigh_List_Type
use mod_math,                 only: dist_sqr
use mod_constants,            only: one,two,zero,pi
use mod_inout,  only: sparse_deinit,full_to_sparse,sparse_distances,matrx2d_init
implicit none
! public subroutines: 
public :: kirchoff_run, kirchoff_crys
! private subroutines 
private :: chkwrt_kirchoff,kirchoff_un, kirchoff_asy

contains

subroutine kirchoff_run(input,asyms,unit_cell,&
                          kirchoff_coor)
use mod_crysbuild, only : mass_weightit,asymlist_to_atomlist
type(inp_par)              ,intent(in out)  :: input
type(asym_List_Type)       ,intent(in)      :: asyms
type(asym_list_Type)       ,intent(in)      :: unit_cell
type(sparse)               ,intent(out)     :: kirchoff_coor    
type(sparse)                             :: mkirchoff_coor
real(dp), allocatable                    :: kirchoff(:,:),distances(:,:),tab_dxyz(:,:,:)
! internal
real(dp)                                 :: dcut 
integer                                  :: i,ii,jj,natom,ndim
real(sp)                                 :: time1,time2

call cpu_time(time1)
ndim=asyms%aunit(1)%natoms
!build kirchoff
if ((trim(input%tpbc) .eq. "pbc") .or. (trim(input%tpbc) .eq. "bvk")) then
  ndim=unit_cell%natoms
  call kirchoff_crys(input,unit_cell, &
                          kirchoff,distances,tab_dxyz)
else if (trim(input%tpbc) .eq. "asy") then
  call kirchoff_asy(input,asyms,&
                          kirchoff,distances,tab_dxyz)
else
  call kirchoff_un(input,asyms%aunit(1),kirchoff,distances,tab_dxyz)
end if

call full_to_sparse(kirchoff,kirchoff_coor)
call sparse_distances(distances,tab_dxyz,kirchoff_coor)

call chkwrt_kirchoff(input,kirchoff_coor)    
  
deallocate(kirchoff,distances,tab_dxyz)
call sparse_deinit(mkirchoff_coor)

if (trim(input%genm) .eq. "gnm") then
!  multiply kirchoff by fcnst
  do i=1,kirchoff_coor%nnzero
    kirchoff_coor%values(i)=kirchoff_coor%values(i)*input%fcnst
  end do

end if

input%gnm_sparsity = kirchoff_coor%sparsity
input%sparsity = input%gnm_sparsity

call cpu_time(time2)

end subroutine

subroutine mat_distances(input,asyms,unit_cell,&
                         distances)
type(inp_par)              ,intent(in out)  :: input
type(asym_List_Type)       ,intent(in)      :: asyms, unit_cell
real(dp), allocatable                       :: distances(:,:)
! internal
real(sp)                                 :: time1,time2

call cpu_time(time1)
!build kirchoff

if ((trim(input%tpbc) .eq. "pbc") .or. (trim(input%tpbc) .eq. "bvk")) then
  call matdist_asy(input,unit_cell, &
                          distances)
else if (trim(input%tpbc) .eq. "asy") then
  call matdist_asy(input,asyms,&
                          distances)
else
  call matdist_atom(input,asyms%aunit(1),&
                          distances)
end if

call cpu_time(time2)

print *, 'time to compute distances', time2-time1


end subroutine

subroutine resres_distances(input,asyms,unit_cell,distances)
use mod_math, only: distance
! compute matrix of min distances between residues
!
type (inp_par), intent(in out)              :: input 
type(asym_List_Type)       ,intent(in)      :: asyms, unit_cell
real(dp), allocatable                       :: distances(:,:)
integer :: i,ii,j,jj,nres,ialloc
real(dp) :: xyzi(3), xyzj(3), dxyz(3),dist

nres = nresidues(asyms%aunit(1))

allocate(distances(nres,nres),stat=ialloc)

distances = 9999.99d0 

do i = 1, nres
  do ii = 1,asyms%aunit(1)%natoms 
    if(asyms%aunit(1)%atom(ii)%ires .eq. i) then
      xyzi = asyms%aunit(1)%atom(ii)%x
      do j = i+1, nres
        do jj = 1,asyms%aunit(1)%natoms 
          if(asyms%aunit(1)%atom(jj)%ires .eq. j) then
            xyzj = asyms%aunit(1)%atom(jj)%x
            dxyz = xyzj-xyzi
            dist=distance(dxyz)
            if (dist .le. distances(i,j)) then
              distances(i,j) = dist
              distances(j,i) = dist
            end if
          end if
        end do
      end do
    end if
  end do
end do
  
end subroutine

subroutine resres_kirch(input,asyms,unit_cell,kirchoff)
use mod_math, only: distance
! compute connectivity for a given cutoff between residues
!  ie. K(i,j) = -4.0 mean there are 4 connections within cutoff
!
type (inp_par), intent(in out)              :: input 
type(asym_List_Type)       ,intent(in)      :: asyms, unit_cell
real(dp), allocatable                       :: kirchoff(:,:)
integer :: i,ii,j,jj,nres,ialloc
real(dp) :: xyzi(3), xyzj(3), dxyz(3),dist

nres = nresidues(asyms%aunit(1))

allocate(kirchoff(nres,nres),stat=ialloc)

kirchoff = zero 

do i = 1, nres
  do ii = 1,asyms%aunit(1)%natoms 
    if(asyms%aunit(1)%atom(ii)%ires .eq. i) then
      xyzi = asyms%aunit(1)%atom(ii)%x
      do j = i+1, nres
        do jj = 1,asyms%aunit(1)%natoms 
          if(asyms%aunit(1)%atom(jj)%ires .eq. j) then
            xyzj = asyms%aunit(1)%atom(jj)%x
            dxyz = xyzj-xyzi
            dist=distance(dxyz)
            if (ii .ne. jj) then  !shouldn't come up
            if (dist .le. input%rcut_start) then
              kirchoff(i,j) = kirchoff(i,j)-one
              kirchoff(j,i) = kirchoff(i,j)
              kirchoff(i,i) = kirchoff(i,i) + one
              kirchoff(j,j) = kirchoff(j,j) + one
            end if
            end if
          end if
        end do
      end do
    end if
  end do
end do
  
end subroutine

integer function nresidues(atoms) 
type(atom_List_Type)       ,intent(in)      :: atoms
integer :: i,j,ires

ires = atoms%atom(1)%ires

j = 1
nresidues = 1
do i = 1, atoms%natoms
  if (ires .ne. atoms%atom(i)%ires) then
    ires = atoms%atom(i)%ires
    nresidues = nresidues+1
  end if
end do

end function

subroutine kirchoff_un (input,molecule, &
                              kirchoff,distances,tab_dxyzs )
use mod_math, only: distance
! inout
type(inp_par)       ,intent(in)   :: input
type(atom_List_Type) ,intent(in)  :: molecule
real(dp),allocatable ,intent(out) :: kirchoff(:,:),distances(:,:),tab_dxyzs(:,:,:)
! internal
real(dp)                          :: dist,dxyz(3),xyzi(3),xyzj(3)
integer                           :: ii,jj,natom,ialloc
real(sp)                          :: time1,time2

call cpu_time(time1)
natom     =molecule%natoms

call matrx2d_init(natom,kirchoff)
call matrx2d_init(natom,distances)
allocate(tab_dxyzs(natom,natom,3), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

do ii=1,natom

  xyzi = molecule%atom(ii)%x

  do jj=ii,natom

      xyzj = molecule%atom(jj)%x 
      dxyz = xyzj-xyzi
      dist = distance(dxyz)

      call kirch_comm(input,ii,jj,dist,dxyz, &
                               kirchoff,distances,tab_dxyzs)
  end do
end do

end subroutine

subroutine matdist_atom (input,molecule, &
                              distances)
use mod_math, only: distance
! inout
type(inp_par)       ,intent(in)   :: input
type(atom_List_Type) ,intent(in)  :: molecule
real(dp),allocatable ,intent(out) :: distances(:,:)
! internal
real(dp)                          :: dist,dxyz(3),xyzi(3),xyzj(3)
integer                           :: ii,jj,natom,ialloc
real(sp)                          :: time1,time2

natom     =molecule%natoms

call matrx2d_init(natom,distances)
distances = 99999.0d0

do ii=1,natom

  xyzi = molecule%atom(ii)%x

  do jj=ii,natom

      xyzj = molecule%atom(jj)%x 
      dxyz = xyzj-xyzi
      dist = distance(dxyz)

      if (dist .le. distances(ii,jj)) then
        distances(ii,jj) = dist
        distances(jj,ii) = dist
      end if
  end do
end do

end subroutine

subroutine matdist_asy (input, asy_uc, &
                            distances)
use mod_math, only: distance

! inout
type(inp_par)        ,intent(in)   :: input
type(asym_List_Type),intent(in)    :: asy_uc
real(dp),allocatable ,intent(out)  :: distances(:,:)
! internal
real(dp)                          :: dist, xyzi(3), xyzj(3),dxyz(3) 
integer                           :: ii,jj,k,j,h,natoms,nneighs,naunits,ialloc
real(sp)                          :: time1, time2
integer :: ja,jb,jc
real(dp),dimension(3) :: avec,bvec,cvec,vec

! mixing crystals here because we have all the neighbors stored in the p1 crystal

nneighs = asy_uc%nneighs
natoms  = asy_uc%aunit(1)%natoms
naunits = asy_uc%naunits

call matrx2d_init(natoms,distances)

! this subroutine computes the shortest distance between atoms under periodic boundary conditions
distances = 99999.0d0


do ii=1,natoms

  xyzi = asy_uc%aunit(input%aunit)%atom(ii)%x
  do j =1, naunits
    do jj=ii,natoms
  ! DMR do for relavent lattice displacements 
      do ja=asy_uc%aneigh(1),asy_uc%aneigh(2)
      avec = real(ja,kind=dp)*asy_uc%avec
      do jb=asy_uc%bneigh(1),asy_uc%bneigh(2)
      bvec = real(jb,kind=dp)*asy_uc%bvec
      do jc=asy_uc%cneigh(1),asy_uc%cneigh(2)
      cvec = real(jc,kind=dp)*asy_uc%cvec
      vec = avec + bvec + cvec

      xyzj = asy_uc%aunit(j)%atom(jj)%x+vec
      dxyz = xyzj-xyzi
      dist = distance(dxyz)

      if (dist .le. distances(ii,jj)) then
        distances(ii,jj) = dist
        distances(jj,ii) = dist
      end if

      end do
      end do
      end do
    end do
  end do
end do

end subroutine

subroutine kirchoff_asy (input, asy_uc, &
                            kirchoff,distances,tab_dxyzs )
use mod_math, only: distance

! inout
type(inp_par)        ,intent(in)   :: input
type(asym_List_Type),intent(in)    :: asy_uc
real(dp),allocatable ,intent(out)  :: kirchoff(:,:),distances(:,:),tab_dxyzs(:,:,:)
! internal
real(dp)                          :: dist, xyzi(3), xyzj(3), dxyz(3) 
integer                           :: ii,jj,k,j,h,natoms,nneighs,naunits,ialloc
real(sp)                          :: time1, time2
integer :: ja,jb,jc
real(dp),dimension(3) :: avec,bvec,cvec,vec

! mixing crystals here because we have all the neighbors stored in the p1 crystal

call cpu_time(time1)
nneighs = asy_uc%nneighs
natoms  = asy_uc%aunit(1)%natoms
naunits = asy_uc%naunits

call matrx2d_init(natoms,kirchoff)
call matrx2d_init(natoms,distances)
allocate(tab_dxyzs(natoms,natoms,3), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

do ii=1,natoms

! input%aunit for something dumb (vcov, I think I have this solved a different way)... 
!     fix dependency at some point.
  xyzi = asy_uc%aunit(input%aunit)%atom(ii)%x
  do j =1, naunits
    do jj=ii,natoms
  ! DMR do for relavent lattice displacements 
      do ja=asy_uc%aneigh(1),asy_uc%aneigh(2)
      avec = real(ja,kind=dp)*asy_uc%avec
      do jb=asy_uc%bneigh(1),asy_uc%bneigh(2)
      bvec = real(jb,kind=dp)*asy_uc%bvec
      do jc=asy_uc%cneigh(1),asy_uc%cneigh(2)
      cvec = real(jc,kind=dp)*asy_uc%cvec
      vec = avec + bvec + cvec

      xyzj = asy_uc%aunit(j)%atom(jj)%x+vec
      dxyz = xyzj-xyzi
      dist = distance(dxyz)

      call kirch_comm(input,ii,jj,dist,dxyz, &
                              kirchoff,distances,tab_dxyzs)
      end do
      end do
      end do
    end do
  end do
end do


call cpu_time(time2)
print *, "Time setting kirchoff_crys matrix:", time2-time1

end subroutine

subroutine kirchoff_crys (input, unit_cell, &
                            kirchoff,distances,tab_dxyzs )
use mod_math, only: distance

! inout
type(inp_par)        ,intent(in)   :: input
type(asym_list_Type),intent(in)    :: unit_cell
real(dp),allocatable ,intent(out)  :: kirchoff(:,:),distances(:,:),tab_dxyzs(:,:,:)
! internal
real(dp)                          :: dist, xyzi(3), xyzj(3), dxyz(3) 
integer                           :: ii,jj,k,h,natoms,nneighs,naunits,ialloc
real(sp)                          :: time1, time2
integer :: ja,jb,jc,nati,nuni,natj,nunj
real(dp),dimension(3) :: avec,bvec,cvec,vec

! mixing crystals here because we have all the neighbors stored in the p1 crystal

call cpu_time(time1)

natoms  = unit_cell%natoms
call matrx2d_init(natoms,kirchoff)
call matrx2d_init(natoms,distances)
allocate(tab_dxyzs(natoms,natoms,3), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

do ii=1,natoms

  xyzi = unit_cell%aunit(1)%atom(ii)%x

  do jj=ii,natoms

    ! DMR do for relavent lattice displacements 
    do ja=unit_cell%aneigh(1),unit_cell%aneigh(2)
      avec = real(ja,kind=dp)*unit_cell%avec
      do jb=unit_cell%bneigh(1),unit_cell%bneigh(2)
        bvec = real(jb,kind=dp)*unit_cell%bvec
        do jc=unit_cell%cneigh(1),unit_cell%cneigh(2)
          cvec = real(jc,kind=dp)*unit_cell%cvec
          vec = avec + bvec + cvec

          xyzj = unit_cell%aunit(1)%atom(jj)%x+vec
          dxyz = xyzj-xyzi
          dist = distance(dxyz)
         ! print *, 'shit',ii,jj,dist
          call kirch_comm(input,ii,jj,dist,dxyz, &
                             kirchoff,distances,tab_dxyzs)
        end do
      end do
    end do
  end do
end do


call cpu_time(time2)
print *, "Time setting in kirchoff_crys matrix:", time2-time1

end subroutine

subroutine kirch_comm(input, ii,jj,dist,dxyz,kirchoff,distances,tab_dxyzs)
type(inp_par),        intent(in)            :: input
integer ,             intent(in)            :: ii,jj
real(dp),             intent(in)            :: dist,dxyz(3)
real(dp),allocatable, intent(in out)        :: kirchoff(:,:),distances(:,:),tab_dxyzs(:,:,:)


if (dist .le. input%rcut_start) then
  if (dist .gt. zero) then

    kirchoff(ii,jj) = kirchoff(ii,jj) - one
! self term
    kirchoff(ii,ii) = kirchoff(ii,ii) + one
    if (ii.ne.jj) kirchoff(jj,jj) = kirchoff(jj,jj) + one

! dealing with multiple interactions...  
! if(multint .eq. .true.) then we will allow -two kirchoff elements
    if (kirchoff(ii,jj) .lt. -one) then
!  we must have an entry for the distances: keep the shortest distance
      if (dist .lt. distances(ii,jj)) then
        tab_dxyzs(ii,jj,:) = dxyz
        distances(ii,jj)   = dist 
      end if
    
      if (input%multint .eqv. .false.) then    
      !if (input%multint) then
            !do nothing...  no unless in fortran??
      !else
        kirchoff(ii,jj) = kirchoff(ii,jj) + one
        kirchoff(ii,ii) = kirchoff(ii,ii) - one
        if(ii.ne.jj) kirchoff(jj,jj) = kirchoff(jj,jj) - one
      end if

    else
      ! kirchoff(ii,jj) = -one  so we just set the dist info
      tab_dxyzs(ii,jj,:) = dxyz
      distances(ii,jj)   = dist
    end if

  end if
end if

end subroutine

subroutine chkwrt_kirchoff(input,kirchoff)
type(inp_par), intent(in)        :: input
type(sparse) ,intent(in out)         :: kirchoff
integer :: ii,jj,ier,natom
real(dp) :: multi

natom=kirchoff%ndim

! write out nonzero entries to kirchoff.txt
if(input%print_level .ge. 3) then
  open(unit=13,file=trim(input%fileroot)//"-kirchoff.txt", &
     status='unknown',action='write', &
     iostat=ier)

  write(13,'(I6,F10.2)')natom,input%rcut_start
end if

multi=zero
do ii = 1, kirchoff%nnzero
  if(input%print_level .ge. 3) then
     write(13,'(2I5,F15.5)') kirchoff%rows(ii), &
                         kirchoff%columns(ii),kirchoff%values(ii)
!     if (kirchoff%rows(ii) .ne. kirchoff%columns(ii)) write(13,'(2I5,F15.5)')  &
!                         kirchoff%columns(ii), kirchoff%rows(ii),kirchoff%values(ii)
  end if
  if ( kirchoff%values(ii) .lt. -one) then
    multi=multi-(kirchoff%values(ii)+one)
    if (input%print_level .ge. 3) then
     print *, 'atoms ',kirchoff%rows(ii),kirchoff%columns(ii),' interact multiple times!'
    end if
  end if
end do

kirchoff%nmulti=nint(multi)
kirchoff%nbonds=kirchoff%nnzero - kirchoff%ndim

close(13)

end subroutine

subroutine distances_kirchoff(input,distances,kirchoff)
! compute kirchoff matrix based on distance matrix
type(inp_par),         intent(in) :: input
real(dp), allocatable, intent(in) :: distances(:,:)
real(dp), allocatable, intent(out) :: kirchoff(:,:)
integer :: i,j,ndim

ndim = size(distances,1)

allocate(kirchoff(ndim,ndim),stat=i)

kirchoff = zero

do i = 1,ndim
  do j = i+1,ndim
    if(distances(i,j) .le. input%rcut_start) then
      kirchoff(i,j) = kirchoff(i,j)-one
      kirchoff(j,i) = kirchoff(j,i)-one
      kirchoff(i,i) = kirchoff(i,i)+one
      kirchoff(j,j) = kirchoff(j,j)+one
    end if
  end do
end do


end subroutine


end module

