module mod_atomvector_store
! contain the subroutine to make bfactor comparisons 
! this code needs to be cleaned up and commented
use cfml_globaldeps,                 only: dp,sp
use mod_constants

implicit none
type :: eigvec
! 3d disp eigvector for atom.. real or complex
complex(kind=8) :: zvec(3)
end type eigvec

type :: byatom_val_vec
! eigen vector stored by atom... eigenvalue also stored
real(dp)                 :: eigenval
type(eigvec),allocatable :: atdisp(:)
end type byatom_val_vec

type :: bybranch_vals_vecs
! all branches (ie. acoustic and optic) for 
integer              :: nfrqs,natoms
real(dp)             :: qvec(3)
type(byatom_val_vec),allocatable :: branch(:)
end type bybranch_vals_vecs

type :: byqvec_branches
! store by wavevector...  could get big.
integer :: natoms,nqvecs
type(bybranch_vals_vecs),allocatable :: wave_vec(:)
end type byqvec_branches

contains

subroutine branches_init(nbranches,natoms,branches)
integer,                      intent(in)    :: nbranches,natoms
type (bybranch_vals_vecs),    intent(out)   :: branches
integer :: i,j,ialloc

branches%natoms = natoms
branches%nfrqs  = nbranches

allocate(branches%branch(branches%nfrqs), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

do i = 1, nbranches

  call byatom_init(natoms,branches%branch(i))

end do

end subroutine

subroutine store_branches(input,vals,zvecs,qvec,branches)
use mod_types,  only: inp_par
type(inp_par),                intent(in)    :: input
real(dp), allocatable,        intent(in)    :: vals(:)
complex(kind=8), allocatable, intent(in)    :: zvecs(:,:)
real(dp),                     intent(in)    :: qvec(3)
type (bybranch_vals_vecs),    intent(out)   :: branches
type (byatom_val_vec)                       :: d_atoms
complex(kind=8), allocatable                :: tzvec(:)
complex(kind=8) :: zvec(3)
real(dp) :: val
integer :: i,j,k,ialloc, natoms,nbranches
if (input%nfrqs .ne. size(vals)) Stop "store_branches> stopping, wrong number of freqs?"
allocate(tzvec(size(zvecs(:,1))),stat = ialloc)
if (ialloc /= 0 ) stop "Not enough memory"

nbranches = input%nfrqs
natoms    = input%natoms

!initialize the branch vects
call branches_init(nbranches,natoms,branches)
branches%qvec = qvec

do i = 1, nbranches
  val = vals(i)
  tzvec(:) = zvecs(:,i)
  call valvecs_byatom(input,val,tzvec,branches%branch(i))
end do

end subroutine


subroutine valvecs_byatom(input,val,zvector,d_atoms)
use mod_types,                only: inp_par

type(inp_par) ,               intent(in)  :: input
real(dp),                     intent(in)  :: val
complex(kind=8), allocatable, intent(in)  :: zvector(:)
type (byatom_val_vec),        intent(out) :: d_atoms
complex(kind=8) :: azvec(3)
integer :: i,j,ialloc

call byatom_init(input%natoms,d_atoms)
d_atoms%eigenval = val

j = 1
do i = 1, input%natoms
  d_atoms%atdisp(i)%zvec(:) = zvector(j:j+2)
  j = j+3
end do

end subroutine

subroutine byatom_init(natoms,d_atoms)
integer,                      intent(in)  :: natoms
type (byatom_val_vec),        intent(out) :: d_atoms
integer :: ialloc

if(allocated(d_atoms%atdisp)) deallocate(d_atoms%atdisp)
allocate(d_atoms%atdisp(natoms), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

end subroutine

end module
