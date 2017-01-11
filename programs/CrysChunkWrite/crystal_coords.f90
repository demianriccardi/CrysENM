program crystal_coords 
! simple pure crysfml crystal writer
!
! invokation:
! ./crystal_coords some.cfl 0 0 0
! ./crystal_coords some.cfl 3 3 3
! crysFML modules
use cfml_atom_typedef,              only: atom_list_type
use cfml_crystallographic_symmetry, only: space_group_type,Write_SpaceGroup, &
                                                             get_orbit,applyso
use cfml_crystal_metrics,           only: Crystal_Cell_Type
use cfml_IO_Formats,                only: Readn_set_Xtal_Structure, file_list_type
use cfml_globaldeps,                only: sp,dp

implicit none

type (file_list_type)       :: fich_cfl
type (space_group_type)     :: SpG
type (Crystal_Cell_Type)    :: Cell
type (atom_list_type)       :: asym_un 
character(len=25)           :: filcfl,filparam,crap,pdbname, cha,chb,chc
integer                     :: narg, ier,ia,ib,ic,stat,i,j,nmultip
real,     allocatable       :: x(:), orb(:,:) ! for get_orbit subroutine 
logical                     :: arggiven
real(sp)                    :: time1,time2
integer                     :: nviab

call cpu_time(time1)
!---- Arguments on the command line ----!
narg=iargc()

!      print *, narg, 'args'
if(narg > 0) then
  call getarg(1,filcfl)
  call getarg(2,cha)
  read(cha,*,iostat=stat) ia
  call getarg(3,chb)
  read(cha,*,iostat=stat) ib
  call getarg(4,chc)
  read(cha,*,iostat=stat) ic
  arggiven=.true.
end if
if(.not. arggiven) then
  print *, "no arguments given"
  stop
end if

! read in the asymetric unit
inquire(file=trim(filcfl)//".cfl",exist=arggiven)
!if(arggiven) then
call Readn_set_Xtal_Structure(trim(filcfl), &
       Cell,SpG,asym_un,Mode="CFL",file_list=fich_cfl)

!  call Write_SpaceGroup(spg,full = .true.)
!else
!  print *, trim(filcfl),'file not found'
!end if

allocate(orb(3,SpG%multip), stat=ier)
allocate(x(3), stat=ier)

! fill up unit cell with appropriate treatment of Wyckoff special positions
! 
do i=1, asym_un%natoms
    call Get_Orbit(asym_un%atom(i)%x(:),SpG,nmultip,orb)
    
    ! get_orbit fills orb with non lattice vector equiv position 
    do j=1,nmultip

        ! switch to cartesian 
        x(:) = matmul( orb(:,j), transpose(Cell%Cr_orth_cel))
        print 2001, asym_un%atom(i)%lab, x(:)

    end do
end do

 2001 format ( A5, 3F10.6 )
 2002 format ( A5, 3F10.6, I4 )

end

