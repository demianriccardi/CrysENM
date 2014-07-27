program chunkwrite
!
! crysFML modules
use cfml_crystallographic_symmetry, only: space_group_type,Write_SpaceGroup
use cfml_Atom_typedef,              only: Atom_List_Type
use cfml_crystal_metrics,           only: Crystal_Cell_Type
use cfml_IO_Formats,                only: Readn_set_Xtal_Structure, file_list_type
use cfml_globaldeps,                only: sp,dp
use mod_crysbuild,                  only: fract_to_cart, cart_to_fract,crys_build_new,p1asym_to_chunk

! external modules
use mod_inout,                only: pdbwriter_new,fileroot,fileroot_set,iso_bc
use mod_constants,            only: one,two,zero,pi
use mod_types,                only: inp_par,asym_list_type
use mod_bnm,                  only: read_blocks_atom
use mod_lattdyn,              only: qvec_fqvec

!use mod_correl_dist
implicit none
!#ifdef _OPENMP 
!   include 'omp_lib.h'  !needed for OMP_GET_NUM_THREADS()
!#endif

type (file_list_type)       :: fich_cfl
type (space_group_type)     :: SpG
type (Crystal_Cell_Type)    :: Cell
type (Atom_list_Type)       :: asym_un 
type(asym_list_type)        :: asyms, unit_cell ,chunk
character(len=25)           :: filcfl,filparam,crap,pdbname
integer                     :: narg, ier,ia,ib,ic
logical                     :: arggiven
real(sp)                    :: time1,time2
integer                     :: nviab
type(inp_par)            :: input

call cpu_time(time1)
!---- Arguments on the command line ----!
narg=iargc()

!      print *, narg, 'args'
if(narg > 0) then
  call getarg(1,filcfl)
  input%fileroot = filcfl
  arggiven=.true.
  if (trim(filcfl) .eq. "help") then
    print *, "to run: ./main 1ruv input"
    print *, "this reads in 1ruv.cfl and the parameters in input"
    print *, ".cfl files are for reading in the asym unit in a way that crysFML likes"
    print *, "input: file with 3 entries on same line... for now"
    print *, "first: rcut"
    print *, "second: gnm,enm,cdyn,pdbwrite,or intensity"
    print *, "ie: if (trim(runtype) .eq. 'gnm')        call gnm"
    print *, "third: 'pbc' if periodic boundaries are desired"
    stop
  end if
end if
if(.not. arggiven) then
  print *, "no arguments given"
  stop
end if

call fileroot_set(filcfl)

if(narg > 1) then
  call getarg(2,filparam)
  open(unit=11,file=filparam, &
     status='old',action='read', position='rewind', &
     iostat=ier)
  if (ier > 0) stop "cannot open param file"
  read(unit=11,fmt=*,iostat=ier) crap, ia,ib,ic
  read(unit=11,fmt=*,iostat=ier) crap, pdbname
  read(unit=11,fmt=*,iostat=ier) crap, input%coord

 
end if

!input%rcut_start = 20.0d0
input%tpbc = "pbc" 
input%genm = "enm"
input%weigh = "mass"
input%fctyp  = "anm"   
input%fcpower=0 
input%fcnst = 1.0d0
input%rcut_end = input%rcut_start
input%eigtyp = "lapack"
input%perc   = 0.05d0
input%first_mode = 0
input%last_mode  = 999999
input%bztyp = "full"
input%qdiv  = 3
input%qmax  = 0.5d0
input%isocorr_gamm = 10.0d0
input%atom_anl        = "all"
input%multint         = .true.
input%runtype         = "isointensity"
input%drcut           = one
input%qdir            = zero 
input%fqvec           = zero
input%hihfkikflilf    = 0
input%aniso           = .true.
input%dosgo           = .true.
input%print_level     = 2 
input%isocorr_cut     = input%rcut_start
input%isocorr_c0      = one
input%animate_mag     = one
input%animate_nframes = 0
input%block_read      = .false.
input%bnm             = 1
input%temperature     = 298.0d0
input%isovcov_scale   = "lsq"

if (iso_bc(input)) then
  input%multint = .false.
else 
  input%multint = .true.
  ! DMR: there is a problem with the multint flag when pbc is used
  !      since i generally don't use it for pbc (min image)
  !      I'm not fixing it for now.
end if

! read in the asymetric unit
inquire(file=trim(filcfl)//".cfl",exist=arggiven)
if(arggiven) then
  call Readn_set_Xtal_Structure(trim(filcfl)//".cfl", &
       Cell,SpG,asym_un,Mode="CFL",file_list=fich_cfl)
!  call Write_SpaceGroup(spg,full = .true.)
else
  print *, trim(filcfl)//".cfl ",'file not found'
end if

if (trim(input%coord) .eq. "fract") then
print *, "shit",SpG%multip
  call fract_to_cart(cell,asym_un)
end if

! is this really optional?
if (input%bnm .eq. 1) then
  call read_blocks_atom(input,asym_un)
end if


input%multiplicity=SpG%multip
call qvec_fqvec(input,cell)
call crys_build_new(input,cell,spg,asym_un,asyms,unit_cell)
call p1asym_to_chunk(ia,ib,ic,unit_cell,chunk)
call pdbwriter_new(chunk,trim(pdbname))

end

