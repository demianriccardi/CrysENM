program crysdyn
!
! crysFML modules
use cfml_crystallographic_symmetry,only: space_group_type,Write_SpaceGroup
use cfml_Atom_typedef,              only: Atom_List_Type
use cfml_crystal_metrics,           only: Crystal_Cell_Type
use cfml_IO_Formats,               only: Readn_set_Xtal_Structure,
file_list_type
use cfml_globaldeps,             only: sp,dp
use mod_crysbuild,            only: fract_to_cart, cart_to_fract,crys_build_new

! external modules
use mod_inout,                only: fileroot,fileroot_set,iso_bc
use mod_constants,            only: one,two,zero,pi,b2u
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
type(asym_list_type)        :: asyms, unit_cell 
character(len=25)           :: filcfl,filparam,crap
integer                     :: narg, ier
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
    print *, ".cfl files are for reading in the asym unit in a way that crysFML
likes"
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
  read(unit=11,fmt=*,iostat=ier) crap,input%atom_anl
  read(unit=11,fmt=*,iostat=ier) crap,input%tpbc,crap,input%genm, &
                                 crap,input%weigh, crap, input%fctyp     
  read(unit=11,fmt=*,iostat=ier) crap,input%fcnst, crap,input%fcpower 
  read(unit=11,fmt=*,iostat=ier) crap, input%rcut_start
  read(unit=11,fmt=*,iostat=ier) input%eigtyp, crap, input%perc
  read(unit=11,fmt=*,iostat=ier) crap, input%first_mode,input%last_mode
  read(unit=11,fmt=*,iostat=ier)
crap,input%bztyp,crap,input%qdiv,crap,input%qmax
  read(unit=11,fmt=*,iostat=ier) crap, input%temperature 
  read(unit=11,fmt=*,iostat=ier) crap, input%dosgo,crap, input%aniso 
  
end if

input%isocorr_gamm    = 10.0d0 
input%rcut_end        = input%rcut_start
input%multint         = .true.
input%runtype         = "isointensity"
input%drcut           = one
input%qdir            = zero 
input%fqvec           = zero
input%hihfkikflilf    = 0
input%print_level     = 2 
input%isocorr_cut     = input%rcut_start
input%isocorr_c0      = one
input%animate_mag     = one
input%animate_nframes = 0
input%coord           = "cart"
input%block_read      = .false.
input%bnm             = 1
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
  call fract_to_cart(cell,asym_un)
end if

! is this really optional?
if (input%bnm .eq. 1) then
  call read_blocks_atom(input,asym_un)
end if


input%multiplicity=SpG%multip
call qvec_fqvec(input,cell)

call crys_build_new(input,cell,spg,asym_un,asyms,unit_cell)
call temperature_factor(input,spg,cell,asyms,unit_cell)

call cpu_time(time2)

write(6,'(a7,1x,a11,1x,I5,1x,a16,1x,I5)') input%weigh, ':: nmodes =', &
                                          input%last_mode-input%first_mode +1,
':: nocc1-bfact =',input%nviab 
write(6,'(a27,1x,I5,1x,A7,1x,I5)') 'vcov constructed from mode:', &
                                          input%first_mode, 'to mode',
input%last_mode
write(6,'(a4,1x,a9,1x,F10.2,1x,a7)') input%fileroot, 'CPU time:', time2-time1,
'seconds'
print *, " "

contains

subroutine temperature_factor(input,spg,cell,asyms,unit_cell)
use mod_types,                only: dos,inp_par,asym_list_type,dos_histgen
!use mod_types,                only:
dos,dos_analysis,inp_par,asym_list_type,dos_histgen
use mod_vcov_store
use mod_bnm,      only: block,read_blocks_atom,blocks_setup
use mod_crysbuild,only: block_assign,asym_shrinker
use mod_math, only: analyt_lsq
use mod_tempfact, only: bfactor_scale,get_expiso
use mod_inout,    only: pdbwriter_new
use mod_kirchoff, only: mat_distances
use mod_inout,    only: isovcov_atshrink
use mod_intensity,            only: rval,reflections, atscatter,  &
                                    reflections_atscatt_setres,&
                                    intensity_zero_iso
use mod_symop, only: anisotempfactor_p1,tempfactor_p1,atcor_p1

type (inp_par)   ,        intent(in out)    :: input
type (space_group_type), intent(in)         :: spg
type (Crystal_Cell_Type), intent(in)        :: Cell
type (asym_list_type),    intent(in out)    :: asyms,unit_cell
type (asym_list_type)                       :: subasy,subunt
type(vcov_store)        :: svcov
type(dos)               :: states
real(dp)                :: correl,fcnst,b,dr
real(dp),   allocatable ::
distances(:,:),atcorr(:,:),subisov(:,:),isovcov(:,:),subiso(:)
type(block),allocatable :: blocks(:)
real(dp),   allocatable :: bfexp(:),bfthr(:)
integer,    allocatable :: resid(:)
type(reflections)                     :: reflect
type(atscatter)                       :: atomic_scatt
real(dp) :: a,rval_iso
integer :: ialloc,i
real(dp), allocatable    :: intns_eiso(:),intns_tiso(:)

dr = 1.0d0

!first we'll calculate the isovcov
 
call vcov_calc(input,cell,asyms%aunit(1),asyms,unit_cell,svcov,states,blocks)
print *, 'done! with vcov_calc'

call asym_shrinker(input, asyms,subasy)
call asym_shrinker(input, unit_cell,subunt)

if(trim(input%tpbc) .eq. "pbc" .or. trim(input%tpbc) .eq. "bvk") then
  call isovcov_atshrink(input,unit_cell%aunit(1),svcov%isovcov,subisov)
else
  call isovcov_atshrink(input,asyms%aunit(1),svcov%isovcov,subisov)
end if


if (trim(input%genm) .eq. "tls" .and. (trim(input%tpbc) .eq. "pbc" .or.
trim(input%tpbc) .eq. "bvk" )) then
  call atcor_p1(spg,subisov,isovcov)
  call atcorr_calc(isovcov,atcorr)
else
  call atcorr_calc(subisov,atcorr)
end if

call mat_distances(input,subasy,subunt,distances)


allocate(subiso(size(subisov,1)),stat = ialloc)

do i = 1, size(subiso)
  subiso(i) = subisov(i,i)
 !  print *, 'shit',subasy%aunit(1)%atom(i)%biso*b2u, subiso(i)
end do

call bfactor_scale(input,subasy%aunit(1),subiso,&
                            bfexp,bfthr,resid,correl,fcnst,b,nviab)

write(6,'(A8,1x,F7.3,1x,A8,1x,E12.5)') &
         'isocor =',correl, ':: fcnst_fit =', one/fcnst

open(unit=22,file=trim(input%fileroot)//"-"//trim(input%genm)//"-"// &
                  trim(input%tpbc)//"-"//trim(input%fctyp) //"-"//   &
                  trim(input%bztyp)//"-isoaniso.txt", &
      status='unknown',action='write', iostat=ier)

 do i=1,size(resid)
    write (22,'(I5,2F12.5)') resid(i),bfexp(i),bfthr(i)
  end do

call list_block_anal(input,dr,atcorr,distances,subasy,subunt)

if (input%genm .ne. "tls") then
  states%nbins = 200
  call dos_histgen(states)
  if (input%dosgo) call dos_analysis(input,states)
end if

open(unit=22,file=trim(input%fileroot)//"-"//trim(input%genm)//"-"// &
                  trim(input%tpbc)//"-"//trim(input%fctyp)//"-freqs.txt", &
      status='unknown',action='write', iostat=ier)


do i = 1, size(states%freqs)
  write(22,'(F10.6)') states%freqs(i)
end do

close(22)

end subroutine

subroutine list_block_anal(input,dr,atcorr,distances,asyms,unit_cell) 
use mod_types,                only: inp_par,asym_list_type
use mod_bnm,                  only:
block,read_blocks_atom,blocks_setup,isocorr_byblock
use mod_crysbuild,            only: block_assign
type (inp_par)   ,        intent(in out)    :: input
real(dp),                intent(in)     :: dr
real(dp),allocatable   , intent(in)     :: atcorr(:,:),distances(:,:)
type (asym_list_type),   intent(in out) :: asyms,unit_cell
type(block),allocatable                 :: blocks(:)
character(len=25)                        :: blockid(100)
integer  :: i,j,nbin,ier,inputstatus,nanalys

open(unit=11,file="blocks.txt", &
   status='old',action='read', position='rewind', &
   iostat=ier)
if (ier>0) stop "blocks.txt input error"

i = 0
do
  i = i +1
  read (11,*,iostat = inputstatus) blockid(i)
  if (inputstatus>0) stop "blocks.txt input error"
  if (inputstatus<0) exit
end do
close(11)

nanalys = i - 1

do i = 1,nanalys
  if (allocated(blocks)) deallocate(blocks)
  call read_blocks_atom(input,asyms%aunit(1),blockid(i))
  call block_assign(asyms%aunit(1), asyms,input%nblocks) 
  input%nblocks = asyms%nblocks
  call block_assign(asyms,unit_cell) 

  if(trim(input%tpbc) .eq. "pbc" .or. trim(input%tpbc) .eq. "bvk") then
    call blocks_setup (input,unit_cell%aunit(1),blocks)
  else
    call blocks_setup (input,asyms%aunit(1),blocks)
  end if

  call isocorr_byblock(input,dr,atcorr,distances,blocks,blockid(i))

end do

end subroutine

subroutine atcorr_calc(isovcov,atcorr)
real(dp), allocatable, intent(in) :: isovcov(:,:)
real(dp), allocatable, intent(out) :: atcorr(:,:)
integer :: i,j,ier

if(allocated(atcorr)) deallocate(atcorr) 
allocate(atcorr(size(isovcov,1),size(isovcov,2)),stat=ier)

do i = 1, size(isovcov,1)
  do j = 1, size(isovcov,2)
    atcorr(i,j) = isovcov(i,j)/dsqrt(isovcov(i,i)*isovcov(j,j))
  end do
end do

end subroutine

subroutine dos_analysis(input,states)
!DMR add in some curve fitting stuff
use mod_curvefit
use mod_constants , only: rkcal, hkcal, convfrq,clight
use mod_types     , only: inp_par, dos
integer, parameter :: max_temp=400
type(inp_par), intent(in) :: input
type(dos),     intent(in) :: states
real(dp), allocatable, dimension(:) :: temps, cv
integer :: i,jtemp,ier
real(dp) :: hfreq,kt,bestpart,hvkt,maxperc,firstperc,deltaperc
type(datafunc) :: adfunc,gdfunc


allocate(temps(max_temp),cv(max_temp), stat=ier)

if (input%print_level .ge. 2) then
      open(unit=22,file=trim(input%fileroot)//"-"//trim(input%genm)//"-"//trim(input%tpbc)//"-"//
&
                        trim(input%fctyp) //"-"//trim(input%bztyp)//"-dos.txt",
&
      status='unknown',action='write', iostat=ier)

      open(unit=23,file=trim(input%fileroot)//"-"//trim(input%genm)//"-"//trim(input%tpbc)//"-"//
&
                        trim(input%fctyp) //"-"//trim(input%bztyp)//"-cv.txt", &
      status='unknown',action='write', iostat=ier)
end if

! write out histogram
if (input%print_level .ge. 2) then
  do i=1,states%nbins
    write (22,'(I4,I7,F10.5,F10.5,F10.5,F15.5)') i,states%histo(i),  &
         states%coor(i),states%gw(i),states%biGw(i) ,states%fgw(i)
  end do
end if

! calculate heat capacity
! something is just not right
    cv=zero
  do jtemp = 1, 400
    temps(jtemp) = real(jtemp,kind=dp)
! start from 2 to avoid the zero
    do i=2,states%nbins !DMR
      hfreq    = hkcal * states%coor(i) * 1.0D12
      !hfreq    = hkcal * real(i)*states%delta_frq * 1.0D12
      kt       = real(jtemp)*rkcal
      bestpart = one/(two*cosh(hfreq/kt) - two )
      cv(jtemp)    = cv(jtemp) + &
                     (one/(kt*real(jtemp)))*states%gw(i)*(hfreq**2)*bestpart*states%delta_frq
    end do
    if (cv(jtemp) .le. 1.0D-16) cv(jtemp)=zero
if (input%print_level .ge. 2) then
    write(23,'(F10.2,E15.5)') real(jtemp), Cv(jtemp)
end if

  end do

close(22)
close(23)

print '(A12)', "BEGINFIT CV>"
print '(A3,2A6,A10,A6,2A10,A6,A10,A6,A10,A6)',
"Cv>","frst%","last%","a","b","msd", &
                             "a_log","b_log","msd_log", "xcent","yshft", "#pnts"

maxperc = 0.1d0
firstperc = 1.0d-02
deltaperc = 5.0d-02

call data_adjust(temps,cv,rkcal,firstperc,maxperc,adfunc)
call minsearch(adfunc)
call loglogfit(adfunc)

print '(A3,2F6.2,E10.3,F6.2,2E10.3,F6.2,E10.3,F6.2,E10.3,I6)', "Cv>",
firstperc,maxperc, &
                                        adfunc%a, adfunc%b,adfunc%msd, &
                                        adfunc%a_log,adfunc%b_log,
adfunc%msd_log, &
                                        adfunc%acoor_offset,adfunc%func_offset,
&
                                        size(adfunc%func)


print '(A10)', "ENDFIT CV>"
print *, " "
print '(A12)', "BEGINFIT GW>"
print '(A3,2A6,A10,A6,2A10,A6,A10,A6,A10,A6)',
"Gw>","frst%","last%","a","b","msd", &
                             "a_log","b_log","msd_log", "xcent","yshft", "#pnts"

maxperc = 0.2d0
firstperc = 3.0d-02

call data_adjust(states%coor,states%biGw,one,firstperc,maxperc,gdfunc)
call minsearch(gdfunc)
call loglogfit(gdfunc)

print '(A3,2F6.2,E10.3,F6.2,2E10.3,F6.2,E10.3,F6.2,E10.3,I6)',
"Gw>",firstperc,maxperc, &
                                          gdfunc%a, gdfunc%b,gdfunc%msd, &
                                          gdfunc%a_log,gdfunc%b_log,
gdfunc%msd_log, &
                                          gdfunc%acoor_offset,gdfunc%func_offset,
&
                                          size(gdfunc%func)


print '(A10)', "ENDFIT GW>"


end subroutine

subroutine input_printer(input)
type(inp_par), intent(in) :: input

print '(A17,1x,A15)', "INPUT> atom_anl :", input%atom_anl
!print '(A15,1x,A15)', "INPUT> isovcov_scale", input%isovcov_scale
print '(A17,1x,A10)', "INPUT> tpbc     :", input%tpbc
print '(A17,1x,A11)', "INPUT> gebnm    :",input%genm
print '(A17,1x,A11)', "INPUT> weigh    :",input%weigh
print '(A17,1x,A11)', "INPUT> fctyp    :",input%fctyp
print '(A17,1x,E12.4)',  "INPUT> fcnst    :",input%fcnst
print '(A17,1x,I1)',     "INPUT> fcpower  :",input%fcpower
print '(A17,1x,F6.2)',   "INPUT> rcuti    :",input%rcut_start
!print '(A15,1x,A15)', "INPUT> rcutf",input%rcut_end
!print '(A15,1x,A15)', "INPUT> hkl rng",input%hihfkikflilf
print '(A17,1x,A15)', "INPUT> eigtyp   :",input%eigtyp
print '(A17,1x,F6.2)',   "INPUT> perc     :",input%perc
print '(A17,1x,I3)',     "INPUT> modei    :",input%first_mode
print '(A17,1x,I5)',     "INPUT> modef    :",input%last_mode
!print '(A15,1x,A15)', "INPUT> bztyp",input%bztyp
!print '(A15,1x,A15)',"INPUT> qdiv",input%qdiv
!print '(A15,1x,A15)',"INPUT> qmax",input%qmax
print '(A17,1x,F8.2)',   "INPUT> temp     :", input%temperature
!print '(A15,1x,A15)',"INPUT> gamma", input%isocorr_gamm
!print '(A15,1x,A15)',"INPUT> multint",input%multint
!print '(A15,1x,A15)',"INPUT> runtype",input%runtype
!print '(A15,1x,A15)',"INPUT> drcut",input%drcut
!print '(A15,1x,A15)',"INPUT> qdir",input%qdir
!print '(A15,1x,A15)',"INPUT> fqvec",input%fqvec
!print '(A15,1x,A15)',"INPUT> aniso",input%aniso
!print '(A15,1x,A15)',"INPUT> dosgo",input%dosgo
print '(A17,1x,I2)',     "INPUT> printlev :",input%print_level
print '(A17,1x,A5)', "INPUT> coord    :",input%coord
!print '(A15,1x,A15)',"INPUT> isocorr_cut",input%isocorr_cut
!print '(A15,1x,A15)',"INPUT> isocorr_c0",input%isocorr_c0
!print '(A15,1x,A15)',"INPUT> animate_mag",input%animate_mag
!print '(A15,1x,A15)',"INPUT> animate_nfrm",input%animate_nframes
!print '(A15,1x,A15)',"INPUT> block_read",input%block_read
!print '(A15,1x,A15)',"INPUT> bnm",input%bnm

end subroutine

end
