module mod_lattdyn
!
! DMR 04-06-2008
!
! module for carrying out crystal dynamics calculations
!
! 1. lay out grid for brillouin zone
! 2. construct dynamical matrix from periodic hessian
! 3. compute eigenvals and vecs for brillouin zone mesh
!
! important to keep factors of 2pi in check
!
!CrysFML modules
use cfml_globaldeps,                  only: dp,sp
use cfml_Atom_typedef,                only: Atom_List_Type
use cfml_crystal_metrics,             only: Crystal_Cell_Type
!external modules
use mod_types,                only: inp_par,sparse,uc_list_type,asym_list_type
use mod_math,                 only: dist_sqr,distance
use mod_constants,            only: one,two,zero,pi
use mod_arpack
use mod_inout

implicit none
integer :: print_level

! public subroutines: 
public ::  dynmatrx, qdispers,speedsnd,qdisp_split
! private subroutines 
private :: wrt_dynm

contains

subroutine dynmatrx (input,unit_cell,kirchoff,qvec, &
                             kirchdyn_coor, hessdyn_coor)
!inout
use mod_inout,     only: zfull_to_zsparse
use mod_hessian,   only: sparsehess_massweigh
type(inp_par)        , intent(in out)     :: input
type(asym_List_Type) , intent(in)         :: unit_cell 
type(sparse)         , intent(in)         :: kirchoff
real(dp)             , intent(in)         :: qvec(3)
type(sparse)   ,       intent(out)        :: kirchdyn_coor,hessdyn_coor
complex(kind=8),allocatable               :: kirchdyn(:,:),hessdyn(:,:)
type(sparse)                              :: tmpdyn
!internal
integer                            :: i,ii,j,jj,iii,jjj,kkk
integer                            :: natom,nneighs,degfred,ialloc, &
                                      chaini,chainj,iaunit,jaunit,ichain,jchain
real(dp)                           :: rcut,dist,distmp
real(dp), dimension(3)             :: avec,bvec,cvec,vec
real(dp), dimension(3)             :: xyzi,xyzj,dxyz
real(sp)                           :: time1, time2
real(dp) :: dottest

call cpu_time(time1)
natom   = unit_cell%natoms
degfred = 3*natom
rcut    = input%rcut_start

if (allocated(kirchdyn)) deallocate(kirchdyn)
if (allocated(hessdyn)) deallocate(hessdyn)
allocate(kirchdyn(natom,natom),hessdyn(degfred,degfred),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

kirchdyn = cmplx(zero,zero,kind=8)
hessdyn  = cmplx(zero,zero,kind=8)

do i=1,kirchoff%nnzero

  ii   = kirchoff%rows(i)
  jj   = kirchoff%columns(i)
  xyzi = unit_cell%aunit(1)%atom(ii)%x

  iaunit = unit_cell%aunit(1)%atom(ii)%iaunit
  jaunit = unit_cell%aunit(1)%atom(jj)%iaunit
  ichain = unit_cell%aunit(1)%atom(ii)%ichain
  jchain = unit_cell%aunit(1)%atom(jj)%ichain

  if(input%multint) then

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

  !  print *, 'multint'
    xyzj = unit_cell%aunit(1)%atom(jj)%x+vec
    dxyz = xyzj-xyzi
    dist  = distance(dxyz)

    if (dist .le. rcut) then
      if (dist .gt. zero) then
        call dyn_common(input, dist,ii,jj,dxyz, qvec,   &
                               kirchdyn, hessdyn,chaini,chainj)
      end if
    end if

    end do
    end do
    end do

  else

  print *, 'dynmat no selfinteractions' 
    distmp = distance(dxyz) 
    
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

      xyzj = unit_cell%aunit(1)%atom(jj)%x+vec
      dxyz = xyzj-xyzi
      dist  = distance(dxyz)
      if (dist .gt. distmp) dist = distmp

    end do
    end do
    end do

    if (dist .le. rcut) then
      if (dist .gt. zero) then
        call dyn_common(input, dist,ii,jj,dxyz, qvec,   &
                               kirchdyn, hessdyn,chaini,chainj)
      end if
    end if

  end if

end do

!convert to sparse
call zfull_to_zsparse(kirchdyn,tmpdyn)
deallocate(kirchdyn)
call  sparse_upper_to_gen(tmpdyn,kirchdyn_coor,"z")
call  sparse_deinit(tmpdyn)
call zfull_to_zsparse(hessdyn,tmpdyn)
deallocate(hessdyn)
call  sparse_upper_to_gen(tmpdyn,hessdyn_coor,"z")
call  sparse_deinit(tmpdyn)

if (trim(input%genm) .eq. "gnm") input%sparsity = kirchdyn_coor%sparsity
if (trim(input%genm) .eq. "enm") input%sparsity = hessdyn_coor%sparsity

call sparsehess_massweigh(input,unit_cell%aunit(1),hessdyn_coor,"z")

!if (trim(input%genm) .eq. "enm") then
!  do i =1, size(hessdyn(1,:))
!    do j =i, size(hessdyn(1,:))
!      if ( abs(real(hessdyn(i,j))) .gt.1.0d-08 ) then
!        write(13,'(2I5,F15.5)') i,j,real(hessdyn(i,j))
!        write(14,'(2I5,2F15.5)') i,j,real(hessdyn(i,j)), aimag(hessdyn(i,j))
!      end if
!    end do
!  end do
!end if


end subroutine

subroutine dynmat_calc (input,unit_cell,kirchoff,qvec, &
                             dynmat_coor)
!dmr added dec 8 2009 to simplify the dispersion calculation
!inout
use mod_inout,     only: zfull_to_zsparse
use mod_crysbuild, only: zmass_weightit
use mod_hessian,   only: sparsehess_massweigh
type(inp_par)        , intent(in out)     :: input
type(asym_List_Type) , intent(in)         :: unit_cell 
type(sparse)         , intent(in)         :: kirchoff
real(dp)             , intent(in)         :: qvec(3)
type(sparse)         , intent(out)         :: dynmat_coor
type(sparse)                               :: tmpdyn
complex(kind=8),allocatable               :: dynmat(:,:)
complex(kind=8),allocatable               :: kirchdyn(:,:),hessdyn(:,:)
!internal
integer                            :: i,ii,j,jj,iii,jjj,kkk
integer                            :: natom,nneighs,degfred,ialloc, &
                                      chaini,chainj,iaunit,jaunit,ichain,jchain
real(dp)                           :: rcut,dist,distmp
real(dp), dimension(3)             :: avec,bvec,cvec,vec
real(dp), dimension(3)             :: xyzi,xyzj,dxyz
real(sp)                           :: time1, time2
real(dp) :: dottest

chaini = 1
chainj = 1


call cpu_time(time1)
natom   = unit_cell%natoms
degfred = 3*natom
rcut    = input%rcut_start

if (allocated(kirchdyn)) deallocate(kirchdyn)
if (allocated(hessdyn)) deallocate(hessdyn)
allocate(kirchdyn(natom,natom),hessdyn(degfred,degfred),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

kirchdyn = cmplx(zero,zero,kind=8)
hessdyn  = cmplx(zero,zero,kind=8)

do i=1,kirchoff%nnzero

  ii   = kirchoff%rows(i)
  jj   = kirchoff%columns(i)
  xyzi = unit_cell%aunit(1)%atom(ii)%x

! for reach inter intra
  iaunit = unit_cell%aunit(1)%atom(ii)%iaunit
  jaunit = unit_cell%aunit(1)%atom(jj)%iaunit
  ichain = unit_cell%aunit(1)%atom(ii)%ichain
  jchain = unit_cell%aunit(1)%atom(jj)%ichain

  if(input%multint) then

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

  !  print *, 'multint'
    xyzj = unit_cell%aunit(1)%atom(jj)%x+vec
    dxyz = xyzj-xyzi
    dist  = distance(dxyz)

    if (dist .le. rcut) then
      if (dist .gt. zero) then
        call dyn_common(input, dist,ii,jj,dxyz, qvec,   &
                               kirchdyn, hessdyn,chaini,chainj)
      end if
    end if

    end do
    end do
    end do

  else

    print *, 'dynmat no selfinteractions' 
    ! no multiple interactions... find the shortest interaction and use it in construct    
    xyzj = unit_cell%aunit(1)%atom(jj)%x
    distmp = distance(dxyz) 
    
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


      xyzj = unit_cell%aunit(1)%atom(jj)%x+vec
      dxyz = xyzj-xyzi
      dist  = distance(dxyz)
      if (dist .gt. distmp) dist = distmp

    end do
    end do
    end do

    if (dist .le. rcut) then
      if (dist .gt. zero) then
        call dyn_common(input, dist,ii,jj,dxyz, qvec,   &
                               kirchdyn, hessdyn,chaini,chainj)
      end if
    end if

  end if

end do

!  call zmass_weightit("enm",unit_cell%aunit(1),hessdyn)
!  call zmass_weightit("gnm",unit_cell%aunit(1),kirchdyn)
!end if

if (trim(input%genm) .eq. "gnm") then
! could use pointer here
  deallocate(hessdyn)
  allocate(dynmat(natom,natom),stat=ialloc)
  if (ialloc /= 0) STOP "*** not enough memory ***"
  dynmat = kirchdyn
  deallocate(kirchdyn)
else
  deallocate(kirchdyn)
  allocate(dynmat(degfred,degfred),stat=ialloc)
  if (ialloc /= 0) STOP "*** not enough memory ***"
  dynmat = hessdyn
  deallocate(hessdyn)
end if  

call zfull_to_zsparse(dynmat,tmpdyn)
deallocate(dynmat)
call  sparse_upper_to_gen(tmpdyn,dynmat_coor,"z")
call  sparse_deinit(tmpdyn)

if(trim(input%genm) .ne. "gnm") call sparsehess_massweigh(input,unit_cell%aunit(1),dynmat_coor,"z")

end subroutine

subroutine bnm_dynmatrx (input,unit_cell,blocks,hessbnm_coor,qvec, &
                             hessdyn_coor)
!inout
use mod_inout,     only: zfull_to_zsparse, sparse_to_full
use mod_crysbuild, only: mass_weightit
use mod_bnm,       only: block,bnm_interact_tf
type(inp_par)           , intent(in out)     :: input
type(asym_List_Type)    , intent(in)         :: unit_cell
type(block),allocatable , intent(in)         :: blocks(:) 
type(sparse)            , intent(in)         :: hessbnm_coor
real(dp)   ,allocatable                      :: hessbnm(:,:)
real(dp)                , intent(in)         :: qvec(3)
type(sparse)   ,          intent(out)        :: hessdyn_coor
complex(kind=8),allocatable                  :: hessdyn(:,:)
type(sparse)                                 :: tmpdyn
!internal
integer                            :: i,ii,j,jj,iii,jjj,kkk
integer                            :: natom,nneighs,degfred,ialloc
real(dp)                           :: rcut,dist,distmp
real(dp), dimension(3)             :: avec,bvec,cvec,vec
real(dp), dimension(3)             :: xyzi,xyzj,dxyz
real(sp)                           :: time1, time2
integer :: nblocks,k,l
complex(kind=8) :: cnum
real(dp) :: q_dot_dr,sin_qdr,cos_qdr
logical :: tfint

print *, 'WRONG! passing in the q=0 and mult by phase factor doesnt work cuz of the implicit sums'
stop

call cpu_time(time1)
call sparse_to_full(hessbnm_coor,hessbnm)

natom   = unit_cell%natoms
degfred = size(hessbnm(1,:)) 
rcut    = input%rcut_start
nblocks = size(blocks)

if (allocated(hessdyn)) deallocate(hessdyn)
allocate(hessdyn(degfred,degfred),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

hessdyn  = cmplx(zero,zero,kind=8)

do i=1,nblocks
  
  xyzi = blocks(i)%cent
  ii = 6*(i-1)+1

  do j = i+1, nblocks
    jj = 6*(j-1)+1
    do iii=unit_cell%aneigh(1),unit_cell%aneigh(2)
    avec = real(iii,kind=dp)*unit_cell%avec
    do jjj=unit_cell%bneigh(1),unit_cell%bneigh(2)
    bvec = real(jjj,kind=dp)*unit_cell%bvec
    do kkk=unit_cell%cneigh(1),unit_cell%cneigh(2)
    cvec = real(kkk,kind=dp)*unit_cell%cvec
    vec = avec + bvec + cvec
! could be slow... this is a tricky thing.
    call bnm_interact_tf(input,unit_cell%aunit(1),unit_cell%aunit(1),vec,blocks(i),blocks(j),tfint)
  
    if (tfint) then
      print '(A30,2I4,3F10.4)', 'blocksij interact with vec:', i,j,vec
      xyzj = blocks(j)%cent+vec
      dxyz = xyzj-xyzi
      dist  = distance(dxyz)
      q_dot_dr = dot_product(qvec,dxyz)
      sin_qdr  = dsin(q_dot_dr)
      cos_qdr  = dcos(q_dot_dr)
      if (dist > zero) then
      ! off diagonal
        do k=0,5
          do l = 0,5    
            cnum = cmplx(hessbnm(ii+k,jj+l)*cos_qdr, &
                     hessbnm(ii+k,jj+l)*sin_qdr, kind=8)
            hessdyn(ii+k,jj+l) = hessdyn(ii+k,jj+l)+cnum
          end do
        end do 

!diagonal
        do k=0,5
          do l = 0,5
                         hessdyn(ii+k,ii+l) = hessdyn(ii+k,ii+l) - hessbnm(ii+k,jj+l)
            if(ii.ne.jj) hessdyn(jj+k,jj+l) = hessdyn(jj+k,jj+l) - hessbnm(ii+k,jj+l)
          end do
        end do
      end if
    end if

    end do
    end do
    end do
  end do
end do


call zfull_to_zsparse(hessdyn,tmpdyn)
!deallocate(hessdyn)
call  sparse_upper_to_gen(tmpdyn,hessdyn_coor,"z")
call  sparse_deinit(tmpdyn)

input%sparsity = hessdyn_coor%sparsity

  do i =1, size(hessdyn(1,:))
    do j =1, size(hessdyn(1,:))
      if ( abs(real(hessdyn(i,j))) .gt.1.0d-08 ) then
        write(13,'(2I5,F15.5)') i,j,real(hessdyn(i,j))
        write(14,'(2I5,2F15.5)') i,j,real(hessdyn(i,j)), aimag(hessdyn(i,j))
      end if
    end do
  end do

end subroutine

subroutine dyn_common(input,dist,ii,jj,dxyz,qvec, &
                          kirchdyn,hessdyn,chaini,chainj)
use mod_hessblock, only: sprnghessblock,reach_hessblock_temperat,reach_hessblock
type(inp_par),        intent(in)            :: input
integer ,             intent(in)            :: ii,jj  
real(dp),             intent(in)            :: dist,dxyz(3),qvec(3)
complex(kind=8),allocatable, intent(in out) :: kirchdyn(:,:),hessdyn(:,:)
integer , optional,   intent(in)            :: chaini,chainj
integer                                     :: i,j,ll,mm
real(dp)                                    :: hessblock(3,3)
complex(kind=8)                             :: dynblock(3,3),cnum
real(dp)                                    :: sin_qdr, cos_qdr, q_dot_dr


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


! factor of two pi included in the bz mesh
q_dot_dr = dot_product(qvec,dxyz)
sin_qdr  = dsin(q_dot_dr)
cos_qdr  = dcos(q_dot_dr)


!kirchdyn calcs 
kirchdyn(ii,jj) = kirchdyn(ii,jj) + cmplx(-one*cos_qdr,-one*sin_qdr,kind=8) ! neg one for interactions
! self term
kirchdyn(ii,ii) = kirchdyn(ii,ii)+one
if (ii.ne.jj) kirchdyn(jj,jj) = kirchdyn(jj,jj)+one 

! hess dyn
! off diagonal blocks
do ll=0,2
  do mm=0,2
    cnum                = cmplx(hessblock(1+ll,1+mm)*cos_qdr, &
                                hessblock(1+ll,1+mm)*sin_qdr, kind=8) 
    hessdyn(i+ll,j+mm)  = hessdyn(i+ll,j+mm) + cnum 
  end do
end do

!diagonal blocks
do ll=0,2
  do mm=0,2
    hessdyn(i+ll,i+mm)                = hessdyn(i+ll,i+mm) - hessblock(1+ll,1+mm)
    if (ii.ne.jj)  hessdyn(j+ll,j+mm) = hessdyn(j+ll,j+mm) - hessblock(1+ll,1+mm)
  end do
end do

end subroutine

subroutine brillouin_mesh(input,cell,linmesh)
! nq will be the parameter that gives the number of q's in each recip direction
!
! modes of fmesh building
! full,oct1,oct2, axes, spec,astr,bstr,cstr,star                    
!         full: all of it
!         axes: scan along recip direction from zero to wherever 
!         oct1: positive octant
!         oct2: negative h octant
!         spec: same as oct1 but with nq=3  combos of 0.0 0.25 0.5 
!         *str: is for a given direction (100, 010, 001) 
!         
use cfml_crystal_metrics,              only: Crystal_Cell_Type
type(inp_par),            intent(in out)  :: input ! bztyp: full, oct1,oct2.. axes, spec
type(crystal_cell_type),  intent(in)      :: cell
real(dp), allocatable,    intent(out)     :: linmesh(:,:)
real(dp), allocatable                     :: fmesh(:,:,:,:)
real(dp), allocatable                     :: axes(:),full(:)        
real(dp)                                  :: dnq,rli    
integer                                   :: ialloc,i,j,k,l,m,nq,ndim,fm_nq                

nq         = input%qdiv
if (nq .lt. 1) stop "empty Brillouin Zone, increase nq"

! set dimension for linearized mesh
! ?star is the direction, a in a* dir and so forth.  Dstr is diagonal
if ( (trim(input%bztyp) .eq. "read")  .or. &
     (trim(input%bztyp) .eq. "astr") .or. &
     (trim(input%bztyp) .eq. "abstr") .or. &
     (trim(input%bztyp) .eq. "bstr")  .or. &
     (trim(input%bztyp) .eq. "cstr")  .or. &
     (trim(input%bztyp) .eq. "dstr") )      then
  ndim = nq
else if (trim(input%bztyp) .eq. "star") then
  ndim = 3*nq-2 ! minus for the two extra 0.0s
else if (trim(input%bztyp) .eq. "full") then
  ndim  = nq**3 + 3*(nq-1)*nq**2 ! number of unique points in full 
else
  ndim = nq*nq*nq
end if

if (allocated(full))  deallocate(full)
if (allocated(linmesh)) deallocate(linmesh)
!allocate(full(4*nq),axes(nq), stat=ialloc)
allocate(axes(nq), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

if (allocated(linmesh)) deallocate(linmesh)
allocate(linmesh(ndim,3), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

linmesh=zero

! DMR mesh gen is not so complicated 
! historical monkhurst and pack
! u_r = (2r-nq-1)/2nq   where r=1,2,3...nq
! 
! adjusted here to just get a mesh that goes out to qmax 
! mesh point spacing will have asymetry of recip unit cell built in

! nq -1 for points not including zero
! if one point is called dnq=1 
if (nq .gt. 1) dnq=real(nq-1,kind=dp)
if (nq .eq. 1) then 
  dnq=real(nq,kind=dp)
  print *, "BZ with one point has been requested qmag=", input%qmax
  if (trim(input%bztyp) .eq. "full" ) stop "full and nq=1 is incompatible"
end if

do i=1,nq
  if (nq .eq. 1) rli=real(i,kind=dp)
  if (nq .gt. 1) rli=real(i-1,kind=dp)

 ! putting factor of two pi in here
  axes(i)=two*pi* rli * input%qmax/dnq
end do

!next we fill in the linmesh for all the different treatments of the BZ sampling

if (trim(input%bztyp) .eq. "read") then
  linmesh(:,1)=axes(:)*input%qdir(1)      
  linmesh(:,2)=axes(:)*input%qdir(2)      
  linmesh(:,3)=axes(:)*input%qdir(3)      
else if (trim(input%bztyp) .eq. "astr") then
  linmesh(:,1)=axes(:)          
else if (trim(input%bztyp) .eq. "abstr") then
  linmesh(:,1)=axes(:)          
  linmesh(:,2)=axes(:)          
else if (trim(input%bztyp) .eq. "bstr") then
  linmesh(:,2)=axes(:)          
else if (trim(input%bztyp) .eq. "cstr") then
  linmesh(:,3)=axes(:)   
else if (trim(input%bztyp) .eq. "dstr") then
  linmesh(:,1)=axes(:)   
  linmesh(:,2)=axes(:)   
  linmesh(:,3)=axes(:) 
else if (trim(input%bztyp) .eq. "star") then
! all the axis together
  linmesh(     1:    nq,1) = axes(:)       
  linmesh(  nq+1:2*nq-1,2) = axes(2:nq)       
  linmesh(  2*nq:3*nq-2,3) = axes(2:nq)       

else
    l=1
    do i=1,nq
      do j=1,nq
        do k=1,nq
! set up to take only unique entries.  the octants share faces
! oct1: + + +
! oct2: - + +
! oct3: + - +
! oct4: + + -
!
          if (trim(input%bztyp) .eq. "full" ) then
            linmesh(l,:) = (/axes(i),axes(j),axes(k)/)
            l=l+1
            if (i.eq.1) then  ! if x=0
              if (j.ne.1) then  ! if y!=0 take oct3
                linmesh(l,:) = (/axes(i),-axes(j),axes(k)/)
                l=l+1
              end if
              if (k.ne.1) then  ! if z!=0 take oct4
                linmesh(l,:) = (/axes(i),axes(j),-axes(k)/)
                l=l+1
              end if
            !d
            else if (j.eq.1) then ! if y=0
              if (i.ne.1) then  ! if x!=0 take oct2
                linmesh(l,:) = (/-axes(i),axes(j),axes(k)/)
                l=l+1
              end if
              if ((i.ne.1) .and. (k.ne.1)) then ! watching out for overcounting
                ! if x and z !=0 take oct4: for nq=3, 4 entries
                linmesh(l,:) = (/axes(i),axes(j),-axes(k)/)
                l=l+1
              end if
            !d
            else if (k.eq.1) then
              if ((i.ne.1) .and. (j.ne.1)) then
                ! if x and y !=0 take oct3 and oct2: for nq=3, 4 entries each 
                linmesh(l,:) = (/-axes(i),axes(j),axes(k)/)
                linmesh(l+1,:) = (/axes(i),-axes(j),axes(k)/)
                l=l+2
              end if
            else ! nonzero block stamped about.  for nq=3, 8*3=24
              linmesh(l,:)   = (/-axes(i),axes(j),axes(k)/)
              linmesh(l+1,:) = (/axes(i),-axes(j),axes(k)/)
              linmesh(l+2,:) = (/axes(i),axes(j),-axes(k)/)
              l=l+3
            end if
          end if

          if (trim(input%bztyp) .eq. "oct1" ) linmesh(l,:) = (/axes(i),axes(j),axes(k)/)
          if (trim(input%bztyp) .eq. "oct2" ) linmesh(l,:) = (/-axes(i),axes(j),axes(k)/)
          if (trim(input%bztyp) .eq. "oct3" ) linmesh(l,:) = (/axes(i),-axes(j),axes(k)/)
          if (trim(input%bztyp) .eq. "oct4" ) linmesh(l,:) = (/axes(i),axes(j),-axes(k)/)
          if (trim(input%bztyp) .eq. "oct5" ) linmesh(l,:) = (/-axes(i),-axes(j),axes(k)/)
          if (trim(input%bztyp) .eq. "oct6" ) linmesh(l,:) = (/-axes(i),axes(j),-axes(k)/)
          if (trim(input%bztyp) .eq. "oct7" ) linmesh(l,:) = (/axes(i),-axes(j),-axes(k)/)
          if (trim(input%bztyp) .eq. "oct8" ) linmesh(l,:) = (/-axes(i),-axes(j),-axes(k)/)
          if (trim(input%bztyp) .ne. "full") l = l+1
        end do
      end do
    end do

end if

!   multiply by transformation matrix that converts cartesian to fractional 
do i=1,size(linmesh(:,1))
!  print '(3F10.2)', linmesh(i,:) 
  linmesh(i,:) = matmul(linmesh(i,:),cell%orth_cr_cel)
end do

end subroutine

subroutine qdisp_split(input,cell,unit_cell,kirchoff,hessian,&
                        vcov,vcov_acoustic,vcov_optic)
! separate the vcov into 3: acoustic modes, optic, all 

use mod_linalg, only : lapack_zeig 
use mod_types,  only : sparse,dos,dos_init,dos_append,dos_histgen
use mod_math,   only : realrecip
type(inp_par)          , intent(in out) :: input
type (Crystal_Cell_Type),intent(in) :: cell
type( asym_list_type), intent(in)  :: unit_cell
type(sparse)          , intent(in)  :: kirchoff,hessian
real(dp),intent(out), allocatable   :: vcov(:,:),vcov_acoustic(:,:),vcov_optic(:,:)
real(dp),             allocatable   :: linmesh(:,:)
type(sparse)                        :: kirchdyn,hessdyn         
complex(kind=8),      allocatable   :: zvals(:),zvecs(:,:),fdynmat(:,:)
!,kirchdyn(:,:),hessdyn(:,:)
real(dp)       ,      allocatable   :: vals(:)
real(dp)                            :: qvec(3),qmag,fqvec(3),test_val 
integer                             :: natoms,i,j,k,ier,znfrq
integer                             :: nq,nfrq

print *, "not quite right, the qvec=000 is overweighted see qdispers"
 
nq = input%qdiv
nfrq = input%nfrqs
if (input%print_level .ge. 2) then
open(unit=32,file=trim(input%fileroot)//"-"//trim(input%fctyp)// &
                  "-"//trim(input%bztyp)//"-dispers.txt", & 
     status='unknown',action='write', &
     iostat=ier)
end if


call brillouin_mesh(input,cell,linmesh)

natoms=unit_cell%natoms
!initialize vcov
if (trim(input%genm) .eq. "enm") then
  call matrx2d_init(natoms*3,vcov)
  call matrx2d_init(natoms*3,vcov_acoustic)
  call matrx2d_init(natoms*3,vcov_optic)
end if
if (trim(input%genm) .eq. "gnm") then
  call matrx2d_init(natoms,vcov)
  call matrx2d_init(natoms,vcov_acoustic)
  call matrx2d_init(natoms,vcov_optic)
end if
vcov=zero ; vcov_acoustic=zero ;  vcov_optic = zero

do i=1,size(linmesh(:,1))
  znfrq=nfrq
  qvec(:) = linmesh(i,:)
  fqvec(:)=matmul(qvec,cell%Cr_orth_cel)
  fqvec=fqvec/(two*pi)
 
  call dynmatrx(input,unit_cell,kirchoff,qvec,kirchdyn,hessdyn)

  if (trim(input%genm) .eq. "enm")  then
    if (trim(input%eigtyp) .eq. "lapack") then
      call zsparse_to_full(hessdyn,fdynmat)
      call lapack_zeig(fdynmat, zvals, zvecs)
    else
      call mod_arpack_zndrv1(input,znfrq,hessdyn,zvals,zvecs)
    end if
  else if(trim(input%genm) .eq. "gnm")  then
    if (trim(input%eigtyp) .eq. "lapack") then
      call zsparse_to_full(kirchdyn,fdynmat)
      call lapack_zeig(fdynmat, zvals, zvecs)
    else
      call mod_arpack_zndrv1(input,znfrq,kirchdyn,zvals,zvecs)
    end if
  end if

  call d_zvcov_split(input%temperature,zvals,zvecs,vcov,vcov_acoustic,vcov_optic)

if (input%print_level .ge. 1) write(32,'(F9.3,$)') fqvec,dsqrt(abs(real(zvals)))
if (input%print_level .ge. 1) write(32, *), ''
end do

vcov=vcov/real(size(linmesh(:,1)),dp)
vcov_acoustic=vcov_acoustic/real(size(linmesh(:,1)),dp)
vcov_optic=vcov_optic/real(size(linmesh(:,1)),dp)

input%nbz = size(linmesh(:,1))
print '(A45,I14,A35)', 'Brillouin zone qvectors sampled:',size(linmesh(:,1))


end subroutine

subroutine qdispers(input,cell,unit_cell,kirchoff,hessian,vcov,states)
!DMR: 01-12-2009     when qmag = 0.0 use real stuff with hessian.
! DMR KEEP NOTE:  qmag is the same as shitdot below
!  fqvec(:)=matmul(qvec,cell%Cr_orth_cel)
!  shit(:)=matmul(cell%GR,fqvec)
!  shitdot=dsqrt(dot_product(fqvec,shit))
! DMR KEEP NOTE: reciprocal lattice does not have factor of 2pi
!  call realrecip(cell%cell,cell%ang,crap,crapp)
!  print '(A2,6F10.4)','D', cell%cell, cell%ang
!  print '(A2,6F10.4)','R', cell%rcell,cell%rang 
!  print '(A2,6F10.4)','R', crap,crapp 
use mod_arpack, only : mod_arpack_zndrv1 
use mod_linalg, only : lapack_zeig 
use mod_types,  only : sparse,dos,dos_init,dos_append,dos_histgen
use mod_math,   only : realrecip
use mod_hessian,   only: fullhess_massweigh,eig_to_THZ
type(inp_par)          , intent(in out) :: input
type (Crystal_Cell_Type),intent(in) :: cell
type( asym_list_type), intent(in)  :: unit_cell 
type(sparse)          , intent(in)  :: kirchoff,hessian
real(dp),intent(out), allocatable   :: vcov(:,:)
type(dos) ,intent(out)              :: states
real(dp),             allocatable   :: linmesh(:,:)
type(sparse)                        :: kirchdyn,hessdyn,hessian_coor         
complex(kind=8),      allocatable   :: zvals(:),zvecs(:,:),fdynmat(:,:),zvcov(:,:)
!,kirchdyn(:,:),hessdyn(:,:)
real(dp)       ,      allocatable   :: vals(:),tmpvcov(:,:),dvals(:),dvecs(:,:),valsthz(:)
real(dp)                            :: vcov_mult,test_val,qvec(3),qmag,fqvec(3) ,sumshit1,sumshit2
integer                             :: natoms,i,j,k,ier,znfrq,zero_sub=0
integer                             :: nq,nfrq,nzeros,prnt_lev

prnt_lev = 0
nq = input%qdiv
nfrq = input%nfrqs

if (input%print_level .ge. prnt_lev) then
open(unit=32,file=trim(input%fileroot)//"-"//trim(input%fctyp)// &
                  "-"//trim(input%bztyp)//"-dispers.txt", & 
     status='unknown',action='write', &
     iostat=ier)
end if

call brillouin_mesh(input,cell,linmesh)

natoms=unit_cell%natoms
!initialize vcov
if (trim(input%genm) .eq. "enm") call matrx2d_init(natoms*3,vcov)
if (trim(input%genm) .eq. "gnm") call matrx2d_init(natoms,vcov)
if (trim(input%genm) .eq. "enm") call matrx2d_init(natoms*3,tmpvcov)
if (trim(input%genm) .eq. "gnm") call matrx2d_init(natoms,tmpvcov)


! use THZ for dos
if (input%dosgo) then
    if (trim(input%eigtyp) .eq. "lapack") then
      call dos_init(75,states)
      !call dos_init(input%nfrqs/10,states)
    else
      call dos_init(input%nfrqs/2,states)
    end if
end if

do i=1,size(linmesh(:,1))
  znfrq=nfrq
  qvec(:) = linmesh(i,:)
  fqvec(:)=matmul(qvec,cell%Cr_orth_cel)
  fqvec=fqvec/(two*pi)

  call dynmat_valvec(input,unit_cell,kirchoff,hessian,qvec,dvals,zvecs)

  test_val = dot_product(fqvec,fqvec)
  if (test_val .eq. zero) then
    zero_sub  = 1
    vcov_mult = one
  else
    vcov_mult = two 
  end if

  ! sum it up right
  tmpvcov=vcov
  vcov=zero
  call d_zvcov(input,dvals,zvecs,vcov)
! get the weighting of the vcov right.  The zero should be counted once
! while nonzeros correspond to 2 qs +-
  vcov = vcov_mult*vcov
  vcov = vcov+tmpvcov
  if (input%dosgo) then
    call eig_to_THz(dvals,vals)
    call dos_append(vals,states)
  end if
  if (input%print_level .ge. prnt_lev) write(32,'(F9.3,$)') fqvec,dsqrt(abs(real(zvals)))
  if (input%print_level .ge. prnt_lev) write(32, *), ''

end do

!compute the appropriate average vcov  2*N_qvec -1 to count qvec = 000 only once
! this way the vcov is always real.
vcov=vcov/real(2*size(linmesh(:,1))-zero_sub,dp)

call fullhess_massweigh(input, unit_cell%aunit(1),vcov)

input%nbz = size(linmesh(:,1))
print '(A45,I14,A35)', 'Brillouin zone qvectors sampled:',size(linmesh(:,1))

if (input%dosgo) call dos_histgen(states)

end subroutine

subroutine dispersion_calc(input,cell,unit_cell,kirchoff,hessian,states,blocks)
use mod_arpack,   only : mod_arpack_zndrv1 
use mod_linalg,   only : lapack_zeig,zheevr_zeig 
use mod_types,    only : sparse,dos,dos_init,dos_append,dos_histgen
use mod_math,     only : realrecip
use mod_hessian,  only : eig_to_THZ
use mod_bnm,      only : block,dynmat_proj
type(inp_par)          ,  intent(in out)        :: input
type (Crystal_Cell_Type), intent(in)            :: cell
type( asym_list_type),    intent(in)            :: unit_cell 
type(sparse)          ,   intent(in)            :: kirchoff,hessian
type(dos) ,               intent(out)           :: states
type(block), allocatable, intent(in), optional  :: blocks(:)
real(dp),             allocatable   :: linmesh(:,:)
type(sparse)                        :: kirchdyn,hessdyn         
type(sparse)   :: dynmat,bnmdyn
!,kirchdyn(:,:),hessdyn(:,:)
real(dp)       ,      allocatable   :: vals(:),valsthz(:)
complex(dp)    ,      allocatable   :: fdynmat(:,:),fbnmdyn(:,:),zval(:),zvec(:,:)
real(dp)                            :: qvec(3),qmag,fqvec(3)
integer                             :: natoms,i,j,k,ier
integer                             :: nq,nfrq,nzeros,prnt_lev
real(sp) :: time1,time2


if (trim(input%genm) .eq. "bnm") then
  if (present(blocks)) then
    ! do nothing, keep going
  else
    print *, "dispersion_calc> BNM called, but blocks not passed"
    stop
  end if
end if

prnt_lev = 2
nq = input%qdiv
nfrq = input%last_mode
if (nfrq .gt. 100) then
  print *, 'nfrq>100, did you set maxfrq correctly?'
  print *, 'setting nfrq = 100'
  print *, 'see mod_lattdyn and change if you want >100'
  nfrq = 100
end if
allocate(vals(nfrq),stat=ier)

if (input%print_level .ge. prnt_lev) then
open(unit=32,file=trim(input%fileroot)//"-"//trim(input%genm)//"-"//trim(input%fctyp)// &
                  "-"//trim(input%bztyp)//"-dispers.txt", & 
     status='unknown',action='write', &
     iostat=ier)
end if

call brillouin_mesh(input,cell,linmesh)

natoms=unit_cell%natoms

! use THZ for dos
if (input%dosgo) then
    if (trim(input%eigtyp) .eq. "lapack") then
      call dos_init(75,states)
    else
      call dos_init(input%nfrqs/2,states)
    end if
end if

do i=1,size(linmesh(:,1))

  qvec(:) = linmesh(i,:)
  fqvec(:)=matmul(qvec,cell%Cr_orth_cel)
  fqvec=fqvec/(two*pi)
  print '(A8,3F10.2)', "wavevect", fqvec
  call dynmat_calc (input,unit_cell,kirchoff,qvec, &
                             dynmat)
  if (trim(input%genm) .eq. "bnm") then
    call dynmat_proj(dynmat,blocks,bnmdyn)
    call zsparse_to_full(bnmdyn,fdynmat)
  else
    call zsparse_to_full(dynmat,fdynmat)
  end if

 ! call zheevr_zeig (1,input%last_mode,fdynmat, &
 !                            vals)
  call lapack_zeig(fdynmat, zval, zvec)

!print *, "about to calc arpack"
!  call cpu_time(time1)
!  if (trim(input%genm) .eq. "bnm") then
!    call mod_arpack_zndrv1 (input,nfrq,bnmdyn,ZVAL,ZVEC,.false.) 
!  else
!    call mod_arpack_zndrv1 (input,nfrq,dynmat,ZVAL,ZVEC,.false.) 
!  end if
  vals(1:nfrq) = real(zval(1:nfrq))
!  call cpu_time(time2)
!  print *, "arpack time:", time2-time1

  if (input%dosgo) then
    call eig_to_THz(vals,valsthz)
    call dos_append(valsthz,states)
  end if
  do j = 1, size(vals)
    if (abs(valsthz(j)) .le. 1.0D-07) valsthz(j) = zero
  end do

  if (input%print_level .ge. prnt_lev) write(32,'(F9.3,$)') fqvec,valsthz
  if (input%print_level .ge. prnt_lev) write(32, *), ''

end do

input%nbz = size(linmesh(:,1))
print '(A45,I14,A35)', 'Brillouin zone qvectors sampled:',size(linmesh(:,1))

!if (input%dosgo) call dos_histgen(states)
call dos_histgen(states)

end subroutine

subroutine bnm_qdispers(input,cell,blocks,unit_cell,kirchoff,hessian,vcov,states)
use mod_arpack,   only : mod_arpack_zndrv1 
use mod_linalg,   only : lapack_zeig 
use mod_types,    only : sparse,dos,dos_init,dos_append,dos_histgen
use mod_math,     only : realrecip
use mod_hessian,  only : fullhess_massweigh,eig_to_THZ
use mod_bnm,      only : block
use mod_inout,  only : sparse_statalloc,sparse_deinit
type(inp_par)          ,  intent(in out)  :: input
type (Crystal_Cell_Type), intent(in)      :: cell
type(block), allocatable, intent(in)      :: blocks(:)
type( asym_list_type),    intent(in)      :: unit_cell 
type(sparse)          ,   intent(in)      :: kirchoff,hessian
!real(dp), allocatable  , intent(in)  :: hessian(:,:)
real(dp),intent(out), allocatable   :: vcov(:,:)
type(dos) ,intent(out)              :: states
real(dp),             allocatable   :: linmesh(:,:)
type(sparse)                        :: kirchdyn,hessdyn,hessian_coor         
complex(kind=8),      allocatable   :: zvals(:),zvecs(:,:),fdynmat(:,:),zvcov(:,:)
!,kirchdyn(:,:),hessdyn(:,:)
real(dp)       ,      allocatable   :: vals(:),tmpvcov(:,:),dvals(:),dvecs(:,:),valsthz(:)
real(dp)                            :: vcov_mult,test_val,qvec(3),qmag,fqvec(3) ,sumshit1,sumshit2
integer                             :: natoms,i,j,k,ier,znfrq,zero_sub=0
integer                             :: nq,nfrq,nzeros,prnt_lev

prnt_lev = 2
nq = input%qdiv
nfrq = input%nfrqs

if (input%print_level .ge. prnt_lev) then
open(unit=32,file=trim(input%fileroot)//"-"//trim(input%fctyp)// &
                  "-"//trim(input%bztyp)//"-dispers.txt", & 
     status='unknown',action='write', &
     iostat=ier)
end if

call brillouin_mesh(input,cell,linmesh)

natoms=unit_cell%natoms
!initialize vcov
call matrx2d_init(6*size(blocks),vcov)
call matrx2d_init(6*size(blocks),tmpvcov)



! use THZ for dos
if (input%dosgo) then
    if (trim(input%eigtyp) .eq. "lapack") then
      call dos_init(75,states)
    else
      call dos_init(input%nfrqs/2,states)
    end if
end if

do i=1,size(linmesh(:,1))
  znfrq=nfrq
  qvec(:) = linmesh(i,:)
  fqvec(:)=matmul(qvec,cell%Cr_orth_cel)
  fqvec=fqvec/(two*pi)

  call bnm_dynmat_valvec(input,blocks,unit_cell,kirchoff,hessian,qvec,dvals,zvecs)

  test_val = dot_product(fqvec,fqvec)
  if (test_val .eq. zero) then
    zero_sub  = 1
    vcov_mult = one
  else
    vcov_mult = two 
  end if

  ! sum it up right
  tmpvcov=vcov
  vcov=zero
  call d_zvcov(input,dvals,zvecs,vcov)
! get the weighting of the vcov right.  The zero should be counted once
! while nonzeros correspond to 2 qs +-
  vcov = vcov_mult*vcov
  vcov = vcov+tmpvcov
  if (input%dosgo) then
    call eig_to_THz(dvals,vals)
    call dos_append(vals,states)
  end if
  if (input%print_level .ge. prnt_lev) write(32,'(F9.3,$)') fqvec,dsqrt(abs(real(zvals)))
  if (input%print_level .ge. prnt_lev) write(32, *), ''

  !call dos_statalloc(states)

end do

!compute the appropriate average vcov  2*N_qvec -1 to count qvec = 000 only once
! this way the vcov is always real.
vcov=vcov/real(2*size(linmesh(:,1))-zero_sub,dp)

input%nbz = size(linmesh(:,1))
print '(A45,I14,A35)', 'Brillouin zone qvectors sampled:',size(linmesh(:,1))

call dos_histgen(states)

end subroutine

subroutine dynmat_valvec(input,unit_cell,kirchoff,hessian,qvec,dvals,zvecs)
! DMR: 01-15-2009  subroutine dynmat_valvec
! constructs the dynamic matrix and computes eigenvalues and eigenvectors 
!   for a give wave vector.
!   input:
!   type(neigh_List_Type)         ,intent(in)   :: crystal
!   type(sparse)                  ,intent(in)   :: kirchoff,hessian
!   real(dp)                      ,intent(in)  :: qvec(3)
!    
! qvec is the fractional wavevector dotted with the cartesian to fractional transformation
!   matrix.  This allows us to use the cartesian coordinates in "crystal"
use mod_arpack, only : mod_arpack_zndrv1 
use mod_linalg, only : lapack_zeig 
use mod_types,  only : sparse
use mod_inout,  only : sparse_statalloc,sparse_deinit
use mod_hessian,only : kirchess_valvec
type(inp_par)        ,         intent(in out)  :: input
type(asym_List_Type),          intent(in)      :: unit_cell
type(sparse)         ,         intent(in)      :: kirchoff,hessian
real(dp)             ,         intent(in)      :: qvec(3)
complex(kind=8), allocatable,  intent(out)     :: zvecs(:,:)
real(dp)       , allocatable,  intent(out)     :: dvals(:) 
complex(kind=8), allocatable                   :: zvals(:)
real(dp)       , allocatable                   :: dvecs(:,:)
type(sparse)                                   :: kirchdyn,hessdyn         
complex(kind=8), allocatable                   :: fdynmat(:,:)
real(dp) :: test_val
integer :: znfrq

  znfrq=input%nfrqs

  test_val = dot_product(qvec,qvec)
  if (test_val .eq. zero) then

    call kirchess_valvec(input,kirchoff,hessian, dvals, dvecs)
    call matrx_real2complex(dvecs,zvecs)

  else

    call dynmatrx(input,unit_cell,kirchoff,qvec,kirchdyn,hessdyn)

    if (trim(input%genm) .eq. "enm")  then
      call sparse_deinit(kirchdyn)
      if (trim(input%eigtyp) .eq. "lapack") then
        call zsparse_to_full(hessdyn,fdynmat)
        call sparse_deinit(hessdyn)
        call lapack_zeig(fdynmat, zvals, zvecs)
        deallocate(fdynmat)
        
      else
        call mod_arpack_zndrv1(input,znfrq,hessdyn,zvals,zvecs)
        call sparse_deinit(hessdyn)
      end if
    else if(trim(input%genm) .eq. "gnm")  then
      if (trim(input%eigtyp) .eq. "lapack") then
        call zsparse_to_full(kirchdyn,fdynmat)
        call sparse_deinit(kirchdyn)
        call lapack_zeig(fdynmat, zvals, zvecs)
        deallocate(fdynmat)
      else
        call mod_arpack_zndrv1(input,znfrq,kirchdyn,zvals,zvecs)
        call sparse_deinit(kirchdyn)
      end if
    end if
    call vector_complex2real(zvals,dvals)

  end if 


end subroutine

subroutine bnm_dynmat_valvec(input,blocks,unit_cell,kirchoff,hessian,qvec,dvals,zvecs)
! DMR: 01-15-2009  subroutine dynmat_valvec
! constructs the dynamic matrix and computes eigenvalues and eigenvectors 
!   for a give wave vector.
!   input:
!   type(neigh_List_Type)         ,intent(in)   :: crystal
!   type(sparse)                  ,intent(in)   :: kirchoff,hessian
!   real(dp)                      ,intent(in)  :: qvec(3)
!    
! qvec is the fractional wavevector dotted with the cartesian to fractional transformation
!   matrix.  This allows us to use the cartesian coordinates in "crystal"
use mod_arpack, only : mod_arpack_zndrv1 
use mod_linalg, only : lapack_zeig 
use mod_types,  only : sparse
use mod_hessian,only : kirchess_valvec
use mod_bnm
use mod_inout,  only : sparse_statalloc,sparse_deinit
type(inp_par)        ,         intent(in out)  :: input
type(block) ,allocatable   ,   intent(in)      :: blocks(:)
type(asym_List_Type),          intent(in)      :: unit_cell
type(sparse)         ,         intent(in)      :: kirchoff,hessian
real(dp)             ,         intent(in)      :: qvec(3)
complex(kind=8), allocatable,  intent(out)     :: zvecs(:,:)
real(dp)       , allocatable,  intent(out)     :: dvals(:) 
complex(kind=8), allocatable                   :: zvals(:)
real(dp)       , allocatable                   :: dvecs(:,:)
type(sparse)                                   :: kirchdyn,hessdyn,bnmhess,bnmdyn         
complex(kind=8), allocatable                   :: fdynmat(:,:)
real(dp) :: test_val
integer :: znfrq,i,ii,j
real(dp) :: time1,time2,time3

  znfrq=input%nfrqs
  if(trim(input%genm) .eq. "gnm")  stop "gnm passed to bnm, which &
                                         should be impossible"

  test_val = dot_product(qvec,qvec)
  if (test_val .eq. zero) then

    call hess_proj(hessian,blocks,bnmhess)
    call bnm_eigenal_valvec(input,bnmhess,dvals,dvecs)
    call matrx_real2complex(dvecs,zvecs)

  else

    call cpu_time(time1)
    print *, "entering dynmatrx"
    call dynmatrx(input,unit_cell,kirchoff,qvec,kirchdyn,hessdyn)
    ! dynmat_proj is the nut buster
    call cpu_time(time2)
    print *, 'time generating dynmat', time2-time1
    print *, "entering dynmat_proj"
    call dynmat_proj(hessdyn,blocks,bnmdyn) 
    call cpu_time(time3)
    print *, 'time projecting it', time3-time2
    call sparse_deinit(hessdyn) ! DMR added 011910 
    call sparse_deinit(kirchdyn) ! DMR added 011910 
    !call sparse_statalloc(hessdyn)

!    do ii = 1, bnmdyn%nnzero
!      write(769,'(2I4,F10.5)') bnmdyn%rows(ii), &
!                         bnmdyn%columns(ii),bnmdyn%zvalues(ii)
!    end do
 
    if (trim(input%eigtyp) .eq. "lapack") then
      call zsparse_to_full(bnmdyn,fdynmat)
      call sparse_deinit(bnmdyn) ! DMR added 011910 
      call lapack_zeig(fdynmat, zvals, zvecs)
      deallocate(fdynmat) ! DMR added 011910
    else
      call mod_arpack_zndrv1(input,znfrq,bnmdyn,zvals,zvecs)
      call sparse_deinit(bnmdyn) ! DMR added 011910 
    end if
    call vector_complex2real(zvals,dvals)
    deallocate(zvals) ! DMR added 011910
  end if 


end subroutine

subroutine matrx_real2complex(matreal,matcomplex)
! creates complex matrix from real one
real(dp), allocatable, intent(in)          :: matreal(:,:)
complex(kind=8), allocatable, intent(out)  :: matcomplex(:,:)
integer :: ialloc

allocate(matcomplex(size(matreal(:,1)),size(matreal(1,:))), stat = ialloc)

matcomplex = (zero,zero)
matcomplex = matreal

end subroutine

subroutine vector_complex2real(vecomplex,vecreal)
! creates real vector from complex one
! also tests that it is indeed real, useful for eigenvals of hermitian matrix
complex(kind=8), allocatable, intent(in)  :: vecomplex(:)
real(dp), allocatable, intent(out)          :: vecreal(:)
real(dp) :: test_complex
integer :: ialloc,i

allocate(vecreal(size(vecomplex)), stat = ialloc)

test_complex = zero

do i = 1,size(vecomplex)
  test_complex = test_complex+ abs(dimag(vecomplex(i)))
  vecreal(i) = dble(vecomplex(i))
end do

if (test_complex .gt. 1.0d-08 ) then 
  print *, "vector_complex2real> ERROR vector has nontrivial complex components."
  print *, "vector_complex2real> sum {|complex_part|} = ", test_complex
end if

end subroutine

!subroutine ordered_sparsematwrite (un,sparsemat)
!integer, intent(in) :: un
!type(sparse)        :: sparsemat

!  do i =1, hessdyn%nnzero
!    write(77,*) hessdyn%rows(i),hessdyn%columns(i),hessdyn%zvalues(i)
!  end do

!end subroutine

subroutine zqdispers(input,cell,unit_Cell,kirchoff,nq,nfrq,zvcov,states)
! DMR 1-6-2009      complex vcov matrix
!  this is no longer used, but will be kept for future referenc
use mod_arpack, only : mod_arpack_zndrv1 
use mod_linalg, only : lapack_zeig 
use mod_types,  only : sparse,dos,dos_init,dos_append,dos_histgen
use mod_math,   only : realrecip
type(inp_par)          , intent(in out) :: input
type (Crystal_Cell_Type),intent(in) :: cell
type( asym_list_type), intent(in)  :: unit_cell
type(sparse)          , intent(in)  :: kirchoff
integer, intent(in)                 :: nq,nfrq
complex(kind=8),intent(out), allocatable   :: zvcov(:,:)
type(dos) ,intent(out)              :: states
real(dp),             allocatable   :: linmesh(:,:)
type(sparse)                        :: kirchdyn,hessdyn         
complex(kind=8),      allocatable   :: zvals(:),zvecs(:,:),fdynmat(:,:),tmpzvcov(:,:)
!,kirchdyn(:,:),hessdyn(:,:)
real(dp)       ,      allocatable   :: vals(:)
real(dp)                            :: qvec(3),qmag,fqvec(3) 
integer                             :: natoms,prntlev,i,j,k,ier,znfrq,zero_sub

!print *, 'called old routine, stopping program'
!stop
prntlev = 2
if (input%print_level .ge. prntlev) then
open(unit=32,file=trim(input%fileroot)//"-"//trim(input%fctyp)// &
                  "-"//trim(input%bztyp)//"-dispers.txt", & 
     status='unknown',action='write', &
     iostat=ier)
end if


call brillouin_mesh(input,cell,linmesh)

natoms=unit_cell%natoms
!initialize vcov
if (trim(input%genm) .eq. "enm") call zmatrx2d_init(natoms*3,zvcov)
if (trim(input%genm) .eq. "gnm") call zmatrx2d_init(natoms,zvcov)
if (trim(input%genm) .eq. "enm") call zmatrx2d_init(natoms*3,tmpzvcov)
if (trim(input%genm) .eq. "gnm") call zmatrx2d_init(natoms,tmpzvcov)

if (.true.) then
    if (trim(input%eigtyp) .eq. "lapack") then
      call dos_init(75,states)
      !call dos_init(input%nfrqs/10,states)
    else
      call dos_init(input%nfrqs/2,states)
    end if
end if

do i=1,size(linmesh(:,1))
  znfrq=nfrq
  qvec(:) = linmesh(i,:)
  fqvec(:)=matmul(qvec,cell%Cr_orth_cel)
  fqvec=fqvec/(two*pi)
  call dynmatrx(input,unit_cell,kirchoff,qvec,kirchdyn,hessdyn)

  if (trim(input%genm) .eq. "enm")  then
    if (trim(input%eigtyp) .eq. "lapack") then
      call zsparse_to_full(hessdyn,fdynmat)
      call lapack_zeig(fdynmat, zvals, zvecs)
    else
      call mod_arpack_zndrv1(input,znfrq,hessdyn,zvals,zvecs)
    end if
  else if(trim(input%genm) .eq. "gnm")  then
    if (trim(input%eigtyp) .eq. "lapack") then
      call zsparse_to_full(kirchdyn,fdynmat)
      call lapack_zeig(fdynmat, zvals, zvecs)
    else
      call mod_arpack_zndrv1(input,znfrq,kirchdyn,zvals,zvecs)
    end if
  end if

  call zeig_to_THz(zvals,vals)
  tmpzvcov=zvcov
  zvcov=(zero,zero)
  !dmr temp
  print *, zvecs(:,input%last_mode)
  call z_zvcov(input%temperature,zvals,zvecs,zvcov)
! get the weighting of the vcov right.  The zero should be counted once
! while nonzeros correspond to 2 qs +-
  if ((qvec(1) .eq. zero) .and. (qvec(2) .eq. zero) .and. (qvec(3).eq.zero)) then
    zero_sub = 1
  else 
    zvcov = two*zvcov
  end if
  zvcov=zvcov+tmpzvcov
  if (.true.) call dos_append(vals,states)

if (input%print_level .ge. prntlev) write(32,'(F9.3,$)') fqvec,dsqrt(abs(real(zvals)))
if (input%print_level .ge. prntlev) write(32, *), ''
end do

!vcov=vcov/real(size(linmesh(:,1)),dp)
!compute the appropriate average vcov
zvcov=zvcov/real(2*size(linmesh(:,1))-zero_sub,dp)

!zvcov=zvcov/real(size(linmesh(:,1)),dp)

!call zmatrx2d_write(66,zvcov)

input%nbz = size(linmesh(:,1))
print '(A45,I14,A35)', 'Brillouin zone qvectors sampled:',size(linmesh(:,1))

if (.true.) call dos_histgen(states)


end subroutine

subroutine speedsnd(input,cell,unit_cell,kirchoff)
use mod_arpack, only : mod_arpack_zndrv1 
use mod_linalg, only : lapack_zeig 
use mod_types,  only : sparse
use mod_math,   only : realrecip
type(inp_par)          , intent(in out) :: input
type (Crystal_Cell_Type),intent(in) :: cell
type (asym_list_Type),  intent(in) :: unit_cell 
type(sparse), intent(in)            :: kirchoff
real(dp),             allocatable   :: linmesh(:,:)
complex(kind=8),      allocatable   :: dynmat(:,:)
complex(kind=8),      allocatable   :: zvals(:),zvecs(:,:)
real(dp)                            :: qvec(3),fqvec(3),qmag
real(dp),dimension(3)               :: spdsnd(3)
integer                             :: i,j,ier,znq=3

input%bztyp="astr"
call brillouin_mesh(input,cell,linmesh)
print '(A27,1x,I3,1x,A21,1x,F6.3)', 'computing speed of sound for', input%qdiv-1, &
                                  'points > 0.000 and <=', input%qmax
call sound_common(input,cell,unit_cell,kirchoff,linmesh,spdsnd)
print '(A12,1x,A5,1x,3F10.3)', "sound speed>",trim(input%bztyp),spdsnd

input%bztyp="bstr"
call brillouin_mesh(input,cell,linmesh)
call sound_common(input,cell,unit_cell,kirchoff,linmesh,spdsnd)
print '(A12,1x,A5,1x,3F10.3)', "sound speed>",trim(input%bztyp),spdsnd

input%bztyp="cstr"
call brillouin_mesh(input,cell,linmesh)
call sound_common(input,cell,unit_cell,kirchoff,linmesh,spdsnd)
print '(A12,1x,A5,1x,3F10.3)', "sound speed>",trim(input%bztyp),spdsnd

input%bztyp="dstr"
call brillouin_mesh(input,cell,linmesh)
call sound_common(input,cell,unit_cell,kirchoff,linmesh,spdsnd)
print '(A12,1x,A5,1x,3F10.3)', "sound speed>",trim(input%bztyp),spdsnd

end subroutine

subroutine sound_common(input,cell,unit_cell,kirchoff,linmesh,avgspd)
use mod_arpack, only : mod_arpack_zndrv1 
use mod_types,  only : sparse
type(inp_par)          , intent(in out)  :: input
type (crystal_Cell_Type),intent(in)  :: cell
type (asym_list_Type),  intent(in)  :: unit_cell
type(sparse),            intent(in)  :: kirchoff
real(dp),   allocatable, intent(in)  :: linmesh(:,:)
real(dp),                intent(out) :: avgspd(3) 
type(sparse)                         :: hessdyn,kirchdyn      
complex(kind=8),      allocatable    :: zvals(:),zvecs(:,:)
real(dp)       ,      allocatable    :: vals(:)
real(dp)                             :: qvec(3),qmag,spdsound(5,3),sumsnd(3),val
integer                              :: i,j,k,ier,istart,istop,itot,znq

sumsnd=zero
!why start at 2? because the first linmesh entry corresponds to the q=0,0,0
istart = 2 
istop  = input%qdiv
if (trim(input%genm) .eq. "enm") then
  itot=3
  znq=3
else if (trim(input%genm) .eq. "gnm") then
  itot=1
  znq=3
end if

k=1
do i=istart,istop

  qvec(:)  = linmesh(i,:)
  qmag     = dsqrt(dot_product(qvec,qvec)) ! right
  
  CALL dynmatrx(input,unit_cell,kirchoff,qvec,kirchdyn,hessdyn)

  if (trim(input%genm) .eq. "enm") call mod_arpack_zndrv1(input,znq,hessdyn,zvals,zvecs) 
  if (trim(input%genm) .eq. "gnm") call mod_arpack_zndrv1(input,znq,kirchdyn,zvals,zvecs) 

  call zeig_to_THz(zvals,vals)

  do j=1,itot
  ! 2pi brings it back to radfreq and 100 is for the 10^12 terahertz vs 10^10 angstrom
  ! final value m/s
    spdsound(k,j) = 100.0d0*two*pi*vals(j)/qmag
    sumsnd(j)     = sumsnd(j) + spdsound(k,j) 
  end do

  k=k+1

end do

avgspd=sumsnd/real(k-1,kind=dp)

end subroutine

subroutine zeig_to_THz(eig,frq)
! convert eigenvalues (square of circ freq) to freqs in THz 
use mod_constants, only : CONVFRQ, CLIGHT
complex(kind=8),        allocatable, intent(in)   :: eig(:)
real(dp),               allocatable, intent(out)  :: frq(:)
integer :: i,ialloc

if (allocated(frq)) deallocate (frq)
allocate(frq(size(eig)),stat=ialloc)

do i=1,size(eig)
  if(real(eig(i)) .lt. 1.0D-07) then
    frq(i)=zero
  else
    frq(i)=dsqrt(real(eig(i)))*CONVFRQ*CLIGHT*1.0D2/1.0D12 
    ! CLIGHT is in m/s 1.0D2 converts CLIGHT to cm/s
    ! this leaves hertz...  /1.0D12 converts to terahertz 
  end if
end do

end subroutine

subroutine d_zvcov_split(temp,zvals,zvecs,vcov,vcov_acoustic,vcov_optic)
use mod_constants, only: Rkcal
! take in the eigenvalues: units rad^2 kcal mol^-1 amu^-1 angs^2
! 
real(dp),                    intent(in)     :: temp
complex(kind=8), allocatable,intent(in)     :: zvals(:),zvecs(:,:)
real(dp),        allocatable,intent(in out) :: vcov(:,:), &
                                          vcov_acoustic(:,:),vcov_optic(:,:)
complex(kind=8) :: tmpval
real(dp)                                    :: kT
integer ::   nvals, ndim, nzeros
integer ::   i,j,k

kT = temp*Rkcal
nvals = size(zvals)
ndim  = size(zvecs(:,1))

do i=1,nvals      
  if (abs(real(zvals(i))).gt.1.0d-08) then
    do j=1,ndim
      do k=1,ndim
        tmpval    = zvecs(k,i)*CONJG(zvecs(j,i))
        vcov(k,j) = vcov(k,j) + real(tmpval)*kT/real(zvals(i))
      end do
    end do
!  print '(A4,4F10.4)','shit', vcov(i,i),real(tmpval),real(zvals(i)) , &
!          4.0d0*(pi**2)*kt/(real(zvals(i)))
  else
    nzeros = nzeros + 1
  end if
end do

do i=4,nvals      
  if (abs(real(zvals(i))).gt.1.0d-08) then
    do j=1,ndim
      do k=1,ndim
        tmpval    = zvecs(k,i)*CONJG(zvecs(j,i))
        vcov_optic(k,j) = vcov_optic(k,j) + real(tmpval)*kT/real(zvals(i))
      end do
    end do
  else
    nzeros = nzeros + 1
  end if
end do

do i=1,3      
  if (abs(real(zvals(i))).gt.1.0d-08) then
    do j=1,ndim
      do k=1,ndim
        tmpval    = zvecs(k,i)*CONJG(zvecs(j,i))
        vcov_acoustic(k,j) = vcov_acoustic(k,j) + real(tmpval)*kT/real(zvals(i))
      end do
    end do
  else
    nzeros = nzeros + 1
  end if
end do

end subroutine

subroutine d_zvcov(input,dvals,zvecs,vcov)
use mod_constants, only: Rkcal
! take in the eigenvalues: units rad^2 kcal mol^-1 amu^-1 angs^2
! 
type(inp_par),               intent(IN)     :: input
real(dp),        allocatable,intent(in)     :: dvals(:)
complex(kind=8), allocatable,intent(in)     :: zvecs(:,:)
real(dp),        allocatable,intent(in out) :: vcov(:,:)
real(dp)                                    :: temp
complex(kind=8)                             :: tmpval
real(dp)                                    :: kT
integer ::   nvals, ndim, nzeros
integer ::   i,j,k

temp = input%temperature
kT = temp*Rkcal
nvals = size(dvals)
ndim  = size(zvecs(:,1))

! this needs to be fixed up like dvcov from mod_hessian

do i=input%first_mode,input%last_mode      
  if (abs(dvals(i)).gt.1.0d-08) then
    do j=1,ndim
      do k=1,ndim
        tmpval    = zvecs(k,i)*CONJG(zvecs(j,i))
        vcov(k,j) = vcov(k,j) + dble(tmpval)*kT/(dvals(i))
      end do
    end do
  else
    nzeros = nzeros + 1
  end if
end do


end subroutine

subroutine z_zvcov(temp,zvals,zvecs,zvcov)
use mod_constants, only: Rkcal
! compute complex vcov matrix from complex eigenvecs
! 
real(dp),                    intent(in)     :: temp
complex(kind=8), allocatable,intent(in)     :: zvals(:),zvecs(:,:)
complex(kind=8), allocatable,intent(in out) :: zvcov(:,:)
complex(kind=8) :: tmpval
real(dp)                                    :: kT
integer ::   nvals, ndim, nzeros
integer ::   i,j,k

kT = temp*Rkcal
nvals = size(zvals)
ndim  = size(zvecs(:,1))

do i=1,nvals      
  if (abs(real(zvals(i))).gt.1.0d-08) then
    do j=1,ndim
      do k=1,ndim
        tmpval    = zvecs(k,i)*CONJG(zvecs(j,i))
        zvcov(k,j) = zvcov(k,j) + tmpval*kT/real(zvals(i))
      end do
    end do
!  print '(A4,4F10.4)','shit', vcov(i,i),real(tmpval),real(zvals(i)) , &
!          4.0d0*(pi**2)*kt/(real(zvals(i)))
  else
    nzeros = nzeros + 1
  end if
end do


end subroutine

subroutine wrt_dynm(rcut,dynmat)
complex(kind=8),allocatable ,intent(in) :: dynmat(:,:)
real(dp), intent(in) :: rcut
integer :: degfred,ii,jj,ier

degfred=size(dynmat(1,:))

! write out nonzero entries to kirchoff.txt
open(unit=13,file=trim(fileroot)//"-cdynmat.txt", &
     status='unknown',action='write', &
     iostat=ier)

write(13,'(I6,F10.2)')degfred,rcut

do ii = 1,degfred
  do jj = 1,degfred
    if (( abs(real(dynmat(ii,jj))) .gt.1.0d-06 ) .or. &
        ( abs(aimag(dynmat(ii,jj))) .gt.1.0d-06 )) then
      write(13,*) ii,jj,real(dynmat(ii,jj)), aimag(dynmat(ii,jj))
    end if
  end do
end do

close(13)

end subroutine

subroutine bnm_vcov(input,cell,asyms,unit_cell,kirchoff,hessian,vcov)
use mod_types,     only : sparse,inp_par,asym_list_type,dos
use mod_hessian,   only : dvcov, fullhess_massweigh
use mod_inout,     only : full_to_sparse,matrx2d_write
use cfml_crystal_metrics, only : Crystal_Cell_Type
use mod_bnm !,       only : block,vcov_proj,hess_proj,blocks_setup,bnm_eigenal_valvec,vect_proj
use mod_inout,     only : vect_atshrink
use mod_crysbuild, only : atom_shrinker
type(inp_par),        intent(in out) :: input
type (Crystal_Cell_Type),intent(in) :: cell
!type(atom_list_type), intent(in out) :: atoms
type(atom_list_type)                 :: subatoms
type(asym_list_type), intent(in out) :: unit_cell,asyms
type (sparse)       , intent(in)     :: kirchoff,hessian
real(dp), allocatable , intent(out)  :: vcov(:,:)
real(dp), allocatable                :: vcov_tmp(:,:),vcov3(:,:)
type (sparse) :: bnmhess,bnmkirch
real(dp), allocatable                :: vals(:),vecs(:,:),atvecs(:,:),tmp_atvecs(:,:)
real(dp), allocatable                :: bigproj(:,:)
type(block),allocatable              :: blocks(:),subblocks(:)
integer :: ialloc,nzeros,nblocks,i,j,natoms
type(sparse) :: bnmhess_coor
type(dos) :: states

!call atom_shrinker(input,atoms,subatoms)

!call blocks_setup (input,atoms,blocks)
if(trim(input%tpbc) .eq. "pbc" .or. trim(input%tpbc) .eq. "bvk") then
  call atom_shrinker(input,unit_cell%aunit(1),subatoms)
  call blocks_setup (input,unit_cell,blocks)
  call blocks_subblocker(input,unit_cell%aunit(1),blocks,subblocks)
  natoms = unit_cell%aunit(1)%natoms
else if(trim(input%tpbc) .eq. "asy") then
  call atom_shrinker(input,asyms%aunit(1),subatoms)
  call blocks_setup (input,asyms,blocks)
  call blocks_subblocker(input,asyms%aunit(1),blocks,subblocks)
  natoms = asyms%aunit(1)%natoms
else
  call atom_shrinker(input,asyms%aunit(1),subatoms)
  call blocks_setup (input,asyms%aunit(1),blocks)
  call blocks_subblocker(input,asyms%aunit(1),blocks,subblocks)
  natoms = asyms%aunit(1)%natoms
end if

if (6*size(blocks) .lt. input%last_mode) input%last_mode = 6*size(blocks)


if (trim(input%tpbc) .ne. "bvk") then
  call hess_proj(hessian,blocks,bnmhess)
  call bnm_eigenal_valvec(input,bnmhess,vals,vecs)

! call in the atshrink again and double check with something that's not just CA
!  call vect_proj(blocks,vecs,tmp_atvecs)
!  call vect_atshrink(input,atoms,tmp_atvecs,atvecs)
!  !allocate(vcov(size(atvecs(:,1)),size(atvecs(:,1))),stat=ialloc)
!  allocate(vcov_tmp(size(tmp_atvecs(:,1)),size(tmp_atvecs(:,1))),stat=ialloc)
!  allocate(vcov(size(atvecs(:,1)),size(atvecs(:,1))),stat=ialloc)
!  vcov=zero
!  vcov_tmp=zero
!  call dvcov(input,vals,tmp_atvecs,vcov_tmp,nzeros)
!  call dvcov(input,vals,atvecs,vcov,nzeros)
!  deallocate(vcov_tmp)
!  deallocate(vcov)

  allocate(vcov(3*subatoms%natoms,3*subatoms%natoms),stat=ialloc)
  vcov=zero
  call dvcov(input,vals,vecs,vcov_tmp,nzeros)
  call vcov_proj(vcov_tmp,subblocks,vcov)

else

  call bnm_qdispers(input,cell,blocks,unit_cell,kirchoff,hessian,vcov_tmp,states)
  call vcov_proj(vcov_tmp,subblocks,vcov)
  
  deallocate(vcov_tmp)
end if

! mass weighting
call fullhess_massweigh(input,unit_cell%aunit(1),vcov)

end subroutine

subroutine bnm_blockvcov(input,cell,blocks,unit_cell,kirchoff,hessian,vcov)
! DMR added Oct 1, 2009
use mod_types,     only : sparse,inp_par,asym_list_type,dos
use mod_hessian,   only : dvcov
use cfml_crystal_metrics, only : Crystal_Cell_Type
use mod_bnm ,       only : block,hess_proj,bnm_eigenal_valvec
type(inp_par),            intent(in out)  :: input
type (Crystal_Cell_Type), intent(in)      :: cell
type(block),allocatable , intent(in)      :: blocks(:)
type(asym_list_type)    , intent(in)      :: unit_cell ! for bvk
type (sparse)          ,  intent(in)      :: kirchoff,hessian
real(dp), allocatable ,   intent(out)     :: vcov(:,:)
type (sparse)                             :: bnmhess
real(dp), allocatable                     :: vals(:),vecs(:,:)
integer :: ialloc,nzeros,nblocks,i,j,natoms
type(dos) :: states

nblocks = size(blocks)

if (6*nblocks .lt. input%last_mode) input%last_mode = 6*nblocks

if (trim(input%tpbc) .ne. "bvk") then
  call hess_proj(hessian,blocks,bnmhess)
  call bnm_eigenal_valvec(input,bnmhess,vals,vecs)
  print *,"first mode: ", input%first_mode, " last mode: ", input%last_mode
  call dvcov(input,vals,vecs,vcov,nzeros)

else

  call bnm_qdispers(input,cell,blocks,unit_cell,kirchoff,hessian,vcov,states)

end if

print *, 'exiting bnm_blockvcov'
end subroutine

subroutine bnm_bigvcov(input,cell,asyms,unit_cell,kirchoff,hessian,vcov,states,blocks_out)
use mod_types,     only : sparse,inp_par,asym_list_type,dos
use mod_hessian,   only : dvcov, fullhess_massweigh,dos_calculator
use mod_inout,     only : full_to_sparse,matrx2d_write
use cfml_crystal_metrics, only : Crystal_Cell_Type
use mod_bnm ,       only : block,vcov_proj,hess_proj,blocks_setup,bnm_eigenal_valvec,vect_proj
use mod_inout,     only : vect_atshrink,sparse_deinit
use mod_crysbuild, only : atom_shrinker
type(inp_par),        intent(in out) :: input
type (Crystal_Cell_Type),intent(in) :: cell
!type(atom_list_type), intent(in out) :: atoms
type(asym_list_type), intent(in out) :: asyms, unit_cell
type (sparse)       , intent(in out)     :: kirchoff,hessian
real(dp), allocatable , intent(out)  :: vcov(:,:)
type(block),allocatable,optional,intent(out)     :: blocks_out(:)
type(dos), intent(out) :: states
real(dp), allocatable                :: vcov_tmp(:,:),vcov3(:,:)
type (sparse) :: bnmhess,bnmkirch
real(dp), allocatable                :: vals(:),vecs(:,:),atvecs(:,:),tmp_atvecs(:,:)
real(dp), allocatable                :: bigproj(:,:)
type(block),allocatable              :: blocks(:)
integer :: ialloc,nzeros,nblocks,i,j,natoms
type(sparse) :: bnmhess_coor


if(trim(input%tpbc) .eq. "pbc" .or. trim(input%tpbc) .eq. "bvk") then
  call blocks_setup (input,unit_cell,blocks)
! dumb way,can we use pointer?
  if(present(blocks_out)) call blocks_setup (input,unit_cell,blocks_out)
  natoms = unit_cell%aunit(1)%natoms
else if(trim(input%tpbc) .eq. "asy") then
  call blocks_setup (input,asyms,blocks)
  if(present(blocks_out)) call blocks_setup (input,asyms,blocks_out)
  natoms = asyms%aunit(1)%natoms
else
  call blocks_setup (input,asyms%aunit(1),blocks)
  if(present(blocks_out)) call blocks_setup (input,asyms%aunit(1),blocks_out)
  natoms = asyms%aunit(1)%natoms
end if

if (6*size(blocks) .lt. input%last_mode) input%last_mode = 6*size(blocks)

if (trim(input%tpbc) .ne. "bvk") then
  call hess_proj(hessian,blocks,bnmhess)
  call sparse_deinit(hessian) 
  call sparse_deinit(kirchoff)
  call bnm_eigenal_valvec(input,bnmhess,vals,vecs)
  call sparse_deinit(bnmhess)
  print *,"shit first mode: ", input%first_mode
  print *,"shit last mode: ", input%last_mode

  allocate(vcov(3*natoms,3*natoms),stat=ialloc)
  vcov=zero
  call dvcov(input,vals,vecs,vcov_tmp,nzeros)
  deallocate(vecs)
  call vcov_proj(vcov_tmp,blocks,vcov)

! density of states
  if (trim(input%eigtyp) .eq. "lapack") then
    call dos_calculator(75,vals,states)
  else
    call dos_calculator(input%nfrqs/2,vals,states)
  end if

else

  call bnm_qdispers(input,cell,blocks,unit_cell,kirchoff,hessian,vcov_tmp,states)
  call vcov_proj(vcov_tmp,blocks,vcov)
  deallocate(vcov_tmp)

end if

call fullhess_massweigh(input,unit_cell%aunit(1),vcov)

end subroutine

subroutine qvec_fqvec(input,cell)
  ! returns qvec times the cartesian 2 fractional transform matrix
type (inp_par) , intent(inout)  :: input
type (Crystal_Cell_Type), intent(in)        :: Cell
input%qvec(:)=two*pi*matmul(input%fqvec,cell%orth_cr_cel) ! right
end subroutine

end module

