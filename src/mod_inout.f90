module mod_inout
!
! DMR-01-11-2008
!
! while it is mod_inout, it's mostly just out.  writes pdbs of single snapshot and modes
!
! pdbwriter (atoms_in,filpdb)
! pdbwriter_cryst (images,filpdb)
! pdbanimate(atoms_in,dispvec,freq,nframe)

! crysFML modules
use cfml_atom_typedef,              only: Atom_List_Type
use cfml_globaldeps,                 only: sp,dp
! external modules
use mod_types,                only: inp_par, uc_list_type,asym_list_type,Neigh_List_Type
use mod_constants,            only: zero,one,two,pi

implicit none
character(len=25)     :: fileroot

!public :: crys_movie, eigen_getter,pdbwriter_new, pdbwriter_cr,&
!          pdbwriter_uc_crys,fileroot_set,psfwriter,pdbwriter_uc, pdbwriter_cryst, pdbanimate

interface pdbwriter_new
  module procedure pdbwriter_atom, pdbwriter_asym, pdbwriter_uc
end interface pdbwriter_new

interface pdbwriter_crys
  module procedure pdbwriter_asym_crys, pdbwriter_uc_crys
end interface pdbwriter_crys

interface xyzwriter_new
  module procedure xyzwriter_atom
end interface xyzwriter_new


contains 


logical function iso_bc(input)
type(inp_par), intent(in) :: input

if ((trim(input%tpbc) .ne. "pbc") .and. &
    (trim(input%tpbc) .ne. "asy") .and. &
    (trim(input%tpbc) .ne. "bvk")) then
  iso_bc = .true.
else
  iso_bc = .false.
end if

end function


logical function oned(input)
type(inp_par), intent(in) :: input
oned = .false.
if ( (trim(input%genm) .eq. "gnm") .or. &
     (trim(input%genm) .eq. "isocorr")) then
  oned = .true.
end if
end function

logical function threed(input)
type(inp_par), intent(in) :: input

threed = .true.
if (oned(input)) threed = .false.

end function

subroutine aniso_atshrink(input,atoms,aniso_in,aniso_out)
! take big aniso and select out parts that we want
type(inp_par),intent(in)           :: input
type(atom_list_type),intent(in)    :: atoms
real(dp), allocatable, intent(in)  :: aniso_in(:,:)
real(dp), allocatable, intent(out) :: aniso_out(:,:)
integer :: nat, i,j,ii,jj,h,hh,k,kk,l,m,ialloc,ndim,mult
! just use one atom type for now


if(trim(input%atom_anl) .eq. "all") then
  allocate (aniso_out(size(aniso_in,1),size(aniso_in,2)), stat=ialloc)
  if (ialloc /= 0) stop "aniso_atshrink> malloc"
  aniso_out = aniso_in
  return
else
  nat = 0
  do i = 1, atoms%natoms
    if(trim(atoms%atom(i)%lab) .eq. trim(input%atom_anl)) then
      nat = nat+1
    end if
  end do


  if (oned(input)) then 
    mult = 1
  else
    mult = 3
  end if
  ndim = mult*nat

!print *, "nat", nat, trim(input%atom_anl),"ndim",ndim

  allocate (aniso_out(ndim,3), stat=ialloc)
  if (ialloc /= 0) stop "aniso_atshrink> malloc"

  aniso_out=zero

  h  = 0
  do i = 1, atoms%natoms
    if (trim(input%atom_anl) .eq. trim(atoms%atom(i)%lab)) then
      h = h + 1
      ii  = mult*(i-1)+1 
      hh  = mult*(h-1)+1 
      aniso_out(hh:hh+2,:) = aniso_in(ii:ii+2,:)
    end if
  end do
end if

!print *, 'exiting aniso_atshrink', atoms%natoms, ndim

end subroutine

subroutine vcov_atshrink(input,atoms,vcov_in,vcov_out)
! take big vcov and select out parts that we want
type(inp_par),intent(in)           :: input
type(atom_list_type),intent(in)    :: atoms
real(dp), allocatable, intent(in)  :: vcov_in(:,:)
real(dp), allocatable, intent(out) :: vcov_out(:,:)
integer :: nat, i,j,ii,jj,h,hh,k,kk,l,m,ialloc,ndim,mult
! just use one atom type for now

if(trim(input%atom_anl) .eq. "all") then
  allocate (vcov_out(size(vcov_in,1),size(vcov_in,2)), stat=ialloc)
  if (ialloc /= 0) stop "vcov_atshrink> malloc"
  vcov_out = vcov_in
  return
else

  nat = 0
  do i = 1, atoms%natoms
    if(trim(atoms%atom(i)%lab) .eq. trim(input%atom_anl)) then
      nat = nat+1
    end if
  end do


  if (oned(input)) then 
    mult = 1
  else
    mult = 3
  end if
    ndim = mult*nat

  print *, "nat", nat, trim(input%atom_anl),"ndim",ndim

  allocate (vcov_out(ndim,ndim), stat=ialloc)
  if (ialloc /= 0) stop "vcov_atshrink> malloc"

  vcov_out=zero

  h  = 0
  do i = 1, atoms%natoms
    if(trim(atoms%atom(i)%lab) .eq. trim(input%atom_anl)) then
      h = h + 1
      ii  = mult*(i-1)+1 
      hh  = mult*(h-1)+1 
      k = 0
      do j = 1, atoms%natoms
        if(trim(atoms%atom(j)%lab) .eq. trim(input%atom_anl)) then
        k = k + 1
        jj = mult*(j-1)+1  
        kk = mult*(k-1)+1
!      print *, 'CACA', h,k,hh,kk

        do l = 1,mult 
        do m = 1,mult 
          vcov_out(hh+l-1,kk+m-1) = vcov_in(ii+l-1,jj+m-1)
        end do
        end do
        end if
      end do
    end if
  end do
end if

print *, 'exiting vcov_atshrink', atoms%natoms, ndim

end subroutine

subroutine isovcov_atshrink(input,atoms,vcov_in,vcov_out)
! take big vcov and select out parts that we want
type(inp_par),intent(in)           :: input
type(atom_list_type),intent(in)    :: atoms
real(dp), allocatable, intent(in)  :: vcov_in(:,:)
real(dp), allocatable, intent(out) :: vcov_out(:,:)
integer :: nat, i,j,ii,jj,h,hh,k,kk,l,m,ialloc,ndim,mult
! just use one atom type for now

if(trim(input%atom_anl) .eq. "all") then
  allocate (vcov_out(size(vcov_in,1),size(vcov_in,2)), stat=ialloc)
  if (ialloc /= 0) stop "vcov_atshrink> malloc"
  vcov_out = vcov_in
  return
else

  nat = 0
  do i = 1, atoms%natoms
    if(trim(atoms%atom(i)%lab) .eq. trim(input%atom_anl)) then
      nat = nat+1
    end if
  end do
  
  ndim = nat
   
  allocate (vcov_out(ndim,ndim), stat=ialloc)
  if (ialloc /= 0) stop "vcov_atshrink> malloc"

  vcov_out=zero

  h  = 0
  do i = 1, atoms%natoms
    if(trim(atoms%atom(i)%lab) .eq. trim(input%atom_anl)) then
      h = h + 1
      k = 0
      do j = 1, atoms%natoms
        if(trim(atoms%atom(j)%lab) .eq. trim(input%atom_anl)) then
          k = k + 1
          vcov_out(h,k) = vcov_in(i,j)
        end if
      end do
    end if
  end do
end if

print *, 'exiting vcov_atshrink', atoms%natoms, ndim

end subroutine

subroutine vect_atshrink(input,atoms,tmp_atvecs,atvecs)
type(inp_par),intent(in)           :: input
type(atom_list_type),intent(in)    :: atoms
real(dp), allocatable, intent(in)  :: tmp_atvecs(:,:)
real(dp), allocatable, intent(out) :: atvecs(:,:)
integer :: nat, i,j,ii,jj,k,ll,ialloc,ndim,mult
! just use one atom type for now

nat = 0
do i = 1, atoms%natoms
  if(trim(atoms%atom(i)%lab) .eq. trim(input%atom_anl)) then
    nat = nat+1
  end if
end do


if (oned(input)) then 
  mult = 1
else
  mult = 3
end if
  ndim = mult*nat

print *, "nat", nat, trim(input%atom_anl),"ndim",ndim

allocate (atvecs(ndim,size(tmp_atvecs,2)), stat=ialloc)
if (ialloc /= 0) stop "vect_atshrink> malloc"

ii = 0
do i = 1, atoms%natoms
  j = mult*(i-1)+1  
  if(trim(atoms%atom(i)%lab) .eq. trim(input%atom_anl)) then
    ii = ii+1
    jj = mult*(ii-1)+1
    do k = 1, size(tmp_atvecs,2) ! copy all columns
      do ll = 1,mult 
       atvecs(jj+ll-1,k)   = tmp_atvecs(j+ll-1,k)
      end do  
    end do
  end if
end do

print *, 'exiting vect_atshrink', atoms%natoms, ndim

end subroutine

subroutine fileroot_set(flnm)
! make something quick to share root across multiple modules
character(len=25)     :: flnm

fileroot = flnm

end subroutine

subroutine eigen_getter(nqvec,dvec,freqs,qvec_cart)
! reads in requested complex eigenvectors stored in two files (imaginary, real) 
!    eigenvectors depend on qvector and are stored accordingly:
!    example line1:  "qvec 1 0.00 192 10  64 real   0.5  0.0  0.0"
!       1st qvector of magnitude 0.00
!       192 lines follow with 10 frequencies per row corresponding to xyz of 64 atoms 
!
!   so vector_getter(nqvec,nfrq,dvec) returns the dvec requested by nqvec and nfrq 
!

integer, intent(in)                      :: nqvec
complex(kind=8),intent(out), allocatable :: dvec(:,:)
real(dp),intent(out),        allocatable :: freqs(:)
integer                                  :: iqvec,tqvec,degfred,ifrq,natom
character*4                              :: cqvc,typd
real(dp)                                 :: qvec_mag, qvec_cart(3)
real(dp),allocatable                     :: rtmp(:),itmp(:),disper(:,:)
integer                                  :: ier,i,ii


open(unit=33,file="dispersion.txt", &
     status='unknown',action='read', &
     iostat=ier)
open(unit=34,file="rzvecs.txt", &
     status='unknown',action='read', &
     iostat=ier)
open(unit=35,file="izvecs.txt", &
     status='unknown',action='read', &
     iostat=ier)

read(34,*) cqvc, iqvec,tqvec, qvec_mag,degfred,ifrq,natom,typd,qvec_cart  
backspace(34)
if(allocated(dvec))deallocate(dvec)
if(allocated(rtmp))deallocate(rtmp)
if(allocated(itmp))deallocate(itmp)
if(allocated(disper))deallocate(disper)
if(allocated(freqs))deallocate(freqs)
allocate(dvec(degfred,ifrq),rtmp(ifrq),itmp(ifrq),freqs(ifrq),disper(tqvec,ifrq+1),stat=ier)

! load up the dispersion.txt

do i=1,tqvec
  read (33,*) disper(i,:)
end do

do i=2,ifrq+1
  freqs(i-1)=disper(nqvec,i)
end do

close(33)
deallocate(disper)

!111 read(35.200) cqvc, iqvec,tqvec, qvec_mag,degfred,ifrq,natom,typd,qvec_cart  
!    read(35,200) cqvc, iqvec,tqvec, qvec_mag,degfred,ifrq,natom,typd,qvec_cart  
111 read(34,*) cqvc, iqvec,tqvec, qvec_mag,degfred,ifrq,natom,typd,qvec_cart  
    read(35,*) cqvc, iqvec,tqvec, qvec_mag,degfred,ifrq,natom,typd,qvec_cart  
    do i=1,degfred
      read (34,300), rtmp(:)
      read (35,300), itmp(:)
      do ii=1,ifrq
        dvec(i,ii)=cmplx(rtmp(ii),itmp(ii),kind=8)
      end do
    end do
    print *, iqvec, nqvec
    if(iqvec .ne. nqvec) goto 111
   
    close(34)
    close(35)

200 format(A4,1x,I3,1x,I3,1x,F5.2,1x,I4,I3,I4,1x,A4,1x,3ES14.5)
300 format(10ES12.3)

end subroutine

subroutine psfwriter (kirchoff,atoms_in,filpdb)
use mod_types, only: sparse
type (sparse),           intent(in) :: kirchoff
type (Atom_list_Type),   intent(in) :: atoms_in
character(*)                        :: filpdb
integer                             :: i,j,ier

open(unit=11,file=filpdb//".psf", &
     status='unknown',action='write', &
     iostat=ier)

write (11, '(A8)')     "PSF CMAP" 
write (11,'(I8,1x,A7)') 6, "!NTITLE" 
write (11,'(A7,1x,A4,1x,A3)') "REMARKS",filpdb,"psf" 
write (11,*) ""
write (11,'(I8,1x,A7)') atoms_in%natoms, "!NATOM" 

    do i=1,atoms_in%natoms
      write(11,'(I8,1x,A4,1x,I4,1x,A2,3x,A2,3x,A2,3x,F11.6,1x,F13.4,I12)') &
        i,filpdb,i,atoms_in%atom(i)%lab,atoms_in%atom(i)%lab,atoms_in%atom(i)%lab, &
        0.0d0,12.0110d0, 0 
    end do

write (11,*) ""
write (11,'(I8,1x,A13)') kirchoff%nbonds, "!NBOND: bonds" 

i=1
  do while (i .lt. kirchoff%nnzero)
    j = 0
    do while ((j .lt. 4) .and. (i.le. kirchoff%nnzero))
      if (kirchoff%values(i) .lt. zero) then
        write (11,'(2I8,$)') kirchoff%rows(i),kirchoff%columns(i)  
        i = i + 1
        j = j + 1
      else  
        i = i + 1
      end if
    end do

        write (11,*) ""

  end do

write (11,*) ""
write (11,'(I8,1x,A15)') 0, "!NTHETA: angles"
write (11,*) ""
write (11,'(I8,1x,A16)') 0, "!NPHI: dihedrals"
write (11,*) ""
write (11,'(I8,1x,A18)') 0, "!NIMPHI: impropers"
write (11,*) ""
write (11,'(I8,1x,A21)') 0, "!NCRTERM: cross-terms"

close(11)

end subroutine

subroutine pdbwriter_atom (atoms_in,filpdb)
type (Atom_list_Type),   intent(in) :: atoms_in
character(*)                        :: filpdb
integer                             :: i,ier

open(unit=11,file=filpdb//".pdb", &
     status='unknown',action='write', &
     iostat=ier)

    do i=1,atoms_in%natoms
      if (len_trim(atoms_in%atom(i)%resn) .eq. 3) then
        write(11,'(A4,1x,I6,1x,A3,2x,A3,I6,2x,F9.2,2F8.2,1x,2F6.2)') "ATOM", i, &
                atoms_in%atom(i)%lab, atoms_in%atom(i)%resn,atoms_in%atom(i)%iblock , &
                atoms_in%atom(i)%x(:),atoms_in%atom(i)%occ,atoms_in%atom(i)%biso
      else
        write(11,'(A4,1x,I6,1x,A3,1x,A4,I6,2x,F9.2,2F8.2,1x,2F6.2)') "ATOM", i, &
                atoms_in%atom(i)%lab,atoms_in%atom(i)%resn,atoms_in%atom(i)%iblock , &
                atoms_in%atom(i)%x(:),atoms_in%atom(i)%occ,atoms_in%atom(i)%biso
      end if
    end do

close(11)

end subroutine

subroutine xyzwriter_atom (atoms_in,filpdb)
type (Atom_list_Type),   intent(in) :: atoms_in
character(*)                        :: filpdb
integer                             :: i,ier

open(unit=11,file=filpdb//".xyz", &
     status='unknown',action='write', &
     iostat=ier)

    do i=1,atoms_in%natoms
        write(11,'(A2,1x,3F10.3)') atoms_in%atom(i)%chemsymb,atoms_in%atom(i)%x(:)
    end do

close(11)

end subroutine

subroutine pdbwriter (atoms_in,filpdb)
type (Atom_list_Type),   intent(in) :: atoms_in
character(*)                        :: filpdb
integer                             :: i,ier

open(unit=11,file=filpdb//".pdb", &
     status='unknown',action='write', &
     iostat=ier)

    do i=1,atoms_in%natoms
      write(11,'(A4,1x,I6,1x,A3,2x,A3,I6,2x,F9.2,2F8.2)') "ATOM", i, &
                atoms_in%atom(i)%lab, atoms_in%atom(i)%lab, atoms_in%atom(i)%iblock , atoms_in%atom(i)%x(:)
    end do

close(11)

end subroutine

subroutine pdbwriter_uc (atoms_in,filpdb)
use mod_types, only: uc_type
type (uc_Type),   intent(in) :: atoms_in
character(*)                        :: filpdb
integer                             :: i,ier

open(unit=11,file=filpdb//".pdb", &
     status='unknown',action='write', &
     iostat=ier)

    do i=1,atoms_in%natoms
      write(11,'(A4,1x,I6,1x,A3,2x,A3,I6,2x,F9.3,2F9.3)') "ATOM", i, &
                atoms_in%atom(i)%lab, atoms_in%atom(i)%lab,atoms_in%atom(i)%iblock  , atoms_in%atom(i)%x(:)
    end do

close(11)

end subroutine

subroutine pdbwriter_asym (uc_asym,filpdb)
type (asym_List_Type),   intent(in) :: uc_asym
character(*)                         :: filpdb
integer                              :: i,j,nat,ier


open(unit=11,file=filpdb//".pdb", &
     status='unknown',action='write', &
     iostat=ier)
  nat=0
  do j=1, uc_asym%naunits
    do i=1,uc_asym%aunit(j)%natoms
      nat=nat+1

      if (len_trim(uc_asym%aunit(j)%atom(i)%resn) .eq. 3) then
        write(11,'(A4,1x,I6,1x,A3,2x,A3,I6,2x,F9.3,2F9.3,2F8.2,A5)') "ATOM", nat, &
                uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
                uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:), &
                uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso, &
                uc_asym%aunit(j)%atom(i)%ChemSymb
      else
        write(11,'(A4,1x,I6,1x,A3,1x,A4,I6,2x,F9.3,2F9.3,2F8.2,A5)') "ATOM", nat, &
                uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
                uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:),&
                uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso,&
                uc_asym%aunit(j)%atom(i)%ChemSymb
      end if

    end do
  end do

!  do i =1, size(uc_asym%block_cents(:,1))
!    nat = nat+1
!            write(11,'(A4,1x,I6,1x,A2,3x,A3,I6,2x,F9.2,2F8.2)') "ATOM", nat, &
!                "PB", "PBB",i, uc_asym%block_cents(i,:)
!  end do

close(11)

end subroutine

subroutine pdbwriter_asym_crys (uc_asym,filpdb)
type (asym_List_Type),   intent(in) :: uc_asym
character(*)                         :: filpdb
integer                              :: i,j,nat,ier
integer :: ii,jj,kk
real(dp), dimension(3) :: avec,bvec,cvec,vec

open(unit=11,file=filpdb//".pdb", &
     status='unknown',action='write', &
     iostat=ier)

! write primary unit cell
nat=0
do j=1, uc_asym%naunits
  do i=1,uc_asym%aunit(j)%natoms
    nat=nat+1
    if (len_trim(uc_asym%aunit(j)%atom(i)%resn) .eq. 3) then
      write(11,'(A4,1x,I6,2x,A3,1x,A3,I6,2x,F10.3,2F8.3,2F6.2,A13)') "ATOM", nat, &
          uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
          uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:)+vec, &
          uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso, &
          uc_asym%aunit(j)%atom(i)%ChemSymb
    else if (len_trim(uc_asym%aunit(j)%atom(i)%resn) .eq. 1) then
      write(11,'(A4,1x,I6,2x,A3,3x,A1,I6,2x,F10.3,2F8.3,2F6.2,A13)') "ATOM", nat, &
          uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
          uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:)+vec, &
          uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso, &
          uc_asym%aunit(j)%atom(i)%ChemSymb
    else if (len_trim(uc_asym%aunit(j)%atom(i)%resn) .eq. 2) then
      write(11,'(A4,1x,I6,2x,A3,2x,A2,I6,2x,F10.3,2F8.3,2F6.2,A13)') "ATOM", nat, &
          uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
          uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:)+vec, &
          uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso, &
          uc_asym%aunit(j)%atom(i)%ChemSymb
    else
      write(11,'(A4,1x,I6,2x,A3,A4,I6,2x,F10.3,2F8.3,2F6.2,A13)') "ATOM", nat, &
      !write(11,'(A4,1x,I6,2x,A3,1x,A4,I6,2x,F9.2,2F8.2,2F8.2,A5)') "ATOM", nat, &
          uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
          uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:)+vec,&
          uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso,&
          uc_asym%aunit(j)%atom(i)%ChemSymb
    end if
  end do
end do

do ii=uc_asym%aneigh(1),uc_asym%aneigh(2)
  avec = real(ii,kind=dp)*uc_asym%avec 
  do jj=uc_asym%bneigh(1),uc_asym%bneigh(2)
    bvec = real(jj,kind=dp)*uc_asym%bvec 
    do kk=uc_asym%cneigh(1),uc_asym%cneigh(2)
      cvec = real(kk,kind=dp)*uc_asym%cvec 
      vec = avec + bvec + cvec
      if ((ii .eq. 0) .and. (jj .eq. 0) .and.(kk .eq. 0)) cycle
      nat=0
      do j=1, uc_asym%naunits
        do i=1,uc_asym%aunit(j)%natoms
          nat=nat+1
    if (len_trim(uc_asym%aunit(j)%atom(i)%resn) .eq. 3) then
      write(11,'(A4,1x,I6,2x,A3,1x,A3,I6,2x,F10.3,2F8.3,2F6.2,A13)') "ATOM", nat, &
          uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
          uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:)+vec, &
          uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso, &
          uc_asym%aunit(j)%atom(i)%ChemSymb
    else if (len_trim(uc_asym%aunit(j)%atom(i)%resn) .eq. 1) then
      write(11,'(A4,1x,I6,2x,A3,3x,A1,I6,2x,F10.3,2F8.3,2F6.2,A13)') "ATOM", nat, &
          uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
          uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:)+vec, &
          uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso, &
          uc_asym%aunit(j)%atom(i)%ChemSymb
    else if (len_trim(uc_asym%aunit(j)%atom(i)%resn) .eq. 2) then
      write(11,'(A4,1x,I6,2x,A3,2x,A2,I6,2x,F10.3,2F8.3,2F6.2,A13)') "ATOM", nat, &
          uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
          uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:)+vec, &
          uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso, &
          uc_asym%aunit(j)%atom(i)%ChemSymb
    else
      write(11,'(A4,1x,I6,2x,A3,A4,I6,2x,F10.3,2F8.3,2F6.2,A13)') "ATOM", nat, &
      !write(11,'(A4,1x,I6,2x,A3,1x,A4,I6,2x,F9.2,2F8.2,2F8.2,A5)') "ATOM", nat, &
          uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
          uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:)+vec,&
          uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso,&
          uc_asym%aunit(j)%atom(i)%ChemSymb
      end if
!          if (len_trim(uc_asym%aunit(j)%atom(i)%resn) .eq. 3) then
!            write(11,'(A4,1x,I6,1x,A3,2x,A3,I6,2x,F9.2,2F8.2,2F8.2,A5)') "ATOM", nat, &
!                uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
!                uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:)+vec, &
!                uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso, &
!                uc_asym%aunit(j)%atom(i)%ChemSymb
!          else if (len_trim(uc_asym%aunit(j)%atom(i)%resn) .eq. 1) then
!            write(11,'(A4,1x,I6,1x,A3,4x,A1,I6,2x,F9.2,2F8.2,2F8.2,A5)') "ATOM", nat, &
!                uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
!                uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:)+vec, &
!                uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso, &
!                uc_asym%aunit(j)%atom(i)%ChemSymb
!          else if (len_trim(uc_asym%aunit(j)%atom(i)%resn) .eq. 2) then
!            write(11,'(A4,1x,I6,1x,A3,3x,A2,I6,2x,F9.2,2F8.2,2F8.2,A5)') "ATOM", nat, &
!                uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
!                uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:)+vec, &
!                uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso, &
!                uc_asym%aunit(j)%atom(i)%ChemSymb
!          else
!             write(11,'(A4,1x,I6,1x,A3,1x,A4,I6,2x,F9.2,2F8.2,2F8.2,A5)') "ATOM", nat, &
!                uc_asym%aunit(j)%atom(i)%lab, uc_asym%aunit(j)%atom(i)%resn, &
!                uc_asym%aunit(j)%atom(i)%iblock , uc_asym%aunit(j)%atom(i)%x(:)+vec,&
!                uc_asym%aunit(j)%atom(i)%occ,uc_asym%aunit(j)%atom(i)%biso,&
!                uc_asym%aunit(j)%atom(i)%ChemSymb
!          end if
        end do
      end do

    end do
  end do
end do  

close(11)

end subroutine

subroutine pdbwriter_uc_crys (unit_cell,filpdb)
use mod_types, only: uc_type
type (uc_Type),   intent(in) :: unit_cell
character(*)                         :: filpdb
integer                              :: i,j,nat,ier
integer :: ii,jj,kk
real(dp), dimension(3) :: avec,bvec,cvec,vec


open(unit=11,file=filpdb//".pdb", &
     status='unknown',action='write', &
     iostat=ier)

! write primary unit cell
nat=0
do j=1, unit_cell%natoms
  write(11,'(A4,1x,I6,1x,A3,2x,A3,I6,2x,F9.2,2F8.2)') "ATOM", j, &
            unit_cell%atom(j)%lab, unit_cell%atom(j)%resn, &
            j, unit_cell%atom(j)%x(:)
end do

do ii=unit_cell%aneigh(1),unit_cell%aneigh(2)
  avec = real(ii,kind=dp)*unit_cell%avec 
  do jj=unit_cell%bneigh(1),unit_cell%bneigh(2)
    bvec = real(jj,kind=dp)*unit_cell%bvec 
    do kk=unit_cell%cneigh(1),unit_cell%cneigh(2)
      cvec = real(kk,kind=dp)*unit_cell%cvec 
      vec = avec + bvec + cvec
      if ((ii .eq. 0) .and. (jj .eq. 0) .and.(kk .eq. 0)) cycle
      do j=1, unit_cell%natoms
        write(11,'(A4,1x,I6,1x,A3,2x,A3,I6,2x,F9.3,2F8.3)') "ATOM", j, &
            unit_cell%atom(j)%lab, unit_cell%atom(j)%resn, &
            j, unit_cell%atom(j)%x(:)+vec
      end do

    end do
  end do
end do  

close(11)

end subroutine

subroutine asym_to_atom (uc_asym,unit_cell)
use cfml_atom_typedef,              only: Atom_List_Type,init_atom_type
type (asym_List_Type),   intent(in)  :: uc_asym
type (Atom_list_Type),   intent(out) :: unit_cell
integer :: i,j,natoms,nat,ier

natoms = uc_asym%naunits*uc_asym%aunit(1)%natoms
unit_cell%natoms=natoms
allocate (unit_cell%atom(unit_cell%natoms),stat=ier)
print *, "number of atoms:", natoms
do i=1,natoms
   call init_atom_type(unit_cell%atom(i))
end do

nat = 0
do j=1, uc_asym%naunits
  do i=1,uc_asym%aunit(j)%natoms
    nat=nat+1
    unit_cell%atom(nat)=uc_asym%aunit(j)%atom(i)
    unit_cell%atom(nat)%x=uc_asym%aunit(j)%atom(i)%x
  end do
end do

end subroutine

subroutine pdbwriter_cr (uc_image,filpdb)
type (uc_list_type),   intent(in) :: uc_image
character(*)                         :: filpdb
integer                              :: i,j,k,nat,ier

open(unit=11,file=filpdb//".pdb", &
     status='unknown',action='write', &
     iostat=ier)
  nat=0
do k=0,uc_image%nneighs
  do j=1, uc_image%neigh(k)%naunits
    do i=1,uc_image%neigh(k)%aunit(j)%natoms
      nat=nat+1
!      print '(A4,1x,I6,1x,A3,2x,A3,I6,2x,F9.2,2F8.2)', "ATOM", nat, &
!                uc_image%neigh(i)%aunit(j)%atom(i)%lab,  &
!                uc_image%neigh(i)%aunit(j)%atom(i)%lab,  &
!                1, uc_image%neigh(i)%aunit(j)%atom(i)%x(:)
     write(11,'(A4,1x,I6,1x,A3,2x,A3,I6,2x,F9.2,2F8.2)') "ATOM", nat, &
                uc_image%neigh(k)%aunit(j)%atom(i)%lab,  &
                uc_image%neigh(k)%aunit(j)%atom(i)%lab,  &
                j, uc_image%neigh(k)%aunit(j)%atom(i)%x(:)
    end do
  end do
end do

close(11)

end subroutine

subroutine pdbwriter_cryst (images,filpdb)
type (Neigh_List_Type),   intent(in) :: images
character(*)                        :: filpdb
integer                             :: i,j,nat,ier


open(unit=11,file=filpdb//".pdb", &
     status='unknown',action='write', &
     iostat=ier)
  nat=0
  do j=0,images%nneighs
    do i=1,images%neigh(j)%natoms
      nat=nat+1
      write(11,'(A4,1x,I6,1x,A3,2x,A3,I2,I5.2x,F10.3,2F8.3,2F6.2,1x,A2)') "ATOM", nat, &
                images%neigh(j)%atom(i)%lab, images%neigh(j)%atom(i)%lab, &
                j+1,nat, images%neigh(j)%atom(i)%x(:),images%neigh(0)%atom(i)%occ, &
                images%neigh(0)%atom(i)%biso,images%neigh(0)%atom(i)%ChemSymb
!      write(11,'(A4,1x,I6,1x,A3,2x,A3,I6,2x,F9.2,2F8.2)') "ATOM", nat, &
!                images%neigh(j)%atom(i)%lab, images%neigh(j)%atom(i)%lab, &
!                1, images%neigh(j)%atom(i)%x(:)
!
    end do
  end do

close(11)
  
end subroutine

subroutine real_animate(input,asyms,unit_cell,vals,vecs)
use mod_types, only: inp_par 
type (inp_par)   ,        intent(in)    :: input
type (asym_list_type),    intent(in)    :: asyms,unit_cell
real(dp) , allocatable,   intent(in)    :: vals(:),vecs(:,:)
character(len=10) :: file_digi,file_name
real(dp) , allocatable                  :: tmpvec(:)
real(dp)                                :: tmpval
integer                  :: nframes
integer :: i,j,ialloc
file_digi = '1234567890'

nframes = input%animate_nframes
allocate(tmpvec(size(vecs(:,1))),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"
tmpvec = zero

do i = input%first_mode,input%last_mode
  file_name = "-vib"//file_digi(i:i)
  tmpvec(:) = vecs(:,i)
  tmpval    = vals(i)
  if (trim(input%tpbc) .eq. "pbc") then
    call pdbanimate(input,unit_cell%aunit(1),tmpvec,tmpval,nframes,trim(file_name)) 
  else
    print *, 'shit',input%first_mode,input%last_mode
    call pdbanimate(input,asyms%aunit(1),& 
                          tmpvec,tmpval,nframes,trim(file_name)) 
    !call pdbanimate(input,asym_un,& 
  end if

end do

end subroutine

subroutine complex_animate(input,asyms,unit_cell,vals,vecs)
use mod_types, only: inp_par 
type (inp_par)   ,        intent(in)    :: input
type (asym_list_type),    intent(in)    :: asyms,unit_cell
real(dp) , allocatable,   intent(in)    :: vals(:)
complex(kind=8) , allocatable,   intent(in)    :: vecs(:,:)
character(len=10) :: file_digi,file_name
complex(kind=8) , allocatable                  :: tmpvec(:)
real(dp)                                :: tmpval
integer :: i,j,ialloc

file_digi = '1234567890'

allocate(tmpvec(size(vecs(:,1))),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"
tmpvec = zero

do i = input%first_mode,input%last_mode
  file_name = "-vib"//file_digi(i:i)
  tmpvec(:) = vecs(:,i)
  tmpval    = vals(i)
  if (trim(input%tpbc) .eq. "bvk") then
    print *, 'shit', tmpval
    call crys_movie(input,unit_Cell,tmpvec,tmpval,trim(file_name)) 
  end if

end do

end subroutine

subroutine pdbanimate(input,atoms_in,dispvec,freq,nframe,filpdb)
! todo: get freqs in
use mod_types, only: inp_par 
use cfml_scattering_chemical_tables, only: get_atomic_mass
type(inp_par)        , intent(in) :: input            
type (Atom_list_Type), intent(in) :: atoms_in
real(dp), allocatable, intent(in) :: dispvec(:)          
real(dp),              intent(in) :: freq
integer,               intent(in) :: nframe
character(*),          intent(in) :: filpdb
integer :: i,j,k,cnt,ier
real(dp) :: dx,dy,dz,x,y,z
real(sp) :: sqrmass,mass

open(unit=22,file=trim(input%fileroot)//"-"//trim(input%tpbc)//filpdb//".pdb", &
     status='unknown',action='write', &
     iostat=ier)

do k=0,nframe-1
  j=1
  write(22,'(A5,I15)') "MODEL ", k
  do i=1,atoms_in%natoms
    if (trim(input%weigh) .eq. "mass") then
      call Get_Atomic_Mass(atoms_in%atom(i)%chemsymb,mass)
      sqrmass = sqrt(mass)
    else
      sqrmass = 1.0d0
    end if
    dx=input%animate_mag*dispvec(j)  *dcos(two*pi*k/nframe)/sqrmass
    dy=input%animate_mag*dispvec(j+1)*dcos(two*pi*k/nframe)/sqrmass
    dz=input%animate_mag*dispvec(j+2)*dcos(two*pi*k/nframe)/sqrmass
    x=atoms_in%atom(i)%x(1)+dx
    y=atoms_in%atom(i)%x(2)+dy
    z=atoms_in%atom(i)%x(3)+dz
    write(22,'(A4,1x,I6,1x,A3,2x,A3,I6,2x,F9.2,2F8.2)') "ATOM", i, &
                atoms_in%atom(i)%lab, atoms_in%atom(i)%lab, atoms_in%atom(i)%iblock, x,y,z
    j=j+3
  end do
  write(22,'(A4)') "ENDMDL "
end do

!300  format(A4,3x,I4,1x,A4,2x,A4,I5.2x, 3F9.3)

end subroutine

subroutine crys_movie(input,unit_cell,dispvec,freq,filpdb)
! todo: get freqs in
use mod_types, only: inp_par 
type(inp_par)        , intent(in) :: input            
type (asym_List_Type),intent(in) :: unit_cell
complex(kind=8),       intent(in) :: dispvec(:)
real(dp),              intent(in) :: freq
character(*),          intent(in) :: filpdb
real(dp)                          :: qvec(3)
integer :: i,j,k,iii,jjj,kkk,cnt,nat,ier
integer :: nframe
real(dp) :: dx,dy,dz,x,y,z,qdotr,sin_qdotr, cos_qdotr, tfrq, sin_tfrq, cos_tfrq
real(dp),dimension(3) :: vec,avec,bvec,cvec

print *, "crys_movie> nframes", input%animate_nframes
nframe = input%animate_nframes
qvec = input%qvec

!do i =1, size(dispvec)
!  if ( (abs(dble(dispvec(i))) .gt. 1.0d-08 ) .or. ( abs(dimag(dispvec(i))) .gt. 1.0d-08 ) ) then
!    print *,i, dispvec(i)
!  end if
!end do
!print *, 'hit'

open(unit=22,file=trim(input%fileroot)//"-"//trim(input%tpbc)//filpdb//".pdb", &
     status='unknown',action='write', &
     iostat=ier)

do k=0,nframe-1
  nat=0
  tfrq=two*pi*(real(k,kind=8)/real(nframe,kind=8))
  sin_tfrq = dsin(tfrq)
  cos_tfrq = dcos(tfrq)
  write(22,'(A5,I15)') "MODEL ", k
  do iii=unit_cell%aneigh(1),unit_cell%aneigh(2)
  avec = real(iii,kind=dp)*unit_cell%avec
  do jjj=unit_cell%bneigh(1),unit_cell%bneigh(2)
  bvec = real(jjj,kind=dp)*unit_cell%bvec
  do kkk=unit_cell%cneigh(1),unit_cell%cneigh(2)
  cvec = real(kkk,kind=dp)*unit_cell%cvec
  vec = avec + bvec + cvec
  do i=1,unit_cell%natoms
    nat=nat+1
! dot qvec and r. qvec must be orthogonal to r
      !qdotr=two*pi*dot_product(qvec,images%neigh(j)%atom(i)%x)
    qdotr=dot_product(qvec,unit_Cell%aunit(1)%atom(i)%x+vec)
    cos_qdotr = dcos(qdotr)
    sin_qdotr = dsin(qdotr)
      !print '(I5,3F10.4,3F10.4,F10.4)', i, qvec, images%neigh(j)%atom(i)%x,qdotr
! compute the displacements
    dx=input%animate_mag*dble(dispvec(i))* &
           (cos_qdotr*cos_tfrq + sin_qdotr*sin_tfrq) + &
            dimag(dispvec(i))  *(cos_qdotr*sin_tfrq - sin_qdotr*cos_tfrq)
    dy=input%animate_mag*dble(dispvec(i+1))* &
           (cos_qdotr*cos_tfrq + sin_qdotr*sin_tfrq) + &
            dimag(dispvec(i+1))*(cos_qdotr*sin_tfrq - sin_qdotr*cos_tfrq)
    dz=input%animate_mag*dble(dispvec(i+2))* &
           (cos_qdotr*cos_tfrq + sin_qdotr*sin_tfrq) + &
            dimag(dispvec(i+2))*(cos_qdotr*sin_tfrq - sin_qdotr*cos_tfrq)
     ! print '(2I5,3F10.4)', k,nat, dble(dispvec(i+1)),dimag(dispvec(i+1))
     ! print '(2I5,3F10.4)', k,nat, dx,dy,dz
     x=unit_Cell%aunit(1)%atom(i)%x(1)+dx + vec(1)
     y=unit_Cell%aunit(1)%atom(i)%x(2)+dy + vec(2)
     z=unit_Cell%aunit(1)%atom(i)%x(3)+dz + vec(3)

      write(22,'(A4,1x,I6,1x,A3,2x,A3,I6,2x,F9.2,2F8.2)') "ATOM", nat, &
                unit_Cell%aunit(1)%atom(i)%lab, unit_Cell%aunit(1)%atom(i)%lab, &
                unit_Cell%aunit(1)%atom(i)%iblock, x,y,z
  end do
  end do
  end do
  end do
  write(22,'(A4)') "ENDMDL "
end do

close(22)
  
end subroutine

subroutine sparse_refinit(refmat,sparsemat,typ)
! initializes sparse matrix from reference.  copies over row and column elements
use mod_types,       only:sparse
character(len=1)             :: typ ! d or z for allocating correct value array
type(sparse), intent(in)     :: refmat
integer :: ndim,nnzero
type(sparse), intent(out)     :: sparsemat
integer :: ialloc

ndim    = refmat%ndim
nnzero  = refmat%nnzero
sparsemat%nnzero= nnzero
sparsemat%ndim  = ndim
sparsemat%sparsity = refmat%sparsity

if (allocated(sparsemat%values))  deallocate (sparsemat%values)
if (allocated(sparsemat%zvalues)) deallocate (sparsemat%zvalues)
if (allocated(sparsemat%distance)) deallocate (sparsemat%distance)
if (allocated(sparsemat%tab_dxyz)) deallocate (sparsemat%tab_dxyz)
if (allocated(sparsemat%columns)) deallocate (sparsemat%columns)
if (allocated(sparsemat%rows))    deallocate (sparsemat%rows)
allocate( sparsemat%columns(nnzero), sparsemat%rows(nnzero),&
          sparsemat%distance(nnzero),sparsemat%tab_dxyz(nnzero,3),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

if (typ .eq. "z") then
  allocate(sparsemat%zvalues(nnzero),stat=ialloc)
  if (ialloc /= 0) STOP "*** not enough memory ***"
  sparsemat%zvalues = cmplx(zero,zero,kind=8)
else if (typ .eq. "d") then
  allocate(sparsemat%values(nnzero), stat=ialloc)
  if (ialloc /= 0) STOP "*** not enough memory ***"
  sparsemat%values=zero
end if

sparsemat%rows     = refmat%rows
sparsemat%columns  = refmat%columns
sparsemat%distance = refmat%distance
sparsemat%tab_dxyz = refmat%tab_dxyz

end subroutine

!subroutine full_to_sparse_rc(refmat,sparsemat)
! get sparse matrix from nxm matrix

!end subroutine

subroutine sparse_upper_to_gen(refmat,sparsemat,typ)
! convert sparse matrix stored upper to general stored sparse matrix
! complex assumed hermitian and double prec assumed symmetric
use mod_types,       only:sparse
character(len=1)             :: typ ! d or z for allocating correct value array
type(sparse), intent(in)     :: refmat
integer :: ndim,nnzero
type(sparse), intent(out)     :: sparsemat
integer :: i,j,ialloc

ndim               = refmat%ndim
nnzero             = 2*refmat%nnzero-ndim
sparsemat%nnzero   = nnzero
sparsemat%ndim     = ndim
sparsemat%rdim     = refmat%rdim
sparsemat%cdim     = refmat%cdim
sparsemat%sparsity = real(nnzero)/real(ndim*ndim)

! we're leaving out the distances and tab_dxyz for now
if (allocated(sparsemat%values))  deallocate (sparsemat%values)
if (allocated(sparsemat%zvalues)) deallocate (sparsemat%zvalues)
if (allocated(sparsemat%columns)) deallocate (sparsemat%columns)
if (allocated(sparsemat%rows))    deallocate (sparsemat%rows)
allocate( sparsemat%columns(nnzero), sparsemat%rows(nnzero),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

! switch from copying the whole thing as commented below as there
! seems to be some limit to the number of entries
do i = 1,refmat%nnzero
  sparsemat%rows(i)     = refmat%rows(i)
  sparsemat%columns(i)  = refmat%columns(i)
  !sparsemat%columns(1:refmat%nnzero)  = refmat%columns
end do


if (typ .eq. "z") then


  allocate(sparsemat%zvalues(nnzero),stat=ialloc)
  if (ialloc /= 0) STOP "*** not enough memory ***"


  !sparsemat%zvalues(1:refmat%nnzero)     = refmat%zvalues(:)

!fill in the rest
  j=1
  do i=1,refmat%nnzero
    sparsemat%zvalues(i) = refmat%zvalues(i)
    if (refmat%rows(i) .ne. refmat%columns(i)) then
      sparsemat%zvalues(refmat%nnzero + j) = CONJG(refmat%zvalues(i))
      sparsemat%rows(refmat%nnzero    + j) = refmat%columns(i)
      sparsemat%columns(refmat%nnzero + j) = refmat%rows(i)
      j=j+1
    end if
  end do

else if (typ .eq. "d") then

  allocate(sparsemat%values(nnzero), stat=ialloc)
  if (ialloc /= 0) STOP "*** not enough memory ***"

  sparsemat%values(1:refmat%nnzero)     = refmat%values

!fill in the rest
  j=1
  do i=1,refmat%nnzero
    if (refmat%rows(i) .ne. refmat%columns(i)) then
      sparsemat%zvalues(refmat%nnzero +j) = refmat%zvalues(i)
      sparsemat%rows(refmat%nnzero+j)     = refmat%columns(i)
      sparsemat%columns(refmat%nnzero+j)  = refmat%rows(i)
      j=j+1
    end if
  end do

end if

end subroutine

subroutine sparse_init(rdim,cdim,nnzero,sparsemat,typ)
! DMR comments added 12/8/2008
! ndim is the dimension of the original matrix
! nnzero is the number of nonzeros above the diagonal 
! typ defines whether it is real or complex.  if complex, we're talking hermitian, not symmetric
use mod_types,       only:sparse
character(len=1)             :: typ ! d or z for allocating correct value array
integer,intent(in)           :: rdim,cdim,nnzero
type(sparse), intent(out)    :: sparsemat
integer :: ialloc

if (rdim .eq. cdim) then
  sparsemat%ndim = rdim
else
  sparsemat%ndim = 0
end if

sparsemat%rdim = rdim
sparsemat%cdim = cdim

sparsemat%nnzero= nnzero

if (allocated(sparsemat%values))   deallocate (sparsemat%values)
if (allocated(sparsemat%zvalues))  deallocate (sparsemat%zvalues)
if (allocated(sparsemat%distance)) deallocate (sparsemat%distance)
if (allocated(sparsemat%tab_dxyz)) deallocate (sparsemat%tab_dxyz)
if (allocated(sparsemat%columns))  deallocate (sparsemat%columns)
if (allocated(sparsemat%rows))     deallocate (sparsemat%rows)
allocate( sparsemat%columns(nnzero), sparsemat%rows(nnzero),&
          sparsemat%distance(nnzero),sparsemat%tab_dxyz(nnzero,3),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

if (typ .eq. "z") then
  allocate(sparsemat%zvalues(nnzero),stat=ialloc)
  if (ialloc /= 0) STOP "*** not enough memory ***"
  sparsemat%zvalues = cmplx(zero,zero,kind=8)
else if (typ .eq. "d") then
  allocate(sparsemat%values(nnzero), stat=ialloc)
  if (ialloc /= 0) STOP "*** not enough memory ***"
  sparsemat%values=zero
end if

  sparsemat%distance = zero
  sparsemat%tab_dxyz = zero
!DMR fix sparsity calc.  12-08-08
if (rdim .eq. cdim) then
  sparsemat%sparsity = real(2*nnzero-rdim)/real(rdim*rdim)
else
  sparsemat%sparsity = real(nnzero)/real(rdim*cdim)
end if

end subroutine

subroutine sparse_statalloc(sparsemat)
! DMR 2-19-2010 size printer
use mod_types,       only:sparse
type(sparse), intent(in)    :: sparsemat

print *, "about to print stats for sparsemat"
print '(A4,1x,I5,A4,1x,I5,A4,1x,I5)', &
           "ndim", sparsemat%ndim, "rdim", sparsemat%rdim, "cdim", sparsemat%cdim

if (allocated(sparsemat%values))   print *, "values:",  size(sparsemat%values) 
if (allocated(sparsemat%zvalues))  print *, "zvalues:", size(sparsemat%zvalues) 
if (allocated(sparsemat%distance)) print *, "distance:", size(sparsemat%distance)
if (allocated(sparsemat%tab_dxyz)) print *, "tab_dxyz:", size(sparsemat%tab_dxyz)
if (allocated(sparsemat%columns)) print *, "columns:", size(sparsemat%columns)
if (allocated(sparsemat%rows))    print *, "rows:", size(sparsemat%rows)

end subroutine


subroutine sparse_deinit(sparsemat)
use mod_types,       only:sparse
type(sparse), intent(in out) :: sparsemat

if (allocated(sparsemat%values))   deallocate (sparsemat%values)
if (allocated(sparsemat%zvalues))  deallocate (sparsemat%zvalues)
if (allocated(sparsemat%columns))  deallocate (sparsemat%columns)
if (allocated(sparsemat%rows))     deallocate (sparsemat%rows)
if (allocated(sparsemat%distance)) deallocate (sparsemat%distance)
if (allocated(sparsemat%tab_dxyz)) deallocate (sparsemat%tab_dxyz)

end subroutine

subroutine full_to_sparse(fullmatrix,sparsemat)
use mod_types,       only:sparse
! sparse storage
! convert upper triangular half of symmetric matrix with full storage to simple
! coor sparse storage:
!   take an NxN matrix with a bunch of zeros and convert it to 3 arrays that collect nonzero elements and 
!   their locations. all arrays have dimensions of number of nonzeros.
!      all diagonals are stored, even if zero
!      values: all values stored in row major order.  ->->->->->
!                                                     v ->->->-> 
!      row   : row(i)    is row location of values(i)    as assoc with A(i,j)
!      column: column(j) is column location of values(j) as assoc with A(i,j)
real(dp), allocatable, intent(in)  :: fullmatrix(:,:)
type(sparse),          intent(out) :: sparsemat
integer                            :: nnzero
integer                            :: i,j,rowcnt,ialloc
real(dp), allocatable              :: tmpval(:)
integer,  allocatable              :: tmpcol(:),tmprow(:)
logical                            :: rowtest
integer                            :: N

N = size(fullmatrix(:,1))
allocate(tmpval(N*N),tmpcol(N*N),tmprow(N*N), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

nnzero=0
do i=1,N
  do j=i,N
    if ((abs(fullmatrix(i,j)) .gt. 1.0d-16) .or. (i .eq. j)) then
      nnzero = nnzero+1
      tmpval(nnzero)=fullmatrix(i,j)
      tmprow(nnzero)=i
      tmpcol(nnzero)=j
    end if
  end do
end do

call sparse_init(N,N,nnzero,sparsemat,"d")

sparsemat%rows       = tmprow(1:nnzero)   ; deallocate(tmprow)
sparsemat%values     = tmpval(1:nnzero)   ; deallocate(tmpval)
sparsemat%columns    = tmpcol(1:nnzero)   ; deallocate(tmpcol)

end subroutine
subroutine genfull_to_sparse(fullmatrix,sparsemat)
use mod_types,       only:sparse
! DMR: 11-19-2009
!      row   : row(i)    is row location of values(i)    as assoc with A(i,j)
!      column: column(j) is column location of values(j) as assoc with A(i,j)
real(dp), allocatable, intent(in)  :: fullmatrix(:,:)
type(sparse),          intent(out) :: sparsemat
integer                            :: nnzero
integer                            :: i,j,nrow,ncol,ialloc
real(dp), allocatable              :: tmpval(:)
integer,  allocatable              :: tmpcol(:),tmprow(:)
logical                            :: rowtest
integer                            :: N

nrow = size(fullmatrix,1)
ncol = size(fullmatrix,2)
allocate(tmpval(nrow*ncol),tmpcol(nrow*ncol),tmprow(nrow*ncol), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

nnzero=0
do i=1,nrow
  do j=1,ncol
    if (abs(fullmatrix(i,j)) .gt. 1.0d-16) then
      nnzero = nnzero+1
      tmpval(nnzero)=fullmatrix(i,j)
      tmprow(nnzero)=i
      tmpcol(nnzero)=j
    end if
  end do
end do

call sparse_init(nrow,ncol,nnzero,sparsemat,"d")

sparsemat%rows       = tmprow(1:nnzero)   ; deallocate(tmprow)
sparsemat%columns    = tmpcol(1:nnzero)   ; deallocate(tmpcol)
sparsemat%values     = tmpval(1:nnzero)   ; deallocate(tmpval)

end subroutine

subroutine genzfull_to_sparse(fullmatrix,sparsemat)
use mod_types,       only:sparse
! DMR: 11-19-2009
!      row   : row(i)    is row location of values(i)    as assoc with A(i,j)
!      column: column(j) is column location of values(j) as assoc with A(i,j)
complex(kind=8), allocatable, intent(in)  :: fullmatrix(:,:)
type(sparse),          intent(out) :: sparsemat
integer                            :: nnzero
integer                            :: i,j,nrow,ncol,ialloc
complex(kind=8), allocatable       :: tmpval(:)
integer,  allocatable              :: tmpcol(:),tmprow(:)
logical                            :: rowtest
integer                            :: N

nrow = size(fullmatrix,1)
ncol = size(fullmatrix,2)
allocate(tmpval(nrow*ncol),tmpcol(nrow*ncol),tmprow(nrow*ncol), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

nnzero=0
do i=1,nrow
  do j=1,ncol
    if ( (abs( real(fullmatrix(i,j))) .gt. 1.0d-16) .or. &
         (abs(aimag(fullmatrix(i,j))) .gt. 1.0d-16) ) then
      nnzero = nnzero+1
      tmpval(nnzero)=fullmatrix(i,j)
      tmprow(nnzero)=i
      tmpcol(nnzero)=j
    end if
  end do
end do

call sparse_init(nrow,ncol,nnzero,sparsemat,"z")

sparsemat%rows       = tmprow(1:nnzero)   ; deallocate(tmprow)
sparsemat%columns    = tmpcol(1:nnzero)   ; deallocate(tmpcol)
sparsemat%values     = tmpval(1:nnzero)   ; deallocate(tmpval)

end subroutine

subroutine full_to_sparse_diagto(fullmatrix,sparsemat)
use mod_types,       only:sparse
! sparse storage
! convert upper triangular half of symmetric matrix with full storage to simple
! coor sparse storage:
!   take an NxN matrix with a bunch of zeros and convert it to 3 arrays that collect nonzero elements and 
!   their locations. all arrays have dimensions of number of nonzeros.
!      values: all values stored in row major order.  ->->->->->
!                                                     v ->->->-> 
!      row   : row(i)    is row location of values(i)    as assoc with A(i,j)
!      column: column(j) is column location of values(j) as assoc with A(i,j)
real(dp), allocatable, intent(in)  :: fullmatrix(:,:)
type(sparse),          intent(out) :: sparsemat
integer                            :: nnzero
integer                            :: i,j,rowcnt,ialloc
real(dp), allocatable              :: tmpval(:)
integer,  allocatable              :: tmpcol(:),tmprow(:)
logical                            :: rowtest
integer                            :: N

N = size(fullmatrix(:,1))
allocate(tmpval(N*N),tmpcol(N*N),tmprow(N*N), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

nnzero=0
do i=1,N
  do j=i,N
    if (abs(fullmatrix(i,j)) .gt. 1.0d-16)  then
      nnzero = nnzero+1
      tmpval(nnzero)=fullmatrix(i,j)
      tmprow(nnzero)=i
      tmpcol(nnzero)=j
    end if
  end do
end do

call sparse_init(N,N,nnzero,sparsemat,"d")

sparsemat%rows       = tmprow(1:nnzero)   ; deallocate(tmprow)
sparsemat%values     = tmpval(1:nnzero)   ; deallocate(tmpval)
sparsemat%columns    = tmpcol(1:nnzero)   ; deallocate(tmpcol)

end subroutine

subroutine zfull_to_zsparse(fullmatrix,sparsemat)
use mod_types,       only:sparse
! sparse storage
! convert upper triangular half of symmetric matrix with full storage to simple
! coor sparse storage:
!   take an NxN matrix with a bunch of zeros and convert it to 3 arrays that collect nonzero elements and 
!   their locations. all arrays have dimensions of number of nonzeros.
!      all diagonals are stored, even if zero
!      values: all values stored in row major order.  ->->->->->
!                                                     v ->->->-> 
!      row   : row(i)    is row location of values(i)    as assoc with A(i,j)
!      column: column(j) is column location of values(j) as assoc with A(i,j)
complex(kind=8), allocatable, intent(in)  :: fullmatrix(:,:)
type(sparse),          intent(out)        :: sparsemat
integer                                   :: nnzero
integer                                   :: n,i,j,rowcnt,ialloc
complex(kind=8), allocatable              :: tmpval(:)
integer,  allocatable                     :: tmpcol(:),tmprow(:)
logical                                   :: rowtest

n=size(fullmatrix(1,:))

allocate(tmpval(n*n),tmpcol(n*n),tmprow(n*n), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

nnzero=0
do i=1,n
  do j=i,n
    if ( &
         (abs( real(fullmatrix(i,j))) .gt. 1.0d-16) .or. &
         (abs(aimag(fullmatrix(i,j))) .gt. 1.0d-16) .or. &
         (i .eq. j) ) then
      nnzero = nnzero+1
      tmpval(nnzero)=fullmatrix(i,j)
      tmprow(nnzero)=i
      tmpcol(nnzero)=j
    end if
  end do
end do

call sparse_init(N,N,nnzero,sparsemat,"z")

sparsemat%rows       = tmprow(1:nnzero)   ; deallocate(tmprow)
sparsemat%columns    = tmpcol(1:nnzero)   ; deallocate(tmpcol)
sparsemat%zvalues    = tmpval(1:nnzero)   ; deallocate(tmpval)


end subroutine

subroutine sparse_distances(distances,tab_dxyz,sparsemat)
use mod_types,       only:sparse
! convert coor sparse storage U to full matrix
real(dp), allocatable, intent(in)     :: distances(:,:),tab_dxyz(:,:,:)
type(sparse),          intent(in out) :: sparsemat
integer                               :: i

!upper part
do i=1,sparsemat%nnzero
  sparsemat%distance(i)   = distances(sparsemat%rows(i),sparsemat%columns(i))
  sparsemat%tab_dxyz(i,:) = tab_dxyz(sparsemat%rows(i),sparsemat%columns(i),:)
end do

end subroutine

subroutine sparse_to_full(sparsemat,fullmatrix)
use mod_types,       only:sparse
! convert coor sparse storage U to full matrix
type(sparse),          intent(in)  :: sparsemat
real(dp), allocatable, intent(out) :: fullmatrix(:,:)
integer                            :: i,j,ialloc
real(dp), parameter :: zero = 0.0d0

allocate(fullmatrix(sparsemat%ndim,sparsemat%ndim), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"


fullmatrix=zero

!upper part
do i=1,sparsemat%nnzero
  fullmatrix(sparsemat%rows(i),sparsemat%columns(i)) = sparsemat%values(i)
end do

!fill in the rest
do i=1,sparsemat%ndim
  do j=i,sparsemat%ndim 
    fullmatrix(j,i)=fullmatrix(i,j)
  end do
end do

end subroutine

subroutine zsparse_to_full(sparsemat,fullmatrix)
use mod_types,       only:sparse
! convert coor sparse storage U to full matrix
type(sparse),          intent(in)        :: sparsemat
complex(kind=8), allocatable, intent(out) :: fullmatrix(:,:)
integer                                   :: i,j,ialloc
real(dp), parameter :: zero = 0.0d0

allocate(fullmatrix(sparsemat%ndim,sparsemat%ndim), stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

fullmatrix=zero

!upper part
do i=1,sparsemat%nnzero
  fullmatrix(sparsemat%rows(i),sparsemat%columns(i)) = sparsemat%zvalues(i)
end do

!fill in the rest
do i=1,sparsemat%ndim
  do j=i,sparsemat%ndim 
    fullmatrix(j,i)=conjg(fullmatrix(i,j))
  end do
end do

end subroutine

subroutine matrx2d_init(ndim,matrx)
integer             , intent(in)  :: ndim
real(dp),allocatable, intent(out) :: matrx(:,:)
integer                           :: ii,jj,ialloc

if (allocated(matrx)) deallocate(matrx)
allocate(matrx(ndim,ndim),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

matrx=zero

end subroutine

subroutine zmatrx2d_init(ndim,matrx)
integer             , intent(in)  :: ndim
complex(kind=8),allocatable, intent(out) :: matrx(:,:)
integer                           :: ii,jj,ialloc

if (allocated(matrx)) deallocate(matrx)
allocate(matrx(ndim,ndim),stat=ialloc)
if (ialloc /= 0) STOP "*** not enough memory ***"

matrx=(zero,zero)

end subroutine

subroutine zmatrx2d_write(nfort,matrx)
integer             , intent(in)        :: nfort
complex(kind=8),allocatable, intent(in) :: matrx(:,:)
integer :: i,j

do i = 1, size(matrx(:,1))
  do j = 1, size(matrx(1,:))
    if (abs(matrx(i,j)).gt.1.0d-06) write(nfort,*) i,j,matrx(i,j)
  end do
end do

end subroutine

subroutine matrx2d_write(nfort,matrx)
integer             , intent(in) :: nfort
real(dp),allocatable, intent(in) :: matrx(:,:)
integer :: i,j

do i = 1, size(matrx(:,1))
  do j = 1, size(matrx(1,:))
    if (abs(matrx(i,j)).gt.1.0d-06) write(nfort,*) i,j,matrx(i,j)
  end do
end do

end subroutine

subroutine z_onebythree_print(zvec)
complex(kind=8), intent(in) :: zvec(3)
character(len=1) :: xchar,ychar,zchar

!return
xchar="+"
ychar="+"
zchar="+"

if(dimag(zvec(1)).lt. zero) xchar = "-"
if(dimag(zvec(2)).lt. zero) ychar = "-"
if(dimag(zvec(3)).lt. zero) zchar = "-"
print '(F10.4,a1,F7.4,a1,F10.4,a1,F7.4,a1,F10.4,a1,F7.4,a1)', &
                 dble(zvec(1)),xchar,abs(dimag(zvec(1))), "i", &
                 dble(zvec(2)),ychar,abs(dimag(zvec(2))), "i", &
                 dble(zvec(3)),zchar,abs(dimag(zvec(3))), "i"

end subroutine

subroutine z_threebythree_print(zmatrx)
complex(kind=8), intent(in) :: zmatrx(3,3)
integer :: i,j
character(len=1) :: xchar,ychar,zchar

do i = 1,3
    xchar="+";ychar="+";zchar="+"
    if(dimag(zmatrx(i,1)).lt. zero) xchar = "-"
    if(dimag(zmatrx(i,2)).lt. zero) ychar = "-"
    if(dimag(zmatrx(i,3)).lt. zero) zchar = "-"
    print '(F10.4,a1,F7.4,a1,F10.4,a1,F7.4,a1,F10.4,a1,F7.4,a1)',&
                     dble(zmatrx(i,1)),xchar,abs(dimag(zmatrx(i,1))), "i", &
                     dble(zmatrx(i,2)),ychar,abs(dimag(zmatrx(i,2))), "i", &
                     dble(zmatrx(i,3)),zchar,abs(dimag(zmatrx(i,3))), "i"
end do

end subroutine


end module



