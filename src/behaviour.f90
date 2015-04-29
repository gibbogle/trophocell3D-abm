! T cell behaviour
module behaviour

use global
use motility

implicit none

integer, parameter :: NORMAL_DIST      = 1
integer, parameter :: LOGNORMAL_DIST   = 2
integer, parameter :: EXPONENTIAL_DIST = 3
integer, parameter :: CONSTANT_DIST    = 4

integer, parameter :: Nbindtime_nc = 10

type(dist_type), allocatable :: stage_dist(:,:,:), life_dist(:), divide_dist(:)

contains

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_normal(p1,p2,kpar)
integer :: kpar
real(REAL_KIND) :: p1,p2
real(REAL_KIND) :: R

R = par_rnor(kpar)
rv_normal = p1+R*p2
end function

!--------------------------------------------------------------------------------------
! When Y is normal N(p1,p2) then X = exp(Y) is lognormal with
!   median = m = exp(p1)
!   shape  = s = exp(p2)
! Also median = m = mean/(s^2/2)
! kpar = parallel process number
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_lognormal(p1,p2,kpar)
integer :: kpar
real(REAL_KIND) :: p1,p2
real(REAL_KIND) :: R,z

R = par_rnor(kpar)
z = p1 + R*p2
rv_lognormal = exp(z)
end function

!--------------------------------------------------------------------------------------
! For testing.
!--------------------------------------------------------------------------------------
real(REAL_KIND) function my_rnor()
real(REAL_KIND) :: sum, R
integer :: k
integer :: kpar=0

sum = 0
do k = 1,12
!    call random_number(R)
    R = par_uni(kpar)
    sum = sum + R
enddo
my_rnor = sum - 6.0
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_exponential(p1)
real(REAL_KIND) :: p1
real(REAL_KIND) :: r
integer :: kpar = 0

r = par_rexp(kpar)
rv_exponential = p1*r
end function

!--------------------------------------------------------------------------------------
! Cumulative probability distribution for a lognormal variate with median m, shape s
! Computes Pr(X < a) where X = exp(Y) and Y is N(p1,p2), p1 = log(m), p2 = log(s)
! Pr(X < a) = Pr(Y < log(a)) = Pr(p1 + p2*R < log(a)) = Pr(R < (log(a)-p1)/p2)
! where R is N(0,1)
!--------------------------------------------------------------------------------------
real(REAL_KIND) function cum_prob_lognormal(a,p1,p2)
real(REAL_KIND) :: a, p1, p2
real(REAL_KIND) :: b, prob

!p1 = log(m)
!p2 = log(s)
b = (log(a) - p1)/p2
prob = 0.5 + 0.5*erf(b/sqrt(2.0))
cum_prob_lognormal = prob
end function

!--------------------------------------------------------------------------------
! The probability histogram for a (supposedly) lognormal variable X is given
! by the probability values p(i), i=1,..,n and the associated interval
! delimiting values x(i), i=1,..,n+1
! where p(i) = Pr{x(i) <= X < x(i+1)}
! The cumulative probability for a lognormal variate with median m, shape s is
! Pr(X < a) where X = exp(Y) and Y is N(p1,p2), p1 = log(m), p2 = log(s)
! Pr(X < a) = Pr(Y < log(a)) = Pr(p1 + p2*R < log(a)) = Pr(R < (log(a)-p1)/p2)
! where R is N(0,1)
!--------------------------------------------------------------------------------
subroutine lognfit(n,x,p,mrange,srange,mean,shape)
integer :: n
real(REAL_KIND) :: x(*), p(*), mrange(2), srange(2)
real(REAL_KIND) :: mean, shape
integer :: i
real(REAL_KIND) :: alo, ahi, p1, p2, c, err

p1 = log(mrange(1))
p2 = log(srange(1))
err = 0
do i = 1,n
	alo = x(i)
	ahi = x(i+1)
	c = cum_prob_lognormal(ahi,p1,p2) - cum_prob_lognormal(alo,p1,p2)
	err = err + (c - p(i))**2
enddo
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine read_cell_params(ok)
logical :: ok
integer :: ncpu_dummy, ntgui, idelay, iwrap, ichemo_1, ichemo_2
!real(REAL_KIND) :: sigma, divide_mean1, divide_shape1, divide_mean2, divide_shape2, real_DCradius, facs_h
!integer :: i, invitro, shownoncog, ncpu_dummy, dcsinjected, ispecial
!integer :: usetraffic, useexitchemo, useDCchemo, cognateonly, useCCL3_0, useCCL3_1, usehev, halveCD69
!character(64) :: specialfile
!character(4) :: logstr
!logical, parameter :: use_chemo = .false.

ok = .false.
write(logmsg,*) 'Read cell parameter file: ',inputfile
call logger(logmsg)
!call logger('Reading cell parameter file')
!call logger(inputfile)
open(nfcell,file=inputfile,status='old')
read(nfcell,*) DELTA_X						! lattice grid spacing
read(nfcell,*) DELTA_T						! time step (min)
read(nfcell,*) BETA							! speed: 0 < beta < 1		(0.65)
read(nfcell,*) RHO							! persistence: 0 < rho < 1	(0.95)
read(nfcell,*) days							! number of days to simulate
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
read(nfcell,*) ncpu_dummy					! # of processors - not used at the moment
read(nfcell,*) nlength						! length of tube (grids)
read(nfcell,*) nradius						! radius of tube (grids)
read(nfcell,*) nplug						! length of plug (grids)
read(nfcell,*) ntgui						! interval between GUI outputs (timesteps)
read(nfcell,*) idelay						! simulation step delay (ms)
read(nfcell,*) ichemo_1
read(nfcell,*) grad_amp(1)					! chemokine gradient amplitude
read(nfcell,*) BG_flow_amp				    ! background velocity amplitude (um/min)
read(nfcell,*) Kdrag
read(nfcell,*) Kadhesion
read(nfcell,*) Kstay
read(nfcell,*) n_cell_positions				! number of cell positions to save each time step
close(nfcell)

call logger('Finished reading cell parameter file')

use_wrapping = .false.
chemo(1)%used = (ichemo_1 == 1)
chemo(2)%used = .false.
PI = 4*atan(1.0d0)
Ve = DELTA_X*DELTA_X*DELTA_X

Nsteps = days*60*24/DELTA_T
Mnodes = 1
call make_outputfilename
open(nfout,file=outputfile,status='replace')

!call setup_dists
ok = .true.
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine make_outputfilename

outputfile = 'tropho.out'
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine save_inputfile(cellfile)
character(LEN=*) :: cellfile
character*(128) :: line

write(nfout,'(a)') '---------------------------------------------------------------------------'
write(nfout,'(a,a)') 'Input file: ',trim(cellfile)
write(nfout,*)
open(nfcell,file=cellfile,status='old')
do
    read(nfcell,'(a)',end=999) line
    write(nfout,'(a)') line
enddo
999 continue
write(nfout,*)
close(nfcell)
end subroutine

!----------------------------------------------------------------------------------------
! Save parameters (values hard-coded, not yet in input file)
!----------------------------------------------------------------------------------------
subroutine save_parameters

write(nfout,'(a)') '---------------------------------------------------------------------------'
write(nfout,'(a)') 'Hard-coded PARAMETERS'
write(nfout,'(a)') '---------------------'
write(nfout,'(a,f8.4)') 'DELTA_X: ',DELTA_X
write(nfout,'(a,f8.2)') 'DELTA_T: ',DELTA_T
write(nfout,*)
end subroutine

!-----------------------------------------------------------------------------------------
! When a T cell dies it is removed from the cell list (%ID -> 0)
! and removed from occupancy()%indx()
! The count of sites to add is decremented, for later adjustment of the blob size.
!-----------------------------------------------------------------------------------------
subroutine Tcell_death(kcell)
integer :: kcell
integer :: k, idc, site(3), indx(2), ctype, stype, region
logical :: cognate

!write(logmsg,*) 'Tcell_death: ',kcell
!call logger(logmsg)
cell_list(kcell)%ID = 0
ctype = cell_list(kcell)%ctype
ngaps = ngaps + 1
if (ngaps > max_ngaps) then
    write(*,*) 'Tcell_death: ngaps > max_ngaps'
    stop
endif
gaplist(ngaps) = kcell

NTcells = NTcells - 1
site = cell_list(kcell)%site
indx = occupancy(site(1),site(2),site(3))%indx
if (indx(1) == kcell) then
    occupancy(site(1),site(2),site(3))%indx(1) = indx(2)
    occupancy(site(1),site(2),site(3))%indx(2) = 0
elseif (indx(2) == kcell) then
    occupancy(site(1),site(2),site(3))%indx(2) = 0
else
    write(logmsg,*) 'ERROR: Tcell_death: cell not at site: ',kcell,site,indx
    call logger(logmsg)
    stop
endif

! Now we need to make a site unavailable
! Remove (make OUTSIDE) a boundary site near to the specified site.
! The criteria for selection of a site to remove are that it is on the blob boundary,
! preferably "sticking out", and that it is vacant.  A site may be made vacant by
! moving a cell to a nearby available site.
! One way to handle this is to maintain a count of the number of sites to be added(removed).
! At regular intervals the counts can be aggregated and if the total is big enough (+ or -)
! they can all be processed within a gather-scatter cycle.
if (use_add_count) then
    nadd_sites = nadd_sites - 1
else
    write(*,*) 'Tcell_death: No site removal code, use add count'
    stop
endif
end subroutine

!-----------------------------------------------------------------------------------------
! A cell divides.  It has already been determined that there is space for the extra cell.
! The other cell is placed in freeslot of site2 (which may be the same site as kcell).
! If the cell is in the periphery the site, slot etc is irrelevant.
! NOT USED CURRENTLY
!-----------------------------------------------------------------------------------------
subroutine cell_division(kcell,site2,freeslot,ok)
integer :: kcell, site2(3), freeslot
logical :: ok
integer :: icnew, ctype, gen, region, site(3), indx(2)
integer :: iseq, tag, kfrom, kto
real(REAL_KIND) :: tnow
!type(cog_type), pointer :: p1, p2
integer :: kpar = 0

ok = .true.
!tnow = istep*DELTA_T
!p1 => cell_list(kcell)%cptr
!cognate = .true.
!gen = get_generation(p1)
!call get_region(p1,region)
!if (gen == TC_MAX_GEN) then
!    write(logmsg,*) 'cell_division: reached maximum generation: ',kcell
!	call logger(logmsg)
!    return
!endif
!ndivided(gen) = ndivided(gen) + 1
!tdivided(gen) = tdivided(gen) + (tnow - p1%dividetime)
!effector = p1%effector
!
!if (ngaps > 0) then
!    icnew = gaplist(ngaps)
!    ngaps = ngaps - 1
!else
!    nlist = nlist + 1
!    if (nlist > max_nlist) then
!		write(logmsg,*) 'Error: cell_division: cell list full: ',nlist
!		call logger(logmsg)
!		ok = .false.
!		return
!	endif
!    icnew = nlist
!endif
!cell_list(kcell)%lastdir = random_int(1,6,kpar)
!gen = gen + 1
!call set_generation(p1,gen)
!ctype = cell_list(kcell)%ctype
!tag = 0
!call create_Tcell(icnew,cell_list(icnew),site2,ctype,cognate,gen,tag,0,region,.true.,ok)
!if (.not.ok) return
!
!cell_list(icnew)%ID = cell_list(kcell)%ID               ! the progeny cell inherits the parent's ID
!cell_list(icnew)%entrytime = cell_list(kcell)%entrytime ! and entrytime
!
!ndivisions = ndivisions + 1
!
!occupancy(site2(1),site2(2),site2(3))%indx(freeslot) = icnew
!if (use_add_count) then
!    nadd_sites = nadd_sites + 1
!else
!    write(logmsg,*) 'cell_division: No site removal code, use add count'
!	call logger(logmsg)
!	ok = .false.
!    return
!endif
!NTcells = NTcells + 1
end subroutine

!-----------------------------------------------------------------------------------------
! Locate a free slot in a site adjacent to site1: site2 (freeslot)
! Returns freeslot = 0 if there is no free space in an adjacent site.
! The occupancy array occ() can be either occupancy() or big_occupancy(), hence
! the need for xlim (= NXX or NX)
!-----------------------------------------------------------------------------------------
subroutine get_free_slot(occ,xlim,site1,site2,freeslot)
type(occupancy_type) :: occ(:,:,:)
integer :: xlim, site1(3), site2(3), freeslot
logical :: Lwall, Rwall
integer :: i, indx2(2)

if (site1(1) == 1) then
    Lwall = .true.
else
    Lwall = .false.
endif
if (site1(1) == xlim) then
    Rwall = .true.
else
    Rwall = .false.
endif

if (site1(1) < 1) then
    write(logmsg,*) 'get_free_slot: bad site1: ',site1
	call logger(logmsg)
    stop
endif
do i = 1,27
    if (i == 14) cycle       ! i = 14 corresponds to the no-jump case
	site2 = site1 + jumpvec(:,i)
    if (Lwall .and. site2(1) < 1) then
        cycle
    endif
    if (Rwall .and. site2(1) > xlim) then
        cycle
    endif
    if (site2(2) < 1 .or. site2(2) > NY .or. site2(3) < 1 .or. site2(3) > NZ) cycle
	indx2 = occ(site2(1),site2(2),site2(3))%indx
	if (indx2(1) >= 0) then         ! not OUTSIDE_TAG or DC
	    if (indx2(1) == 0) then     ! slot 1 is free
	        freeslot = 1
            return
	    elseif (indx2(2) == 0) then ! slot 2 is free
	        freeslot = 2
            return
        endif
    endif
enddo
freeslot = 0
end subroutine

!-----------------------------------------------------------------------------------------
! Create a new T cell.
! We need to give a progeny cell (the result of cell division) the same ID as its parent.
! This implies that there needs to be a flag to indicate that the cell results from division.
!-----------------------------------------------------------------------------------------
subroutine create_Tcell(kcell,cell,site,ctype,gen,tag,region,dividing,ok)
type(cell_type) :: cell
integer :: kcell, site(3), ctype, gen, tag, region
logical :: dividing, ok
real(REAL_KIND) :: tnow
integer :: kpar = 0

ok = .true.
tnow = istep*DELTA_T
cell%entrytime = tnow
if (dividing) then
    cell%ID = 0
else
    lastID = lastID + 1
    cell%ID = lastID
endif
cell%site = site
cell%ctype = ctype
cell%tag = tag
cell%step = 0
cell%lastdir = random_int(1,6,kpar)
cell%dtotal = 0

end subroutine

!--------------------------------------------------------------------------------
! Add a cell (kcell) with characteristics (ctype, gen, stage) at site.
!--------------------------------------------------------------------------------
subroutine add_Tcell(site,ctype,cognate,gen,tag,stage,region,kcell,ok)
integer :: site(3), ctype, gen, tag, stage, region, kcell
logical :: cognate, ok
integer :: indx(2)

ok = .true.
if (ngaps > 0) then
    kcell = gaplist(ngaps)
    ngaps = ngaps - 1
else
    nlist = nlist + 1
    if (nlist > max_nlist) then
		write(logmsg,*) 'Error: add_Tcell: cell list full: ',nlist
		call logger(logmsg)
		ok = .false.
		return
	endif
    kcell = nlist
endif
if (dbug) then
    write(*,'(a,9i7,L2)') 'add_Tcell: ',istep,kcell,site,ctype,gen,stage,region,cognate
endif
call create_Tcell(kcell,cell_list(kcell),site,ctype,gen,tag,region,.false.,ok)
if (.not.ok) return

indx = occupancy(site(1),site(2),site(3))%indx
if (indx(1) == 0) then
    indx(1) = kcell
elseif (indx(2) == 0) then
    indx(2) = kcell
else
    write(logmsg,*) 'ERROR: add_Tcell: no free slot: ',site,indx
    call logger(logmsg)
    ok = .false.
    return
endif
occupancy(site(1),site(2),site(3))%indx = indx

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine add_site_local(site)
integer :: site(3)

if (use_add_count) then
    nadd_sites = nadd_sites - 1
else
    write(*,*) 'add_site_local: no code to add site, use add count'
    stop
endif
end subroutine

!--------------------------------------------------------------------------------
! Add a vacant site at a boundary to account for the T cell added at site
! In the case of an end node (me = 0 or me = Mnodes-1) the added site can be
! at a different x value, but for internal nodes the added site  must go at
! (or near) the same x value.
! For the spherical blob case, try to place the site at a point near the bdry
! on a line drawn through site from the centre.
! NOT USED
!--------------------------------------------------------------------------------
subroutine add_vacant_site(site,kpar)
integer :: site(3),kpar
integer :: k, site0(3),newsite(3)
real(REAL_KIND) :: R
real(REAL_KIND) :: dxyz(3)
logical :: redo

if (dbug) write(*,'(a,4i6)') 'add_vacant_site: ',site
site0 = site
dxyz = real(site0) - Centre
do k = 2,3
!    call random_number(R)
    R = par_uni(kpar)
    dxyz(k) = dxyz(k) + (R - 0.5)
enddo
call normalize(dxyz)
redo = .false.
k = 0
do
    k = k+1
    newsite = site0 + k*0.5*dxyz
    if (newsite(1) < 1 .or. newsite(1) > NX) then
        site0 = site0 + (k-1)*0.5*dxyz
        dxyz(1) = 0
        call normalize(dxyz)
        redo = .true.
        exit
    endif
    if (newsite(2) < 1 .or. newsite(2) > NY .or. newsite(3) < 1 .or. newsite(3) > NZ) then
        write(*,*) 'ERROR: add_vacant_site: reached grid limits (a): ',k,site,dxyz
        stop
    endif
enddo
if (redo) then
    if (dbug) write(*,*) 'redo: ', site0
    k = 0
    do
        k = k+1
        newsite = site0 + k*0.5*dxyz
        if (newsite(2) < 1 .or. newsite(2) > NY .or. newsite(3) < 1 .or. newsite(3) > NZ) then
            write(*,*) 'ERROR: add_vacant_site: reached grid limits (b): ',k,site,dxyz
            newsite = site0 + (k-1)*0.5*dxyz
            write(*,*) newsite,occupancy(newsite(1),newsite(2),newsite(3))%indx(1)
            stop
        endif
    enddo
endif
if (dbug) write(*,'(a,4i6)') 'newsite: ',newsite
occupancy(newsite(1),newsite(2),newsite(3))%indx = 0
if (dbug) write(*,'(a,7i6)') 'site, vacant site: ',site,newsite

end subroutine

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine test_exp
integer :: i, N=100
real(REAL_KIND) :: R, sum
integer :: kpar = 0

sum = 0
write(*,*) 'start'
do i = 1,N
	R = par_rexp(kpar)
	write(*,'(i6,f12.8)') i,R
	sum = sum + R
enddo
write(*,*) 'mean = ',sum/N
sum = 0
write(*,*) 'start'
do i = 1,N
!	call random_number(R)
	R = par_uni(kpar)
	R = -log(R)
	write(*,'(i6,f12.8)') i,R
	sum = sum + R
enddo
write(*,*) 'mean = ',sum/N
end subroutine

!-----------------------------------------------------------------------------------------
! Initial cell position data is loaded into occupancy() and cell_list().
! The 
!-----------------------------------------------------------------------------------------
subroutine PlaceCells(ok)
logical :: ok
integer :: x, y, z, kcell, site(3)
integer :: gen, tag, region, ctype
integer :: x0
real(REAL_KIND) :: R
integer :: kpar = 0

ok = .false.
write(nflog,*) nlength,nradius,nplug
write(nflog,*) NX,NY,NZ
istep = 0
gen = 1
ctype = 1
region = 1
tag = 0
x0 = NX/2 - nplug/2
kcell = 0
do x = x0,x0+nplug-1
	do y = 1,NY
		do z = 1,NZ
			if (occupancy(x,y,z)%indx(1) == 0) then ! vacant site
				kcell = kcell + 1
				site = (/x,y,z/)
				gen = 1
				region = 0	! don't know about regions yet
				ctype = TROPHO_CELL
				tag = 0
				call create_Tcell(kcell,cell_list(kcell),site,ctype,gen,tag,region,.false.,ok)
				if (.not.ok) return
				occupancy(x,y,z)%indx(1) = kcell
			endif
		enddo
	enddo
enddo
NTcells0 = kcell
nlist = NTcells0
NTcells = NTcells0

end subroutine

end module
