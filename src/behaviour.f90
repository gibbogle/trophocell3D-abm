! T cell behaviour
module behaviour

use global
use packer
use chemokine
!use motility

implicit none

integer, parameter :: NORMAL_DIST      = 1
integer, parameter :: LOGNORMAL_DIST   = 2
integer, parameter :: EXPONENTIAL_DIST = 3
integer, parameter :: CONSTANT_DIST    = 4

integer, parameter :: Nbindtime_nc = 10

type(dist_type), allocatable :: stage_dist(:,:,:), life_dist(:)	!, divide_dist(:)

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

ok = .false.
write(logmsg,*) 'Read cell parameter file: ',inputfile
call logger(logmsg)
!call logger('Reading cell parameter file')
!call logger(inputfile)
open(nfcell,file=inputfile,status='old')
read(nfcell,*) DELTA_T						! time step (min)
read(nfcell,*) BETA							! speed: 0 < beta < 1		(0.65)
read(nfcell,*) RHO							! persistence: 0 < rho < 1	(0.95)
read(nfcell,*) days							! number of days to simulate
read(nfcell,*) settle_hrs					! number of hours for settling
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
read(nfcell,*) ncpu_dummy					! # of processors - not used at the moment
read(nfcell,*) tube_length					! length of tube (um)
read(nfcell,*) tube_radius					! radius of tube (um)
read(nfcell,*) plug_zmin					! plug z from (um)
read(nfcell,*) plug_zmax					! plug z to (um)
read(nfcell,*) plug_hmax					! plug height (um)
read(nfcell,*) Raverage						! average cell radius (um)
read(nfcell,*) ntgui						! interval between GUI outputs (timesteps)
read(nfcell,*) idelay						! simulation step delay (ms)
read(nfcell,*) ichemo_1
read(nfcell,*) grad_amp(1)					! chemokine gradient amplitude
read(nfcell,*) BG_flow_amp				    ! background velocity amplitude (um/min)
!read(nfcell,*) Kdrag
!read(nfcell,*) Kadhesion
!read(nfcell,*) Kstay
read(nfcell,*) a_separation
read(nfcell,*) a_force
read(nfcell,*) c_force
read(nfcell,*) x0_force
read(nfcell,*) x1_force
read(nfcell,*) kdrag
read(nfcell,*) frandom
read(nfcell,*) n_cell_positions				! number of cell positions to save each time step
close(nfcell)

call logger('Finished reading cell parameter file')

chemo(1)%used = (ichemo_1 == 1)
chemo(2)%used = .false.
PI = 4*atan(1.0d0)

DELTA_T = 60*DELTA_T			! min -> sec
Nsteps = days*3600*24/DELTA_T
!Mnodes = 1
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
! The count of sites to add is decremented, for later adjustment of the blob size.
!-----------------------------------------------------------------------------------------
subroutine Tcell_death(kcell)
integer :: kcell
integer :: k, idc, site(3), indx(2), ctype, stype, region
logical :: cognate

!write(logmsg,*) 'Tcell_death: ',kcell
!call logger(logmsg)
cell_list(kcell)%ID = 0
cell_list(kcell)%state = DEAD
ctype = cell_list(kcell)%ctype
ngaps = ngaps + 1
if (ngaps > max_ngaps) then
    write(*,*) 'Tcell_death: ngaps > max_ngaps'
    stop
endif
gaplist(ngaps) = kcell

Ncells = Ncells - 1
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
end subroutine

!-----------------------------------------------------------------------------------------
! Create a new T cell.
! We need to give a progeny cell (the result of cell division) the same ID as its parent.
! This implies that there needs to be a flag to indicate that the cell results from division.
!-----------------------------------------------------------------------------------------
subroutine create_Tcell(kcell,rsite,ctype,gen,tag,region,dividing,ok)
type(cell_type), pointer :: cp
integer :: kcell, ctype, gen, tag, region
real(REAL_KIND) :: rsite(3)
logical :: dividing, ok
real(REAL_KIND) :: tnow
integer :: kpar = 0

ok = .true.
tnow = istep*DELTA_T
cp => cell_list(kcell)
!cell%entrytime = tnow
cp%birthtime = tnow
if (dividing) then	!???????????????
    cp%ID = 0
else
    lastID = lastID + 1
    cp%ID = lastID
endif
cp%state = ALIVE
cp%nspheres = 1
cp%centre(:,1) = rsite
cp%radius = Raverage
cp%ctype = ctype
cp%tag = tag
cp%step = 0
!cell%lastdir = random_int(1,6,kpar)
cp%dtotal = 0
cp%mitosis = 0
cp%Iphase = .true.

end subroutine

!--------------------------------------------------------------------------------
! Add a cell (kcell) with characteristics (ctype, gen, stage) at site.


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
! Possible cell locations are in a close-packed rectangular block.  Locations outside the
! cylindrical vessel are eliminated, then cells in the assumed pre-existing gap are
! eliminated.
!
! The origin is at the centre of the tube, at the placenta end.
! The Z axis is along the tube centre-line.
!
! Tube dimensions:
!	radius = tube_radius (um)
!	length = tube_length (um)
! Plug dimensions:
!	from   = plug_zmin (um)
!	to     = plug_zmax (um)
!	height = plug_hmax (um)
!-----------------------------------------------------------------------------------------
subroutine PlaceCells(ok)
logical :: ok
integer :: kcell, i, ix, iy, iz, nxx, nyy, nzz
integer :: gen, region, ctype, tag
real(REAL_KIND) :: cell_radius, xrng, zrng, z0, dz, u, h
real(REAL_KIND) :: centre(3), rsite(3), block_centre(3), plug_centre(3), x, y, z, r

z0 = (plug_zmax + plug_zmin)/2
xrng = tube_radius
zrng = plug_zmax - plug_zmin
dz = zrng/2
plug_centre = [0.d0, 0.d0, z0]
cell_radius = Raverage
kcell = 0
if (use_packing) then
	nxx = 1.2*xrng/cell_radius
	nzz = 1.5*dz/cell_radius
	nyy = 1.3*nxx
	call SelectCellLocations(nxx, nyy, nzz, cell_radius, block_centre)

	write(nflog,*) 'nxx,nyy,nzz: ',nxx,nyy,nzz
	write(nflog,'(a,3f8.1)') 'plug_centre: ',plug_centre
	write(nflog,'(a,3f8.1)') 'block_centre: ',block_centre
	do i = 1,nxx*nyy*nzz
		if (cdist(i)+cell_radius > tube_radius) cycle
		ix = xyz_lookup(i)%x
		iy = xyz_lookup(i)%y
		iz = xyz_lookup(i)%z
		centre = cloc(ix,iy,iz)%centre - block_centre + plug_centre
		if (centre(3) < plug_zmin .or. centre(3) > plug_zmax) cycle
		u = centre(3) - z0
		h = (1 + cos(u*PI/dz))*plug_hmax/2
		if (cdist(i) < tube_radius - h) cycle
		kcell = kcell + 1
		rsite = centre
	!	write(nflog,'(2i6,3i4,3f8.1)') kcell,i,ix,iy,iz,rsite
		gen = 1
		region = 0	! don't know about regions yet
		ctype = TROPHO_CELL
		tag = 0
		call create_Tcell(kcell,rsite,ctype,gen,tag,region,.false.,ok)
		if (.not.ok) return
	enddo
	call FreeCellLocations
else
	nxx = tube_radius/(2*cell_radius) + 1
	nzz = dz/(2*cell_radius) + 1
	do ix = -nxx,nxx
		do iy = -nxx,nxx
			do iz = -nzz,nzz
				x = ix*2*cell_radius	! (x,y,z) is the offset from the centre of the plug
				y = iy*2*cell_radius
				z = iz*2*cell_radius
				r = sqrt(x*x + y*y)
				if (r + cell_radius > tube_radius) cycle
				if (z < -dz .or. z > dz) cycle
				h = (1 + cos(z*PI/dz))*plug_hmax/2
				if (r < tube_radius - h) cycle
				kcell = kcell + 1
				rsite = [x, y, z0+z]
				write(nflog,'(i6,3i4,3f8.1)') kcell,ix,iy,iz,rsite
				gen = 1
				region = 0	! don't know about regions yet
				ctype = TROPHO_CELL
				tag = 0
				call create_Tcell(kcell,rsite,ctype,gen,tag,region,.false.,ok)
				if (.not.ok) return
			enddo
		enddo
	enddo
endif
nlist = kcell
ncells = kcell


end subroutine



end module
