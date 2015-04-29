! To solve 3D diffusion-decay eqtn as an ODE system
! This formulation is based on chemokine secretion into sites adjacent to DCs

module ode_diffuse_secretion

use rkf45
use behaviour

implicit none

!integer, parameter :: REAL_KIND = 4

integer, allocatable :: ivar(:,:,:)
integer, allocatable :: varsite(:,:)
real(REAL_KIND), allocatable :: coef(:,:)
integer, allocatable :: icoef(:,:)
!integer :: NX, NY, NZ
integer :: nvars, mid
real (REAL_KIND) :: DX
real(REAL_KIND), parameter :: Kdiff = 0.01
real(REAL_KIND), parameter :: Kdecay = 0.002
real(REAL_KIND), parameter :: BASE_CHEMOKINE_SECRETION = 0.012
!real(REAL_KIND), parameter :: secretion = 1.0

!type occupancy_type
!	integer :: indx(2)
!end type
!type(occupancy_type) :: occupancy(50,50,50)
!
!type DC_type
!	integer :: site(3)
!	real :: secretion
!end type
!type(DC_type) :: DClist(100)
!integer :: NDC


contains

!----------------------------------------------------------------------------------
! Solve for a test case with solute flux into the gridcells with x = 1.
! For each active site (x,y,z) mapping to i = ivar(x,y,z) there are
! 7 coefficients of the associated variables, i.e. those at:
! 1 <-- (x,y,z)
! 2 <-- (x-1,y,z)
! 3 <-- (x+1,y,z)
! 4 <-- (x,y-1,z)
! 5 <-- (x,y+1,z)
! 6 <-- (x,y,z-1)
! 7 <-- (x,y,z+1)
! In each case, provided ivar(:,:,:) > 0
!----------------------------------------------------------------------------------
subroutine setup_diffusion_secretion
integer :: x, y, z, i, site(3), idc
real(REAL_KIND) :: DX2, c(7)
real :: secretion

!NX = 41
!NY = 41
!NZ = 41
DX = 1.0
DX2 = DX*DX
mid = (NX+1)/2
allocate(ivar(NX,NY,NZ))
allocate(varsite(NX*NY*NZ,3))
allocate(coef(NX*NY*NZ,7))
allocate(icoef(NX*NY*NZ,7))
ivar = 0
i = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			if (occupancy(x,y,z)%indx(1) >= 0) then
				i = i+1
				ivar(x,y,z) = i
				varsite(i,:) = (/x,y,z/)
			endif
		enddo
	enddo
enddo
nvars = i

do i = 1,nvars
	c = 0
	site = varsite(i,:)
	x = site(1)
	y = site(2)
	z = site(3)
	c(1) = -Kdecay
	icoef(i,1) = i
	if (x==1) then 
		c(2) = 0
		c(3) = 2*Kdiff/DX2
		icoef(i,3) = ivar(x+1,y,z)
	elseif (ivar(x-1,y,z) == 0) then
		c(2) = 0
		c(3) = 2*Kdiff/DX2
		icoef(i,3) = ivar(x+1,y,z)
	elseif (x==NX) then 
		c(2) = 2*Kdiff/DX2
		c(3) = 0
		icoef(i,2) = ivar(x-1,y,z)
	elseif (ivar(x+1,y,z) == 0) then
		c(2) = 2*Kdiff/DX2
		c(3) = 0
		icoef(i,2) = ivar(x-1,y,z)
	else
		c(2) = Kdiff/DX2
		c(3) = Kdiff/DX2
		icoef(i,2) = ivar(x-1,y,z)
		icoef(i,3) = ivar(x+1,y,z)
	endif
	if (y==1) then 
		c(4) = 0
		c(5) = 2*Kdiff/DX2
		icoef(i,5) = ivar(x,y+1,z)
	elseif (ivar(x,y-1,z) == 0) then
		c(4) = 0
		c(5) = 2*Kdiff/DX2
		icoef(i,5) = ivar(x,y+1,z)
	elseif (y==NY) then 
		c(4) = 2*Kdiff/DX2
		c(5) = 0
		icoef(i,4) = ivar(x,y-1,z)
	elseif (ivar(x,y+1,z) == 0) then
		c(4) = 2*Kdiff/DX2
		c(5) = 0
		icoef(i,4) = ivar(x,y-1,z)
	else
		c(4) = Kdiff/DX2
		c(5) = Kdiff/DX2
		icoef(i,4) = ivar(x,y-1,z)
		icoef(i,5) = ivar(x,y+1,z)
	endif
	if (z==1) then 
		c(6) = 0
		c(7) = 2*Kdiff/DX2
		icoef(i,7) = ivar(x,y,z+1)
	elseif (ivar(x,y,z-1) == 0) then
		c(6) = 0
		c(7) = 2*Kdiff/DX2
		icoef(i,7) = ivar(x,y,z+1)
	elseif (z==NZ) then 
		c(6) = 2*Kdiff/DX2
		c(7) = 0
		icoef(i,6) = ivar(x,y,z-1)
	elseif (ivar(x,y,z+1) == 0) then
		c(6) = 2*Kdiff/DX2
		c(7) = 0
		icoef(i,6) = ivar(x,y,z-1)
	else
		c(6) = Kdiff/DX2
		c(7) = Kdiff/DX2
		icoef(i,6) = ivar(x,y,z-1)
		icoef(i,7) = ivar(x,y,z+1)
	endif
	c(1) = c(1) - 6*Kdiff/DX2
	coef(i,:) = c
enddo
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!subroutine add_DC(site,secretion)
!integer :: site(3)
!real :: secretion
!integer :: del
!
!NDC = NDC + 1
!DClist(NDC)%site = site
!DClist(NDC)%secretion = secretion
!occupancy(site(1),site(2),site(3))%indx(1) = -NDC
!do del = -1,1,2
!	occupancy(site(1)+del,site(2),site(3))%indx(1) = -NDC
!	occupancy(site(1),site(2)+del,site(3))%indx(1) = -NDC
!	occupancy(site(1),site(2),site(3)+del)%indx(1) = -NDC
!enddo
!end subroutine

!----------------------------------------------------------------------------------
! The concentration field is created by secretion at sites near DCs.  The sites
! receiving secretion are those with r^2 <= r2lim that are not DC sites.
!----------------------------------------------------------------------------------
subroutine deriv_secretion(t,v,dv)
real(REAL_KIND) :: t, v(*), dv(*)
integer :: i, k, n, x, y, z, idbug, idc, site(3), dx, dy, dz
real(REAL_KIND) :: sum, s, r2lim = 4.1
logical :: dbug

idbug = ivar(mid+1,mid,mid)
n = nvars
do i = 1,n
	if (i == -idbug) then
		dbug = .true.
	else
		dbug = .false.
	endif
	sum = 0
	do k = 1,7
		if (coef(i,k) /= 0) then
			sum = sum + coef(i,k)*v(icoef(i,k))
			if (dbug) then
				write(*,'(i4,3f8.4)') k,coef(i,k),v(icoef(i,k)),sum
			endif
		endif
	enddo
	dv(i) = sum
	if (dbug) then
		write(*,'(a,2f8.5)') 'dv: ',t,dv(i)
	endif
enddo
!i = ivar(mid,mid,mid)
!dv(i) = dv(i) + secretion
!do idc = 1,NDC
!	site = DClist(idc)%site
!	do dx = -2,2
!		x = site(1)+dx
!		do dy = -2,2
!			y = site(2)+dy
!			do dz = -2,2
!				z = site(3)+dz
!				if (dx*dx+dy*dy+dz*dz > r2lim) cycle
!				if (occupancy(x,y,z)%indx(1) >= 0) then
!					i = ivar(x,y,z)
!					dv(i) = dv(i) + DClist(idc)%secretion
!				endif
!			enddo
!		enddo
!	enddo
!enddo

end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine diffuse_steadystate_secretion
integer flag, i, j, k
real(REAL_KIND) :: tstart, tend, relerr, abserr
real(REAL_KIND), allocatable :: state(:), statep(:)
real(REAL_KIND), parameter :: dt = 20.0
real(REAL_KIND) :: res(10)
integer :: x, y, z, dx, dy, dz
real(REAL_KIND) :: grad(3), g(3), gamp, gmax
real :: amp, r

allocate(state(nvars))
allocate(statep(nvars))
state = 0
statep = 0
!do x = 1,mid
!	do y = 1,NY
!		do z = 1,NZ
!			i = ivar(x,y,z)
!			state(i) = 1.0
!		enddo
!	enddo
!enddo
tstart = 0
call deriv_secretion(tstart,state,statep)

abserr = sqrt ( epsilon ( abserr ) )
relerr = sqrt ( epsilon ( relerr ) )

flag = 1
do k = 1,200
	tstart = (k-1)*dt
	tend = tstart + dt
	if (REAL_KIND == 4) then	
		call r4_rkf45 ( deriv_secretion, nvars, state, statep, tstart, tend, relerr, abserr, flag )
	else
!		call r8_rkf45 ( deriv_secretion, nvars, state, statep, tstart, tend, relerr, abserr, flag )
	endif
	if (flag /= 2) then
		write(*,*) 'Bad flag: ',flag
	endif
	flag = 2
	do j = 1,8
		i = ivar(mid+j-4,mid,mid)
		if (i > 0) then
			res(j) = state(i)
		else
			res(j) = 0
		endif
	enddo
	write(*,'(i6,10f8.2)') k,tend,(res(j),j=1,8)
enddo
write(*,*)
do k = 1,NX
	i = ivar(k,mid,mid)
	if (i > 0) then
		write(*,'(i4,f8.2)') k,state(i)
	endif
enddo
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			if (occupancy(x,y,z)%indx(1) < 0) cycle
			i = ivar(x,y,z)
!			occupancy(x,y,z)%chemo_conc = state(i)
!			call compute_gradient(state,x,y,z,grad)
			chemo(1)%conc(x,y,z) = state(i)
			call compute_gradient(state,x,y,z,chemo(1)%grad(:,x,y,z))
!			occupancy(x,y,z)%chemo_grad = grad
		enddo
	enddo
enddo
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			g = chemo(1)%grad(:,x,y,z)
			gamp = sqrt(dot_product(g,g))
			gmax = max(gamp,gmax)
		enddo
	enddo
enddo
write(logmsg,'(a,f8.3)') 'Max gradient: ',gmax
call logger(logmsg)
!if (use_single_DC) then
!	write(*,*) 'Single DC site: ',DClist(1)%site
!	do dx = -7,7
!		x = DClist(1)%site(1) + dx
!		do dy = -7,7
!			y = DClist(1)%site(2) + dy
!			do dz = -7,7
!				z = DClist(1)%site(3) + dz
!!				amp = norm(real(occupancy(x,y,z)%chemo_grad))
!				amp = norm(real(chemo(1)%grad(:,x,y,z)))
!				if (amp > 0) then
!					r = sqrt(real(dx*dx+dy*dy+dz*dz))
!					write(nfout,'(2f8.3)') r,amp
!				endif
!			enddo
!		enddo
!	enddo
!endif
deallocate(state)
deallocate(statep)
end subroutine

!----------------------------------------------------------------------------------
! Compute the gradient vector for chemokine concentration at (x,y,z).
! The gradient determination is 2-sided if possible, 1-sided if only one adjacent
! value is available, and the gradient is set to 0 if neither is possible.
!----------------------------------------------------------------------------------
subroutine compute_gradient(conc,x,y,z,grad)
real(REAL_KIND) :: conc(:)
integer :: x, y, z
real(REAL_KIND) :: grad(3)
integer :: i0, i1, i2
real(REAL_KIND) :: c0, c1, c2, del

i0 = ivar(x,y,z)
c0 = conc(i0)
del = 2
if (x > 1) then
	i1 = ivar(x-1,y,z)
else
	i1 = 0
endif
if (x < NX) then
	i2 = ivar(x+1,y,z)
else
	i2 = 0
endif
if (i1 > 0) then
	c1 = conc(i1)
else
	c1 = c0
	del = del - 1
endif			
if (i2 > 0) then
	c2 = conc(i2)
else
	c2 = c0
	del = del - 1
endif
if (del > 0) then
	grad(1) = (c2 - c1)/del
else
	grad(1) = 0
endif			

del = 2
if (y > 1) then
	i1 = ivar(x,y-1,z)
else
	i1 = 0
endif
if (y < NY) then
	i2 = ivar(x,y+1,z)
else
	i2 = 0
endif
if (i1 > 0) then
	c1 = conc(i1)
else
	c1 = c0
	del = del - 1
endif			
if (i2 > 0) then
	c2 = conc(i2)
else
	c2 = c0
	del = del - 1
endif
if (del > 0) then
	grad(2) = (c2 - c1)/del
else
	grad(2) = 0
endif

del = 2
if (z > 1) then
	i1 = ivar(x,y,z-1)
else
	i1 = 0
endif
if (z < NZ) then
	i2 = ivar(x,y,z+1)
else
	i2 = 0
endif
if (i1 > 0) then
	c1 = conc(i1)
else
	c1 = c0
	del = del - 1
endif			
if (i2 > 0) then
	c2 = conc(i2)
else
	c2 = c0
	del = del - 1
endif
if (del > 0) then
	grad(3) = (c2 - c1)/del
else
	grad(3) = 0
endif
end subroutine

end module

