! Solving 3D diffusion-decay eqtn. 
! 
! The boundary conditions are either:
! (a) Chemokine concentrations at specified sites, or
! (b) Chemokine secretion rates at specified sites.
!
! For each chemokine there is a list of sites with a specified concentration or secretion rate at each site.
! The lists may need to be updated when the boundary changes, or when the chemokine-producing cells
! move and/or the number changes.

module ode_diffuse_general

use chemokine
use rkf45

implicit none

integer :: ivdbug

contains

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine SetupVars
integer :: i, x, y, z

ODEdiff%ivar = 0
i = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			if (occupancy(x,y,z)%indx(1) >= 0) then
				i = i+1
				ODEdiff%ivar(x,y,z) = i
				ODEdiff%varsite(i,:) = (/x,y,z/)
			endif
		enddo
	enddo
enddo
ODEdiff%nvars = i
end subroutine

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
!
! If diffusion in an axis direction is suppressed (there is a boundary on one or
! the other neighbour sites), the current approximation used is to set the
! first derivative in that direction to zero, and compute the second derivative
! in the usual way.
! E.g. if we are considering the point (x,y,z) and ivar(x-1,y,z)=0,
! then the Laplacian d2C/dx2 + d2C/dy2 + d2C/dz2 that determines diffusion
! is reduced to d2C/dy2 + d2C/dz2.
! This is a crude simplification - in fact we should use a one-sided expression
! for d2C/dx2, derived by polynomial fitting.  The approach is as follows:
! Consider that x=0 is a boundary (dC/dx = 0), and let C(x) = a + bx + cx^2
! Then using the function values at the boundary and two grid points in,
! C0 = a
! C1 = a + b.dx + c.dx^2
! C2 = a + 2b.dx + 4c.dx^2
! solving these equations yields
! a = C0
! b = (-3C0 + 4C1 - C2)/(2dx)
! c = (C1 - C0 - b.dx)/(dx^2) = (C0 - 2C1 + C2)/(2dx^2)
! => d2C/dx2 = 2c = (C0 - 2C1 + C2)/dx^2
! In other words, the second derivative can be approximated by the second derivative
! at the neighbouring point (x+1,y,z), but this requires knowledge of concentration
! at more points.
!----------------------------------------------------------------------------------
subroutine SetupODEDiffusion
integer :: x, y, z, i, site(3), ifdc, ichemo, k, nc
real(REAL_KIND) :: DX, DX2, c(MAX_CHEMO,7)
real(REAL_KIND) :: secretion
logical :: left, right

call logger('Set diffusion parameters')
DX = 1.0
DX2 = DX*DX

if (.not.use_ODE_diffusion) return

do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	if (allocated(chemo(ichemo)%coef)) deallocate(chemo(ichemo)%coef)
	allocate(chemo(ichemo)%coef(ODEdiff%nvars,7))
enddo

do i = 1,ODEdiff%nvars
	c = 0
	nc = 0
	site = ODEdiff%varsite(i,:)
	x = site(1)
	y = site(2)
	z = site(3)
	do ichemo = 1,MAX_CHEMO
		c(ichemo,1) = -chemo(ichemo)%decay_rate
	enddo
	ODEdiff%icoef(i,1) = i
	
	left = .true.
	right = .true.
	if (x==1) then
		left = .false.
	elseif (ODEdiff%ivar(x-1,y,z) == 0) then
		left = .false.
	endif
	if (x==NX) then 
		right = .false.
	elseif (ODEdiff%ivar(x+1,y,z) == 0) then
		right = .false.
	endif
	if (left) then
		nc = nc + 1
		ODEdiff%icoef(i,2) = ODEdiff%ivar(x-1,y,z)
		do ichemo = 1,MAX_CHEMO
			c(ichemo,2) = chemo(ichemo)%diff_coef/DX2
		enddo
	endif
	if (right) then
		nc = nc + 1
		ODEdiff%icoef(i,3) = ODEdiff%ivar(x+1,y,z)
		do ichemo = 1,MAX_CHEMO
			c(ichemo,3) = chemo(ichemo)%diff_coef/DX2
		enddo
	endif
	left = .true.
	right = .true.
	if (y==1) then
		left = .false.
	elseif (ODEdiff%ivar(x,y-1,z) == 0) then
		left = .false.
	endif
	if (y==NY) then 
		right = .false.
	elseif (ODEdiff%ivar(x,y+1,z) == 0) then
		right = .false.
	endif
	if (left) then
		nc = nc + 1
		ODEdiff%icoef(i,4) = ODEdiff%ivar(x,y-1,z)
		do ichemo = 1,MAX_CHEMO
			c(ichemo,4) = chemo(ichemo)%diff_coef/DX2
		enddo
	endif
	if (right) then
		nc = nc + 1
		ODEdiff%icoef(i,5) = ODEdiff%ivar(x,y+1,z)
		do ichemo = 1,MAX_CHEMO
			c(ichemo,5) = chemo(ichemo)%diff_coef/DX2
		enddo
	endif
	left = .true.
	right = .true.
	if (z==1) then
		left = .false.
	elseif (ODEdiff%ivar(x,y,z-1) == 0) then
		left = .false.
	endif
	if (z==NZ) then 
		right = .false.
	elseif (ODEdiff%ivar(x,y,z+1) == 0) then
		right = .false.
	endif
	if (left) then
		nc = nc + 1
		ODEdiff%icoef(i,6) = ODEdiff%ivar(x,y,z-1)
		do ichemo = 1,MAX_CHEMO
			c(ichemo,6) = chemo(ichemo)%diff_coef/DX2
		enddo
	endif
	if (right) then
		nc = nc + 1
		ODEdiff%icoef(i,7) = ODEdiff%ivar(x,y,z+1)
		do ichemo = 1,MAX_CHEMO
			c(ichemo,7) = chemo(ichemo)%diff_coef/DX2
		enddo
	endif
	
	do ichemo = 1,MAX_CHEMO
		c(ichemo,1) = c(ichemo,1) - nc*chemo(ichemo)%diff_coef/DX2
		if (chemo(ichemo)%used) then
			chemo(ichemo)%coef(i,:) = c(ichemo,:) 
		endif
	enddo
!	if (i == ivdbug) then
!		write(*,*) 'ivdbug: ',nc, ODEdiff%varsite(ivdbug,:)
!		write(*,*) ODEdiff%ivar(x-1,y,z),ODEdiff%ivar(x+1,y,z)
!		write(*,*) ODEdiff%ivar(x,y-1,z),ODEdiff%ivar(x,y+1,z)
!		write(*,*) ODEdiff%ivar(x,y,z-1),ODEdiff%ivar(x,y,z+1)
!		write(*,'(a,7i6)') 'icoef: ',ODEdiff%icoef(i,:)
!		write(*,'(a,7f8.4)') 'coef: ',chemo(1)%coef(i,:)
!	endif
enddo
end subroutine


!----------------------------------------------------------------------------------
! The sources of chemokine are now being recorded in:
! occupancy%bdry (for those chemokines that enter at the follicle boundary), and 
! occupancy%DC_nbdry
!----------------------------------------------------------------------------------
subroutine deriv(t,v,dv)
real(REAL_KIND) :: t, v(*), dv(*)
integer :: i, k, n, x, y, z, idbug, ifdc, site(3), dx, dy, dz, ic, nf_DC
real(REAL_KIND) :: sum, s, vtemp, ctemp
logical :: dbug

n = ODEdiff%nvars
do i = 1,n
	sum = 0
	do k = 1,7
		if (chemo(ODEdiff%ichemo)%coef(i,k) /= 0) then
			sum = sum + chemo(ODEdiff%ichemo)%coef(i,k)*v(ODEdiff%icoef(i,k))
		endif
	enddo
	dv(i) = sum
	site = ODEdiff%varsite(i,:)
!	nf_DC = occupancy(site(1),site(2),site(3))%DC_nbdry 
!	if (nf_DC > 0) then
!		if (chemo(ODEdiff%ichemo)%use_secretion) then
!			dv(i) = dv(i) + nf_DC*chemo(ODEdiff%ichemo)%bdry_rate
!		else
!			dv(i) = 0
!		endif
!	endif
enddo

end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine SolveSteadystate_B
integer flag, i, j, k, ichemo
real(REAL_KIND) :: tstart, tend, relerr, abserr
real(REAL_KIND), allocatable :: state(:), prev_state(:), statep(:)
real(REAL_KIND), parameter :: dt = 100.0
real(REAL_KIND) :: res(10)
integer :: x, y, z, dx, dy, dz, xmid, ymid, zmid, ibnd
real(REAL_KIND) :: grad(3)
real(REAL_KIND) :: amp, r, ctemp(10)
integer :: nt(2) = (/20,100/)	! Number of rkf45 iterations with bdry concentration (1) and bdry secretion (2)
logical :: ok

write(logmsg,*) 'SolveSteadystate_B: nvars: ',ODEdiff%nvars
call logger(logmsg)
allocate(state(ODEdiff%nvars))
allocate(prev_state(ODEdiff%nvars))
allocate(statep(ODEdiff%nvars))
xmid = NX/2
ymid = NY/2
zmid = NZ/2
do ichemo = 1,MAX_CHEMO
	ODEdiff%ichemo = ichemo
	if (.not.chemo(ichemo)%used) cycle
	if (chemo(ichemo)%use_secretion) then
		ibnd = 2
	else
		ibnd = 1
	endif
	call InitState(ichemo,state)
	statep = 0
	tstart = 0
	call deriv(tstart,state,statep)

	abserr = sqrt ( epsilon ( abserr ) )
	relerr = sqrt ( epsilon ( relerr ) )

	prev_state = 0
	flag = 1
	do k = 1,nt(ibnd)
		tstart = (k-1)*dt
		tend = tstart + dt
		if (REAL_KIND == 4) then	
			call r4_rkf45 ( deriv, ODEdiff%nvars, state, statep, tstart, tend, relerr, abserr, flag )
		else
	!		call r8_rkf45 ( deriv, ODEdiff%nvars, state, statep, tstart, tend, relerr, abserr, flag )
		endif
		if (flag /= 2) then
			write(logmsg,*) 'Bad flag: ',flag
			call logger(logmsg)
		endif
		flag = 2
		call CheckConvergence(state,prev_state,ok)
		if (ok) exit
		prev_state = state
	enddo
	do x = 1,NX
		do y = 1,NY
			do z = 1,NZ
				if (occupancy(x,y,z)%indx(1) < 0) cycle
				i = ODEdiff%ivar(x,y,z)
				chemo(ichemo)%conc(x,y,z) = state(i)
				call compute_gradient(state,x,y,z,grad)
				chemo(ichemo)%grad(:,x,y,z) = grad
			enddo
		enddo
	enddo
enddo
deallocate(state)
deallocate(prev_state)
deallocate(statep)
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine CheckConvergence(s1,s2,ok)
real(REAL_KIND) :: s1(:), s2(:)
logical :: ok
integer :: i, imax
real(REAL_KIND) :: ds, df, dfmax, dsmax
real(REAL_KIND), parameter :: tol = 0.001

dfmax = 0
dsmax = 0
do i = 1,ODEdiff%nvars
	if (s1(i) > 0) then
		ds = (s1(i) - s2(i))
		df = ds/s1(i)
		if (abs(df) > abs(dfmax)) then
			imax = i
			dfmax = df
			dsmax = ds
		endif
	endif
enddo
if (abs(dfmax) < tol) then
	ok = .true.
else
	ok = .false.
endif
end subroutine

!----------------------------------------------------------------------------------
! Initialise the state vector to the current concentrations
!----------------------------------------------------------------------------------
subroutine InitState(ichemo,state)
integer :: ichemo
real(REAL_KIND) :: state(:)
integer :: x, y, z, i, site(3), nz
real(REAL_KIND) :: smin, smax

write(logmsg,*) 'InitState: ',chemo(ichemo)%name
call logger(logmsg)
smin = 1.0e10
smax = -smin
do i = 1,ODEdiff%nvars
	site = ODEdiff%varsite(i,:)
	state(i) = chemo(ichemo)%conc(site(1),site(2),site(3))
	if (state(i) < smin) smin = state(i)
	if (state(i) > smax) smax = state(i)
enddo
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

i0 = ODEdiff%ivar(x,y,z)
c0 = conc(i0)
del = 2
if (x > 1) then
	i1 = ODEdiff%ivar(x-1,y,z)
else
	i1 = 0
endif
if (x < NX) then
	i2 = ODEdiff%ivar(x+1,y,z)
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
	i1 = ODEdiff%ivar(x,y-1,z)
else
	i1 = 0
endif
if (y < NY) then
	i2 = ODEdiff%ivar(x,y+1,z)
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
	i1 = ODEdiff%ivar(x,y,z-1)
else
	i1 = 0
endif
if (z < NZ) then
	i2 = ODEdiff%ivar(x,y,z+1)
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

