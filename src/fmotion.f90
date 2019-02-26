! Motion based on forces
module fmotion

use global
use behaviour

implicit none

contains

!-----------------------------------------------------------------------------------------
! The time step size dt_move is set dynamically here.
!-----------------------------------------------------------------------------------------
subroutine fmover(dt_move, done, ok)
real(REAL_KIND) :: dt_move
logical :: done, ok
integer :: k1, kcell, kpar, nd, nr, nc(0:8), kfrom(0:8), kto(0:8), tMnodes
real(REAL_KIND), allocatable :: force(:,:,:)
real(REAL_KIND) :: fmax, dx1(3), dx2(3)
type(cell_type), pointer :: cp

done = .false.
ok = .true.
dt_move = ndt*delta_tmove
!write(*,*) 'fmover: ',ndt,delta_tmove,dt_move
if (ncells <= 1 .and. cell_list(1)%Iphase) then
	done = .true.
	return
endif
allocate(force(3,ncells,2))
call forces(force,fmax,ok)

if (fmax > 0) then
	do
		if (dt_move*fmax/kdrag > delta_max) then
			ndt = ndt/2
			if (ndt <= 1) then
				ndt = 1
				dt_move = ndt*delta_tmove
				exit
			endif
		elseif (dt_move*fmax/kdrag < delta_min) then
			ndt = ndt + 1
			if (ndt >= ndt_max) then
				ndt = ndt_max
				dt_move = ndt*delta_tmove
				exit
			endif
		endif
		exit
	enddo
endif
if (dt_move + t_fmover > DELTA_T) then
	dt_move = DELTA_T - t_fmover
	done = .true.
endif

if (ncells < 500) then
	tMnodes = 1
else
	tMnodes = Mnodes
endif
if (ncells <= tMnodes) then
    do kpar = 0,tMnodes-1
        if (kpar < ncells) then
            kfrom(kpar) = kpar+1
            kto(kpar) = kpar+1
        else
            kfrom(kpar) = 1
            kto(kpar) = 0
        endif
    enddo
else
    nd = ncells/tMnodes
    nr = ncells - nd*tMnodes
    do kpar = 0,tMnodes-1
        kfrom(kpar) = kpar*nd + 1
        kto(kpar) = (kpar+1)*nd
        if (kpar == tMnodes-1) kto(kpar) = ncells
    enddo
endif

call omp_set_num_threads(tMnodes)
!$omp parallel do private(k1,kcell,cp,dx1,dx2)
do kpar = 0,tMnodes-1
    do k1 = kfrom(kpar),kto(kpar)
		kcell = perm_index(k1)
		cp => cell_list(kcell)
		if (cp%state == GONE_BACK .or. cp%state == GONE_THROUGH) cycle
		if (cp%Iphase) then
			dx1 = dt_move*force(:,k1,1)/kdrag
			cp%centre(:,1) = cp%centre(:,1) + dx1
		else
			dx1 = dt_move*force(:,k1,1)/kdrag
			dx2 = dt_move*force(:,k1,2)/kdrag
			cp%centre(:,1) = cp%centre(:,1) + dx1
			cp%centre(:,2) = cp%centre(:,2) + dx2
		endif
		if (cp%centre(3,1)>tube_length) then
         !   call cell_death(kcell)
            cp%state=GONE_BACK
        elseif (cp%centre(3,1)<0) then
            cp%state=GONE_THROUGH
        endif
	enddo
enddo
!omp end parallel do
deallocate(force)
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! Contact forces are cell-cell or cell-wall
!-----------------------------------------------------------------------------------------
subroutine forces(force,fmax,ok)
real(REAL_KIND) :: fmax, force(:,:,:)
logical :: ok
integer :: k1, kcell, kpar, nd, nr, nc(0:8), kfrom(0:8), kto(0:8), tMnodes, k2
real(REAL_KIND) :: F(3,2), r(3), c(3), radius, amp, fsum, dFwall(3), dFflow(3), dFchemo(3)
real(REAL_KIND),allocatable :: cell_fmax(:)
real(REAL_KIND) :: fflow = 0.005

if (ncells < 100) then
	tMnodes = 1
else
	tMnodes = Mnodes
endif
if (ncells <= tMnodes) then
    do kpar = 0,tMnodes-1
        if (kpar < ncells) then
            kfrom(kpar) = kpar+1
            kto(kpar) = kpar+1
        else
            kfrom(kpar) = 1
            kto(kpar) = 0
        endif
    enddo
else
    nd = ncells/tMnodes
    nr = ncells - nd*tMnodes
    do kpar = 0,tMnodes-1
        kfrom(kpar) = kpar*nd + 1
        kto(kpar) = (kpar+1)*nd
        if (kpar == tMnodes-1) kto(kpar) = ncells
    enddo
endif
!write(*,'(19i6)') ncells,Mnodes,tmnodes,kfrom(0:tMnodes-1),kto(0:tMnodes-1)
allocate(cell_fmax(ncells))
cell_fmax = 0
fsum = 0
call omp_set_num_threads(tMnodes)

!write(*,'(a)') 'total force'

!$omp parallel do private(k1,kcell,c,radius,r,dFflow,dFwall,dFchemo,amp,F,ok)
do kpar = 0,tMnodes-1
    do k1 = kfrom(kpar),kto(kpar)
	    kcell = perm_index(k1)
	    !write(*,'(a,1i3)') 'cell num: ',kcell
	    if (cell_list(kcell)%state == GONE_BACK .or. cell_list(kcell)%state == GONE_THROUGH) cycle
	    c = cell_list(kcell)%centre(:,1)
	    radius = cell_list(kcell)%radius(1)
	    call get_cell_force(kcell,F,ok)
	    !write(*,'(a,3f8.1)') 'Cell-cell force: ',F
	    call get_random_dr(r)
	    F(:,1) = F(:,1) + frandom*r
        call get_flow_force(kcell,c,radius,dFflow)
		F(:,1) = F(:,1) + dFflow
        call get_chemo_force(c,radius,dFchemo)
        F(:,1) = F(:,1) + dFchemo
		kcell_debug = kcell
	    call get_wall_force(c,radius,dFwall,ok)
	    if (.not.ok) then
			write(*,*) 'Bad wall force'
			stop
		endif
	    F(:,1) = F(:,1) + dFwall
	    force(:,k1,1) = F(:,1)
        !write(*, '(3f10.1,$)') force(:,k1,1)
	    amp = sqrt(dot_product(F(:,1),F(:,1)))
	    cell_fmax(k1) = max(cell_fmax(k1),amp)
	    if (.not.cell_list(kcell)%Iphase) then
			call get_random_dr(r)
			F(:,2) = F(:,2) + frandom*r
		    force(:,k1,2) = F(:,2)
		    cell_fmax(k1) = max(cell_fmax(k1),sqrt(dot_product(F(:,2),F(:,2))))
        endif
    enddo
enddo
!write(nflog,*)
!omp end parallel do
fmax = maxval(cell_fmax(:))
deallocate(cell_fmax)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_cell_force(kcell,F,ok)
integer :: kcell
real(REAL_KIND) :: F(3,2)
logical :: ok
integer :: k2, knbr
integer :: isphere1, isphere2, nspheres1, nspheres2
real(REAL_KIND) :: R1, c1(3), R2, c2(3), v(3), d, d_hat, dF, fmult
logical :: incontact, Mphase
type(cell_type), pointer :: cp1, cp2
!integer :: nvar	! nv = ncells + # of Mphase cells

ok = .true.
F = 0
cp1 => cell_list(kcell)
Mphase = .not.cp1%Iphase
if (Mphase) then
	nspheres1 = 2
else
	nspheres1 = 1
endif
F = 0
do k2 = 1,cp1%nbrs
	knbr = cp1%nbrlist(k2)%indx
	cp2 => cell_list(knbr)
	fmult = 1
	if (cp2%state == WALL) then
		fmult = 2
	endif
	if (cp2%Iphase) then
		nspheres2 = 1
	else
		nspheres2 = 2
	endif
	do isphere1 = 1,nspheres1
		R1 = cp1%radius(isphere1)
		c1 = cp1%centre(:,isphere1)
!		if (isphere1 == 1) then
!			iv = cp1%varindex
!		else
!			iv = cp1%varindex + 1
!		endif
		do isphere2 = 1,nspheres2
			R2 = cp2%radius(isphere2)
			c2 = cp2%centre(:,isphere2)
			v = c1 - c2
			d = sqrt(dot_product(v,v))
			v = v/d
!			if (cp2%state == WALL) then
!				write(*,'(a,2i6,f8.2)') 'WALL nbr d: ',kcell,knbr,d
!			endif
			incontact = cp1%nbrlist(k2)%contact(isphere1,isphere2)
			dF = get_force(R1,R2,d,incontact,ok)	! returns magnitude and sign of force
            if (.not.ok) then
                write(*,*) kcell,knbr,isphere1,isphere2
                write(*,'(a,i6,3e12.3)') 'nbrlist: ',kcell,c1
                write(*,'(10i6)') cp1%nbrlist(1:cp1%nbrs)%indx
                write(*,'(a,i6,3e12.3)') 'nbrlist: ',knbr,c2
                write(*,'(10i6)') cp2%nbrlist(1:cp2%nbrs)%indx
                stop
            endif
            dF = fmult*dF
			if (isphere1 == 1) F(:,isphere1) = F(:,isphere1) + dF*v
			if (Mphase .and. isphere1 == 2) F(:,isphere1) = F(:,isphere1) + dF*v
		enddo
	enddo
enddo
if (Mphase) then	! compute force between separating spheres, based on deviation from desired d = cp1%d
	c1 = cp1%centre(:,1)
	c2 = cp1%centre(:,2)
	v = c1 - c2
	d = sqrt(dot_product(v,v))
	v = v/d
	d_hat = max(cp1%mitosis*cp1%d_divide,small_d)
	dF = get_separating_force(d,d_hat)
!		write(*,'(a,3e12.3)') 'separating dF: ',d,d_hat,dF
	F(:,1) = F(:,1) + dF*v	! dF acts in the direction of v on sphere #1 when dF > 0
	F(:,2) = F(:,2) - dF*v	! dF acts in the direction of -v on sphere #2 when dF > 0
endif

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_wall_force(c,radius,F,ok)
real(REAL_KIND) :: c(3), radius, F(3)
logical :: ok
real(REAL_KIND) :: dF(3)
real(REAL_KIND) :: R, d, x, Fset, Famp, b_force, dx, delta
real(REAL_KIND) :: fwall = 0.005
real(REAL_KIND) :: settling_factor = 1

dx = x1_force - x0_force
!b_force = -c_force*110 - 4*a_force/dx**2
b_force = -c_force_wall - 4*a_force/dx**2
delta = dx**2 + 4*a_force/b_force
xcross2_force = (x0_force+x1_force)/2 + 0.5*sqrt(delta)

ok = .true.
if (c(3) < 0 .or. c(3) > tube_length) then
	F = 0
	return
endif
c(3) = 0
R = sqrt(dot_product(c,c))
d = tube_radius - R	! distance of the cell centre from the wall
x = d/radius
if (x < x0_force) then
	write(logmsg,'(a,4e12.3)') 'Error: get_wall_force: x < x0: ',R,d,x,x0_force
	call logger(logmsg)
	F = 0
	ok = .false.
	return
endif
if (x > xcross2_force) then
	Famp = 0
else
	Famp = a_force/((x-x0_force)*(x1_force-x))+ b_force
	if (abs(Famp) > 10) then
        write(*,'(a,i6)') 'No of Cell pushing on wall', kcell_debug
		write(*,'(a,6f8.2)') 'Big wall force: R, tube_radius, d, radius, x, Famp: ',R, tube_radius, d, radius, x, Famp
		write(*,'(a,4e12.3)') 'a_force,x0_force,x1_force,b_force: ',a_force,x0_force,x1_force,b_force
		write(nflog,'(a,6f8.2)') 'Big wall force: R, tube_radius, d, radius, x, Famp: ',R, tube_radius, d, radius, x, Famp
		write(nflog,'(a,4e12.3)') 'a_force,x0_force,x1_force,b_force: ',a_force,x0_force,x1_force,b_force
		stop
	endif
endif
! Famp is +ve for repulsion
Famp = -Famp	! to convert to direction of c(:)
if (settling) then
	Famp = settling_factor*Famp
	if (d > radius) then
		Fset = fwall*(radius/d)**2
		Famp = Famp + Fset
	endif
endif
F = Famp*c/R
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_flow_force(cellNo,c,radius,F)
real(REAL_KIND) :: c(3), radius, F(3)
real(REAL_KIND) :: r, v(3), k1, shear
real(REAL_KIND) :: vmax = 0		! no idea what this should be
real(REAL_KIND) :: kvdrag = 0.01	! convert velocity to drag force
real(REAL_KIND) :: mu = 0.003
real(REAL_KIND) :: Delta_P = 80
integer :: cellNo

!c(3) = 0
!r = sqrt(dot_product(c,c))
!vmax = tube_radius**2/(4*mu) * Delta_P/tube_length
!v = vmax*(1-(r/tube_radius)**2) *[0, 0, -1]!cos((r/tube_radius)*(PI/2))*[0, 0, -1]

k1 = -0.01
if (calibration_run) then
	!parametrisation:
	shear = 0.6
	F = k1*shear*[0, 0, 1]
else
	F(1) = k1*(0.0002*u_cell_x(cellNo)-0.0002)
	F(2) = k1*(0.0002*u_cell_y(cellNo)-0.0002)
	F(3) = k1*(0.0002*u_cell_z(cellNo)-0.0002)
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_chemo_force(c,radius,F)
real(REAL_KIND) :: c(3), radius, F(3)
real(REAL_KIND) :: r, d
real(REAL_KIND) :: fchemo1, fchemo2, theta

fchemo1 = grad_amp(1)
fchemo2 = grad_amp(2)
d = c(3)	!tube_length+c(3)
r = sqrt(c(1)*c(1) + c(2)*c(2))
theta = ATAN2(c(2),c(1))
if (theta < 0) then
    theta = theta+2*PI
endif
F = 0
if (chemo(1)%used) then
	F = F + fchemo1 *[0, 0, 1]*d/(tube_length)*EXP(chemo_coef1*r/tube_radius)
endif
if (chemo(2)%used) then
	F = F + fchemo2 *[1, 0, 0]*r*cos(theta)*exp(chemo_coef2*r/tube_radius) + &
		fchemo2 *[0, 1, 0]*r*sin(theta)*exp(chemo_coef2*r/tube_radius)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! This force combines cell-cell repulsion and cell-cell attraction.
! F > 0 = repulsion
! F < 0 = attraction
!-----------------------------------------------------------------------------------------
function get_force(R1,R2,d,incontact,ok) result(F)
real(REAL_KIND) :: R1, R2, d
logical :: incontact, ok
real(REAL_KIND) :: F
real(REAL_KIND) :: x, b_force, dx, delta

dx = x1_force - x0_force
!b_force = -c_force*13 - 4*a_force/dx**2
b_force = -c_force_cell - 4*a_force/dx**2
delta = dx**2 + 4*a_force/b_force
xcross2_force = (x0_force+x1_force)/2 + 0.5*sqrt(delta)

ok = .true.
x = d/(R1+R2)

if (x > xcross2_force) then
	F = 0
	return
endif
if (x < x0_force) then
	write(logmsg,'(a,3e12.3,2f8.3,a,L2)') 'Error: get_force: x < x0: ',R1,R2,d,x,x0_force,' incontact: ',incontact
	call logger(logmsg)
	F = 0
	ok = .false.
!	stop
else
	F = a_force/((x-x0_force)*(x1_force-x)) + b_force
endif
end function

!-----------------------------------------------------------------------------------------
! This force tries to maintain the desired separation distance d^ of the two parts of a
! dividing cell.
! The force F > 0, tending to push the two "spheres" apart, when d < d^, and
! the force F < 0, tending to push them together, when d > d^.
! The original setup used distances in um, so to keep forces the same convert separation
! distance to um.
!-----------------------------------------------------------------------------------------
function get_separating_force(d,d_hat) result(F)
real(REAL_KIND) :: d, d_hat
real(REAL_KIND) :: F
real(REAL_KIND) :: dum

!dum = (d_hat - d)*1.0e4		! cm -> um
dum = (d_hat - d)		! um
F = a_separation*dum**3
end function

!-----------------------------------------------------------------------------------------
! The force is expressed as a function of x = d/(R1+R2)
! where d is the centre-centre distance, R1 and R2 are the two cell radii.
! The cell-cell force parameters are:
! a_force		scales the basic force function 1/((x-x0)(x1-x))
! c_force		magnitude of maximum attractive force, i.e. minimum of function at x = (x0 + x1)/2
! x0_force		location of lower asymptote
! x1_force		location of upper asymptote
! xhat_force	location of zero crossing
!
! F(x) = a/((x-x0)(x1-x)) - c - 4a/(x1-x0)^2
! the minimum value -c occurs at x = (x0+x1)/2
! the relationship between the parameters is:
! xhat = (x0+x1)/2 + (x1-x0)/2.(c(x1-x0)^2/(c(x1-x0)^2 + 4a))^(1/2)
!-----------------------------------------------------------------------------------------
subroutine setup_force_parameters
integer :: i
real(REAL_KIND) :: b_force, dx, delta
real(REAL_KIND) :: epsilon, es_e, shift, sqr_es_e, k_v

write(logmsg,'(a,f6.2)') 'force parameters: a_separation: ',a_separation
call logger(logmsg)
write(logmsg,'(a,4f6.2)') 'force parameters: a_force,c_force_cell,x0_force,x1_force: ',a_force,c_force_cell,x0_force,x1_force
call logger(logmsg)

dx = x1_force - x0_force
b_force = -c_force_cell - 4*a_force/dx**2
delta = dx**2 + 4*a_force/b_force
xcross1_force = (x0_force+x1_force)/2 - 0.5*sqrt(delta)
xcross2_force = (x0_force+x1_force)/2 + 0.5*sqrt(delta)
write(logmsg,*) 'cell-cell force parameters: xcross:',xcross1_force,xcross2_force
call logger(logmsg)

dt_min = 0.01
!delta_min = 0.02e-4		! um -> cm
!delta_max = 0.30e-4		! um -> cm
delta_min = 0.02		! um
delta_max = 0.10		! um
delta_tmove = dt_min
ndt = 5

! For cell-cell force hysteresis
!alpha_v = 0.25
!epsilon = 7.5
!es_e = 1
!shift = -6
!sqr_es_e = sqrt(es_e)
!k_v = 2/alpha_v - sqr_es_e + sqrt(es_e - shift/epsilon)
!k_detach = k_v*alpha_v/2

end subroutine

!--------------------------------------------------------------------------------
! Assume a parabolic velocity distribution = f(r)
! f(0) = Vmax
! f(R) = 0
! df(0)/dr = 0
! ==> f(r) = Vmax(1 - (r/R)^2)
! NOT USED
!--------------------------------------------------------------------------------
real(REAL_KIND) function FlowVelocity(site)
integer :: site(3)
real(REAL_KIND) :: r
integer :: ndist = 2

FlowVelocity = 0
!r = sqrt((site(2)-y0)**2 + (site(3)-z0)**2)
!FlowVelocity = BG_flow_amp*(1 - (r/nradius)**ndist)
end function

!--------------------------------------------------------------------------------
! NOT USED
!--------------------------------------------------------------------------------
subroutine make_wall
real(REAL_KIND) :: R, dtheta, dz, theta, x, y, z
integer :: ncirc, nlong, ic, il
type(cell_type), pointer :: cp

R = tube_radius + Raverage
ncirc = 1.5*2*PI*R/Raverage
nlong = 1.5*tube_length/Raverage + 1
dtheta = 2*PI/ncirc
dz = tube_length/(nlong-1)
nwallcells = 0
do il = 1,nlong
	z = (il-1)*dz
	do ic = 1,ncirc
		theta = (ic-1)*dtheta + mod(il,2)*dtheta/2
		x = R*cos(theta)
		y = R*sin(theta)
		nwallcells = nwallcells + 1
		nlist = nlist+1
		cp => cell_list(nlist)
		cp%centre(:,1) = [x,y,z]
		cp%radius = Raverage
		cp%state = WALL
		cp%Iphase = .true.
		cp%nspheres = 1
	enddo
enddo
ncells = ncells + nwallcells
end subroutine

end module


