module motility

use global
use chemokine

implicit none

contains

#if 0
!-----------------------------------------------------------------------------------------
! This is the version from bcell-abm
!-----------------------------------------------------------------------------------------
subroutine chemo_jumper(kcell,indx1,kslot1,go,kpar)
integer :: kpar,kcell,indx1(2),kslot1
logical :: go
type (cell_type), pointer :: cell
integer :: fullslots1,fullslots2,site1(3),site2(3),kslot2
integer :: irel,dir1,lastdir1,indx2(2),k,kr,rv(3),id, ichemo, nfull, nrest, nout
integer :: savesite2a(3,MAXRELDIR+1), saveslots2a(MAXRELDIR+1)
real(REAL_KIND) :: p(MAXRELDIR+1),psum, R, pR, psumm, stay_prob,  psave(MAXRELDIR+1)
real(REAL_KIND) :: tnow, v(3), vsum(3), f
logical :: ischemo, cognate

tnow = istep*DELTA_T
cell => cell_list(kcell)

id = cell%id
site1 = cell%site
if (site1(1) < 1) then
    write(logmsg,*) 'chemo_jumper: bad site1: ',site1
    call logger(logmsg)
    stop
endif
fullslots1 = 0
do k = 1,2
    if (indx1(k) > 0) then
        fullslots1 = fullslots1 + k
    endif
enddo
stay_prob = dirprob(0)

ischemo = .false.
do kr = 1,MAX_RECEPTOR
    if (cell%receptor_saturation_time(kr) /= 0) then
        if (tnow > cell%receptor_saturation_time(kr) + receptor(kr)%refractory_time) then
            cell%receptor_saturation_time(kr) = 0
!			call logger('Receptor resensitized')
        endif
    endif
	f = cell%receptor_level(kr)*receptor(kr)%strength
    if (receptor(kr)%used .and. (f > 0) .and. (cell%receptor_saturation_time(kr) == 0)) then
        ischemo = .true.
        exit
    endif
enddo

vsum = 0
if (ischemo) then
    do kr = 1,MAX_RECEPTOR
        if (receptor(kr)%used .and. (cell%receptor_saturation_time(kr) == 0)) then
			ichemo = receptor(kr)%chemokine
			f = receptor(kr)%sign*cell%receptor_level(kr)*receptor(kr)%strength
			v = chemo(ichemo)%grad(:,site1(1),site1(2),site1(3))
			if (receptor_saturation(kr,site1,f,v)) then
			    cell%receptor_saturation_time(kr) = tnow
			    f = f/2
!			    call logger('Receptor saturation')
			endif
	    	vsum = vsum + f*v
	    endif
	enddo
	! For exit chemotaxis:
	! Need to create estimate of v() that corresponds to the direction of vsum,
	! nearest discrete location on the 3D lattice (for chemo_p(x,y,z))
	! This is an approximation to increase speed - it enables a table lookup.
	! Note that we use only the direction of v (magnitude is insignificant at this stage,
	! since it is accounted for in f)
    if (norm(vsum) > 0) then
		f = min(1.0,norm(vsum))     ! Note: f is in (0-1)
		v = vsum/norm(vsum)
		rv = chemo_N*v
	else
		f = 0
		ischemo = .false.
	endif
    stay_prob = dirprob(0)
    stay_prob = (1-f)*stay_prob
else
    stay_prob = dirprob(0)
endif

if (fullslots1 /= BOTH) then
    R = par_uni(kpar)
    if (R <= stay_prob) then    ! case of no jump
	    go = .false.
        return
    endif
endif
! Now we must jump (if possible)
! Compute jump probabilities in the absence of chemotaxis
site1 = cell%site
lastdir1 = cell%lastdir
p = 0
savesite2a = 0
saveslots2a = 0
nfull = 0
nrest = 0
nout = 0
do irel = 1,nreldir
	dir1 = reldir(lastdir1,irel)
	site2 = site1 + jumpvec(:,dir1)
	if (inside_xyz(site2)) then
	    indx2 = occupancy(site2(1),site2(2),site2(3))%indx
		if (indx2(1) >= 0) then     ! not OUTSIDE_TAG or DC
            fullslots2 = 0
            do k = 1,2
                if (indx2(k) > 0) then
                    fullslots2 = fullslots2 + k
                endif
            enddo
            if (fullslots2 == BOTH) then
                nfull = nfull + 1
                cycle
            elseif (fullslots2 /= 0) then
                nrest = nrest + 1
                ! Try this!
                p(dir1) = dirprob(irel)*GAMMA
            else
                nrest = nrest + 1
                p(dir1) = dirprob(irel)
            endif
            saveslots2a(dir1) = fullslots2
        else
	        nout = nout + 1
		endif
	endif
	savesite2a(:,dir1) = site2
enddo
if (sum(p) == 0) then
    go = .false.
    return
endif

if (ischemo) then
	psave = p
	call chemo_probs_pre(p,rv,f)     ! this is the precomputed version
endif
psum = sum(p)

if (psum == 0) then
	go = .false.
	return
else
    go = .true.
endif

! Now choose a direction on the basis of these probs p()
R = par_uni(kpar)
pR = psum*R
psumm = 0
do dir1 = 1,njumpdirs
   	psumm = psumm + p(dir1)
   	if (pR <= psumm) then
   		exit
   	endif
enddo
if (dir1 > njumpdirs) then
    dir1 = 0
    do k = 1,njumpdirs
        if (p(k) > 0) then
            dir1 = k
            exit
        endif
    enddo
endif
site2 = savesite2a(:,dir1)

fullslots2 = saveslots2a(dir1)
! new code
if (diagonal_jumps) then
	dir1 = fix_lastdir(dir1,kpar)
elseif (dir1 == 0) then
	dir1 = random_int(1,6,kpar)
endif

if (fullslots2 == 0) then       ! randomly select a slot
    R = par_uni(kpar)
    if (R <= 0.5) then
        kslot2 = SLOT_NUM1
    else
        kslot2 = SLOT_NUM2
    endif
elseif (fullslots2 == SLOT_NUM1) then
    kslot2 = SLOT_NUM2
elseif (fullslots2 == SLOT_NUM2) then
    kslot2 = SLOT_NUM1
else
    write(logmsg,*) 'ERROR in jumper: jump to crowded site'
	call logger(logmsg)
    stop
endif
cell%site = site2
cell%lastdir = dir1
occupancy(site2(1),site2(2),site2(3))%indx(kslot2) = kcell
occupancy(site1(1),site1(2),site1(3))%indx(kslot1) = 0
end subroutine

#endif

!-----------------------------------------------------------------------------------------
! Determine whether a receptor at site1(:) is saturated.
! Currently the decision is made from:
! f = total strength, which is cell%receptor_level(kr)*receptor(kr)%strength
! v = chemokine gradient
!-----------------------------------------------------------------------------------------
logical function receptor_saturation(kr,site,f,v)
integer :: kr, site(3)
real(REAL_KIND) :: f, v(3), a(3)
real(REAL_KIND) :: magnitude
!real(REAL_KIND) :: saturation_threshold = 0.8

! For now, turn off receptor saturation
receptor_saturation = .false.
return

a = f*v
magnitude = norm(a)
if (magnitude > receptor(kr)%saturation_threshold) then
    receptor_saturation = .true.
else
    receptor_saturation = .false.
endif
end function

!-----------------------------------------------------------------------------------------
! Not parallel
!-----------------------------------------------------------------------------------------
subroutine mover(ok)
logical :: ok
type (cell_type), pointer :: cell
integer :: kcell	!, site(3),indx(2),slot
integer :: kpar=0

ok = .true.
!return
! Randomise sequence!!!!!
do kcell = 1,nlist
    cell => cell_list(kcell)
    if (cell%ID == 0) cycle             ! skip gaps in the list
	call jumper(kcell,kpar)
enddo

end subroutine

!--------------------------------------------------------------------------------
! A cell that can move, i.e. has a vacant neighbour site, is subject to four
! influences:
! 1. Drag effect from the flow velocity at that location.
! 2. Chemotactic attraction.
! 3. Cell-cell adhesion
! 4. Cell-vessel wall adhesion
!--------------------------------------------------------------------------------
subroutine jumper(kcell,kpar)
integer :: kpar, kcell
logical :: go
integer :: site1(3)
integer :: k, k2, site2(3), site3(3)
real(REAL_KIND) :: vx, chemo_attraction, upvalue, downvalue, q(27), p(27), qsum, psum, R
real(REAL_KIND) :: qsum0, factor
logical :: wall

chemo_attraction = grad_amp(1)
site1 = cell_list(kcell)%site
if (site1(1) == 1 .or. site1(1) == NX) then		! lose the cell
	cell_list(kcell)%ID = 0
	occupancy(site1(1),site1(2),site1(3))%indx(1) = 0
	NTcells = NTcells - 1
	return
endif
vx = FlowVelocity(site1)
downvalue = max((kdrag*vx - chemo_attraction), 0.0)
upvalue = max((-kdrag*vx + chemo_attraction), 0.0)
wall = .false.
qsum0 = 0
do k = 1,njumpdirs
	site2 = site1 + jumpvec(:,k)
	if (k == 14) cycle
	if (site2(1) < 1 .or. site2(1) > NX) cycle
	if (occupancy(site2(1),site2(2),site2(3))%indx(1) == 0) cycle	! no cell or wall
	if (occupancy(site2(1),site2(2),site2(3))%indx(1) == OUTSIDE_TAG) then
		wall = .true.
	endif
	if (occupancy(site2(1),site2(2),site2(3))%indx(1) == OUTSIDE_TAG) then
		factor = 2*Kadhesion
	else
		factor = Kadhesion
	endif
	qsum0 = qsum0 + 1/jumpdist(k)**2
enddo
q = 0;
do k = 1,njumpdirs
	site2 = site1 + jumpvec(:,k)
	if (occupancy(site2(1),site2(2),site2(3))%indx(1) /= 0) cycle
	! Determine proximity attractiveness of destination site
	! by looking at the neighbours.  The contribution from a site
	! that contains a cell or is wall (OUTSIDE TAG) is the inverse
	! of the distance.  We need to exclude the current cell site,
	! since this will be vacant if the cell moves.
	qsum = 0
	do k2 = 1,njumpdirs
		site3 = site2 + jumpvec(:,k2)
		if (k2 == 14) cycle
		if (site1(1)==site3(1) .and. site1(2)==site3(2) .and. site1(3)==site3(3)) cycle
		if (site3(1) < 1 .or. site3(1) > NX) cycle
		if (occupancy(site3(1),site3(2),site3(3))%indx(1) == 0) cycle	! no cell or wall
		if (occupancy(site3(1),site3(2),site3(3))%indx(1) == OUTSIDE_TAG) then
			factor = 2*Kadhesion
		else
			factor = Kadhesion
		endif
		qsum = qsum + 1/jumpdist(k2)**2
	enddo
	if (k == 14) then	
		cycle
	else
		q(k) = q(k) + qsum
	endif
	! Determine attractiveness of upstream or downstream jump
	if (jumpvec(1,k) > 0) then
		q(k) = q(k) + downvalue
	elseif (jumpvec(1,k) < 0) then
		q(k) = q(k) + upvalue
	endif
	q(k) = max(q(k) - Kstay*qsum0, 0.0)
	if (wall .and. jumpvec(1,k) < 0) q(k) = 1.2*q(k)
enddo
qsum = sum(q)
if (qsum == 0) return

p = q/qsum
R = par_uni(kpar)
!write(*,'(10f7.3)') p
!write(*,*) 'R: sum(p) ',R,sum(p)
psum = 0
do k = 1,njumpdirs
   	psum = psum + p(k)
   	if (R <= psum) then
   		exit
   	endif
enddo
!if (.not.wall .and. jumpvec(1,k) < 0) return
site2 = site1 + jumpvec(:,k)
cell_list(kcell)%site = site2
occupancy(site1(1),site1(2),site1(3))%indx(1) = 0
occupancy(site2(1),site2(2),site2(3))%indx(1) = kcell
end subroutine

!--------------------------------------------------------------------------------
! Assume a parabolic velocity distribution = f(r)
! f(0) = Vmax
! f(R) = 0
! df(0)/dr = 0
! ==> f(r) = Vmax(1 - (r/R)^2)
!--------------------------------------------------------------------------------
real(REAL_KIND) function FlowVelocity(site)
integer :: site(3)
real(REAL_KIND) :: r
integer :: ndist = 2

r = sqrt((site(2)-y0)**2 + (site(3)-z0)**2)
FlowVelocity = BG_flow_amp*(1 - (r/nradius)**ndist)
end function

!--------------------------------------------------------------------------------
! 2D motility
! In this case only one cell can occupy a site - no slot 2 - and no passing.
! This is valid because cells are not crammed together in the petri dish.
! We now wish to relax this constraint.  A site might be crowded, as a result
! of cell division, but in jumper2D cells can move only to vacant sites.
! Only near a DC are there sites available with z > 1.
!--------------------------------------------------------------------------------
subroutine jumper2D(kcell,indx1,kslot1,go,kpar)
integer :: kpar,kcell,indx1(2),kslot1
logical :: go
type (cell_type), pointer :: cell
integer :: site1(3),site2(3)
integer :: irel,dir1,lastdir1,z,irel_z,nslots
integer :: savesite2(3,26)
real(REAL_KIND) :: psum, p(8), R, pfac
!real(REAL_KIND), parameter :: upfactor(2) = (/0.5, 0.3/)		! reducing prob of jumps up by 1, 2 sites
!real(REAL_KIND), parameter :: downfactor(2) = (/1.3, 1.5/)		! increasing prob of jumps down by 1, 2 sites

cell => cell_list(kcell)
site1 = cell%site
if (indx1(1) > 0 .and. indx1(2) > 0) then
	nslots = 2
else
	nslots = 1
endif
if (nslots == 1) then
	R = par_uni(kpar)
	if (R <= dirprob2D(0)) then    ! case of no jump
		go = .false.
		return
	endif
endif

! Now we must jump (if possible)

lastdir1 = cell%lastdir
irel_z = nreldir2D
p = 0
psum = 0
!if (site1(3) > 1) then
!	z = site1(3) - 1
!	if (occupancy(site1(1),site1(2),z)%indx(1) == 0) then
!		if (zrange2D(site1(1),site1(2),1) <= z) then
!!			write(*,*) 'in space: ',kcell,site1
!			irel_z = irel_z + 1		! drop cell into vacant site below with high probability
!			p(irel_z) = 999
!			psum = psum + p(irel_z)
!			savesite2(:,irel_z) = (/ site1(1),site1(2),z /)
!		endif
!	endif
!endif
!do irel = 1,nreldir2D
!    p(irel) = 0
!	dir1 = reldir2D(lastdir1,irel)
!	site2 = site1 + jumpvec2D(:,dir1)
!	if (inside_xy(site2)) then
!		! Look for the lowest free site
!		if (zrange2D(site2(1),site2(2),1) == 0) cycle	! outside site
!		do z = zrange2D(site2(1),site2(2),1),zrange2D(site2(1),site2(2),2)
!			if (occupancy(site2(1),site2(2),z)%indx(1) == 0) then
!				pfac = 1.0
!				if (z > site1(3)) then
!					pfac = upfactor(z-site1(3))
!				elseif (z < site1(3)) then
!					pfac = downfactor(site1(3)-z)
!				endif
!				p(irel) = pfac*dirprob2D(irel)
!				psum = psum + p(irel)
!				savesite2(:,irel) = (/ site2(1),site2(2),z /)
!				exit
!			endif
!		enddo
!	endif
!enddo

!if (kcell == 1) then
!	write(*,*) 'jumper2D: lastdir1: ',lastdir1
!endif
do irel = 1,nreldir2D
    p(irel) = 0
	dir1 = reldir2D(lastdir1,irel)
!	if (kcell == 1) write(*,*) irel,dir1
	site2 = site1 + jumpvec2D(:,dir1)
	if (inside_xy(site2)) then
		if (occupancy(site2(1),site2(2),1)%indx(1) == 0) then
			p(irel) = dirprob2D(irel)
			psum = psum + p(irel)
			savesite2(:,irel) = [ site2(1), site2(2), 1 ]
		endif
	endif
enddo

if (psum == 0) then
	go = .false.
	cell%lastdir = random_int(1,8,kpar)
	return
else
    go = .true.
endif

!if (kcell == 1) write(*,'(a,8f8.4)') 'p: ',p

! Now choose a direction on the basis of these probs p()
R = par_uni(kpar)
R = psum*R
psum = 0
do irel = 1,irel_z
   	psum = psum + p(irel)
   	if (R <= psum) then
   		exit
   	endif
enddo
if (irel > irel_z) then
	if (irel_z == nreldir2D) then
		go = .false.
		return
	else
		irel = irel_z
	endif
endif
site2 = savesite2(:,irel)
if (irel <= nreldir2D) then
	dir1 = reldir2D(lastdir1,irel)
else	! note that the direction of a stack jump is lost here.
	dir1 = 0
endif
if (dir1 == 0) then
	dir1 = random_int(1,8,kpar)
endif

cell%site = site2
cell%lastdir = dir1
occupancy(site2(1),site2(2),site2(3))%indx(1) = kcell
if (nslots == 1) then
	if (indx1(2) > 0) then
		write(*,*) 'Error: jumper2D: cell in slot 2 only: ',site1,indx1
		stop
	endif
	occupancy(site1(1),site1(2),site1(3))%indx(1) = 0
elseif (indx1(2) == kcell) then
	occupancy(site1(1),site1(2),site1(3))%indx(2) = 0
else
	occupancy(site1(1),site1(2),site1(3))%indx(1) = indx1(2)
	occupancy(site1(1),site1(2),site1(3))%indx(2) = 0
endif
end subroutine

!-----------------------------------------------------------------------------------------
! This is the version from bcell-abm
!-----------------------------------------------------------------------------------------
subroutine chemo_jumper2D(kcell,indx1,kslot1,go,kpar)
integer :: kpar,kcell,indx1(2),kslot1
logical :: go
type (cell_type), pointer :: cell
integer :: fullslots1,fullslots2,site1(3),site2(3),kslot2
integer :: irel,dir1,lastdir1,indx2(2),k,kr,rv(3),id, ichemo, nfull, nrest, nout
integer :: savesite2a(3,MAXRELDIR2D+1), saveslots2a(MAXRELDIR2D+1)
real(REAL_KIND) :: p(MAXRELDIR2D+1),psum, R, pR, psumm, stay_prob,  psave(MAXRELDIR2D+1)
real(REAL_KIND) :: tnow, v(3), vsum(3), f
logical :: ischemo

GAMMA = 0.5
tnow = istep*DELTA_T
cell => cell_list(kcell)

id = cell%id
site1 = cell%site
if (site1(1) < 1) then
    write(logmsg,*) 'chemo_jumper: bad site1: ',site1
    call logger(logmsg)
    stop
endif
fullslots1 = 0
do k = 1,2
    if (indx1(k) > 0) then
        fullslots1 = fullslots1 + k
    endif
enddo
stay_prob = dirprob2D(0)

ischemo = .false.
do kr = 1,MAX_RECEPTOR
    if (cell%receptor_saturation_time(kr) /= 0) then
        if (tnow > cell%receptor_saturation_time(kr) + receptor(kr)%refractory_time) then
            cell%receptor_saturation_time(kr) = 0
!			call logger('Receptor resensitized')
        endif
    endif
	f = cell%receptor_level(kr)*receptor(kr)%strength
    if (receptor(kr)%used .and. (f > 0) .and. (cell%receptor_saturation_time(kr) == 0)) then
        ischemo = .true.
        exit
    endif
enddo

vsum = 0
if (ischemo) then
    do kr = 1,MAX_RECEPTOR
        if (receptor(kr)%used .and. (cell%receptor_saturation_time(kr) == 0)) then
			ichemo = receptor(kr)%chemokine
			f = receptor(kr)%sign*cell%receptor_level(kr)*receptor(kr)%strength

			v = chemo(ichemo)%grad(:,site1(1),site1(2),site1(3))
!			v = [1,0,0]	! need the chemo field gradient

			if (receptor_saturation(kr,site1,f,v)) then
			    cell%receptor_saturation_time(kr) = tnow
			    f = f/2
!			    call logger('Receptor saturation')
			endif
	    	vsum = vsum + f*v
!	    	write(*,'(2i4,7f8.4)') kcell,kr,f,v,vsum
	    endif
	enddo
	! Note that we use only the direction of v (magnitude is insignificant at this stage,
	! since it is accounted for in f)
    if (norm(vsum) > 0) then
		f = min(1.0,norm(vsum))     ! Note: f is in (0-1)
		v = vsum/norm(vsum)
		rv = chemo_N*v
!		if (kcell == 1) write(*,'(a,4f8.3,3i4)') 'f,v,rv: ',f,v,rv
	else
		f = 0
		ischemo = .false.
	endif
    stay_prob = dirprob2D(0)
    stay_prob = (1-f)*stay_prob
else
    stay_prob = dirprob2D(0)
endif

if (fullslots1 /= BOTH) then
    R = par_uni(kpar)
    if (R <= stay_prob) then    ! case of no jump
	    go = .false.
        return
    endif
endif

! Now we must jump (if possible)
! Compute jump probabilities in the absence of chemotaxis
site1 = cell%site
lastdir1 = cell%lastdir
!write(*,*) 'kcell, lastdir1: ',kcell,lastdir1
p = 0
savesite2a = 0
saveslots2a = 0
nfull = 0
nrest = 0
nout = 0
do irel = 1,nreldir2D
	dir1 = reldir2D(lastdir1,irel)
	site2 = site1 + jumpvec2D(:,dir1)
	if (inside_xy(site2)) then
	    indx2 = occupancy(site2(1),site2(2),site2(3))%indx
		if (indx2(1) >= 0) then     ! not OUTSIDE_TAG
            fullslots2 = 0
            do k = 1,2
                if (indx2(k) > 0) then
                    fullslots2 = fullslots2 + k
                endif
            enddo
            if (fullslots2 == BOTH) then
                nfull = nfull + 1
                cycle
            elseif (fullslots2 /= 0) then
                nrest = nrest + 1
                p(dir1) = dirprob2D(irel)*GAMMA
            else
                nrest = nrest + 1
                p(dir1) = dirprob2D(irel)
            endif
            saveslots2a(dir1) = fullslots2
        else
	        nout = nout + 1
		endif
	endif
	savesite2a(:,dir1) = site2
enddo
if (sum(p) == 0) then
    go = .false.
    return
endif
!write(*,'(a,8f8.4)') 'before: ',p/sum(p)

if (ischemo) then
	psave = p
	call chemo_probs_pre2D(p,rv,f)     ! this is the precomputed version
endif
psum = sum(p)
!write(*,'(a,8f8.4)') 'after: ',p
!stop

if (psum == 0) then
	go = .false.
	return
else
    go = .true.
endif

! Now choose a direction on the basis of these probs p()
R = par_uni(kpar)
pR = psum*R
psumm = 0
do dir1 = 1,njumpdirs2D
   	psumm = psumm + p(dir1)
   	if (pR <= psumm) then
   		exit
   	endif
enddo
if (dir1 > njumpdirs2D) then
    dir1 = 0
    do k = 1,njumpdirs2D
        if (p(k) > 0) then
            dir1 = k
            exit
        endif
    enddo
endif
site2 = savesite2a(:,dir1)
fullslots2 = saveslots2a(dir1)

if (dir1 == 0) then
	dir1 = random_int(1,8,kpar)
endif

if (fullslots2 == 0) then       ! randomly select a slot
    R = par_uni(kpar)
    if (R <= 0.5) then
        kslot2 = SLOT_NUM1
    else
        kslot2 = SLOT_NUM2
    endif
elseif (fullslots2 == SLOT_NUM1) then
    kslot2 = SLOT_NUM2
elseif (fullslots2 == SLOT_NUM2) then
    kslot2 = SLOT_NUM1
else
    write(logmsg,*) 'ERROR in jumper2D: jump to crowded site'
	call logger(logmsg)
    stop
endif
cell%site = site2
cell%lastdir = dir1
occupancy(site2(1),site2(2),site2(3))%indx(kslot2) = kcell
occupancy(site1(1),site1(2),site1(3))%indx(kslot1) = 0
end subroutine

!--------------------------------------------------------------------------------
! This code follows the flow chart.
! (note that there was a mistake in the chart - the P0 case goes to "Choose free slot at S1")
!--------------------------------------------------------------------------------
subroutine flow_jumper2D(kcell,indx1,kslot1,go,kpar)
integer :: kcell, indx1(2), kslot1, kpar
logical :: go
integer :: indx(2), site1(3), site2(3), jrel
integer :: kslot2, irel, j, k, nv(3)
logical :: in, try_S1
real(REAL_KIND) :: vx, vy, vz, v(3), d(3), Arad, u1(3), u2(3), p(0:2), R, psum
type (cell_type), pointer :: cell

cell => cell_list(kcell)

Arad = 0
vx = BG_flow_amp*cos(Arad)/DELTA_X
vy = BG_flow_amp*sin(Arad)/DELTA_X
vz = 0
v = (/vx,vy,vz/)
! At this point v is the velocity in gridcells/min
! to get the displacement in gridcells, need to multiply by DELTA_T
d = v*DELTA_T	! note that Fortran allows operations on vectors (more generally on arrays)
! then find the integer jumps
nv = int(d)	
d = d - nv		! this is the residual displacement

! Determine intermediate site S1.  Note that we do not jump here unless S2 is not available
site1 = cell%site + nv
if (nv(1) /= 0 .or. nv(2) /= 0) then
	in = inside_xy(site1)	! This is only to adjust site1 if wrapping is on.  At this stage we do not need anything else.
endif

! Determine bracketting jump directions for residual d
Arad = atan2(d(2),d(1))
if (Arad < 0) Arad = Arad + 2*PI	! the range of atan2 is (-PI, PI)

jrel = floor(1 + 4*Arad/PI)

if (jrel == 8) then
    u1 = jumpvec2D(:,jrel)
    u2 = jumpvec2D(:,1)
else
    u1 = jumpvec2D(:,jrel)
    u2 = jumpvec2D(:,jrel+1)
endif

p(1) = (d(2)*u2(1)- d(1)*u2(2))/(u1(2)*u2(1)-u1(1)*u2(2))
p(2) = (d(2)*u1(1)- d(1)*u1(2))/(u1(1)*u2(2)-u1(2)*u2(1))
p(0) = 1 - (p(1)+p(2))

R = par_uni(kpar)
psum = 0
do irel = 0,2
   	psum = psum + p(irel)
   	if (R <= psum) then
   		exit
   	endif
enddo

try_S1 = .false.
if (irel == 0) then
	try_S1 = .true.
elseif (irel == 1) then
    site2 = site1 + u1
else
    site2 = site1 + u2
endif

if (.not.try_S1) then
	in = inside_xy(site2)
	if (in) then
		! choose a free slot at S2
	    indx = occupancy(site2(1),site2(2),site2(3))%indx
		if (indx(1) == 0) then
			kslot2 = 1
		elseif (indx(2) == 0) then
			kslot2 = 2
		else
			try_S1 = .true.		! no free slot
		endif
	else
		try_S1 = .true.
	endif
endif

if (try_S1) then
	! choose a free slot at S1
    indx = occupancy(site1(1),site1(2),site1(3))%indx
	if (indx(1) == 0) then
		kslot2 = 1
	elseif (indx(2) == 0) then
		kslot2 = 2
	else
		go = .false.		! no free slot, stay at S0
		return
	endif
	site2 = site1
else
	cell%dtotal(1:2) = cell%dtotal(1:2) + jumpvec2D(1:2,irel) 
endif
cell%dtotal = cell%dtotal + nv

! move to kslot2 at site2
cell%site = site2
occupancy(site2(1),site2(2),site2(3))%indx(kslot2) = kcell
occupancy(site1(1),site1(2),site1(3))%indx(kslot1) = 0
end subroutine


!--------------------------------------------------------------------------------
! Crude check for a site inside, just looking at the allowable ranges of x, y and z.
!--------------------------------------------------------------------------------
logical function inside_xyz(site)
integer :: site(3)

if (site(1) < 1 .or. site(1) > NX .or. site(2) < 1 .or. site(2) > NY .or. site(3) < 1 .or. site(3) > NZ) then
    inside_xyz = .false.
else
    inside_xyz = .true.
endif
end function

!--------------------------------------------------------------------------------
! Crude check for a site inside 2D, just looking at the allowable ranges of x, y.
! If use_wrapping, site() is modified to be inside, and .true. is returned.
!--------------------------------------------------------------------------------
logical function inside_xy(site)
integer :: site(3)

if (use_wrapping) then
	if (site(1) < 1) then
		site(1) = site(1) + NX
	elseif (site(1) > NX) then
		site(1) = site(1) - NX
	endif
	if (site(2) < 1) then
		site(2) = site(2) + NY
	elseif (site(2) > NY) then
		site(2) = site(2) - NY
	endif
	inside_xy = .true.
else
	if (site(1) < 1 .or. site(1) > NX .or. site(2) < 1 .or. site(2) > NY) then
		inside_xy = .false.
	else
		inside_xy = .true.
	endif
endif
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checker2D
integer, allocatable :: xcnt(:), ycnt(:)
integer :: kcell, site(3)

allocate(xcnt(NX))
allocate(ycnt(NY))
xcnt = 0
ycnt = 0
do kcell = 1,nlist
	if (cell_list(kcell)%ID == 0) cycle
	site = cell_list(kcell)%site
	xcnt(site(1)) = xcnt(site(1)) + 1
	ycnt(site(2)) = ycnt(site(2)) + 1
enddo
write(*,*) 'xcnt'
write(*,'(20i4)') xcnt
write(*,*) 'ycnt'
write(*,'(20i4)') ycnt
deallocate(xcnt)
deallocate(ycnt)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine compute_Cm(tagsite,ntagged,nvar0,nvar,dt,Cm)
integer :: tagsite(3,ntagged,0:nvar)
integer :: ntagged,nvar0,nvar
real(REAL_KIND) :: dt,Cm
integer :: k,j,d(3),d2,d2sum
real(REAL_KIND), allocatable :: r2mean(:)

allocate(r2mean(nvar))

do k = 1,nvar
    d2sum = 0
    do j = 1,ntagged
        d = tagsite(:,j,k) - tagsite(:,j,0)
        d2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
        d2sum = d2sum + d2
    enddo
    r2mean(k) = d2sum/ntagged
enddo
write(*,'(10f6.1)') r2mean
call bestfit(r2mean,nvar0,nvar,dt,Cm)
Cm = Cm*DELTA_X*DELTA_X
!write(*,*) 'Cm: ',Cm

deallocate(r2mean)

end subroutine

!--------------------------------------------------------------------------------------
! Use least squares to find the best straight line fit to r2mean() from n1 to n2.
! In fact we might need to reduce the range of points.
! We need to use a range within which the slope (dr2/dt) doesn't vary too much from
! an "average" value.  Note that the larger the persistence parameter rho the longer
! it takes for the r2 plot to become linear.
!--------------------------------------------------------------------------------------
subroutine bestfit(r2mean,n1,n2,dt,Cm)
real(REAL_KIND) :: r2mean(:),dt,Cm
integer :: n1,n2
integer :: n,i
real(REAL_KIND) :: sumt,sumt2,sumy,sumyt,t,y,a,b

n = n2 - n1 + 1		! number of data points

sumt = 0
sumt2 = 0
sumy = 0
sumyt = 0
do i = n1,n2
	t = i*dt
	y = r2mean(i)
	sumt = sumt + t
	sumt2 = sumt2 + t*t
	sumy = sumy + y
	sumyt = sumyt + y*t
enddo
a = (n*sumyt - sumy*sumt)/(n*sumt2 - sumt*sumt)		! slope of line
b = (sumy - a*sumt)/n								! intercept
if (SIMULATE_2D) then
	Cm = a/4
else
	Cm = a/6
endif
write(*,*) 'a,b: ',a,b

end subroutine

!---------------------------------------------------------------------
! Need to precompute array reldir(:,:)
! FOR NRELDIR = 6 (diagonal_jumps = .false.)
! The directions 1 - 6 are defined according to the neumann array,
! i.e. 1 is -x, 2 is +x, 3 is -y, 4 is +y, 5 is -z, 6 is +z
! The relative directions are numbered as follows:
! irel = 1 gives the same direction as lastdir
! irel = 6 gives the opposite direction to lastdir
! irel = 2,3,4,5 cover the remaining directions - order not important
! at the moment since all directions normal to lastdir are equally likely
! FOR NRELDIR = 17 (diagonal_jumps = .true.)
! The reldir(:,:) array uses the same set of 6 previous jump directions,
! restricted to the axes.  Diagonal previous jumps are accommodated by
! making a random selection of an axis direction from either 2 or 3
! possibilities (depending on whether it was a D2 or D3 jump).
! The possible jumps, now including diagonal moves, are stored in jumpvec(:)
! The jumpvec(:) entries are ordered with z varying fastest, then y, then x,
! each ranging -1, 0, +1.
! For dirprob(:) calculation:
! There are 3 groups of jumps, with sets of probs
! Group 1: Q(1) = P(1), D1 jump in the same as last direction
! Group 2: Q(2) = P(2)+P(3)+P(4)+P(5) (D2 jumps) + P(6)+P(7)+P(8)+P(9) (D3 jumps)
! Group 3: Q(3) = P(10)+P(11)+P(12)+P(13) (D1 jumps) + P(14)+P(15)+P(16)+P(17) (D2 jumps)
! i.e. Q(1) = prob of no direction change, Q(2) = prob of direction change < 90 deg,
! Q(3) = prob of direction change = 90 deg.
! Setting D1 jumps (von Neumann) to 1, the D2 jumps have length L2 = sqrt(2)
! and the D3 jumps have length L3 = sqrt(3)
! Therefore jump probs must be scaled appropriately to give symmetry within a set
! In other words:
!		P(i)/P(j) = L2/L3 for j=6,7,8,9 j=2,3,4,5
! and	P(i)/P(j) = 1/L2  for i=14,15,16,17 j=10,11,12,13
! Then choosing P(2) to represent D2 jump prob in Group 2, and P(10) to
! represent D1 jump prob in Group 3:
!		Q(2) = 4(1 + L2/L3).P(2)
!		Q(3) = 4(1 + 1/L2).P(10)
! It is always the case that Q(1) + Q(2) + Q(3) = BETA
!
! When RHO = 1, Q(1) = BETA = p1, Q(2) + Q(3) = 0
!
! When RHO = 0, P(10) = P(1), and P(2) = P(1)/L2
! therefore:	Q(3) = 4(1 + 1/L2).P(1)
! and			Q(2) = 4(1 + L2/L3).P(1)/L2
! giving:		P(1).[1 + 4(1 + L2/L3)/L2 + 4(1 + 1/L2)] = BETA
!				P(1) = BETA/[1 + 4(1 + L2/L3)/L2 + 4(1 + 1/L2)] = p0
! and:			Q(3)/Q(2) = (1 + L2)/(1 + L2/L3) = q320
!
! Now for 0 < RHO < 1
! first set Q(1) = P(1) to vary linearly with RHO from
! the value at 0 (p0) to the value at 1 (p1)
!	P(1) = p0 + RHO*(p1-p0)
! Now we know that Q(2) + Q(3) = BETA - Q(1)
! and the third parameter ALPHA is used to determine how Q(3)/Q(2)
! varies with RHO, by making it change linearly from q23 at RHO=0
! to ALPHA*q23 at RHO = 1
! Now we have all the info to solve for Q(1), Q(2), Q(3) and all P(:)
!
! Revised Moore model.
! Now allow reverse directions, but disallow D3 diagonal jumps, giving
! 18 possible directions.
! In the preferred direction the surface is a prolate spheroid, with
! a = 1, b = 1-RHO.
! In the reverse direction the surface is either a sphere (a = b) or
! an oblate spheroid (a = b^2).
!
! SIMULATE_2D case - 2D motion
! This case is simple, because the 8 jump directions can be ordered
! anti-clockwise (i.e. increasing theta), with the x-axis dir = 1.
! The relative directions are then found by starting at a different
! point in the sequence 1,2,..,8.
!---------------------------------------------------------------------
subroutine make_reldir
integer :: lastdir,irel,k,ix,iy,iz,i,site(3)
integer :: reldir18(6,18),reldir26(6,26)

if (SIMULATE_2D) then
	njumpdirs2D = 8
	nreldir2D = 8
!	do k = 1,njumpdirs2D
!		write(*,'(a,i4,3f8.4)') 'jumpvec2D: ',k,jumpvec2D(:,k)
!	enddo
	do lastdir = 1,nreldir2D
		do k = 1,njumpdirs2D
			irel = lastdir + k-1
			if (irel > 8) irel = irel - 8
			reldir2D(lastdir,k) = irel
		enddo
	enddo
elseif (MODEL == NEUMANN_MODEL) then
	diagonal_jumps = .false.
	njumpdirs = 6
	nreldir = 6
	reldir = 0
	do lastdir = 1,6
		reldir(lastdir,1) = lastdir
		if (mod(lastdir,2) == 0) then
			reldir(lastdir,6) = lastdir - 1
		else
			reldir(lastdir,6) = lastdir + 1
		endif
		irel = 2
		do k = 1,6
			if (k /= reldir(lastdir,1) .and. k /= reldir(lastdir,6)) then
				reldir(lastdir,irel) = k
				irel = irel + 1
			endif
		enddo
	enddo
	write(*,*) 'reldir'
    do lastdir = 1,6
	    write(*,'(6i6)') (reldir(lastdir,k),k=1,6)
    enddo
	jumpvec(:,1:6) = neumann(:,1:6)
	jumpdist(1:6) = 1
else
	if (MODEL == MOORE18_MODEL) then
    	nreldir = 18
	elseif (MODEL == MOORE26_MODEL) then
    	nreldir = 26
    endif
    njumpdirs = 27
	diagonal_jumps = .true.
	k = 0
	do ix = -1,1
		do iy = -1,1
			do iz = -1,1
				k = k+1
				jumpvec(:,k) = (/ix,iy,iz/)
				jumpdist(k) = sqrt(real(ix**2 + iy**2 + iz**2))
			enddo
		enddo
	enddo

! Data for revised Moore18 model.  D3 jumps are excluded, but reverse directions are allowed.
	reldir18(1,:) = (/  5,  2, 4, 6, 8, 11,13,15,17, 10,12,16,18, 22,24,26,20, 23 /)	! -x
	reldir18(2,:) = (/ 23, 20,22,24,26, 11,13,15,17, 10,12,16,18,  2, 4, 6, 8,  5 /)	! +x
	reldir18(3,:) = (/ 11,  2,10,12,20,  5,13,15,23,  4, 6,22,24,  8,16,18,26, 17 /)	! -y
	reldir18(4,:) = (/ 17,  8,16,18,26,  5,13,15,23,  4, 6,22,24,  2,10,12,20, 11/)	! +y
	reldir18(5,:) = (/ 13,  4,10,16,22,  5,11,17,23,  2, 8,20,26,  6,12,18,24, 15 /)	! -z
	reldir18(6,:) = (/ 15,  6,12,18,24,  5,11,17,23,  2, 8,20,26,  4,10,16,22, 13 /)	! +z

! Data for revised Moore26 model.  D3 jumps are excluded, but reverse directions are allowed.
	reldir26(1,:) = (/  5,  2, 4, 6, 8,  1, 3, 7, 9, 11,13,15,17, 10,12,16,18, 19,21,27,25, 20,22,24,26, 23 /)	! -x
	reldir26(2,:) = (/ 23, 20,22,24,26, 19,21,27,25, 11,13,15,17, 10,12,16,18,  1, 3, 7, 9,  2, 4, 6, 8,  5 /)	! +x
	reldir26(3,:) = (/ 11,  2,10,12,20,  1, 3,19,21,  5,13,15,23,  4, 6,22,24,  7, 9,25,27,  8,16,18,26, 17 /)	! -y
	reldir26(4,:) = (/ 17,  8,16,18,26,  7, 9,25,27,  5,13,15,23,  4, 6,22,24,  1, 3,19,21,  2,10,12,20, 11 /)	! +y
	reldir26(5,:) = (/ 13,  4,10,16,22,  1, 7,19,25,  5,11,17,23,  2, 8,20,26,  3, 9,21,27,  6,12,18,24, 15 /)	! -z
	reldir26(6,:) = (/ 15,  6,12,18,24,  3, 9,21,27,  5,11,17,23,  2, 8,20,26,  1, 7,19,25,  4,10,16,22, 13 /)	! +z

    if (MODEL == MOORE18_MODEL) then
        reldir(1:6,1:18) = reldir18
    elseif (MODEL == MOORE26_MODEL) then
        reldir = reldir26
    endif

endif
!call compute_dirprobs
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine compute_dirprobs
integer :: k
real(REAL_KIND) :: L2,L3,b,c,e,d(MAXRELDIR),R(MAXRELDIR)
real(REAL_KIND) :: theta, p(8), psum

if (SIMULATE_2D) then
	! Consider an ellipse with eccentricity = e, with major radius = a, minor radius = b,
	! and e = sqrt((a^2 - b^2)/a^2)  The distances from the negative focus at -ea to
	! the ellipse, in the jump directions, are proportional to the jump probabilities,
	! suitably weighted to account for the jump distances (1 or sqrt(2)).
	! Let p(i) = prob. of jump i
	!     d(i) = length of jump i
	!     r(i) = distance from focus to the ellipse in direction i = a(1-e^2)/(1-e.cos(theta))
	! then p(i) = k.r(i)/d(i)

	dirprob2D(0) = 1 - BETA	! this is the prob of no jump
	e = RHO
	psum = 0
	do k = 1,8
		theta = (k-1)*PI/4
		p(k) = (1 - e*e)/(1 - e*cos(theta))
		if (mod(k,2) == 0) then
			p(k) = p(k)/sqrt(2.0d0)
		endif
		psum = psum + p(k)
	enddo
	dirprob2D(1:8) = BETA*p/psum
!	write(*,'(a,2f8.3)') 'compute_dirprobs: BETA, RHO: ',BETA,RHO
!	write(*,'(a,9f7.4)') 'dirprob2D: ',dirprob2D
	return
endif
if (MODEL == NEUMANN_MODEL) then
	diagonal_jumps = .false.
	nreldir = 6
! New Neumann model (with reverses)
	b = 1-RHO
	b = b*b
	e = BETA/(1 + 4*b + b**2)
	dirprob(0) = 1 - BETA	! this is the prob of no jump
	dirprob(1) = e
	dirprob(2:5) = e*b
	dirprob(6) = e*b**2
else
    ! New prolate spheroid approach
	diagonal_jumps = .true.
    if (MODEL == MOORE18_MODEL) then
	    nreldir = 18
	    L2 = sqrt(2.0d0)
	    b = 1 - RHO
	    b = b*b		! revised version
	    c = b*b		! oblate spheroid on reverse side
	    d(1) = 1
	    R(1) = 1
	    d(2:5) = L2
	    R(2:5) = L2*b/sqrt(b*b+1)
	    d(6:9) = 1
	    R(6:9) = b
	    d(10:13) = L2
	    R(10:13) = b
	    d(14:17) = L2
	    R(14:17) = L2*b*c/sqrt(b*b+c*c)
	    d(18) = 1
	    R(18) = c
	    e = BETA/(1 + 4*b/sqrt(b*b+1) + 4*b + 4*b/L2 + 4*b*c/sqrt(b*b+c*c) + c)
	    dirprob(0) = 1 - BETA
	    dirprob(1:nreldir) = e*R(1:nreldir)/d(1:nreldir)
    elseif (MODEL == MOORE26_MODEL) then
    	nreldir = 26
	    L2 = sqrt(2.0d0)
	    L3 = sqrt(3.0d0)
	    b = 1 - RHO
	    b = b*b		! revised version
	    c = b*b		! oblate spheroid on reverse side
	    d(1) = 1
	    R(1) = 1
	    d(2:5) = L2
	    R(2:5) = L2*b/sqrt(b*b+1)
	    d(6:9) = L3
	    R(6:9) = L2*b/sqrt(b*b+1)
	    d(10:13) = 1
	    R(10:13) = b
	    d(14:17) = L2
	    R(14:17) = b
	    d(18:21) = L3
	    R(18:21) = L2*b*c/sqrt(b*b+c*c)
	    d(22:25) = L2
	    R(22:25) = L2*b*c/sqrt(b*b+c*c)
	    d(26) = 1
	    R(26) = c
	    e = BETA/(1 + 4*(1+L2/L3)*b/sqrt(b*b+1) + 4*b + 4*b/L2 + 4*(1+L2/L3)*b*c/sqrt(b*b+c*c) + c)
	    dirprob(0) = 1 - BETA
	    dirprob(1:nreldir) = e*R(1:nreldir)/d(1:nreldir)
    endif
endif
end subroutine

!--------------------------------------------------------------------------------
! The motivation for precomputing the array chemo_p() is to avoid having to
! compute the jump direction relative probabilities as a result of chemotaxis
! every time a cell moves.
! The array chemo_p() holds the chemotaxis-only relative probabilities of jumps in
! the set of possible jump directions, for a comprehensive set of possible
! chemokine gradient directions.  For a large enough value of chemo_N (e.g. 5),
! the set of points (x,y,z) in the cube given by:
!    x: -chemo_N,..,chemo_N
!    y: -chemo_N,..,chemo_N
!    z: -chemo_N,..,chemo_N
! generates a sufficiently extensive set of directions from (0,0,0).
! The vector v is the chemokine concentration gradient scaled by chemo_N and
! converted to integer.  Effectively it is a discretised approximation to
! the direction of the concentration gradient.
! Note that the weight given to a jump direction is found from the cosine^2 of
! the angle between the jump direction and the chemokine gradient vector.
! The weights are then normalised.
!--------------------------------------------------------------------------------
subroutine chemo_p_setup
integer :: x, y, z, k, r(3), s(3), z1, z2, ndirs
real(REAL_KIND) :: r2, s2, rmod,smod,cosa
real(REAL_KIND), allocatable :: w(:)

write(logmsg,*) 'chemo_p_setup'
call logger(logmsg)
if (SIMULATE_2D) then
	ndirs = njumpdirs2D
	allocate(chemo_p(-chemo_N:chemo_N,-chemo_N:chemo_N,0:0,ndirs))
	z1 = 0
	z2 = 0
else
	ndirs = njumpdirs
	allocate(chemo_p(-chemo_N:chemo_N,-chemo_N:chemo_N,-chemo_N:chemo_N,ndirs))
	z1 = -chemo_N
	z2 = chemo_N
endif
allocate(w(ndirs))

do x = -chemo_N,chemo_N
    do y = -chemo_N,chemo_N
        do z = z1,z2
            r = (/x,y,z/)
            r2 = dot_product(r,r)
            rmod = sqrt(r2)
            w = 0
            do k = 1,ndirs
                if (SIMULATE_2D) then
	                s = jumpvec2D(:,k)
	            else
	                if (k == 14) cycle
		            s = jumpvec(:,k)
		        endif
                s2 = dot_product(s,s)
                smod = sqrt(s2)
                cosa = dot_product(r,s)/(rmod*smod)
                if (cosa > 0) then
                    w(k) = cosa*cosa/smod
                endif
            enddo
            w = w/sum(w)
            chemo_p(x,y,z,:) = w	! Note that these probs sum to 1 - they are relative probabilities
        enddo
    enddo
enddo
deallocate(w)
end subroutine


!--------------------------------------------------------------------------------
! Precomputes the jump probabilities (absolute directions) accounting for chemotaxis
! On input:
!   p(:) holds the jump probabilities not accounting for  chemotaxis
!   v(:) was the site offset relative to the exit, used only to approximate direction
!        by allowing a discrete number of directions (chemo_N determines these, in fact
!        only a small subset of the array positions are used - those roughly falling
!        on a sphere of radius chemo_N)
!   f is the amount of chemotactic influence
! Note that f incorporates both the magnitude of chemokine gradients and the cell's
! susceptibility to chemotaxis.
! On return p(:) holds the modified jump probabilities.
! Note: when njumpdirs = 27, jump 14 corresponds to (0,0,0) - unused.
! Note: code modifications now base v and f on net chemotactic attraction of
! multiple attractors (exits and DCs) by summing the chemokine gradient vectors.
!--------------------------------------------------------------------------------
subroutine chemo_probs_pre(p,v,f)
real(REAL_KIND) :: p(:)
integer :: v(:)
real(REAL_KIND) :: f
integer :: k
real(REAL_KIND) :: pc(MAXRELDIR+1)

if (f == 0) then
    return
endif
if (f > 1) then
	write(logmsg,*) 'ERROR: chemo_probs_pre: f > 0: ',f
	call logger(logmsg)
	return
endif
p = p/sum(p)
pc(1:njumpdirs) = chemo_p(v(1),v(2),v(3),:)
do k = 1,njumpdirs
    if (p(k) > 0) then      ! prevents jumps in disallowed directions
        p(k) = (1-f)*p(k) + f*pc(k)
    endif
enddo
end subroutine

subroutine chemo_probs_pre2D(p,v,f)
real(REAL_KIND) :: p(:)
integer :: v(:)
real(REAL_KIND) :: f
integer :: k
real(REAL_KIND) :: pc(MAXRELDIR2D+1)

if (f == 0) then
    return
endif
if (f > 1) then
	write(logmsg,*) 'ERROR: chemo_probs_pre: f > 0: ',f
	call logger(logmsg)
	return
endif
p = p/sum(p)
pc(1:njumpdirs2D) = chemo_p(v(1),v(2),v(3),:)
!write(*,'(a,3i4)') 'v: ',v
!write(*,'(a,8f8.4)') 'pc: ',pc(1:njumpdirs2D)
do k = 1,njumpdirs2D
    if (p(k) > 0) then      ! prevents jumps in disallowed directions
        p(k) = (1-f)*p(k) + f*pc(k)
    endif
enddo
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine make_probvectors
integer :: icase, i, k, dir, jump(3)
real(REAL_KIND) :: vec(3), veclen, scale
real(REAL_KIND) :: xang = 0.3, yang = -0.15, zang = -0.5  ! radians
real(REAL_KIND) :: scale_M18 = 10  ! for M18
real(REAL_KIND) :: scale_N = 6     ! for N
real(REAL_KIND) :: arrow_head = 0.07

!write(*,*) 'make_probvectors'
call make_reldir
if (MODEL == NEUMANN_MODEL) then
    scale = scale_N
elseif (MODEL == MOORE18_MODEL) then
    scale = scale_M18
endif
open(nfvec,file='d:\immunesystem\lattice\motilitypaper\jump_vectors\cmgui\prob.exnode',status='replace')
beta = 1.0
k = 0
do icase = 0,4
    rho = (icase-1)*0.2
    call compute_dirprobs
    write(nfvec,'(a,i1)') 'Group name : probability',icase
    write(nfvec,'(a)') ' #Fields=1'
    write(nfvec,'(a)') ' 1) vector, field, rectangular cartesian, #Components=3'
    write(nfvec,'(a)') ' x.  Value index= 1, #Derivatives= 0'
    write(nfvec,'(a)') ' y.  Value index= 2, #Derivatives= 0'
    write(nfvec,'(a)') ' z.  Value index= 3, #Derivatives= 0'
    if (icase > 0) then
        do i = 1,nreldir
            k = k+1
            dir = reldir(2,i)
            jump = jumpvec(:,dir)
            vec = jump*dirprob(i)
            call rotate(vec,xang,yang,zang)
            vec = scale*vec
            veclen = norm(vec)
            if (veclen > arrow_head) then
                vec = vec*(veclen-arrow_head)/veclen
            else
                vec = 0
            endif
            write(nfvec,'(a,i3)') 'Node: ',k
            write(nfvec,'(3f8.4)') vec
        enddo
    else
        k = k+1
        vec = (/0.1,0.0,0.0/)
        call rotate(vec,xang,yang,zang)
        vec = scale*vec
        veclen = norm(vec)
        if (veclen > arrow_head) then
            vec = vec*(veclen-arrow_head)/veclen
        else
            vec = 0
        endif
         write(nfvec,'(a,i3)') 'Node: ',k
        write(nfvec,'(3f8.4)') vec
    endif
enddo
close(nfvec)
end subroutine

!---------------------------------------------------------------------
! Rotate the vector about the X-axis by xang, then about the Y-axis by yang
!---------------------------------------------------------------------
subroutine rotate(vec,xang,yang,zang)
real(REAL_KIND) :: vec(3),xang,yang,zang
real(REAL_KIND) :: new(3), cosa, sina

cosa = cos(xang)
sina = sin(xang)
new(1) = vec(1)
new(2) = cosa*vec(2) - sina*vec(3)
new(3) = sina*vec(2) + cosa*vec(3)
vec = new
cosa = cos(zang)
sina = sin(zang)
new(3) = vec(3)
new(1) = cosa*vec(1) - sina*vec(2)
new(2) = sina*vec(1) + cosa*vec(2)
vec = new
cosa = cos(yang)
sina = sin(yang)
new(2) = vec(2)
new(1) = cosa*vec(1) - sina*vec(3)
new(3) = sina*vec(1) + cosa*vec(3)
vec = new

end subroutine

!---------------------------------------------------------------------
! Choose one of the 6 axis directions to allocate to the direction of
! previous step given by jump, which takes values from 0 - 27.
! This could be already on an axis, or could be a D2 or D3 diagonal,
! in which case the choice is random.
!---------------------------------------------------------------------
integer function fix_lastdir(jump,kpar)
integer :: jump, kpar
integer :: k,nax
!                            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
integer :: naxes(0:27) = (/  0, 3, 2, 3, 2, 1, 2, 3, 2, 3, 2, 1, 2, 1, 0, 1, 2, 1, 2, 3, 2, 3, 2, 1, 2, 3, 2, 3 /)
integer :: axes(3,27) = reshape( (/ &
!   1      2      3      4      5      6      7      8      9     10     11     12     13     14
1,3,5, 1,3,0, 1,3,6, 1,5,0, 1,0,0, 1,6,0, 1,4,5, 1,4,0, 1,4,6, 3,5,0, 3,0,0, 3,6,0, 5,0,0, 0,0,0, &
!  15     16     17     18     19     20     21     22     23     24     25     26     27
6,0,0, 4,5,0, 4,0,0, 4,6,0, 2,3,5, 2,3,0, 2,3,6, 2,5,0, 2,0,0, 2,6,0, 2,4,5, 2,4,0, 2,4,6 /), (/3,27/) )

if (diagonal_jumps) then
	nax = naxes(jump)
	if (nax == 0) then
	    write(*,*) 'Should not get here: fix_lastdir: nax=0'
	    stop
		fix_lastdir = random_int(1,6,kpar)
	elseif (nax == 1) then
		fix_lastdir = axes(1,jump)
	else
		k = random_int(1,nax,kpar)
		fix_lastdir = axes(k,jump)
	endif
else
	write(*,*) 'fix_lastdir: Not MOORE model: ',jump
	stop
endif
end function

end module
