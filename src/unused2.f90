!--------------------------------------------------------------------------------------
! Note: _nc refers to non-cognate T cells
!--------------------------------------------------------------------------------------
subroutine make_bind_probs
integer :: i
real :: psum, ptsum
real, allocatable :: pt(:)

allocate(pt(Nbindtime_nc))
psum = 0
ptsum = 0
do i = 1,Nbindtime_nc
	psum = psum + dc_bindprob_nc(i)
	pt(i) = dc_bindprob_nc(i)*dc_bindtime_nc(i)
	ptsum = ptsum + pt(i)
enddo
do i = 1,Nbindtime_nc
	dc_bindprob_nc(i) = dc_bindprob_nc(i)/psum
	pt(i) = pt(i)/ptsum
enddo
psum = 0
ptsum = 0
do i = 1,Nbindtime_nc
	psum = psum + dc_bindprob_nc(i)
	dc_cummul_prob_nc(i) = psum
	ptsum = ptsum + pt(i)
	dc_initial_cummul_prob_nc(i) = ptsum
enddo
!write(*,*) 'dc_cummul_prob_nc:'
!write(*,'(10f5.2)') dc_cummul_prob_nc
!write(*,*) 'dc_initial_cummul_prob_nc:'
!write(*,'(10f5.2)') dc_initial_cummul_prob_nc 
!write(*,*)

deallocate(pt)

end subroutine

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function DClifetime(kpar)
integer :: kpar
real :: p1,p2

p1 = log(DC_LIFETIME_MEDIAN*24*60)
p2 = log(DC_LIFETIME_SHAPE)
DClifetime = rv_lognormal(p1,p2,kpar)
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function DCdensity(kpar)
integer :: kpar
real :: p1,p2

p1 = log(DC_ANTIGEN_MEDIAN)
p2 = log(DC_ANTIGEN_SHAPE)
DCdensity = rv_lognormal(p1,p2,kpar)
end function

!--------------------------------------------------------------------------------------
! If a cognate cell is NAIVE it is long-lived - effectively it does not die.
! When a cognate cell has contact with a cognate DC its lifetime is immediately limited.
!--------------------------------------------------------------------------------------
real function TClifetime(ptr)
type (cog_type), pointer :: ptr
integer :: gen, stage, region
real :: p1, p2
integer :: kpar = 0

TClifetime = 0
!stage = get_stage(ptr)
call get_stage(ptr,stage,region)
if (stage == NAIVE) then
    TClifetime = BIG_TIME
    return
endif
gen = get_generation(ptr)
if (gen < 1 .or. gen > TC_MAX_GEN) then
    write(*,*) 'TClifetime: bad gen: ',gen
!    stop
endif
p1 = life_dist(gen)%p1
p2 = life_dist(gen)%p2
select case (life_dist(gen)%class)
case (NORMAL_DIST)
	TClifetime = rv_normal(p1,p2,kpar)
case (LOGNORMAL_DIST)
	TClifetime = rv_lognormal(p1,p2,kpar)
case (CONSTANT_DIST)
	TClifetime = p1
end select
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function dividetime(gen,ictype)
integer :: gen,ictype
real :: p1, p2
integer :: kpar = 0

dividetime = 0
p1 = divide_dist(gen)%p1
p2 = divide_dist(gen)%p2
!write(*,*) 'dividetime: ',istep,gen,divide_dist(gen)%p1,divide_dist(gen)%p2
select case (divide_dist(gen)%class)
case (NORMAL_DIST)
	dividetime = rv_normal(p1,p2,kpar)
case (LOGNORMAL_DIST)
!    write(*,*) 'dividetime: ',istep,gen,p1,p2
	dividetime = rv_lognormal(p1,p2,kpar)
!	write(*,*) 'dividetime: ',dividetime
case (CONSTANT_DIST)
	dividetime = p1
end select
!if (ictype == CD8_CELL) then
!	dividetime = dividetime*CD8_DIVTIME_FACTOR
!endif
end function

!--------------------------------------------------------------------------------
! Place n DCs
! Note that we may have NDCsites > 1 !!!
! The procedure for placing DCs is as follows:
! The central DC is placed on a site that satisfies these conditions:
! (1) The site is not OUTSIDE_TAG or a DC site
! (2) There are not two T cells at this site
! (3) The site is not too near another DC, i.e. it is more than DC_DCprox*DC_RADIUS
!     from any DC.
! (4) The site is not too near an exit, i.e. it is more than exit_DCprox from any exit.
! (5) The site is not too near the blob boundary, i.e. it is more than
!     bdry_DCprox*DC_RADIUS from the sphere with radius = Radius
! (6) Either the site is free, or it is occupied by a T cell that can be moved
!     to a neighbouring site.
! A site meeting these requirements is selected for the DC, then as many as
! possible of the neighbouring sites are also allocated to the NDCsites-1 other
! sites occupied by a DC, by subroutine addDCsite().  The count of DC sites
! allocated is stored in DC%nsites.
!--------------------------------------------------------------------------------
subroutine place_DCs(n,nadded)
integer :: n, nadded
integer :: i, x, y, z, site1(3), site2(3), freeslot, err
integer :: indx(2), jslot, idc, k, kcell
integer :: xmin, xmax, ymin, ymax, zmin, zmax
integer :: kpar = 0
real(DP) :: R
real :: tnow, dist, rvec(3), prox, tmins
logical :: OK
type(DC_type),pointer :: DC
integer, parameter :: kmax = 10000

!ndbug = DClist(idbug)%ncogbound
!write(*,*) 'place_DCs (1): ',DClist(idbug)%ncogbound
tnow = istep*DELTA_T
xmin = x0 - Radius
xmax = x0 + Radius + 1
ymin = y0 - Radius
ymax = y0 + Radius + 1
zmin = z0 - Radius
zmax = z0 + Radius + 1
xmin = max(1,xmin)
xmax = min(NX,xmax)
ymin = max(1,ymin)
ymax = min(NY,ymax)
zmin = max(1,zmin)
zmax = min(NZ,zmax)
nadded = 0
do i = 1,n
	OK = .false.
	k = 0
    do
		k = k+1
		if (k > kmax) then
			write(logmsg,*) 'Error: place_DCs: unable to find space for a DC'
			call logger(logmsg)
			OK = .false.
			stop
		endif
        R = par_uni(kpar)
        x = xmin + R*(xmax-xmin)
        R = par_uni(kpar)
        y = ymin + R*(ymax-ymin)
        R = par_uni(kpar)
        z = zmin + R*(zmax-zmin)
        site1 = (/x,y,z/)   ! global location
        indx = occupancy(x,y,z)%indx
        if (indx(1) < 0) cycle                          ! OUTSIDE_TAG or DC
        if (indx(1) /= 0 .and. indx(2) /= 0) cycle      ! two T cells at this site
!        prox = proximity_limit(DC_SITE,DC_SITE)
!        if (tooNearDC(site1,prox)) cycle
!        prox = proximity_limit(DC_SITE,EXIT_SITE)
!        if (tooNearExit(site1,prox)) cycle
!!        prox = 0.5*DC_DCprox*DC_RADIUS
!!        prox = bdry_DCprox*DC_RADIUS
!        prox = proximity_limit(DC_SITE,BDRY_SITE)
!        if (.not.tooNearBdry(site1,prox)) then
		call CheckSite(DC_SITE,site1,ok)
		if (ok) then
            jslot = 0
            if (indx(1) /= 0) then
                jslot = 1
            elseif (indx(2) /= 0) then
                jslot = 2
            endif
            if (jslot == 0) then ! free site, use it
                exit
            else  ! one T cell here in jslot, it must be moved
                ! Can the T cell be bumped to a neighbour site?
                call get_free_slot(occupancy,NX,site1,site2,freeslot)
                if (freeslot == 0) cycle    ! cannot be bumped, try again
               ! Move the cell in site1/jslot to site2/freeslot
                kcell = indx(jslot)
!                write(*,*) 'place_DCs: bump cell: ',kcell,site1,site2
                occupancy(x,y,z)%indx = 0
                occupancy(site2(1),site2(2),site2(3))%indx(freeslot) = kcell
                cellist(kcell)%site = site2
                ! Now site1 is free to use
                exit
            endif
        endif
    enddo
!    if (DClist(idbug)%ncogbound /= ndbug) then
!        write(*,*) 'place_DCs (a): ndbug changed: ',i,DClist(idbug)%ncogbound
!        stop
!    endif
    if (.not.OK) exit
    idc = 0
    if (reuse_DC_index) then
	    do k = 1,NDC
		    if (.not.DClist(k)%alive) then
		        idc = k
		        exit
		    endif
		enddo
	endif
    if (idc == 0) then    ! If there isn't a free spot in DClist()
        NDC = NDC + 1
        if (NDC > max_DC) then
			write(logmsg,'(a,i6)') 'Error: place_DCs: number of DCs exceeds limit: ',max_DC
			call logger(logmsg)
			ok = .false.
			return
		endif
        idc = NDC
    endif
    
    DC => DClist(idc)   ! revised
    
    DC%ID = idc
    DC%alive = .true.
    DC%capable = .true.
    DC%site = site1
    DC%nsites = 1
    DC%stimulation = 0
    DC%nbound = 0
    DC%ncogbound = 0
    allocate(DC%cogbound(MAX_COG_BIND))
    DC%cogbound = 0
    DC%density = DCdensity(kpar)
    if (DC_INJECTION) then
        ! If the DCs were injected into an experimental animal, we need to account for the antigen decay since
        ! the time of injection, which was at T_DC_INJECTION hours
        tmins = -T_DC_INJECTION*60 + istep*DELTA_T
        DC%density = DC%density*exp(-DCdecayrate*tmins)
    endif
    DC%dietime = tnow + DClifetime(kpar)
    occupancy(site1(1),site1(2),site1(3))%indx(1) = -idc
    do k = 2,NDCsites
        call addDCsite(idc,site1,k,err)
!        if (DClist(idbug)%ncogbound /= ndbug) then
!            write(*,*) 'place_DCs (b): ndbug changed: ',i,k,DClist(idbug)%ncogbound
!            stop
!        endif
        if (err == 0) then
            DC%nsites = DC%nsites + 1
        else
!            write(*,*) 'addDCsite: idc,k,err: ',idc,k,err
        endif
    enddo
!    DClist(idc) = DC
    nadded = nadded + 1
!    write(*,*) 'Added DC at: ',site1,' with nsites: ',DC%nsites
    ! now the DC proximity data in occupancy()%DC must be updated
    ! this is done by reassign_DC() called from balancer()
!    if (DClist(idbug)%ncogbound /= ndbug) then
!        write(*,*) 'place_DCs (c): ndbug changed: ',i,DClist(idbug)%ncogbound
!        stop
!    endif
enddo
NDCtotal = NDCtotal + nadded

! Now check the DC locations wrt the blob perimeter
do k = 1,NDC
    if (DClist(k)%alive) then
        rvec = DClist(k)%site - (/x0,y0,z0/)
        dist = norm(rvec)
        if (dist > Radius) then
            write(logmsg,*) 'Place_DCs: warning: DC distance: ',k,dist,Radius
			call logger(logmsg)
        endif
    endif
enddo
!write(*,*) 'place_DCs (2): ',DClist(idbug)%ncogbound
end subroutine

!
! Note: Need to account for DCs with cognate and non-cognate antigen.
! In the case of DC_INJECTION, there needs to be a distinction between the DCs initially
! resident in the LN, and injected DCs.
!
NDC = 0
NDCalive = 0
if (TC_TO_DC > 0) then   ! DC placement needs to be checked to account for SOI overlap
!if (use_DC) then   ! DC placement needs to be checked to account for SOI overlap
    NDCrequired = DC_FACTOR*nlist/TC_TO_DC
    if (use_single_DC) then
		NDCrequired = 1
	endif
    NDCtotal = NDCrequired
    if (NDCrequired > MAX_DC) then
        write(logmsg,'(a,i6)') 'Error: PlaceCells: NDC exceeds MAX_DC: ',MAX_DC
        call logger(logmsg)
		ok = .false.
        return
    endif
    p1 = log(DC_ANTIGEN_MEDIAN)
    p2 = log(DC_ANTIGEN_SHAPE)
    do idc = 1,NDCrequired
        do
            R = par_uni(kpar)
            x = 1 + R*NX
            R = par_uni(kpar)
            y = 1 + R*NY
			if (IN_VITRO) then
				z = 1
			else
	            R = par_uni(kpar)
		        z = 1 + R*NZ
		    endif
		    if (use_single_DC) then
				x = NX/2
				y = NY/2
				z = NZ/2
			endif
            site = (/x,y,z/)
            if (occupancy(x,y,z)%indx(1) < 0) cycle     ! OUTSIDE_TAG or DC
			call CheckSite(DC_SITE,site,ok)
			if (ok) exit
        enddo
        DClist(idc)%ID = idc
        DClist(idc)%site = site
        DClist(idc)%nsites = 1
        NDC = NDC + 1
		NDCalive = NDC
! Changed the treatment of injected DCs.  Now assume that all DCs originally in the paracortex are non-cognate,
! while the schedule of arriving DCs with peptide is read from an input file.
!        if (DC_INJECTION) then
!            ! If the DCs were injected into an experimental animal, we need to account for the antigen decay since
!            ! the time of injection, which was at T_DC_INJECTION hours
!            tmins = -T_DC_INJECTION*60
!            DClist(idc)%density = DClist(idc)%density*exp(-DCdecayrate*tmins)
!        endif
	    DClist(idc)%dietime = tnow + DClifetime(kpar)
	    DClist(idc)%alive = .true.
	    DClist(idc)%stimulation = 0
	    DClist(idc)%nbound = 0
	    DClist(idc)%ncogbound = 0
	    allocate(DClist(idc)%cogbound(MAX_COG_BIND))
	    DClist(idc)%cogbound = 0
		DClist(idc)%cognate = .true.
        if (DC_INJECTION) then
		    DClist(idc)%capable = .false.
	        DClist(idc)%density = 0
        else
		    DClist(idc)%capable = .true.
	        DClist(idc)%density = DCdensity(kpar)
		endif
    enddo
    do idc = 1,NDC
		call assignDCsites(idc,nassigned,ok)
		if (.not.ok) return
        DClist(idc)%nsites = nassigned
        nlist = nlist - nassigned
		site = DClist(idc)%site
		xdc = site(1)
		ydc = site(2)
		zdc = site(3)
		xmin = xdc - DC_RADIUS
        xmax = xdc + DC_RADIUS
        ymin = ydc - DC_RADIUS
        ymax = ydc + DC_RADIUS
        zmin = zdc - DC_RADIUS
        zmax = zdc + DC_RADIUS
        xmin = max(1,xmin)
        xmax = min(NX,xmax)
        ymin = max(1,ymin)
        ymax = min(NY,ymax)
        zmin = max(1,zmin)
        zmax = min(NZ,zmax)
        do x = xmin,xmax
	        x2 = (x-xdc)*(x-xdc)
	        do y = ymin,ymax
		        y2 = (y-ydc)*(y-ydc)
		        do z = zmin,zmax
			        z2 = (z-zdc)*(z-zdc)
			        d2 = x2 + y2 + z2
			        added = .false.
					! The following procedure is correct only if DC_RADIUS <= chemo_radius
			        if (d2 <= DC_RADIUS*DC_RADIUS) then
			            kdc = occupancy(x,y,z)%DC(0)
				        if (kdc < 0) cycle   ! don't touch a DC site
			            if (kdc < DCDIM-1) then     ! can add to the list
			                kdc = kdc + 1
			                occupancy(x,y,z)%DC(0) = kdc
			                occupancy(x,y,z)%DC(kdc) = idc
			                added = .true.
				        endif
				    endif
			        if (.not.added .and. (d2 <= chemo_radius*chemo_radius)) then
			            kdc = occupancy(x,y,z)%cDC(0)
			            if (kdc == cDCDIM-1) cycle     ! can add no more to the list
		                kdc = kdc + 1
		                occupancy(x,y,z)%cDC(0) = kdc
		                occupancy(x,y,z)%cDC(kdc) = idc
			        endif
		        enddo
	        enddo
        enddo
        if (IN_VITRO) then
			call setDCzrange(xdc,ydc)
		endif
    enddo
    call check_DCproximity
endif
write(logmsg,*) 'Number of DCs: ',NDC
call logger(logmsg)
if (USE_DC_COGNATE) then
	ncogDC = DC_COGNATE_FRACTION*NDC
	do idc = ncogDC+1,NDC
		DClist(idc)%cognate = .false.
	enddo
	write(logmsg,*) 'Number of DCs bearing cognate antigen: ',ncogDC
	call logger(logmsg)
endif

!--------------------------------------------------------------------------------
! Note: USE_STAGETIME(icstage) means use stagetime for the transition from icstage.
! p%stagetime = time that the cell is expected to make transition to the next stage.
! In this simplified version, the sequence is (revised_staging = .true.):
! NAIVE -> TRANSIENT -> CLUSTERS -> SWARMS -> DIVIDING -> SWARMS
! On cell division, the stage is set to SWARMS, and the time for the next stage
! transition (to DIVIDING, provided candivide()) is the division time.
!--------------------------------------------------------------------------------
subroutine updatestage(kcell,tnow,divide_flag)
integer :: kcell
logical :: divide_flag
real :: tnow, t
integer :: stage, region, gen, ctype
real :: stagetime
type(cog_type), pointer :: p

divide_flag = .false.
p => cellist(kcell)%cptr

! In the UNSTAGED case, cells have a simplified stage development:
! NAIVE -> SWARMS -> ...  
! Division depends on two conditions being met: 
!   time since first signalling must exceed UNSTAGED_MIN_DIVIDE_T
!   stimulation must exceed the threshold level for division at this generation 

call get_stage(p,stage,region)
!if (stage >= CLUSTERS) then
!    if (.not.cansurvive(p)) then
!        write(logmsg,*) 'cell IL2 store too low: ',kcell,p%cogID
!        call logger(logmsg)
!        p%dietime = tnow
!        stage = FINISHED
!        call set_stage(p,stage)
!    endif
!endif
if (stage == FINISHED) return
ctype = cellist(kcell)%ctype
stagetime = p%stagetime
if (tnow > stagetime) then		! time constraint to move to next stage is met
	! May be possible to make the transition from stage
	select case(stage)
	case (NAIVE)        ! possible transition from NAIVE to TRANSIENT
	    if (p%stimulation > 0) then
	        gen = get_generation(p)
            p%firstDCtime = tnow
            t = tnow - cellist(kcell)%entrytime
            call log_count(DCfirstcontact_count,t)
            p%dietime = tnow + TClifetime(p)
	        if (activation_mode == STAGED_MODE) then
                call set_stage(p,TRANSIENT)
		        if (USE_STAGETIME(TRANSIENT)) then
			        p%stagetime = tnow + get_stagetime(p,ctype)
		        else
			        p%stagetime = 0		! time is not criterion for next transition
                endif
            else
                call set_stage(p,SWARMS)
		        p%stagetime = 0		! time is not criterion for next transition
            endif
        endif
	case (TRANSIENT)    ! possible transition from TRANSIENT to CLUSTERS
	    if (reached_IL2_threshold(p)) then
	        nIL2thresh = nIL2thresh + 1
	        tIL2thresh = tIL2thresh + (tnow - cellist(kcell)%entrytime)
!	        write(*,'(a,i6,2f8.1)') '========= Reached IL2 threshold: ',kcell,(tnow - cellist(kcell)%entrytime),p%stimulation
            call set_stage(p,CLUSTERS)
		    if (USE_STAGETIME(CLUSTERS)) then
			    p%stagetime = tnow + get_stagetime(p,ctype)
		    else
			    p%stagetime = 0		! time is not criterion for next transition
		    endif
        endif
	case (CLUSTERS)     ! possible transition from CLUSTERS to SWARMS
	    if (reached_act_threshold(p)) then
            call set_stage(p,SWARMS)
            gen = get_generation(p)
            ctype = cellist(kcell)%ctype
		    if (USE_STAGETIME(SWARMS)) then
				p%stagetime = tnow + dividetime(gen,ctype)
		    else
			    p%stagetime = 0		! time is not criterion for next transition
		    endif
        endif
	case (SWARMS)       ! possible transition from SWARMS to DIVIDING
        gen = get_generation(p)
        ctype = cellist(kcell)%ctype
		if (gen == TC_MAX_GEN) then
            call set_stage(p,FINISHED)
			p%stagetime = BIG_TIME
		elseif (candivide(p,ctype)) then
            call set_stage(p,DIVIDING)
		    if (USE_STAGETIME(DIVIDING)) then
!				p%stagetime = tnow + CYTOKINESIS_TIME
				p%stagetime = tnow + get_stagetime(p,ctype)
		    else
			    p%stagetime = 0		! time is not criterion for next transition
		    endif
		endif
	case (DIVIDING)
	    divide_flag = .true.
	case (FINISHED)

    end select
endif
end subroutine

!--------------------------------------------------------------------------------------
! The bind time depends on the cell stage, and possibly on cell type (CD4/CD8).
! Rules currently apply only during initial activation, Phases 1 & 2 (TRANSIENT, CLUSTERS)
! CT_CONSTANT
! Just use mean bind times for cognate cells.
! CT_HILL
! Now try making bind time depend on the TCR stimulation rate of the binding
! CT_HENRICKSON
! The bind time is drawn from a lognormal distn, fixed shape, median increases with
! the level of TCR stimulation, up to a limiting value.
! HYPOTHESIS: bind time should also increase with level of stimulation rate,
! i.e. CT_HILL should be combined with CT_HENRICKSON - but how?
!--------------------------------------------------------------------------------------
real function get_bindtime(p,signal,ctype,pMHC,stimrate_norm,kpar)
type(cog_type), pointer :: p
logical :: signal
integer :: ctype, kpar
real :: pMHC, stimrate_norm
integer :: stage, region, i
real(DP) :: R
real :: stimrate_hill, btime, p1, median, h, a, d
real, parameter :: CT_shape = 2.0, CT_max_factor = 4.0
real, parameter :: p2 = log(CT_shape)

get_bindtime = 0
if (signal) then
    if (activation_mode == STAGED_MODE) then
	    call get_stage(p,stage,region)
	    btime = dc_mean_bindtime_c(stage,ctype)
	    stimrate_hill = stimulation_rate_hill(pMHC,p%avidity)
!	    write(logmsg,'(a,a,i2,a,f5.3,a,f6.1)') 'STAGED: ',' stage: ',stage,' stimrate_hill: ',stimrate_hill,' T: ',btime
!	    call logger(logmsg)
	    if (stage > NAIVE .and. stage < SWARMS) then
	        select case (STAGED_CONTACT_RULE)
	        case (CT_CONSTANT)
	            get_bindtime = btime
	        case (CT_HILL)
	            stimrate_hill = stimulation_rate_hill(pMHC,p%avidity)
	            get_bindtime = btime*(0.10 + stimrate_hill)
	        case (CT_HENRICKSON)
	            median = CT_median(p,pMHC)
	            p1 = log(median)
	            get_bindtime = min(CT_max_factor*median,rv_lognormal(p1,p2,kpar))
	        end select
	    else
	        get_bindtime = btime
        endif
    else
        a = min(1.0,p%avidity/MAXIMUM_AVIDITY)
        d = min(1.0,pMHC/MAXIMUM_ANTIGEN)
        h = stimrate_norm**BINDTIME_HILL_N/(stimrate_norm**BINDTIME_HILL_N + BINDTIME_HILL_C**BINDTIME_HILL_N)
        h = (1 + BINDTIME_HILL_C**BINDTIME_HILL_N)*h
        get_bindtime = h*BINDTIME_MAX
!	    write(logmsg,'(a,a,f5.3,a,f5.3,a,f5.3,a,f5.3,a,f6.1)') &
!	        'UNSTAGED: ','A: ',a,' D: ',d,' dS/dt: ',stimrate_norm,' H: ',h,' T: ',get_bindtime
!	    call logger(logmsg)
    endif
else
    R = par_uni(kpar)
	do i = 1,Nbindtime_nc
		if (R < dc_cummul_prob_nc(i)) then
			get_bindtime = dc_bindtime_nc(i)
			return
		endif
	enddo
	get_bindtime = dc_bindtime_nc(Nbindtime_nc)
endif
end function

!----------------------------------------------------------------------------------------
! Estimates the median DC contact time on the basis of the level of TCR stimulation.
! The time varies linearly with stimulation, from a minimum value at zero stimulation,
! bounded by a maximum value.
! Use a ramp function for now.
! Min = 5, max = 70
! The slope is set by requiring that the median when S = IL2_THRESHOLD, i.e. the
! Phase 1/Phase 2 transition, takes a specified value T12_median
! (Currently pMHC is not used)
!----------------------------------------------------------------------------------------
real function CT_median(p,pMHC)
type(cog_type), pointer :: p
real :: pMHC
real, parameter :: min_median = 5, max_median = 70, T12_median = 25		!<------ hard-wired
real :: slope

slope = (T12_median - min_median)/IL2_THRESHOLD
CT_median = min(max_median,min_median + slope*p%stimulation)
end function

    if (.not.associated(cell%cptr)) then
        allocate(cell%cptr)
    endif
    param1 = log(TC_AVIDITY_MEDIAN)
    param2 = log(TC_AVIDITY_SHAPE)
    cell%cptr%status = 0
    call set_activation(cell%cptr,0)
    call set_generation(cell%cptr,gen)
    call set_stage_region(cell%cptr,stage,region)
    if (fix_avidity) then
        i = mod(navid,avidity_nlevels)
        navid = navid + 1
        cell%cptr%avidity = avidity_level(i+1)
!        if (avidity_logscale) then
!            cell%cptr%avidity = 10**(avidity_min + i*avidity_step)
!        else
!            cell%cptr%avidity = avidity_min + i*avidity_step
!        endif
!        write(nfout,'(a,4i7,2f8.4)') 'create_Tcell: add: ',i,navid,kcell,lastcogID + 1, &
!            cell%cptr%avidity,cellist(kcell)%cptr%avidity
!        if (i == 3 .and. istep >= 4500) then
!            avid_debug = .true.
!        endif
    else
        cell%cptr%avidity = rv_lognormal(param1,param2,kpar)
    endif
    cell%cptr%stimulation = 0
    cell%cptr%stimrate = 0
    cell%cptr%effector = .false.
!    cell%cptr%IL_state = 0
!    cell%cptr%IL_statep = 0
    if (use_cytokines) then
        call IL2_init_state(cell%cptr%IL_state,cell%cptr%IL_statep)
    endif
    ! What should the initial CD69 level be?  If 0, can a cognate cell exit immediately?
    ! We would prefer not to impose a time or generation constraint on the exit of
    ! cognate T cells, but otherwise if CD69 is initially 0, a cell will be susceptible
    ! to chemotaxis and exit until it has received enough TCR signal to drive CD69 high.
    cell%cptr%CD69 = 0
    cell%cptr%S1PR1 = 0
    cell%cptr%CFSE = generate_CFSE(1.0)
!	cell%cptr%DCchemo = BASE_DCchemo
	cell%cptr%firstDCtime = 0
    cell%cptr%dietime = tnow + TClifetime(cell%cptr)
    cell%cptr%dividetime = tnow
    cell%cptr%stagetime = BIG_TIME
    cell%cptr%totalDCtime = 0
	cell%cptr%cnt = 0

! Maintain cognate_list at start or if we are running on a single node
! Otherwise cogID and cognate_list is maintained by make_cognate_list
!    if (istep == 0 .or. Mnodes == 1) then
        lastcogID = lastcogID + 1
        if (lastcogID > MAX_COG) then
            write(logmsg,'(a,i6)') 'Error: create_Tcell: cognate_list dimension exceeded: ',MAX_COG
            call logger(logmsg)
            ok = .false.
            return
        endif
        cogID = lastcogID
        cell%cptr%cogID = cogID
        cognate_list(cogID) = kcell
!    else
!        cell%cptr%cogID = 0
!    endif

!--------------------------------------------------------------------------------
! Note: USE_STAGETIME(icstage) means use stagetime for the transition from icstage.
! p%stagetime = time that the cell is expected to make transition to the next stage.
!--------------------------------------------------------------------------------
subroutine updatestage1(kcell,tnow,divide_flag)
integer :: kcell
logical :: divide_flag
real :: tnow
integer :: stage, region, gen, ctype
real :: stagetime
type(cog_type), pointer :: p

divide_flag = .false.
p => cellist(kcell)%cptr
call get_stage(p,stage,region)
if (stage >= CLUSTERS) then
    if (.not.cansurvive(p)) then
        write(logmsg,*) 'cell IL2 store too low: ',kcell,p%cogID
        call logger(logmsg)
        p%dietime = tnow
        stage = FINISHED
        call set_stage(p,stage)
    endif
endif
if (stage == FINISHED) return
ctype = cellist(kcell)%ctype
stagetime = p%stagetime
if (tnow > stagetime) then		! time constraint to move to next stage is met
	! May be possible to make the transition from stage
	select case(stage)
	case (NAIVE)        ! possible transition from NAIVE to TRANSIENT
	    if (p%stimulation > 0) then
	        gen = get_generation(p)
            call set_stage(p,TRANSIENT)
            p%dietime = tnow + TClifetime(p)
		    if (USE_STAGETIME(TRANSIENT)) then
			    p%stagetime = tnow + get_stagetime(p,ctype)
		    else
			    p%stagetime = 0		! time is not criterion for next transition
            endif
        endif
	case (TRANSIENT)    ! possible transition from TRANSIENT to CLUSTERS
	    if (reached_IL2_threshold(p)) then
	        nIL2thresh = nIL2thresh + 1
	        tIL2thresh = tIL2thresh + (tnow - cellist(kcell)%entrytime)
!	        write(*,'(a,i6,2f8.1)') '========= Reached IL2 threshold: ',kcell,(tnow - cellist(kcell)%entrytime),p%stimulation
            call set_stage(p,CLUSTERS)
		    if (USE_STAGETIME(CLUSTERS)) then
			    p%stagetime = tnow + get_stagetime(p,ctype)
		    else
			    p%stagetime = 0		! time is not criterion for next transition
		    endif
        endif
	case (CLUSTERS)     ! possible transition from CLUSTERS to SWARMS
	    if (reached_act_threshold(p)) then
            call set_stage(p,SWARMS)
		    if (USE_STAGETIME(SWARMS)) then     ! Must use stagetime(SWARMS) to get to ACTIVATED
			    p%stagetime = tnow + get_stagetime(p,ctype)
		    else
			    p%stagetime = 0
		    endif
        endif
! Code commented because we are now using a revised staging without ACTIVATED
!	case (SWARMS)       ! possible transition from SWARMS to ACTIVATED
!        call set_stage(p,ACTIVATED)
!		if (USE_STAGETIME(ACTIVATED)) then      ! Must use stagetime(ACTIVATED) to get to ACTIVATED
!            gen = get_generation(p)
!            ctype = cellist(kcell)%ctype
!			if (gen == 1) then
!				p%stagetime = tnow
!			else
!				p%stagetime = tnow + dividetime(gen,ctype)
!			endif
!		else
!			p%stagetime = 0
!		endif
!	case (ACTIVATED)     ! possible transition from ACTIVATED to DIVIDING
!        gen = get_generation(p)
!        ctype = cellist(kcell)%ctype
!		if (gen == TC_MAX_GEN) then
!            call set_stage(p,FINISHED)
!			p%stagetime = BIG_TIME
!		elseif (candivide(p,ctype)) then
!            call set_stage(p,DIVIDING)
!			p%stagetime = tnow + CYTOKINESIS_TIME
!		endif
	case (DIVIDING)
	    divide_flag = .true.
	case (FINISHED)

    end select
endif
end subroutine

!--------------------------------------------------------------------------------------
! The stage time is the time that a cell will spend in a stage (after naive stage).
! May possibly on cell type (CD4/CD8).  Applies only to cognate cells.
! Note that the concept of time in a stage is just a surrogate for accurate biological
! mechanisms that determine stage transitions.
! We might be able to override time in stage 4 (SWARMS) by use of activation threshold.
!--------------------------------------------------------------------------------------
real function get_stagetime(p,ctype)
type(cog_type), pointer :: p
integer :: ctype
integer :: stage, region, gen

!stage = get_stage(p)
call get_stage(p,stage,region)
!cd4_8 = ctype - 1       ! converts ctype -> 1=CD4, 2=CD8
gen = get_generation(p)
if (gen == 1) then
	if (stage == SWARMS) then	! use 1st division time (only if revised_stagetime = .false.)
		get_stagetime = dividetime(gen,ctype)
	else
		get_stagetime = 60*mean_stagetime_1(stage,ctype)		! hr -> min
	endif
else
	get_stagetime = 60*mean_stagetime_n(stage,ctype)		! hr -> min
endif
end function

!--------------------------------------------------------------------------------
! Returns .true. if the cognate cell is permitted to bind to a DC.
!--------------------------------------------------------------------------------
logical function bindable(p)
type(cog_type), pointer :: p
integer :: stage, region, kpar=0

bindable = .true.
call get_stage(p,stage,region)
if (stage == FINISHED) then
    bindable = .false.
endif
end function

!-----------------------------------------------------------------------------------------
! Processes T cell binding and unbinding from DC.
! T cells are allowed to bind to MAX_DC_BIND DC, (MAX_DC_BIND = 1 or 2)
! Note that when a cell is able to bind (less than MAX_DC_BIND bindings) the time of the
! most recent unbinding is stored in %unbindtime(2).
! Note that when there is a single bound DC it is always in slot 1.
! Any DC that has become not "alive" is unbound.
! We do not record a list of T cells bound to each DC, therefore it is not possible to
! handle unbinding by traversing the DC list.
! NOTE: Check that non-cognate cells can make contact with only one DC at a time.
!-----------------------------------------------------------------------------------------
subroutine binder(ok)
logical :: ok
integer :: kcell, nbnd, i, k, site(3), ctype, stage, region, neardc(0:DCDIM-1), idc
integer :: nadd, nsub, nu, nb, nc, nbt, nbc, nt
real :: tnow, bindtime, pMHC, stimrate, t_travel, ttravel
logical :: unbound, cognate, bound
type(cell_type), pointer :: cell
type(cog_type), pointer :: cog_ptr
real, parameter :: fast_bind_prob = 0.25	! TESTING THIS -------------------------------
real, parameter :: killing_time = 10        ! time taken for a DC to die
integer :: kpar = 0

nu = 0
nb = 0
nc = 0
nbt = 0
nbc = 0
nt = 0
nadd = 0
nsub = 0
tnow = istep*DELTA_T
do kcell = 1,nlist
    cell => cellist(kcell)
    if (cell%ID == 0) cycle
    nt = nt+1
    if (evaluate_residence_time .or. track_DCvisits) then
        cognate = .false.
    else
        cog_ptr => cell%cptr
        cognate = (associated(cog_ptr))
    endif
    if (cognate) then
		nc = nc+1
		call get_region(cog_ptr,region)
		if (region /= LYMPHNODE) cycle
	endif

    ! Handle unbinding first
    nbnd = 0
    unbound = .false.
    do k = 1,MAX_DC_BIND
        idc = cell%DCbound(k)
        if (idc == 0) cycle
        nbnd = nbnd + 1
        if (tnow >= cell%unbindtime(k) .or. .not.DClist(idc)%capable) then   ! unbind
            DClist(idc)%nbound = DClist(idc)%nbound - 1
            if (DClist(idc)%nbound < 0) then
                write(logmsg,'(a,5i6,3f8.2)') 'Error: binder: DClist(idc)%nbound < 0: ', &
					kcell,nbnd,idc,cell%DCbound,cell%unbindtime,tnow
                call logger(logmsg)
                ok = .false.
                return
            endif
            if (cognate) then
                DClist(idc)%ncogbound = DClist(idc)%ncogbound - 1
                call RemoveCogbound(idc,kcell)
                if (DClist(idc)%ncogbound < 0) then
                    write(logmsg,'(a,5i6,3f8.2)') 'Error: binder: DClist(idc)%ncogbound < 0: ', &
						kcell,nbnd,idc,cell%DCbound,cell%unbindtime,tnow
	                call logger(logmsg)
		            ok = .false.
			        return
                endif
            endif
            cell%DCbound(k) = 0
            nbnd = nbnd - 1
            unbound = .true.
            cell%signalling = .false.
!           nsub = nsub + 1
        endif
     enddo 
     if (unbound) then  ! ensure that a remaining binding is stored in slot 1
        if (cell%DCbound(1) == 0 .and. cell%DCbound(2) /= 0) then
            cell%DCbound(1) = cell%DCbound(2)
            cell%DCbound(2) = 0
            cell%unbindtime(1) = cell%unbindtime(2)
        endif
        cell%unbindtime(2) = tnow
        if (cognate) nu = nu + 1
    endif

    ! Handle new binding
    bound = .false.
    if (nbnd < MAX_DC_BIND .and. tnow >= cell%unbindtime(2) + DC_BIND_DELAY) then
        if (cognate) then
            if (.not.bindable(cog_ptr)) cycle
		    call get_stage(cog_ptr,stage,region)
            if (stage == NAIVE) then    ! on first binding a cognate cell is ready for next stage
                cog_ptr%stagetime = tnow
            endif
        endif
        site = cell%site
        neardc = occupancy(site(1),site(2),site(3))%DC      ! list of DC near this site
        if (neardc(0) /= 0) then
            do k = 1,neardc(0)
                idc = neardc(k)
                if (.not.DClist(idc)%capable) cycle
                if (idc == cell%DCbound(1)) cycle      ! valid for MAX_DC_BIND <= 2
                if (DClist(idc)%ncogbound > MAX_COG_BIND) then      ! ERROR CHECKING
                    write(logmsg,'(a,i6,i16)') 'Error: binder: ncogbound: ',idc,DClist(idc)%ncogbound
			        call logger(logmsg)
		            ok = .false.
	                return
                endif
                bound = bindDC(idc,kpar)    ! This just means the all-cell limit of the DC has not been reached
                cell%signalling = .false.
                if (FAST .and. (par_uni(kpar) > fast_bind_prob)) then
					bound = .false.
				endif
                if (bound .and. cognate) then
                    if (cog_ptr%effector) then
                        ! DC is killed
                        DClist(idc)%dietime = tnow + killing_time
                        bound = .false.
                        call logger('DC is killed')
                        cycle
                    endif
                    cell%signalling = .true.
                    if (activation_mode == UNSTAGED_MODE) then
                        stimrate = stimulation_rate_norm(DClist(idc)%density,cog_ptr%avidity)
                        if (stimrate < BINDTIME_HILL_THRESHOLD) then
                            cell%signalling = .false.
                            stimrate = 0
!                            cognate = .false.	! the idea is that if the signal strength is < threshold, 
												! it is like a non-cognate interaction.  Use cell%signalling instead
                        endif
                    else
                        stimrate = 0    ! not used
                    endif
                endif
                if (cognate .and. DClist(idc)%ncogbound == MAX_COG_BIND) cycle
                if (bound) then
					if (log_firstDCcontact .and. cell%unbindtime(1) < 0) then
						call logfirstDCcontact(cell,idc)
					endif
	                ttravel = tnow - cell%unbindtime(2)		! this is the time of the previous unbinding
                    nbnd = nbnd + 1
                    if (cell%DCbound(nbnd) /= 0) then
                        write(logmsg,'(a,i6)') 'Error: binder: DCbound(nbnd) /= 0: ',cell%DCbound(nbnd)
		                call logger(logmsg)
				        ok = .false.
						return
                    endif
                    cell%DCbound(nbnd) = idc
                    DClist(idc)%nbound = DClist(idc)%nbound + 1
                    ctype = cell%ctype
                    pMHC = DClist(idc)%density
                    bindtime = get_bindtime(cog_ptr,cell%signalling,ctype,pMHC,stimrate,kpar)
                    cell%unbindtime(nbnd) = tnow + bindtime
                    nadd = nadd + 1
                    if (cognate) then
                        DClist(idc)%ncogbound = DClist(idc)%ncogbound + 1
	                    call AddCogbound(idc,kcell)
	                    call log_count(DCbindtime_count,bindtime)
!	                    write(nflog,'(a,f8.1)') 'bindtime: ',bindtime
	                    call log_count(DCtraveltime_count,ttravel)
                    endif
                    if (cell%DCbound(1) == 0 .and. cell%DCbound(2) /= 0) then
                        write(logmsg,'(a,3i6)') 'Error: binder: DCbound order: ',kcell,cell%DCbound
		                call logger(logmsg)
				        ok = .false.
						return
                    endif
					if (track_DCvisits .and. cell%tag == TAGGED_CELL) then
						call logDCvisit(kcell,idc)
					endif
                    if (nbnd == MAX_DC_BIND) exit
                endif
            enddo
        endif
    endif
    if (bound .and. cognate) nb = nb+1
    if (cell%dcbound(1) /= 0) then
		nbt = nbt + 1
		if (cognate) nbc = nbc + 1
	endif
enddo
!write(*,*) 'me, nadd, nsub: ',me,nadd,nsub,nadd-nsub
!write(*,'(f6.2,i4,2f6.3)') tnow/60,nc,real(nbc)/nc,real(nbt)/nt
!write(nflog,'(f8.2,i8,2f8.3)') tnow/60,nc,real(nbc)/nc,real(nbt)/nt
ok = .true.
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
logical function revisit(ic,idc)
integer :: ic
integer :: idc
integer :: i,n

revisit = .false.
n = cellist(ic)%ndclist
if (n == 0) return
do i = 1,n
	if (cellist(ic)%dclist(i) == idc) then
		revisit = .true.
		return
	endif
enddo
end function

!-----------------------------------------------------------------------------------------
! The DC is either chemoattracting (secreting chemokine) or not.
! The T cell is susceptible to chemotaxis at the LO or the HI level
!-----------------------------------------------------------------------------------------
subroutine logDCvisit(kcell,idc)
integer :: kcell, idc
integer :: iTCchemo, iDCchemo
type(cell_type), pointer :: cell
real :: tnow

cell => cellist(kcell)
if (retain_tagged_cells) then
	tnow = istep*DELTA_T
	if (tnow - cell%entrytime > t_log_DCvisits) return		! log DC visits over 2 days
endif
if (USE_DC_COGNATE) then		! in this case iDCchemo indicates DC chemokine secretion status
	if (DClist(idc)%cognate) then
		iDCchemo = 1
	else
		iDCchemo = 2
	endif
else
	iDCchemo = 1
endif
!if (DC_CHEMO_NOTRAFFIC) then	! also for use_traffic case
!	if (cell%DCchemo < (LO_CHEMO + HI_CHEMO)/2) then
	if (cell%receptor_level(CCR1) < (LO_CHEMO + HI_CHEMO)/2) then
		iTCchemo = 2	! LOW DCchemo
	else
		iTCchemo = 1	! HI DCchemo
	endif
!else
!	iTCchemo = 1
!endif
!if (revisit(kcell,idc)) then
!	Nrevisits(iTCchemo,iDCchemo) = Nrevisits(iTCchemo,iDCchemo) + 1
!	cell%revisits(iDCchemo) = cell%revisits(iDCchemo) + 1
!else
!	Nvisits(iTCchemo,iDCchemo) = Nvisits(iTCchemo,iDCchemo) + 1
!	cell%visits(iDCchemo) = cell%visits(iDCchemo) + 1		! distinct DC visits 
!	cell%ndclist = cell%ndclist + 1
!	cell%dclist(cell%ndclist) = idc
!endif
!write(logmsg,*) 'DC visit: ',kcell,idc
!call logger(logmsg)
!write(*,*) "No good!"
end subroutine

!-----------------------------------------------------------------------------------------
! Now record visits to cognate (bearing cognate antigen, secreting chemokine) and 
! non-cognate (no antigen, no chemokine secretion) DCs separately.
! Also the T cell has either HI (iTCchemo = 1) or LO (iTCchemo = 2) susceptibility.
! Also the T cell is either CD4 or CD8.
!-----------------------------------------------------------------------------------------
subroutine recordDCvisits(kcell)
integer :: kcell
integer :: nvtot, iTCchemo, iDCchemo, ctype, nv
!real :: transit_time, ave_dt, totave_dt

ctype = cellist(kcell)%ctype
ntagged_left = ntagged_left + 1
!transit_time = tnow - cellist(kcell)%entrytime
!if (cellist(kcell)%DCchemo < (LO_CHEMO + HI_CHEMO)/2) then
if (cellist(kcell)%receptor_level(CCR1) < (LO_CHEMO + HI_CHEMO)/2) then
	iTCchemo = 2
else
	iTCchemo = 1
endif
do iDCchemo = 1,2
	nv = cellist(kcell)%visits(iDCchemo)
	DCvisits(nv,iTCchemo,iDCchemo,ctype) = DCvisits(nv,iTCchemo,iDCchemo,ctype) + 1	! distinct visits
	nvtot = cellist(kcell)%visits(iDCchemo) + cellist(kcell)%revisits(iDCchemo)
	DCtotvisits(nvtot,iTCchemo,iDCchemo,ctype) = DCtotvisits(nvtot,iTCchemo,iDCchemo,ctype) + 1
enddo
!ave_dt = transit_time/cellist(kcell)%visits		! mean time between distinct DC visits
!totave_dt = transit_time/nvtot		! mean time between DC visits
end subroutine

!-----------------------------------------------------------------------------------------
! iTC = 1  DCchemo = HI_CHEMO
!     = 2  DCchemo = LO_CHEMO
! iDC = 1  chemoattracting DC
!     = 2  non-chemoattracting DC
!-----------------------------------------------------------------------------------------
subroutine write_DCvisit_dist
integer :: i, nd, nt, ndt, ndd(2,2), ntt(2,2), iTC, iDC, kcell, nTCcases, nDCcases, ctype
real :: total_distinct, total, ave_distinct(2,2), ave_total(2,2), tnow
real :: visit_dist(0:1000,2,2), totvisit_dist(0:1000,2,2)

! First, need to log DC visit data for tagged cells that remain
tnow = istep*DELTA_T
do kcell = 1,nlist
	if (cellist(kcell)%ID == 0) cycle
	if (cellist(kcell)%tag /= TAGGED_CELL) cycle
	if (.not.DC_CHEMO_NOTRAFFIC .and. tnow - cellist(kcell)%entrytime < t_log_DCvisits) cycle
	call recordDCvisits(kcell)
enddo
write(nfout,*)
write(nfout,'(a)') 'DC visits'
write(nfout,'(a,i4)') 'Total DCs: ',NDC
write(nfout,'(a,L2,f4.1)') 'USE_DC_COGNATE, DC_COGNATE_FRACTION: ',USE_DC_COGNATE,DC_COGNATE_FRACTION
write(nfout,'(a,L2,f4.1)') 'DC_CHEMO_NOTRAFFIC, DC_CHEMO_FRACTION: ',DC_CHEMO_NOTRAFFIC,DC_CHEMO_FRACTION
write(nfout,'(a,3f4.1)') 'HI_CHEMO_FRACTION, LO_CHEMO, HI_CHEMO: ', HI_CHEMO_FRACTION, LO_CHEMO, HI_CHEMO

!if (DC_CHEMO_NOTRAFFIC) then
	nTCcases = 2
!else
!	nTCcases = 1
!endif
if (USE_DC_COGNATE .and. DC_COGNATE_FRACTION > 0) then
	nDCcases = 2
else
	nDCcases = 1
endif

write(nfout,'(a,i6)') 'Number of tagged cells that contributed: ',ntagged_left

do ctype = 1,2
	visit_dist = 0
	totvisit_dist = 0
	ave_distinct = 0
	ave_total = 0
	if (CTYPE_FRACTION(ctype) == 0) cycle
	if (ctype == CD4) then
		write(nfout,*) '-----------------------------------------------------------------------------'
		write(nfout,'(a)') 'CD4 cells'
		write(nfout,*) '-----------------------------------------------------------------------------'
	else
		write(nfout,*) '-----------------------------------------------------------------------------'
		write(nfout,'(a)') 'CD8 cells'
		write(nfout,*) '-----------------------------------------------------------------------------'
	endif
	do iTC = 1,nTCcases
!		write(nfout,*)
!		if (iTC == 1) then
!			write(nfout,'(a,f4.1)') 'High chemo-susceptibility: ',HI_CHEMO
!		else
!			if (.not.DC_CHEMO_NOTRAFFIC) exit
!			write(nfout,'(a,f4.1)') 'Low chemo-susceptibility: ',LO_CHEMO
!		endif
		do iDC = 1,nDCcases
!			if (iDC == 1) then
!				write(nfout,'(a)') 'Visits to cognate (attracting) DCs'
!			else
	!			if (.not.USE_DC_COGNATE) exit
	!			write(nfout,*)
!				write(nfout,'(a)') 'Visits to non-cognate (non-attracting) DCs'
!			endif
			total_distinct = 0
			do i = NDC,0,-1
				write(*,*) i,iTC,iDC,ctype
				total_distinct = total_distinct + DCvisits(i,iTC,iDC,ctype)
			enddo
			nd = 0
			if (total_distinct > 0) then
				do i = NDC,0,-1
					if (nd == 0 .and. DCvisits(i,iTC,iDC,ctype)/total_distinct >= 0.0001) then
						nd = i
					endif
				enddo
	!			write(nfout,'(a)') 'DC distinct visits distribution'
				do i = 0,nd
	!				write(nfout,'(2i6,f8.4)') i,DCvisits(i,iTC,iDC),DCvisits(i,iTC,iDC)/total_distinct
					visit_dist(i,iTC,iDC) = DCvisits(i,iTC,iDC,ctype)/total_distinct
					ave_distinct(iTC,iDC) = ave_distinct(iTC,iDC) + i*DCvisits(i,iTC,iDC,ctype)/total_distinct
				enddo
			endif
	!		write(nfout,'(a,f8.2)') 'Average number of distinct DC contacts: ',ave_distinct
	!		write(nfout,*)
			total = 0
			do i = 0,1000
				total = total + DCtotvisits(i,iTC,iDC,ctype)
			enddo
			nt = 0
			if (total > 0) then
				do i = 1000,0,-1
					if (nt == 0 .and. DCtotvisits(i,iTC,iDC,ctype)/total >= 0.0001) then
						nt = i
					endif
				enddo
	!			write(nfout,'(a)') 'DC total visits distribution'
				do i = 0,nt
	!				write(nfout,'(2i6,f8.4)') i,DCtotvisits(i,iTC,iDC),DCtotvisits(i,iTC,iDC)/total
					totvisit_dist(i,iTC,iDC) = DCtotvisits(i,iTC,iDC,ctype)/total
					write(*,*) i,iTC,iDC,ctype
					ave_total(iTC,iDC) = ave_total(iTC,iDC) + i*DCtotvisits(i,iTC,iDC,ctype)/total
				enddo
			endif
	!		write(nfout,'(a,f8.2)') 'Average number of total DC contacts: ',ave_total
	!		write(nfout,*)
			ndd(iTC,iDC) = nd
			ntt(iTC,iDC) = nt
			write(*,*) 'iTC, iDC, nd, nt: ',iTC, iDC, nd, nt
			write(nfout,*) 'iTC, iDC, nd, nt: ',iTC, iDC, nd, nt
		enddo
	enddo
	ndt = 0
	do iTC = 1,2
		do iDC = 1,2
			ndt = max(ndt,ndd(iTC,iDC))
			ndt = max(ndt,ntt(iTC,iDC))
		enddo
	enddo
	write(nfout,'(a,4f6.1)') 'Average distinct visits: ',((ave_distinct(iTC,iDC),iTC=1,2),iDC=1,2)
	write(nfout,'(a,4f6.1)') 'Average total visits:    ',((ave_total(iTC,iDC),iTC=1,2),iDC=1,2)
	write(nfout,'(a)') 'Distributions:'
	do i = 0,ndt
		write(nfout,'(i4,8f8.4)') i,((visit_dist(i,iTC,iDC),iTC=1,2),iDC=1,2),((totvisit_dist(i,iTC,iDC),iTC=1,2),iDC=1,2)
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine write_Tres_dist
integer :: k, kmax
real(8) :: Tres

do k = int(days*24),1,-1
    if (Tres_dist(k) /= 0) then
        kmax = k
        exit
    endif
enddo
Tres = 0
do k = 1,kmax
    Tres = Tres + (k-0.5)*Tres_dist(k)/noutflow_tag
enddo
write(nfout,'(a)') 'Transit time distribution'
write(nfout,'(a)') 'Parameters: '
if (TAGGED_EXIT_CHEMOTAXIS) then
	write(nfout,'(a,L)') '  TAGGED_EXIT_CHEMOTAXIS:     ',TAGGED_EXIT_CHEMOTAXIS
	write(nfout,'(a,f6.3)') '  TAGGED_CHEMO_FRACTION: ',TAGGED_CHEMO_FRACTION 
	write(nfout,'(a,f6.3)') '  TAGGED_CHEMO_ACTIVITY: ',TAGGED_CHEMO_ACTIVITY
endif
write(nfout,'(a,f6.3)') '  K1_S1PR1:               ',K1_S1PR1
write(nfout,'(a,3i8)') 'Results: noutflow_tag, ninflow_tag,kmax: ',noutflow_tag,ninflow_tag,kmax
write(nfout,'(a,f6.1)') 'Residence time: ',Tres
!write(nfout,'(a,f6.1)') 'Residence time from restime_tot: ',restime_tot/60.
write(nfout,'(a)') ' Hour   Probability'
do k = 1,kmax
    write(nfout,'(f6.1,e12.4)') k-0.5,Tres_dist(k)/noutflow_tag
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine write_firstDC_dist
integer :: k, kmax
integer :: iDCchemo

write(nfout,'(a)') 'First DC contacts'
write(nfout,'(a)') '1 = HI_CHEMO, 2 = LO_CHEMO'
do iDCchemo = 1,2
	if (iDCchemo == 1) then
		write(nfout,'(a)') 'With cognate (attracting) DC'
		if (DC_COGNATE_FRACTION == 0) then
			write(nfout,*) 'No cognate DCs'
			cycle
		endif
	else
		if (.not.USE_DC_COGNATE) exit
		write(nfout,'(a)') 'With non-cognate (non-attracting) DC'
		if (DC_COGNATE_FRACTION == 1) then
			write(nfout,*) 'No non-cognate DCs'
			cycle
		endif
	endif
	write(nfout,'(a,2i8)') 'Number: ',firstDC_n(:,iDCchemo)
	write(nfout,'(a,2f6.1)') 'Mean time: ',firstDC_tot(1,iDCchemo)/firstDC_n(1,iDCchemo), &
										   firstDC_tot(2,iDCchemo)/firstDC_n(2,iDCchemo)
	do k = 10000,1,-1
		if (firstDC_dist(1,iDCchemo,k) /= 0 .or. firstDC_dist(2,iDCchemo,k) /= 0) then
			kmax = k
			exit
		endif
	enddo
	write(nfout,'(a)') 'Min  Pr(1)   Pr(2)'
	do k = 1,kmax
		write(nfout,'(i4,2f8.3)') k,firstDC_dist(1,iDCchemo,k)/firstDC_n(1,iDCchemo), &
									firstDC_dist(2,iDCchemo,k)/firstDC_n(2,iDCchemo)
	enddo
	write(nfout,*)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! The distributions gendist, tcrdist and aviddist are all for the current DCU population.
! avid_count is for the recent efferent population.
!-----------------------------------------------------------------------------------------
subroutine write_results
integer :: kcell, ctype, stype, ncog, ntot, i
integer :: gen
real :: tcr, avid, dtcr, hour
type (cog_type), pointer :: p
integer :: gendist(TC_MAX_GEN),aviddist(MAX_AVID_LEVELS),tcrdist(tcr_nlevels)
character*(60) :: fmtstr = '(f6.2,2i8,4x,15f7.4,4x,10f7.4,4x,10f7.4,4x,10i7)'

write(fmtstr(14:15),'(i2)') TC_MAX_GEN
write(fmtstr(24:25),'(i2)') tcr_nlevels
write(fmtstr(34:35),'(i2)') avidity_nlevels
write(fmtstr(44:45),'(i2)') avidity_nlevels
hour = istep*DELTA_T/60
dtcr = TCR_limit/TCR_nlevels
gendist = 0
aviddist = 0
tcrdist = 0
ntot = 0
ncog = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    p => cellist(kcell)%cptr
    ntot = ntot + 1
    ctype = cellist(kcell)%ctype
!    stype = struct_type(ctype)
	if (associated(cellist(kcell)%cptr)) then
		stype = COG_TYPE_TAG
	else
		stype = NONCOG_TYPE_TAG
	endif
    if (stype == COG_TYPE_TAG) then
        ncog = ncog + 1
        ! TCR stimulation distribution
        tcr = p%stimulation
        i = tcr/dtcr + 1
        i = min(i,TCR_nlevels)
        tcrdist(i) = tcrdist(i) + 1
        ! T cell generation distribution
        gen = get_generation(p)
        gendist(gen) = gendist(gen) + 1
        ! T cell avidity distribution
        avid = p%avidity
        if (avidity_logscale) then
            avid = log10(avid)
        endif
        if (avidity_nlevels == 1) then
            i = 1
        else
            i = (avid-avidity_min)/avidity_step + 1.5
!           write(nfout,'(a,2f8.4,3i7)') 'Count: ',p%avidity,avid,i,kcell,p%cogID
            i = max(i,1)
            i = min(i,avidity_nlevels)
        endif
        aviddist(i) = aviddist(i) + 1
    endif
enddo
if (fix_avidity) then
    write(nfres,fmtstr) hour,ntot,ncog,real(gendist)/ncog,real(tcrdist)/ncog,real(aviddist)/ncog, &
        avid_count%bincount
else
    write(nfres,fmtstr) hour,ntot,ncog,real(gendist)/ncog,real(tcrdist)/ncog
endif
!avid_count_total%bincount = avid_count_total%bincount + avid_count%bincount
avid_count%bincount = 0
!write(nfout,'(8i6)') aviddist
end subroutine

!--------------------------------------------------------------------------------------
! Motility behaviour
! There are potentially 4 motility states:
! (1) Naive cell, very short non-cognate DC contacts 3 min (motility level 1)
! (2) Short cognate DC contacts 11 min (motility level 2)
! (3) Long cognate DC contacts > 1 hr clusters (motility level 3)
! (4) Short cognate DC contacts 18 min  swarms (motility level 4)
! It isn't clear what level of motility occurs at each stage. Possibly
! the contact duration is related to the motility, since reduced motility
! leads to a higher probability of rebinding.
! A simple way to vary motility is to keep the persistence parameter rho fixed
! and vary alpha (plongjump), the probability of 2 steps, and possibly beta,
! the probability of moving at all.  In any case the jump parameters for each case
! are stored as pjump(2) and dirprob(0:6).
! In principle: speed & Cm -> rho, beta, delta -> dirprob(), pjump().
!
! Mark Miller says that their motility measurements didn't exclude periods when
! T cells were in contact with DC, therefore we have no info about motility in
! various stages of activation.  For now it is safest to use the same motility
! parameters for all stages.
!--------------------------------------------------------------------------------------

!---------------------------------------------------------------------
! TCR stimulation rate depends on the DC's pMHC count and the T cell's
! TCR avidity for the antigen.
! New formulation of TCR stimulation and bind-time.
! The DC antigen density is a measure of the number of pMHC complexes on the cell surface
! For now, this is set initially proportional to the concentration of peptide (pM) used to "pulse" DCs
! Henrickson2008 measured 2.5x10^4 = 25000 pMHC/DC after pulsing 5T33 DCs for 3 hr with 10 uM M- or C-peptide
! The half-life of the pMHC complex was estimated at 6.01 hrs for M-peptide and 2.36 hrs for C-peptide
! Therefore the initial pMHC counts after 18 hours were 25000/8 = 3125 for M-peptide and 130 for C-peptide
! Since 100 pM of M-peptide induced only about 10% proliferation, while 200 pM induced about 80%, it seems
! that the threshold for TCR signalling is about 30 pMHC (assuming that the initial pMHC count is linear
! with antigen concentration).
! THRESHOLD CHANGED.  It is now a stimulation rate threshold, i.e. on level of pMHC*avidity
! Base stimulation rate r is a Hill function of a = pMHC*avidity, max value 1
! Actual rate of change of S is this rate scaled by TC_STIM_RATE_CONSTANT
! THIS HAS BEEN CHANGED.  TC_STIM_RATE_CONSTANT is now effectively 1.0, and default threshold values have 
! been scaled to give the same results as if TC_STIM_RATE_CONSTANT = 5.
! Duration of binding is also determined from r.  Currently the specified bind time is
! treated as the maximum, and the actual value depends linearly on r up to this max.
! This formulation is used in the STAGED_MODE simulations
! NOTE:
! It makes sense to work with normalized pMHC and avidity in this case too 
! (as in UNSTAGED case with stimulation_rate_norm).  This will make it easier to maintain consistency
! with the activation threshold values.
!----------------------------------------------------------------------------------------------------------
real function stimulation_rate_hill(pMHC,avidity)
real :: pMHC, avidity
real :: a, d, x
!a = max(pMHC - STIM_HILL_THRESHOLD,0.0)*avidity
a = min(1.0,avidity/MAXIMUM_AVIDITY)
d = min(1.0,pMHC/MAXIMUM_ANTIGEN)
x = a*d
if (x > STIM_HILL_THRESHOLD) then
    stimulation_rate_hill = (1 + STIM_HILL_C**STIM_HILL_N)*x**STIM_HILL_N/(x**STIM_HILL_N + STIM_HILL_C**STIM_HILL_N)
else
    stimulation_rate_hill = 0
endif
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
real function stimulation_rate_norm(pMHC,avidity)
real :: pMHC, avidity
real :: a, d
a = min(1.0,avidity/MAXIMUM_AVIDITY)
d = min(1.0,pMHC/MAXIMUM_ANTIGEN)
stimulation_rate_norm = a*d
end function

!---------------------------------------------------------------------
! Updates T cell state, and DC %stimulation (for use in modifying
! %density in update_DCstate).
!---------------------------------------------------------------------
subroutine updater(ok)
logical :: ok
integer :: kcell, ctype, stype, region, iseq, tag, kfrom, kto, k, ncog, ntot
integer :: site(3), site2(3), freeslot, indx(2), status, DC(2), idc
real :: C(N_CYT), mrate(N_CYT), tnow, dstim, S, cyt_conc, mols_pM, Ctemp, dstimrate, stimrate
logical :: divide_flag, producing, first, dbg, unbound, flag, flag1
!logical, save :: first = .true.
type (cog_type), pointer :: p

!write(*,*) 'updater: ',me
ok = .true.
dbg = .false.
flag = .false.
flag1 = .false.
ntot = 0
ncog = 0
tnow = istep*DELTA_T
! Scaling factor to convert total number of molecules in the region to conc in pM
mols_pM = L_um3*M_pM/(NTcells*Vc*Navo)

do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    ntot = ntot + 1
    if (dbg) write(*,*) 'kcell: ',kcell
    ctype = cellist(kcell)%ctype
!    stype = struct_type(ctype)
	if (associated(cellist(kcell)%cptr)) then
		stype = COG_TYPE_TAG
	else
		stype = NONCOG_TYPE_TAG
	endif
    if (stype /= COG_TYPE_TAG) cycle
    ! Only cognate cells considered
    ncog = ncog + 1
    p => cellist(kcell)%cptr
    if (.not.associated(p)) then
        write(*,*) 'ERROR: updater: p not associated: ',kcell
        stop
    endif
    call get_region(p,region)
	! Cell death
    if (tnow > p%dietime) then
        call Tcell_death(kcell)
        cycle
    endif
	! TCR stimulation decay
    p%stimulation = p%stimulation*(1 - TCRdecayrate*DELTA_T)

	if (region == LYMPHNODE) then
	! TCR stimulation
		unbound = .false.
		DC = cellist(kcell)%DCbound
		stimrate = 0
		do k = 1,MAX_DC_BIND
			idc = DC(k)
	        dstimrate = 0
			if (idc /= 0) then
				if (DClist(idc)%capable) then
	!               dstimrate = TC_STIM_RATE_CONSTANT*DClist(idc)%density*p%avidity
	!				dstimrate = TC_STIM_RATE_CONSTANT*stimulation_rate(DClist(idc)%density,p%avidity)
	                if (activation_mode == STAGED_MODE) then
    					dstimrate = stimulation_rate_hill(DClist(idc)%density,p%avidity) ! now use THRESHOLD_FACTOR
    				elseif (cellist(kcell)%signalling) then
    				    dstimrate = stimulation_rate_norm(DClist(idc)%density,p%avidity)
    				    if (dstimrate < BINDTIME_HILL_THRESHOLD) then
    				        dstimrate = 0
    				    endif
    				endif
					stimrate = stimrate + dstimrate
					dstim = dstimrate*DELTA_T
					p%stimulation = p%stimulation + dstim
!					if (p%stimulation > FIRST_DIVISION_THRESHOLD(1)) then
!					    write(logmsg,*) 'Reached FIRST_DIVISION_THRESHOLD: ',kcell,tnow
!					    call logger(logmsg)
!					endif
					p%stimulation = min(p%stimulation, STIMULATION_LIMIT)
					DClist(idc)%stimulation = DClist(idc)%stimulation + dstim
					p%totalDCtime = p%totalDCtime + DELTA_T
				else    ! unbind T cell from incapable DC
					idc = cellist(kcell)%DCbound(k)
					DClist(idc)%nbound = DClist(idc)%nbound - 1
					DClist(idc)%ncogbound = DClist(idc)%ncogbound - 1
                    call RemoveCogbound(idc,kcell)
					cellist(kcell)%DCbound(k) = 0
					unbound = .true.
				endif
			endif
		enddo
		p%stimrate = stimrate
		if (unbound .and. cellist(kcell)%DCbound(1) == 0 .and. cellist(kcell)%DCbound(2) /= 0) then
			cellist(kcell)%DCbound(1) = cellist(kcell)%DCbound(2)
			cellist(kcell)%DCbound(2) = 0
			cellist(kcell)%unbindtime(1) = cellist(kcell)%unbindtime(2)
			cellist(kcell)%unbindtime(2) = tnow
		endif

		site = cellist(kcell)%site
		if (use_cytokines) then
			! IL receptor stimulation
			status = p%status
			if (use_diffusion) then
				C = cyt(site(1),site(2),site(3),:)
			endif
			S = p%stimulation
			do iseq = 1,Ncytokines
				tag = cyt_tag(iseq)
				kfrom = NP_offset(iseq)+1
				kto = NP_offset(iseq+1)
				if (use_diffusion) then
					cyt_conc = C(iseq)
				else
					cyt_conc = cyt_mols(iseq)*mols_pM   ! -> conc in pM
				endif
				ctemp = cyt_conc
				select case(tag)
				case(IL2_TAG)
					producing = IL2_production_status(p,tnow)
					if (p%IL_state(kfrom) == 0) then    ! temporary measure to detect first call of IL2_update
						first = .true.
					else
						first = .false.
					endif
					call IL2_update(p%cogID,ctype,tnow,S,p%IL_state(kfrom:kto),p%IL_statep(kfrom:kto), &
						first,producing,cyt_conc,Vc,DELTA_T,mrate(iseq),flag1)
				case(IL4_TAG)
					call IL4_update(p%IL_state(kfrom:kto))
				case(IL7_TAG)
					call IL7_update(p%IL_state(kfrom:kto))
				case(IL9_TAG)
					call IL9_update(p%IL_state(kfrom:kto))
				case(IL15_TAG)
					call IL15_update(p%IL_state(kfrom:kto))
				case(IL21_TAG)
					call IL21_update(p%IL_state(kfrom:kto))
				end select

				if (use_diffusion) then
	! Concentration units
	! mrate = mass rate of flow in molecules/min
	! mrate*L_um3/Vc = molecules/L/min
	! mrate*L_um3*DELTA_T/Vc = molecules/L
	! mrate*L_um3*DELTA_T/Vc/Navo = moles/L = M
	! mrate*L_um3*DELTA_T*M_pM/Vc/Navo = pM
	! To convert total number of molecules in the region to conc in pM
	! mols * L_um3*M_pM/(NTcells*Vc*Navo)
					C(iseq) = cyt_conc
	!                C(iseq) = C(iseq) + (mrate(iseq)/Vc)*DELTA_T*M_pM*L_um3/Navo    ! Vc/L_um3 ->free vol in L
					if (C(iseq) < 0) then
						write(*,'(a,6i6,4f8.4)') 'WARNING: cyt < 0: ',kcell,p%cogID,iseq,site,Ctemp,C(iseq)
						C(iseq) = 0
					endif
				else
					! increment total number of molecules of this cytokine
					dcyt_mols(iseq) = dcyt_mols(iseq) + (cyt_conc - ctemp)/(mols_pM*NTcells)
	!                write(*,*) cyt_mols(iseq),dcyt_mols(iseq),ctemp,cyt_conc
				endif
			enddo
			if (use_diffusion) then
				cyt(site(1),site(2),site(3),:) = C
			endif
		endif

		if (.not.L_selectin) then
			! Note that stimrate is normalized with TC_STIM_RATE_CONSTANT for consistency
			! with the calibration of S1PR1/CD69 dynamics, which was carried out with
			! normalized stimulation of approx. 1  This is easier than adjusting K1_CD69
			call S1PR1_update(p%CD69,p%S1PR1,p%stimrate/TC_STIM_RATE_CONSTANT,DELTA_T)
	!        call S1PR1_update(p%CD69,p%S1PR1,p%stimrate,DELTA_T)
		endif
	endif

! Stage transition
    call updatestage(kcell, tnow, divide_flag)

! Cell division
    if (divide_flag) then
		if (region == LYMPHNODE) then
			indx = occupancy(site(1),site(2),site(3))%indx
			freeslot = 0
			if (indx(1) == 0) then
				site2 = site
				freeslot = 1
			elseif (indx(2) == 0) then
				site2 = site
				freeslot = 2
			else
				call get_free_slot(occupancy,NX,site,site2,freeslot)
			endif
		else
			site2 = 0
			freeslot = -1
!			call logger("T cell divides in the periphery")
		endif
        if (freeslot /= 0) then     ! there is a free slot at site2 (which may be = site)
            call cell_division(kcell,site2,freeslot,ok)
            if (.not.ok) return
        endif
    endif

    if (flag) then
        call show_cognate_cell(kcell)
    endif

enddo
call UpdateHelp
end subroutine

!--------------------------------------------------------------------------------
! The idea is to record the time spent by a cognate CD8 cell bound to a DC at the
! same time as a cognate activated CD4 cell.
! At this time it is not clear how or when help is effective.  For a start we can
! record total timesteps.
!--------------------------------------------------------------------------------
subroutine UpdateHelp
integer :: idc, ncb, k1, k2, n1, n2, kcell1, kcell2
integer :: ctype1, ctype2
logical :: help

do idc = 1,NDC
	ncb = DClist(idc)%ncogbound
	if (ncb == 0) cycle
	n1 = 0
	do k1 = 1,MAX_COG_BIND
		kcell1 = DClist(idc)%cogbound(k1)
		if (kcell1 == 0) cycle
		if (.not.associated(cellist(kcell1)%cptr)) then
			write(*,*) 'Error: UpdateHelp: cptr not associated: ',kcell1
			stop
		endif
		n1 = n1+1
		ctype1 = cellist(kcell1)%ctype
		help = .false.
		n2 = 0
		do k2 = 1,MAX_COG_BIND
			kcell2 = DClist(idc)%cogbound(k2)
			if (kcell2 == 0) cycle
			if (.not.associated(cellist(kcell2)%cptr)) then
				write(*,*) 'Error: UpdateHelp: cptr not associated: ',kcell2
				stop
			endif
			n2 = n2+1
			ctype2 = cellist(kcell2)%ctype
			if (ctype1 == CD4 .and. ctype2 == CD8) then
				help = .true.
				exit
			elseif (ctype1 == CD8 .and. ctype2 == CD4) then
				help = .true.
				exit
			endif
			if (n2 == ncb) exit
		enddo
		if (help) then
!			write(*,*) 'help! ',idc,k1,k2,kcell1,kcell2
			cellist(kcell1)%cptr%cnt(1) = cellist(kcell1)%cptr%cnt(1) + 1
		endif
		if (n1 == ncb) exit
	enddo
enddo
end subroutine





!----------------------------------------------------------------------------------------
! A cell gets permission to divide when it is in the ACTIVATED stage and the activation
! (weighted sum of TCR stimulation and CD25/IL2 signal) exceeds a threshold.
! Note that when optionA = 1, the CD25 signal alone must also exceed a threshold for division,
! and if CD25_SWITCH is true and gen = 1 failure to reach the thresholds cancels division
! permanently.
! In this revised version, the use of a weighted sum of signals for stimulation (with the
! parameter TC_STIM_WEIGHT) has been separated from the optionA cases.
!----------------------------------------------------------------------------------------
logical function candivide(p,ctype)
type(cog_type), pointer :: p
integer :: ctype
integer :: gen
real :: tnow, div_thresh, stim, CD25signal
!
! NOTE: Need to better account for first division time, and need to check CD4/CD8
!
tnow = istep*DELTA_T
candivide = .false.
gen = get_generation(p)
if (gen == 1) then			! undivided cell
	div_thresh = FIRST_DIVISION_THRESHOLD(ctype)
    if (p%stimulation > div_thresh) then
		call set_activation(p,1)
	endif
else											! clone, lower threshold for division
	div_thresh = DIVISION_THRESHOLD(ctype)
endif
if (activation_mode == UNSTAGED_MODE) then
    if (gen == 1 .and. tnow - p%firstDCtime < UNSTAGED_MIN_DIVIDE_T) return
    if (p%stimulation < div_thresh) return
    candivide = .true.
    return 
endif

stim = get_stimulation(p)     ! weighted sum of TCR signal and CD25/IL2 signal
if (optionA == 1) then
    CD25signal = get_IL2store(p)
    if (stim > div_thresh .and. CD25signal > CD25_DIVISION_THRESHOLD) then ! allow division
	    if (gen > TC_MAX_GEN) then
		    write(*,*) 'candivide: bad gen: ',gen
		    stop
	    endif
	    candivide = .true.
    elseif (gen == 1 .and. CD25_SWITCH) then
    !	Tcell(icell)%cog%IL2_pass = .false.		! not enough IL2, division is cancelled
    !	Tcell(icell)%cog%IL2status = ibclr(Tcell(icell)%cog%IL2status,IL2_FLAG3_BIT)
    !	write(*,*) 'Division cancelled for: ',kcell
        call set_stage(p,FINISHED)
	    p%stagetime = BIG_TIME
	endif
elseif (optionA == 2) then
    if (stim > div_thresh) then ! allow division
	    if (gen > TC_MAX_GEN) then
		    write(*,*) 'candivide: bad gen: ',gen
		    stop
	    endif
	    candivide = .true.
	endif
endif
end function

!----------------------------------------------------------------------------------------
! For now base this decision on the cytokine production threshold.
! This must happen only once for a given cell.(?)
!----------------------------------------------------------------------------------------
logical function reached_IL2_threshold(p)
type(cog_type), pointer :: p

!if (p%stimulation > cytokine(IL2_CYT)%cyt_threshold) then
if (p%stimulation > IL2_THRESHOLD) then
	reached_IL2_threshold = .true.
else
	reached_IL2_threshold = .false.
endif
end function

!----------------------------------------------------------------------------------------
! This must happen only once for a given cell.(?)
!----------------------------------------------------------------------------------------
logical function reached_act_threshold(p)
type(cog_type), pointer :: p
real :: stimulation

stimulation = get_stimulation(p)
if (stimulation > ACTIVATION_THRESHOLD) then
	reached_act_threshold = .true.
else
	reached_act_threshold = .false.
endif
end function

!----------------------------------------------------------------------------------------
! T cell activation is determined as the weighted sum of TCR stimulation and IL2_Store,
! which is the current level of integrated CD25 signal.
!----------------------------------------------------------------------------------------
real function get_stimulation(p)
type(cog_type), pointer :: p

get_stimulation = TC_STIM_WEIGHT*p%stimulation + (1-TC_STIM_WEIGHT)*get_IL2store(p)
end function

!----------------------------------------------------------------------------------------
! IL2_production_status is .true. if IL-2 and CD25 are being produced, else .false.
! OptionB:
! = 1  IL2 is produced for a maximum period of IL2_PRODUCTION_TIME in gen = 1 only
! = 2  IL2 is produced in gen = 1 only
! = 3  IL2 is produced always
!----------------------------------------------------------------------------------------
logical function IL2_production_status(p,t)
type(cog_type), pointer :: p
real :: t
integer :: gen, stage, region

IL2_production_status = .false.
!stage = get_stage(p)
call get_stage(p,stage,region)
if (stage == FINISHED) then
    return
endif
gen = get_generation(p)
select case (optionB)
case (1)
    if (t < IL2_PRODUCTION_TIME*60 .and. gen == 1) then
        IL2_production_status = .true.
    else
        IL2_production_status = .false.
    endif
case (2)
    if (gen == 1) then
        IL2_production_status = .true.
    else
        IL2_production_status = .false.
    endif
case (3)
    IL2_production_status = .true.
end select
end function

!----------------------------------------------------------------------------------------
! If optionC = 1, an activated cell cannot survive if IL2 store drops below
! CD25_SURVIVAL_THRESHOLD
!----------------------------------------------------------------------------------------
logical function cansurvive(p)
type(cog_type), pointer :: p
real :: CD25signal

cansurvive = .true.
if (optionC == 1) then
    CD25signal = get_IL2Store(p)
    if (CD25signal < CD25_SURVIVAL_THRESHOLD) then
        cansurvive = .false.
    endif
endif
end function

!---------------------------------------------------------------------
! Remove the last nremex exits in the list.
! Perhaps it would make more sense to remove exits that are close to 
! others, or to adjust the locations of remaining exits.
! NOT USED
!---------------------------------------------------------------------
subroutine removeExits1(nremex)
integer :: nremex
integer :: Nex, k, iexit, site(3)

!write(*,*) 'removeExits: ',nremex
Nex = lastexit
k = 0
do iexit = Nex,1,-1
    if (exitlist(iexit)%ID == 0) cycle
    k = k+1
    site = exitlist(iexit)%site
!    write(logmsg,'(a,4i6)') 'removeExits: removeExitPortal: ',iexit,site
!    call logger(logmsg)
    call RemoveExitPortal(site)
!    write(logmsg,'(a,4i6)') 'did removeExitPortal'
!    call logger(logmsg)
    if (k == nremex) exit
enddo
!call checkExits("in removeExits")
end subroutine

!-----------------------------------------------------------------------------------------
! Determine whether a cognate T cell is licensed to exit the DCU.
! Possible exit rules:
! (0) gen >= NGEN_EXIT
! (1) gen > 1 and act < EXIT_THRESHOLD
! (2) S1PR1 > S1PR1_EXIT_THRESHOLD
! (3) Allow exit to any cognate cell that gets there 
!-----------------------------------------------------------------------------------------
logical function exitOK(p,ctype)
type(cog_type), pointer :: p
integer :: ctype
integer :: gen
real :: act
!integer :: kpar = 0

gen = get_generation(p)

exitOK = .false.
if (exit_rule == EXIT_GEN_THRESHOLD) then
    if (gen >= NGEN_EXIT) then
        exitOK = .true.
    endif
elseif (exit_rule == EXIT_STIM_THRESHOLD) then
    if (gen > 1) then
        act = get_stimulation(p)
!        cd4_8 = ctype - 1
        if (act < EXIT_THRESHOLD(ctype)) then
            exitOK = .true.
        endif
    endif
elseif (exit_rule == EXIT_S1PR1_THRESHOLD) then
    if (p%S1PR1 > S1PR1_EXIT_THRESHOLD) then
        exitOK = .true.
    else
        exitOK = .false.
    endif
elseif (exit_rule == EXIT_UNLIMITED) then
    exitOK = .true.
endif
end function

!-----------------------------------------------------------------------------------------
! If RELAX_INLET_EXIT_PROXIMITY ingress locations are restricted to a layer near the
! blob surface, of thickness INLET_LAYER_THICKNESS.  The distance from an exit portal
! must not be less than INLET_EXIT_LIMIT.
!-----------------------------------------------------------------------------------------
logical function inletOK(site)
integer :: site(3)
real :: d

d = cdistance(site)
if (radius - d > INLET_LAYER_THICKNESS) then
	inletOK = .false.
	return
endif
if (tooNearExit(site,INLET_EXIT_LIMIT)) then
	inletOK = .false.
	return
endif
inletOK = .true.
end function

!-----------------------------------------------------------------------------------------
! This version aims to determine outflow on the basis of the number of exit sites (which
! depends on NTcells) and on how many cells arrive at the exits, which also involves the
! calibration parameter chemo_K_exit and possibly another.
! The number of exits is greater than 1/2 the desired outflow, and the calibration
! parameter exit_prob is used to adjust the mean outflow.
! PROBLEM:
! It seems that CHEMO_K_EXIT has no effect on the exit rate
! For some reason, CHEMO_K_EXIT = 0 is NOT the same as USE_EXIT_CHEMOTAXIS = false,
! but it SHOULD be.
! SOLVED.  The reason was that different coefficients were being used to compute
! requiredExitPortals.
!
! Enhanced egress
! ---------------
! To cover the possibility of activated cells having an enhanced probability of egress,
! the exit region can be expanded to some neighbourhood of the portal site, e.g.
! the Moore27 neighbourhood.  Any cell within the expanded region that meets a
! specified criterion (e.g. generation > 4) has the possibility of egress.
!
! The number of exit portals and the base exit_prob value (0.02) are based on the residence
! time for CD4 cells.  For CD8 cells the exit probability is reduced, so that the exit probs
! are in inverse ratio to the residence times.
!   Pe(1)/Pe(2) = Tres(2)/Tres(1)
!-----------------------------------------------------------------------------------------
subroutine portal_traffic(ok)
logical :: ok
integer :: iexit, esite(3), esite0(3),site(3), i, j, k, slot, indx(2), kcell, ne, ipermex, ihr, nv, ihev, it
integer :: x, y, z, ctype, gen, region, node_inflow, node_outflow, net_inflow, iloop, nloops, nbrs, tag, site0(3)
integer :: ncogin, noncogin
real(DP) :: R, df, cogin, p_CD4
logical :: left, cognate, isbound, hit
logical :: central, egress_possible, use_full_neighbourhood, tagcell
integer, allocatable :: permex(:)
integer :: kpar = 0
real :: tnow, ssfract, exfract
real :: exit_prob(2)

ok = .true.
region = LYMPHNODE
node_inflow = InflowTotal
node_outflow = OutflowTotal
df = InflowTotal - node_inflow
R = par_uni(kpar)
if (R < df) then
	node_inflow = node_inflow + 1
endif

ncogin = 0
noncogin = 0
if (FAST) then     !inflow of cognate cells is accounted for explicitly, noncognate inflow just increases NTcells
    cogin = InflowTotal*(CTYPE_FRACTION(CD4)*TC_COGNATE_FRACTION(CD4) + CTYPE_FRACTION(CD8)*TC_COGNATE_FRACTION(CD8))
    ncogin = cogin
    df = cogin - ncogin
    R = par_uni(kpar)
    if (R < df) then
	    ncogin = ncogin + 1
    endif
    noncogin = node_inflow - ncogin
    node_inflow = ncogin
endif

if (L_selectin) then
	node_inflow = 0
endif

if (steadystate) then
	ssfract = real(NTcells)/NTcells0
	if (ssfract > 1.01) then
		node_inflow = node_inflow - 1
	elseif (ssfract < 0.99) then
		node_inflow = node_inflow + 1
	endif
!		exit_prob = ssfract*exit_prob
else
	df = OutflowTotal - node_outflow
	R = par_uni(kpar)
	if (R < df) then
		node_outflow = node_outflow + 1
	endif
endif
if (computed_outflow) then
	nloops = 5
else
	nloops = 1
endif

tnow = istep*DELTA_T
exfract = 1
if (transient_egress_suppression) then
	exfract = egressFraction(tnow)
endif

! Inflow
if (exit_region == EXIT_LOWERHALF) then    
	write(*,*) 'EXIT_LOWERHALF not simulated with chemotaxis'
	stop
endif
gen = 1
k = 0
it = 0
do while (k < node_inflow)
	it = it + 1
	if (it > 10000) then
		write(*,*) 'portal_traffic: unable to place cell'
		ok = .false.
		return
	endif
	if (use_HEV_portals) then
		ihev = random_int(1,NHEV,kpar)
		if (ihev < 0) write(*,*) 'k,NHEV,ihev: ',k,NHEV,ihev
		site0 = HEVlist(ihev)%site
		indx = occupancy(site0(1),site0(2),site0(3))%indx
!		write(*,*) 'site,indx: ',site,indx
	else
		R = par_uni(kpar)
		x = 1 + R*NX
		R = par_uni(kpar)
		y = 1 + R*NY
		R = par_uni(kpar)
		z = 1 + R*NZ
		site0 = (/x,y,z/)
		indx = occupancy(x,y,z)%indx
		if (indx(1) < 0) cycle      ! OUTSIDE_TAG or DC
		if (RELAX_INLET_EXIT_PROXIMITY) then
			if (.not.inletOK(site0)) cycle
		else
			if (cdistance(site0) > INLET_R_FRACTION*radius) cycle
		endif
	endif
	hit = .false.
	do i = 1,27
		site = site0 + jumpvec(:,i)
		indx = occupancy(site(1),site(2),site(3))%indx
		if (indx(1) < 0) cycle      ! OUTSIDE_TAG or DC
		if (indx(1) == 0 .or. indx(2) == 0) then
			cognate = .false.
			if (evaluate_residence_time) then
				! This is for determining the transit time distribution
				if (TAGGED_EXIT_CHEMOTAXIS) then		! tag a fraction of incoming cells
					R = par_uni(kpar)
					if (istep*DELTA_T < 48*60 .and. R < TAGGED_CHEMO_FRACTION) then
						tag = RES_TAGGED_CELL
						ninflow_tag = ninflow_tag + 1
					else
						tag = 0
					endif
				else							! tag all cells entering over a specified interval
					if (istep > istep_res1 .and. istep <= istep_res2) then
						tag = RES_TAGGED_CELL   
						ninflow_tag = ninflow_tag + 1
					else
						tag = 0
					endif
				endif
			elseif (track_DCvisits .and. istep > istep_DCvisits) then	! This is the normal procedure for computing the distinct DC visit distribution
				if (TAGGED_DC_CHEMOTAXIS) then	! This is to quantify the effect of DC chemotaxis on DC visits
					tagcell = (par_uni(kpar) < TAGGED_CHEMO_FRACTION)
				else
					tagcell = .true.
				endif
				if (tagcell .and. ntagged < ntaglimit) then
					tag = TAGGED_CELL
					if (ntagged == ntaglimit-1) then	! All candidate cells will be tagged
						t_taglimit = tnow + t_log_DCvisits
					endif
				else
					tag = 0
				endif
			else
				tag = 0
			endif
			if (FAST) then
			    cognate = .true.
			    p_CD4 = CTYPE_FRACTION(CD4)*TC_COGNATE_FRACTION(CD4)/(CTYPE_FRACTION(CD4)*TC_COGNATE_FRACTION(CD4) + CTYPE_FRACTION(CD8)*TC_COGNATE_FRACTION(CD8))
                if (par_uni(kpar) < p_CD4) then
                    ctype = CD4
                else
                    ctype = CD8
                endif
!			    ctype = select_CD4_CD8()
            else
    			call select_cell_type(ctype,cognate,kpar)
    			if (cognate) then
    			    ncogin = ncogin + 1
    			else
    			    noncogin = noncogin + 1
    			endif
            endif
			if (cognate) then
				ncogseed(ctype) = ncogseed(ctype) + 1
			endif

			call add_Tcell(site,ctype,cognate,gen,tag,NAIVE,region,kcell,ok)
			if (.not.ok) return
			if (debug_DCchemotaxis) then
				if (tag == TAGGED_CELL .and. idbug == 0) then
					idbug = kcell
					open(nfchemo,file= 'chemo.log',status='replace')
					write(nfchemo,*) 'Logging DC chemotaxis for tagged T cell: ',idbug
				endif
			endif

			hit = .true.
			exit
		endif
	enddo
	if (hit) then
		k = k+1
		cycle
	endif
enddo
check_inflow = check_inflow + node_inflow

! Outflow
! Need to preserve (for now) the code that was used to simulate the adoptive transfer expt,
! i.e. using L_selectin.  In this case:
!    exit_prob = 1
!    use_exit_chemotaxis = false
!    computed_outflow = false
!    egress possible for central site only

if (L_selectin) then
	exit_prob = 1
	use_full_neighbourhood = .false.
	nbrs = 1
else
	exit_prob(1) = 0.02		! TESTING------------------ (0.02 works when chemo_K_exit = 0)
	exit_prob(2) = exit_prob(1)*residence_time(CD4)/residence_time(CD8)
	use_full_neighbourhood = .true.
	nbrs = 26
endif

ne = 0
if (lastexit > max_exits) then
	write(logmsg,'(a,2i4)') 'Error: portal_traffic: lastexit > max_exits: ',lastexit, max_exits
	call logger(logmsg)
	ok = .false.
	return
endif
allocate(permex(lastexit))
do k = 1,lastexit
	permex(k) = k
enddo
call permute(permex,lastexit,kpar)
do iloop = 1,nloops
do ipermex = 1,lastexit
	iexit = permex(ipermex)
	if (exitlist(iexit)%ID == 0) cycle
	if (exfract /= 1) then
		R = par_uni(kpar)
		if (R > exfract) cycle
	endif
	esite0 = exitlist(iexit)%site
	
!	do j = 0,1	! was 0,1 ERROR
!		if (j == 0) then
!	do j = 1,27
!		if (j == 14) then
		
	do j = 1,nbrs+1
		if (use_full_neighbourhood) then
			if (j == 14) then
				esite = esite0
				central = .true.
			else
				esite = esite0 + jumpvec(:,j)
				central = .false.
			endif
		else
			if (j == 1) then
				esite = esite0
				central = .true.
			else
				esite = esite0 + neumann(:,1)
!				esite = esite0 + jumpvec(:,1)
				central = .false.
			endif
		endif
	
		indx = occupancy(esite(1),esite(2),esite(3))%indx
		do slot = 2,1,-1
			kcell = indx(slot)
			if (kcell > 0) then
				isbound = .false.
				do i = 1,MAX_DC_BIND
					if (cellist(kcell)%DCbound(1) /= 0) isbound = .true.
				enddo
				if (isbound) cycle
				! Determine egress_possible from kcell, will depend on kcell - activated cognate cell is allowed
				! The idea is that this is used if use_exit_chemotaxis = .false.
!				write(*,*) kcell,cellist(kcell)%ctype
				if (par_uni(kpar) < exit_prob(cellist(kcell)%ctype)) then
					egress_possible = .true.
				else
					egress_possible = .false.
				endif				
				if (L_selectin) then	! overrides preceding code, no egress for non-cognate cells
					if (associated(cellist(kcell)%cptr)) then	! cognate cell, exit is possible
						egress_possible = .true.
						if (nloops == 2 .and. .not.central .and. par_uni(kpar) < 0.5) then
							egress_possible = .false.
						endif
					else
						egress_possible = .false.
					endif
				endif		
				if (egress_possible) then	
					call cell_exit(kcell,slot,esite,left)
					if (left) then
						ne = ne + 1
						check_egress(iexit) = check_egress(iexit) + 1
						if (evaluate_residence_time) then
							if (cellist(kcell)%tag == RES_TAGGED_CELL) then
								noutflow_tag = noutflow_tag + 1
								restime_tot = restime_tot + tnow - cellist(kcell)%entrytime
								ihr = (tnow - cellist(kcell)%entrytime)/60. + 1
								Tres_dist(ihr) = Tres_dist(ihr) + 1
		!                        write(*,*) 'entry: ',cellist(kcell)%entrytime
							endif
						endif
						if (track_DCvisits .and. cellist(kcell)%tag == TAGGED_CELL) then
							call recordDCvisits(kcell)
							if (debug_DCchemotaxis .and. kcell == idbug) then
								write(nfchemo,*) 'Tagged T cell has left'
								write(logmsg,*) 'Tagged T cell has left'
								call logger(logmsg)
							endif
						endif
						if (TAGGED_LOG_PATHS) then
							tag = cellist(kcell)%tag
							if (tag == TAGGED_CELL) then
								call end_log_path(kcell,1)
							elseif (tag == CHEMO_TAGGED_CELL) then
								call end_log_path(kcell,2)
							endif
						endif
						if (ne == node_outflow .and. computed_outflow) exit
					endif
				endif
			endif
		enddo
	enddo
	
	if (ne == node_outflow .and. computed_outflow) exit
enddo
if (ne == node_outflow .and. computed_outflow) exit
enddo
deallocate(permex)

!write(logmsg,*) 'portal_traffic: ',InflowTotal,noncogin,ncogin,node_outflow 
!call logger(logmsg)

if (FAST) then
    NTcells = NTcells + noncogin - node_outflow
    Radius = (NTcells*3/(4*PI))**0.33333
    return
endif
   
net_inflow = node_inflow - ne
nadd_sites = nadd_sites + net_inflow
total_in = total_in + node_inflow
total_out = total_out + ne
NTcells = NTcells + net_inflow
NTcellsPer = NTcellsPer + node_outflow
Radius = (NTcells*3/(4*PI))**0.33333
	!write(*,'(a,4i8)') 'portal_traffic: ',node_inflow,ne,net_inflow,NTcells 
end subroutine

!-----------------------------------------------------------------------------------------
! The cell kcell in slot at esite(:) is a candidate to exit.  If it meets the criteria,
! left = .true.
! Currently, when a T cell exits the LN its ID is set to 0, and a gap is recorded in cellist(:).
! To keep track of a T cell in the periphery, we need to keep it in cellist(:), but somehow
! flag it to indicate that it is not in the LN.  Use region.
!-----------------------------------------------------------------------------------------
subroutine cell_exit(kcell,slot,esite,left)
integer :: kcell, slot, esite(3)
logical :: left
integer :: x, y, z, ctype, gen, stage, region
real :: tnow
logical :: cognate, activated
type(cog_type), pointer :: p

left = .false.
tnow = istep*DELTA_T
if (retain_tagged_cells .and. cellist(kcell)%tag == TAGGED_CELL) then
	if (tnow - cellist(kcell)%entrytime < t_log_DCvisits) return
endif
if (cellist(kcell)%DCbound(1) /= 0) return     ! MUST NOT BE BOUND TO A DC!!!!!!!!!!!!!!!
if (evaluate_residence_time .or. track_DCvisits) then
    cognate = .false.
elseif (associated(cellist(kcell)%cptr)) then
    cognate = .true.
    p => cellist(kcell)%cptr
	call get_stage(p,stage,region)
    gen = get_generation(p)
    ctype = cellist(kcell)%ctype
    if (.not.exitOK(p,ctype)) then
        return
    endif
else
    cognate = .false.
endif
! For initial testing, remove cells that leave the LN
cellist(kcell)%DCbound = 0
x = esite(1)
y = esite(2)
z = esite(3)
if (slot == 2) then
    occupancy(x,y,z)%indx(2) = 0
else
    if (occupancy(x,y,z)%indx(2) == 0) then
        occupancy(x,y,z)%indx(1) = 0
    else    ! shift cell on indx(2) to indx(1)
        occupancy(x,y,z)%indx(1) = occupancy(x,y,z)%indx(2)
        occupancy(x,y,z)%indx(2) = 0
    endif
endif
if (cognate) then
    if (.not.evaluate_residence_time) then
		call efferent(p,ctype)
	endif
!	if (SIMULATE_PERIPHERY) then ! Useful to keep a record of egressing cognate cells, even if they are not touched again
		region = PERIPHERY
		call set_stage_region(p,stage,region)
!	else
!		ngaps = ngaps + 1
!		gaplist(ngaps) = kcell
!		cellist(kcell)%ID = 0
!	    cognate_list(p%cogID) = 0
	    ncogleft = ncogleft + 1
!	endif
else
	ngaps = ngaps + 1
	gaplist(ngaps) = kcell
	cellist(kcell)%ID = 0
endif
left = .true.
nleft = nleft + 1
end subroutine

!--------------------------------------------------------------------------------
!-------------------------------------------------------------------------------- 
subroutine initialise_vascularity
! VEGF_MODEL = 1
!	VEGF_beta = 4.0e-8
	VEGF_baserate = VEGF_beta*NTcells0
!    VEGF_decayrate = 0.002		! delta_G
!    vasc_maxrate = 0.0006		! alpha_V
    VEGF = VEGF_baserate/VEGF_decayrate    ! steady-state VEGF level M_G0
    c_vegf_0 = VEGF/NTcells0	! taking K_V = 1.0
!    vasc_beta = 1.5				! beta_V
    vasc_decayrate = vasc_maxrate*hill(c_vegf_0,vasc_beta*c_vegf_0,vasc_n)	! delta_V
 !   write(*,*) 'Vascularity parameters:'
 !   write(*,*) 'alpha_G,beta_G,delta_G: ',VEGF_alpha, VEGF_beta, VEGF_decayrate
 !   write(*,*) 'alpha_V,beta_V,delta_V: ',vasc_maxrate,vasc_beta,vasc_decayrate
 !   write(*,*) 'c_vegf_0,VEGF0,VEGF_baserate: ',c_vegf_0,VEGF,VEGF_baserate
 !   write(*,*) 'vasc_beta*c_vegf_0: ',vasc_beta*c_vegf_0
vascularity = 1.00
!write(*,*) 'VEGF_MODEL, VEGF_baserate: ',VEGF_MODEL, VEGF_baserate
end subroutine

!-----------------------------------------------------------------------------------------
! Vascularity responds to the VEGF level.  The rate of production of VEGF is proportional to
! either:
! (VEGF_MODEL_1) a specified inflammation signal
! or
! (VEGF_MODEL_2) the total DC activity level (i.e. antigen load).
!
! In Model 1
!	VEGFsignal depends on the inflammation level
! In Model 2
!   VEGFsignal depends on total DC activity = total antigen density (normalized)
!
!   VEGF_baserate = constitutive rate of VEGF production (per LN volume, i.e. per T cell)
!   dVEGFdt = current rate of VEGF secretion (artificial units)
!   VEGF = current total mass of VEGF (artificial units)
!   c_vegf = current VEGF concentration (artificial units)
!   vasc_beta, vasc_n = Hill function parameters for dependence of vascularity growth
!                    on VEGF concentration
!   vasc_maxrate = maximum rate constant for growth of relative vascularity
!   vasc_decayrate = rate constant for decline of relative vascularity
!   Vascularity = current relative vascularity
!   Note: if inflammation = 0, i.e. VEGFsignal = 0, there should be no change in vascularity,
!   i.e. we should have dVdt = 0.
!-----------------------------------------------------------------------------------------
subroutine vascular
real :: VEGFsignal=0, dVEGFdt, c_vegf, dVdt, Nfactor
real :: Cck1 = 1.5, Cck2 = 0.1

!write(*,*) 'vascular'
if (.not.vary_vascularity) then
    Vascularity = 1.0
    return
endif
Nfactor = real(NTcells)/NTcells0
!if (VEGF_MODEL == 1) then
    VEGFsignal = get_inflammation() ! Rate of secretion of VEGF is proportional to inflammation  
!elseif (VEGF_MODEL == 2) then
!    VEGFsignal = get_DCactivity()   ! Rate of secretion of VEGF is proportional to total DC antigen activity
!endif
dVEGFdt = VEGFsignal*VEGF_alpha + VEGF_baserate - VEGF_decayrate*VEGF
! Mass of VEGF is augmented by rate, and is subject to decay
VEGF = VEGF + dVEGFdt*DELTA_T
c_vegf = VEGF/NTcells   ! concentration (proportional to, anyway) 
c_vegf = c_vegf
if (VEGF_MODEL == 2) then ! not used
    dVdt = vasc_maxrate*hill(c_vegf,vasc_beta,vasc_n)*Vascularity - vasc_decayrate*(Vascularity - 1)
else	! VEGF_MODEL = 1
! WRONG
!    dVdt = vasc_maxrate*hill(c_vegf,vasc_beta*c_vegf_0,vasc_n)*Vascularity - vasc_decayrate*Vascularity
! Use this
!	vasc_decayrate = 0
    dVdt = vasc_maxrate*hill(c_vegf,vasc_beta*c_vegf_0,vasc_n)*Vascularity - vasc_decayrate*Vascularity	!  this works!
! Try
!	dVdt = vasc_maxrate*Ck*(VEGF/(VEGF_baserate/VEGF_decayrate) - Vascularity)	! no good
!	dVdt = vasc_maxrate*Cck1*((c_vegf/c_vegf_0 - Vascularity) + Cck2*(1-Nfactor))				!  this works!
endif
Vascularity = max(Vascularity + dVdt*DELTA_T, 1.0)
dVdt = dVdt
!if (mod(istep,240) == 0) then
!	write(logmsg,'(a,i6,e12.3,5f10.3,i5)') 'vasc: ',istep/240,VEGFsignal, &
!		VEGF/(VEGF_baserate/VEGF_decayrate), &	! M/M0
!		c_vegf/c_vegf_0, &									! C/C0
!		Nfactor, &											! N/N0
!		dVdt, &												! dV/dt
!		Vascularity, &							! V
!		Nexits									! Nexits
!	write(logmsg,'(i6,4f12.6)') NTcells,VEGF,c_vegf,hill(c_vegf,vasc_beta*c_vegf_0,vasc_n)
!	call logger(logmsg)
!endif
!write(*,*) 'dVEGFdt, c_vegf, dVdt: ',dVEGFdt, c_vegf, dVdt, Vascularity
end subroutine

!-----------------------------------------------------------------------------------------
! The total number of T cells added to the blob since the last balancing is nadd_sites.
! (Note that this can be negative, indicating a net loss of T cells). 
! A balancing is triggered either when this count exceeds a limit nadd_limit, which is a
! specified fraction of the original T cell population NTcells0, or when the time since
! the last balancing exceeds BALANCER_INTERVAL, and an adjustment to site count is needed.
! Sites needed or made available by a change to the DC population are accounted for.
! Either new sites are added (made available) or existing sites are removed (made OUTSIDE).
!-----------------------------------------------------------------------------------------
subroutine Balancer(ok)
logical :: ok
integer :: nadd_total, nadd_limit, n, nadded
integer :: k, idc, naddDC, naddex, nremex, Nexits0, dexit
real :: tnow
!real, save :: lasttime = 0
integer :: kpar = 0
logical :: blob_changed

ok = .true.
!dbug = .false.
if (FAST) then
    nadd_sites = NTcells + NDCalive*NDCsites - Nsites
!    write(logmsg,'(a,5i10)') 'balancer: ', &
!		NTcells,NDCalive,NTcells+NDCalive*NDCsites,Nsites,nadd_sites
!    call logger(logmsg)
endif

if (IN_VITRO) then
	if (ndeadDC > 0) then
		do k = 1,ndeadDC
			idc = DCdeadlist(k)
			call clearDC(idc)
		enddo
		ndeadDC = 0
	endif
	return
endif
blob_changed = .false.
tnow = istep*DELTA_T
nadd_limit = 0.01*NTcells0
nadd_total = nadd_sites
if ((abs(nadd_total) > nadd_limit) .or. (tnow > lastbalancetime + BALANCER_INTERVAL)) then
!    write(*,*) 'balancer: ',istep
    if (dbug) write(nflog,*) 'balancer: nadd_total: ',nadd_total,nadd_limit,lastbalancetime,BALANCER_INTERVAL
    if (dbug) write(nflog,*) 'call squeezer'
	call squeezer(.false.)
	if (dbug) write(nflog,*) 'did squeezer'
    ! adding/removal of sites on occupancy
    if (ndeadDC > 0) then
		if (dbug) write(nflog,*) 'remove DCs'
        do k = 1,ndeadDC
            idc = DCdeadlist(k)
!            write(*,*) 'balancer: remove DC: ',k,ndeadDC,idc
            call clearDC(idc)
        enddo
        nadd_total = nadd_total - ndeadDC*NDCsites
        ndeadDC = 0
        if (dbug) write(nflog,*) 'balancer: removed DCs: nadd_total: ',nadd_total
    endif
    naddDC = 0
    if (use_DC .and. use_DCflux) then    ! need to account for incoming DCs
		naddDC = DCinflux(lastbalancetime,tnow,kpar)
		if (naddDC > 0) then
			call place_DCs(naddDC,nadded)
			nadd_total = nadd_total + nadded*NDCsites
			if (dbug) write(nflog,*) 'balancer: added DCs: NDCtotal: ',NDCtotal
		endif
    endif
    if (FAST) then
        nadd_total = NTcells + (NDCalive+naddDC)*NDCsites - Nsites
    endif
    if (dbug) write(nflog,*) 'nadd_total: ',nadd_total
    if (nadd_total > 0) then
        n = nadd_total
	    if (dbug) write(nflog,*) 'call addSites: ',n
        call addSites(n,ok)
        if (.not.ok) return
	    if (dbug) write(nflog,*) 'did addSites'
        blob_changed = .true.
    elseif (nadd_total < 0) then
        n = -nadd_total
	    if (dbug) write(nflog,*) 'call removeSites: ',n
        call removeSites(n,ok)
        if (.not.ok) return
	    if (dbug) write(nflog,*) 'did removeSites'
!	    write(logmsg,*) 'did removeSites: exit #10: ',exitlist(10)%site,exitlist(5)%site
!	    call logger(logmsg) 
!	    call checkexits("after removeSites")
        blob_changed = .true.
    else
        n = 0
    endif
    if (n /= 0) then
        write(logmsg,'(a,4i6)') 'Error: balancer: remaining sites: ',n, &
			NTcells,NDCalive,Nsites
		call logger(logmsg)
        ok = .false.
    stop
        return
    endif
!    if (DC_motion) then
!        call test_moveDC
!    endif
   if (naddDC > 0) then
		if (dbug) write(nflog,*) 'call reassign_DC'
        call reassign_DC(kpar,ok)
        if (.not.ok) return
		if (dbug) write(nflog,*) 'did reassign_DC'
    endif
    call growDC
! The cognate list is maintained at the time that a cell arrives or leaves
    if (NTcells+NDCalive*NDCsites /= Nsites) then
!	    write(logmsg,'(a,4i10)') 'Error: balancer: cells /= sites: ', &
!			NTcells,NDCalive,NTcells+NDCalive*NDCsites,Nsites
!	    call logger(logmsg)
	    ok = .false.
	    return
	endif
    lastbalancetime = tnow
    nadd_sites = 0
    if (blob_changed) then
		if (dbug) write(nflog,*) 'call make_split'
        call make_split(.false.)
        if (use_diffusion) then
            call setup_minmax
        endif
        if (dbug) write(nflog,'(a,2i6)') 'balancer: nadd_total, radius: ',nadd_total,int(radius)
    endif
    if (USE_PORTAL_EGRESS) then
		call AdjustExitPortals
	endif
else
    call set_globalvar
endif
! Do not add or remove portals more than once per hour.  To prevent repetitive adding and removing.
if (use_portal_egress .and. use_traffic .and. tnow - last_portal_update_time > 60) then		
	if (NTcells < NTcells0 .and. .not.SURFACE_PORTALS) then
		! Set Nexits = Nexits0 for steady-state maintenance.  To wrap up the end of the response.
!        Nexits0 = exit_fraction*NTcells0
		Nexits0 = requiredExitPortals(NTcells0)
		if (Nexits < Nexits0) then
			last_portal_update_time = tnow
			do k = 1,Nexits0 - Nexits
			    if (dbug) then
			        write(nflog,*) 'AddExitPortal: ',k
			    endif
				call AddExitPortal()
			enddo
!			write(*,*) '-------------------------------------'
!			write(*,*) 'Set Nexits to base steady-state level'
!			write(*,*) '--------------------------'
!			write(*,*) 'Nexits: ',Nexits
!			write(*,*) '--------------------------'
		elseif (Nexits > Nexits0) then
			last_portal_update_time = tnow
		    if (dbug) then
		        write(nflog,*) 'RemoveExitPortals: ',Nexits,Nexits0,Nexits - Nexits0
		    endif
			call removeExitPortals(Nexits - Nexits0)
!			write(*,*) '-------------------------------------'
!			write(*,*) 'Set Nexits to base steady-state level'
!			write(*,*) '--------------------------'
!			write(*,*) 'Nexits: ',Nexits
!			write(*,*) '--------------------------'
		endif
		return
	endif
    dexit = requiredExitPortals(NTcells) - Nexits
    if (dexit > 0) then
        naddex = dexit
		last_portal_update_time = tnow
!       write(*,*) '--------------------------'
!       write(*,*) 'Nexits: ',Nexits
!       write(*,*) '--------------------------'
!       write(nflog,*) 'added: Nexits: ',naddex, Nexits
!        write(logmsg,*) 'add: Nexits: ',naddex, NTcells, Nexits, requiredExitPortals(NTcells),istep
!        call logger(logmsg) 
        if (dbug) then
            write(nflog,*) 'dexit > 0: ',dexit
        endif
        do k = 1,naddex
            call AddExitPortal()
        enddo
    elseif (dexit < 0) then
        nremex = -dexit
        if (dbug) then
            write(nflog,*) 'dexit < 0: ',dexit
        endif
		last_portal_update_time = tnow
        call removeExitPortals(nremex)
    endif
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Sites are added to occupancy().  The general idea is to preserve the shape of the 
! T cell zone, i.e. a spherical blob remains roughly spherical.
!-----------------------------------------------------------------------------------------
subroutine addSites(n,ok)
integer :: n
logical :: ok
integer :: maxblist,x,y,z,i,k,nb,nadd,idc,kdc,site0(3),site(3)
real :: x2, y2, z2, r2, r2min, r2max, r(3)
logical :: bdryflag
integer, allocatable :: t(:), bdrylist(:,:)
real, allocatable :: r2list(:)

ok = .true.
!write(logmsg,'(a,i5,i7,f8.2)') 'addSites: ',n,Nsites,Radius
!call logger(logmsg)
r2 = Radius*Radius
maxblist = 4*PI*r2*0.1*Radius
allocate(t(maxblist))
allocate(bdrylist(3,maxblist))
allocate(r2list(maxblist))

r2max = 1.3*r2
r2min = 0.8*r2

!write(*,*) 'addSites: r2min,r2max: ',r2min,r2max
k = 0
do z = 2,NZ-1
    z2 = (z-z0)*(z-z0)
    if (z2 > r2max) cycle
    do y = 2,NY-1
        y2 = (y-y0)*(y-y0)
        if (y2+z2 > r2max) cycle
        do x = 2,NX-1
            x2 = (x-x0)*(x-x0)
            r2 = x2 + y2 + z2
            if (r2 > r2min .and. r2 < r2max) then
                if (occupancy(x,y,z)%indx(1) == OUTSIDE_TAG) then
                    bdryflag = .false.
                    site0 = (/x,y,z/)
                    do i = 1,6
                        site = site0 + neumann(:,i)
                        if (occupancy(site(1),site(2),site(3))%indx(1) >= 0) then
                            bdryflag = .true.
                            exit
                        endif
                    enddo
                    if (bdryflag .and. k < maxblist) then
                        k = k+1
                        bdrylist(:,k) = site0
                        r2list(k) = r2
                    endif
                endif
            endif
        enddo
    enddo
enddo
nb = k
!write(*,*) 'bdry sites: ',nb,maxblist
if (nb < n) then
    write(logmsg,'(a,2i8,a)') 'Error: addSites: insufficient candidate sites: ',nb,n,' Increase lattice size NX'
    call logger(logmsg)
    ok = .false.
    return
endif
do i = 1,nb
    t(i) = i
enddo
call qsort(r2list,nb,t)     ! sort in increasing order
nadd = min(n,nb)
do i = 1,nadd
    site = bdrylist(:,t(i))
    occupancy(site(1),site(2),site(3))%indx = 0
    occupancy(site(1),site(2),site(3))%DC = 0
    if (use_cytokines .and. use_diffusion) then
        cyt(site(1),site(2),site(3),:) = 0
    endif
!    call smear_cyt(site)
    ! Now need to determine %DC(0:DCDIM-1)     ! DC(0) = number of near DCs, DC(k) = ID of kth near DC
    if (use_DC) then
        do idc = 1,NDC
            r = DClist(idc)%site - site
            if (norm(r) <= DC_RADIUS) then
	            kdc = occupancy(site(1),site(2),site(3))%DC(0)
		        if (kdc < 0) cycle   ! don't touch a DC site
	            if (kdc < DCDIM-1) then     ! can add to the list
	                kdc = kdc + 1
	                occupancy(site(1),site(2),site(3))%DC(0) = kdc
	                occupancy(site(1),site(2),site(3))%DC(kdc) = idc
		        endif
            elseif (norm(r) <= chemo_radius) then
	            kdc = occupancy(site(1),site(2),site(3))%cDC(0)
		        if (kdc < 0) cycle   ! don't touch a DC site
	            if (kdc < cDCDIM-1) then     ! can add to the list
	                kdc = kdc + 1
	                occupancy(site(1),site(2),site(3))%cDC(0) = kdc
	                occupancy(site(1),site(2),site(3))%cDC(kdc) = idc
		        endif
	        endif
        enddo
    endif
    ! Check for adjacent exit site (portal), and if necessary move it.
!    call adjustExit(site)
enddo
n = n - nadd
Nsites = Nsites + nadd
deallocate(t)
deallocate(bdrylist)
deallocate(r2list)
end subroutine

!-----------------------------------------------------------------------------------------
! Moves exits both out towards the blob boundary - when the blob is growing -
! and back towards the centre - when the blob is shrinking.
!-----------------------------------------------------------------------------------------
subroutine AdjustExitPortals
integer :: iexit, site1(3), site2(3), site(3), k, kmin, kmax
real :: u(3), jump(3), proj, pmin, pmax
integer :: nin1, nin2
logical :: ok

if (TAGGED_LOG_PATHS) return
!write(*,*) 'AdjustExitPortals'
do iexit = 1,lastexit
    if (exitlist(iexit)%ID == 0) cycle
	site1 = exitlist(iexit)%site
	u = site1 - Centre
	u = u/norm(u)	! outward-pointing unit vector
	nin1 = neighbourhoodCount(site1)
	if (nin1 < 27) then
		! need to move the portal inwards, choose the site with jump closest to -u direction
		kmin = 0
		pmin = 0
		do k = 1,27
			if (k == 14) cycle
			jump = jumpvec(:,k)
			site = site1 + jump
			if (.not.portalOK(site)) cycle
			proj = dot_product(u,jump)/norm(jump)
			if (proj < pmin) then
				pmin = proj
				kmin = k
			endif
		enddo
		if (kmin == 0) cycle
		site2 = site1 + jumpvec(:,kmin)
	else
		! may need to move portal outwards - try a new site
		pmax = 0
		do k = 1,27
			if (k == 14) cycle
			jump = jumpvec(:,k)
			site = site1 + jump
			if (.not.portalOK(site)) cycle
			proj = dot_product(u,jump)/norm(jump)
			if (proj > pmax) then
				pmax = proj
				kmax = k
			endif
		enddo
		site2 = site1 + jumpvec(:,kmax)
		nin2 = neighbourhoodCount(site2)
		if (nin2 < 27) cycle	! don't move here
	endif
	call RemoveExitPortal(site1)
	
	! Note: since we are reusing an existing exit index, no need to increment %lastexit
	Nexits = Nexits + 1
	if (Nexits > max_exits) then
		write(*,*) 'Error: AdjustExitPortals: too many exits: need to increase max_exits: ',max_exits
		stop
	endif
	call PlaceExitPortal(iexit,site2)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Try to remove n sites from the boundary of the blob (occupancy).
! A boundary site is one with at least 3 Neumann neighbours outside the blob.
! NOTE: The code for moving an exit should be changed.  Rather than just removing the exit
! and letting balancer() restore it wherever it wants, we should move the exit.  This
! becomes significant when exit chemotaxis is operating.
!-----------------------------------------------------------------------------------------
subroutine removeSites(n,ok)
integer :: n
logical :: ok
integer :: maxblist,x,y,z,i,nout,k,nb,nr,count,site0(3),site(3),indx(2),kcell
integer :: freeslot
real :: x2, y2, z2, r2, r2min, r2max
logical :: moved
integer, allocatable :: t(:), bdrylist(:,:)
real, allocatable :: r2list(:)

!write(logmag,'(a,i6)') 'removeSites: ',n
!call logger(logmsg)
ok = .true.
r2 = Radius*Radius
maxblist = 4*PI*r2*0.1*Radius
allocate(t(maxblist))
allocate(bdrylist(3,maxblist))
allocate(r2list(maxblist))
r2max = 1.3*r2
r2min = 0.8*r2

! Make a list of boundary sites.
k = 0
do z = 1,NZ
    z2 = (z-z0)*(z-z0)
    if (z2 > r2max) cycle
    do y = 1,NY
        y2 = (y-y0)*(y-y0)
        if (y2+z2 > r2max) cycle
        do x = 1,NX
            x2 = (x-x0)*(x-x0)
            r2 = x2 + y2 + z2
            if (r2 > r2min .and. r2 < r2max) then
                indx = occupancy(x,y,z)%indx
                if (indx(1) >= 0) then
                    if (indx(1) > 0 .and. indx(2) > 0) cycle    ! consider only sites with 0 or 1 cell
                    site0 = (/x,y,z/)
                    nout = 0
                    do i = 1,6
                        site = site0 + neumann(:,i)
                        if (site(1) < 1 .or. site(1) > NX) cycle
                        if (site(2) < 1 .or. site(2) > NY) cycle
                        if (site(3) < 1 .or. site(3) > NZ) cycle
                        if (occupancy(site(1),site(2),site(3))%indx(1) == OUTSIDE_TAG) then
                            nout = nout + 1
                        endif
                    enddo
                    if (nout >= 1 .and. k < maxblist) then
                        k = k+1
                        bdrylist(:,k) = site0
                        if (r2 > r2max) then
                            write(logmsg,'(a,4i6,f8.2)') 'Error: removeSites: bad r2: ',k,x,y,z,r2
                            call logger(logmsg)
                            ok = .false.
                            return
                        endif
                        r2list(k) = r2
                    endif
                endif
            endif
        enddo
    enddo
enddo
nb = k
!write(*,*) 'bdry sites: ',nb,maxblist,r2max
do i = 1,nb
    t(i) = i
enddo
call qsort(r2list,nb,t)     ! sort in increasing order 
nr = min(nb,n)
count = 0
do i = nb,1,-1
    if (count == nr) exit
    site0 = bdrylist(:,t(i))
    indx = occupancy(site0(1),site0(2),site0(3))%indx
    if (indx(1) > 0 .and. indx(2) > 0) cycle    ! consider only sites with 0 or 1 cell
    if (indx(1) < 0) then
        write(logmsg,'(a,i6)') 'Error: removeSites: site has a DC: ',indx(1)
        call logger(logmsg)
        ok = .false.
        return
    endif
    if (indx(1) == 0 .and. indx(2) ==  0) then
        occupancy(site0(1),site0(2),site0(3))%indx = OUTSIDE_TAG
        occupancy(site0(1),site0(2),site0(3))%DC = 0
        Nsites = Nsites - 1
        count = count + 1
        cycle
    endif
    moved = .false.
    do k = 1,2
        kcell = indx(k)
        if (kcell > 0) then
            ! This cell needs to be moved to a nearby site
            call get_free_slot(occupancy,NX,site0,site,freeslot)
            if (freeslot == 0) cycle
            occupancy(site(1),site(2),site(3))%indx(freeslot) = kcell
            cellist(kcell)%site = site
            occupancy(site0(1),site0(2),site0(3))%indx = OUTSIDE_TAG
            occupancy(site0(1),site0(2),site0(3))%DC = 0
            if (occupancy(site0(1),site0(2),site0(3))%exitnum < 0) then
            ! This is an exit site that must be moved
				write(logmsg,*) 'removeSites:  need to move exit: ',occupancy(site0(1),site0(2),site0(3))%exitnum,site0
				call logger(logmsg)
!                call removeExitPortal(site0)	! Let the balancer add an exit portal back again (if needed) 
                call moveExitPortalInwards(site0)
!                write(logmsg,'(a,4i6)') 'removeSites: did moveExitPortalInwards: ',site0
!                call logger(logmsg)
!                call checkExits("after moveExitPortalInwards")
            endif
            Nsites = Nsites - 1
            count = count + 1
            moved = .true.
        elseif (kcell < 0) then
            write(logmsg,'(a,5i6)') 'Error: removeSites: negative indx: ',site0,indx
            call logger(logmsg)
            ok = .false.
            return
        endif
        if (moved) exit
    enddo
enddo
n = n - count
deallocate(t)
deallocate(bdrylist)
deallocate(r2list)
end subroutine

!--------------------------------------------------------------------------------
! When a new empty site is added, the cytokines in the neighbouring sites are
! mixed in.
! NOT USED
!--------------------------------------------------------------------------------
subroutine smear_cyt(site0)
integer :: site0(3)
real :: csum(MAX_CYT)
integer :: k, n, site(3), jlist(MAXRELDIR)

n = 0
csum(1:Ncytokines) = 0
do k = 1,27
    if (k == 14) cycle
    site = site0 + jumpvec(:,k)
    if (occupancy(site(1),site(2),site(3))%indx(1) >= 0) then
        n = n+1
        jlist(n) = k
        csum(1:Ncytokines) = csum(1:Ncytokines) + cyt(site(1),site(2),site(3),1:Ncytokines)
    endif
enddo
csum = csum/(n+1)
do k = 1,n
    site = site0 + jumpvec(:,jlist(k))
    cyt(site(1),site(2),site(3),1:Ncytokines) = csum(1:Ncytokines)
enddo
cyt(site0(1),site0(2),site0(3),:) = csum(1:Ncytokines)

end subroutine

!--------------------------------------------------------------------------------
! Determines the average radius of the blob, which is needed for exitter().
! Note: this code only applies to the spherical blob case.
! Get range by averaging over a small number of grid sites near the centre.
! NOT USED
!--------------------------------------------------------------------------------
subroutine sizer(nb,r2list)
integer :: nb
real :: r2list(:)
integer :: k
real :: rmean

rmean = 0
do k = 1,nb
    rmean = rmean + sqrt(r2list(k))
enddo
rmean = rmean/nb
write(*,*) 'sizer: nb,rmean: ',nb,rmean
end subroutine

!--------------------------------------------------------------------------------
! Counts efferent cognate cells, and records their generation distribution.
! Only activated cells (stage >= CLUSTERS) are counted
!--------------------------------------------------------------------------------
subroutine efferent(p,ctype)
type (cog_type), pointer :: p
integer :: ctype, gen, region, i
real :: avid

!call get_stage(p,stage,region)
if (is_activated(p)) then
	gen = get_generation(p)
else
	gen = 0
endif
if (ctype > NCTYPES) then
    write(*,*) 'efferent: bad cell type:', ctype
    stop
endif
!localres%dN_EffCogTC(ctype)  = localres%dN_EffCogTC(ctype) + 1
!localres%dN_EffCogTCGen(gen) = localres%dN_EffCogTCGen(gen) + 1
!localres%N_EffCogTC(ctype)   = localres%N_EffCogTC(ctype) + 1
!localres%N_EffCogTCGen(gen)  = localres%N_EffCogTCGen(gen) + 1
totalres%dN_EffCogTC(ctype)  = totalres%dN_EffCogTC(ctype) + 1
totalres%dN_EffCogTCGen(gen) = totalres%dN_EffCogTCGen(gen) + 1
totalres%N_EffCogTC(ctype)   = totalres%N_EffCogTC(ctype) + 1
totalres%N_EffCogTCGen(gen)  = totalres%N_EffCogTCGen(gen) + 1
!write(*,'(a,2i3,a,2i6)') 'Cognate egress: ctype, gen :',ctype,gen,' totals: ',totalres%N_EffCogTCGen(0),totalres%N_EffCogTCGen(1)

! Debugging early egress
!write(nflog,'(a,2i6,3f8.3)') 'efferent: ',istep,gen,p%stimulation,p%CD69,p%S1PR1
if (log_results) then
    ! Record avidity statistics for exiting cells
    avid = p%avidity
    call log_count(avid_count,avid)
!    if (avid_count%logscale) then
!        avid = log10(avid)
!    endif
!    if (avid_count%nbins == 1) then
!        i = 1
!    else
!        !i = (avid-avidity_min)*1.01/avidity_step + 1
!        i = (avid-avid_count%binmin)/avid_count%binstep + 1.5
!        i = max(i,1)
!        !i = min(i,avidity_nlevels)
!        i = min(i,avid_count%nbins)
!    endif
!    avid_count%bincount(i) = avid_count%bincount(i) + 1
endif
end subroutine

!--------------------------------------------------------------------------------------
! Compute parameters for probability distributions of lifetime and time-to-divide
! In the case of NORMAL distributions, p1 = mean, p2 = std dev
! For X LOGNORMAL distributed, p1 and p2 are mean and std dev of Y, and X = exp(Y)
! Lifetime and dividetime parameters are given in hours, then converted to minutes
! when used.  Program time is in minutes.
!--------------------------------------------------------------------------------------
subroutine setup_dists
!real :: divide_median1, divide_median2
!real :: divide_shape1, divide_shape2
!real :: life_tc1 = 15, life_tc2 = 18
!real :: life_median1 = 48, life_median2 = 24
!real :: life_median1 = 196, life_median2 = 196      !<------------- Note: hard-coded values
!real :: life_shape1 = 1.5, life_shape2 = 1.4
integer :: i

!life_dist(1)%class = EXPONENTIAL_DIST
!life_dist(i)%p1 = life_tc1
!do i = 1,TC_MAX_GEN
!	life_dist(i)%class = EXPONENTIAL_DIST
!	life_dist(i)%p1 = life_tc2
!enddo

allocate(life_dist(TC_MAX_GEN))        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(divide_dist(TC_MAX_GEN))      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

life_dist(1)%class = LOGNORMAL_DIST
life_dist(1)%p1 = log(60*TC_life_median1)
life_dist(1)%p2 = log(TC_life_shape)
do i = 2,TC_MAX_GEN
	life_dist(i)%class = LOGNORMAL_DIST
	life_dist(i)%p1 = log(60*TC_life_median2)
	life_dist(i)%p2 = log(TC_life_shape)
enddo

!life_dist(1)%class = NORMAL_DIST
!life_dist(1)%p1 = log(60*life_median1)
!life_dist(1)%p2 = log(life_shape1)
!do i = 2,TC_MAX_GEN
!	life_dist(i)%class = NORMAL_DIST
!	life_dist(i)%p1 = log(60*life_median1)
!	life_dist(i)%p2 = log(life_shape1)
!enddo


!divide_dist(1)%class = LOGNORMAL_DIST
!divide_dist(1)%class = CONSTANT_DIST
!if (divide_dist(1)%class == LOGNORMAL_DIST) then
!	divide_dist(1)%p1 = log(60*divide_median1)
!	divide_dist(1)%p2 = log(divide_shape1)
!	do i = 2,TC_MAX_GEN
!		divide_dist(i)%class = LOGNORMAL_DIST
!		divide_dist(i)%p1 = log(60*divide_median2)
!		divide_dist(i)%p2 = log(divide_shape2)
!	enddo
!elseif (divide_dist(1)%class == CONSTANT_DIST) then
!	do i = 1,TC_MAX_GEN
!		divide_dist(i)%class = CONSTANT_DIST
!		divide_dist(i)%p1 = 60*divide_median2
!		divide_dist(i)%p2 = 0
!	enddo
!endif


divide_dist(1)%class = divide_dist1%class
divide_dist(1)%p1 = divide_dist1%p1
divide_dist(1)%p2 = divide_dist1%p2
divide_dist(2)%class = divide_dist2%class
divide_dist(2)%p1 = divide_dist2%p1
divide_dist(2)%p2 = divide_dist2%p2
do i = 3,TC_MAX_GEN
	divide_dist(i)%class = divide_dist(2)%class
	divide_dist(i)%p1 = divide_dist(2)%p1
	divide_dist(i)%p2 = divide_dist(2)%p2
enddo

!do i = 1,TC_MAX_GEN
!	if (divide_dist(i)%class == NORMAL_DIST) then
!		divide_dist(i)%p1 = 60*divide_dist(i)%p1
!		divide_dist(i)%p2 = 60*divide_dist(i)%p2
!	elseif (divide_dist(i)%class == LOGNORMAL_DIST) then
!		divide_dist(i)%p1 = log(60*divide_dist(i)%p1)
!		divide_dist(i)%p2 = log(divide_dist(i)%p2)
!	elseif (divide_dist(i)%class == CONSTANT_DIST) then
!		divide_dist(i)%p1 = 60*divide_dist(i)%p1
!		divide_dist(i)%p2 = 0
!	endif
!enddo

call make_bind_probs

end subroutine

!-----------------------------------------------------------------------------------------
! To compare with the Matlab code: Lauff_a.m
!-----------------------------------------------------------------------------------------
subroutine testupdater
real :: state(IL2_NP), statep(IL2_NP)
integer :: ctype = CD4, ID = 1
real :: t, dt, th, S_TCR, C_IL2, mrate, crate, dm, dc
logical :: producing = .true., dbgflag = .false.
logical :: first, ok
integer :: i, nsub = 1, nhours = 24
real ::  C0 = 100

call read_cell_params(ok)
first = .true.
!state = 0  ! for Lauff_a.m comparison
state = (/1200,1300,300,1400,4000,90000/)
t = 0
crate = 0.01
S_TCR = 5000
C_IL2 = C0
dt = DELTA_T/nsub
!write(*,*)
!write(*,'(6x,9a8)') 'Min ','CD25_S','Comp_S','CD25_I','Comp_I','IL2_I','S_Rec','IL-2','rate'
do i = 1,nhours*60*4*nsub
    t = t + dt
    th = t/60

    C_IL2 = C0

!    C_IL2 = 20.0*(1 - 0.5*th/nhours)
!    if (th > 24) producing = .false.
!    call IL2_update(ctype,t,S_TCR,state,producing,C_IL2,dt,mrate,dbgflag) ! old IL2_update, no xcase
    call IL2_update(ID,ctype,t,S_TCR,state,statep,first,producing,C_IL2,Vc,dt,mrate,dbgflag)
    dm = (mrate/Vc)*DELTA_T*M_pM*L_um3/Navo
    dc = crate*DELTA_T*M_pM*L_um3/Navo
!    if (mod(i,60*4*nsub) == 0) then
!        write(*,'(i6,7f8.0,2f8.1)') i,t,state,C_IL2,mrate
!    endif
    first = .false.
enddo

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine vascular_test
integer :: nsteps = 10.0*24*60/DELTA_T
real :: inflow0, act, tnow, exfract, nsum, Fin0, Fin, Fout, Tres

Tres = 24
NTcells0 = 100000
use_exit_chemotaxis = .false.
call initialise_vascularity

NDC = DC_FACTOR*NTcells0/TC_TO_DC
NTcells = NTcells0
Fin0 = NTcells/(Tres*60)
write(*,*) TC_TO_DC,NDC,nsteps,NTcells

nsum = 0
do istep = 1,nsteps
	tnow = istep*DELTA_T
    call vascular
!    call generate_traffic(inflow0)
	Fin = Fin0*vascularity*DELTA_T
	Fout = NTcells*DELTA_T/(Tres*60)
    if (transient_egress_suppression) then
		exfract = egressFraction(tnow)
		Fout = exfract*Fout
	endif
    if (mod(istep,60) == 0) then
        write(nfout,'(i6,f8.3,5e14.5,i8)') istep,tnow, &
			VEGF, &
            dVdt, &
            vascularity, &
            Fin, Fout, &
            NTcells
    endif
    nsum = nsum + (Fin - Fout)
    NTcells = NTcells0 + nsum
enddo

end subroutine

!-----------------------------------------------------------------------------------------
! The number of cells that leave exactly matches the number that enter.
! ZIN_FRACTION = fraction of radius over which inflow traffic occurs
!-----------------------------------------------------------------------------------------
subroutine traffic(ok)
logical :: ok
integer :: x, y, z, k, kcell, indx(2), ctype, gen, stage, region, site(3), n, slot
integer :: zin_min, zin_max, node_inflow, node_outflow, add(3), net_inflow, ihr, tag
real(DP) :: R, df, prob
real :: tnow, exfract
logical :: left, cognate
integer :: kpar=0

if (.not.use_blob) then
    write(*,*) 'ERROR: traffic: only for blobs'
    stop
endif
ok = .true.
!write(*,*) 'traffic'
tnow = istep*DELTA_T
region = LYMPHNODE
node_inflow = InflowTotal
node_outflow = OutflowTotal
df = InflowTotal - node_inflow
R = par_uni(kpar)
if (dbug) write(nfres,'(a,i2,2f10.6)') 'traffic: in: ',k,df,R
if (R < df) then
	node_inflow = node_inflow + 1
endif
exfract = 1
if (transient_egress_suppression) then
	exfract = egressFraction(tnow)
endif
node_outflow = exfract*node_outflow

if (steadystate) then
	node_outflow = node_inflow
else
	df = OutflowTotal - node_outflow
	R = par_uni(kpar)
	if (dbug) write(nfres,'(a,i2,2f10.6)') 'traffic: out: ',k,df,R
	if (R < df) then
		node_outflow = node_outflow + 1
	endif
endif
net_inflow = node_inflow - node_outflow
if (dbug) write(nfres,*) 'traffic: in,out: ',node_inflow,node_outflow
gen = 1
add = 0
zin_max = min(z0 + Radius,real(NZ))
zin_min = max(zin_max - ZIN_FRACTION*Radius,1.0)

! Inflow
k = 0
do while (k < node_inflow)
	R = par_uni(kpar)
	x = 1 + R*NX
	if (dbug) write(nfres,'(a,i4,f15.9)') 'in x R: ',x,R
	R = par_uni(kpar)
	y = 1 + R*NY
	if (dbug) write(nfres,'(a,i4,f15.9)') 'in y R: ',y,R
	if (exit_region == EXIT_EVERYWHERE) then    ! currently => entry everywhere
		R = par_uni(kpar)
		z = 1 + R*NZ
	if (dbug) write(nfres,'(a,i4,f15.9)') 'in z R: ',z,R
	else
		z = random_int(zin_min,zin_max,kpar)
	endif
	indx = occupancy(x,y,z)%indx
	if (dbug) write(nfres,*) 'site,indx: ',x,y,z,indx
	if (indx(1) < 0) cycle      ! OUTSIDE_TAG or DC
	if (indx(1) == 0 .or. indx(2) == 0) then
		site = (/x,y,z/)
		cognate = .false.
		if (evaluate_residence_time) then
			! This is for measuring residence time
			if (istep > istep_res1 .and. istep <= istep_res2) then
				tag = RES_TAGGED_CELL   
				ninflow_tag = ninflow_tag + 1
			else
				tag = 0
			endif
		elseif (track_DCvisits .and. istep > istep_DCvisits) then
			if ((par_uni(kpar) < TAGGED_CHEMO_FRACTION) .and. (ntagged < ntaglimit)) then
				tag = TAGGED_CELL
			else
				tag = 0
			endif
		else
			tag = 0
		endif
		call select_cell_type(ctype,cognate,kpar)
		if (cognate) then
			ncogseed(ctype) = ncogseed(ctype) + 1
		endif
		call add_Tcell(site,ctype,cognate,gen,tag,NAIVE,region,kcell,ok)
		if (dbug) then
			write(nfres,'(a,5i4,i6)') 'added cell: ',k,site,ctype,kcell
		 endif
		if (.not.ok) return
		k = k+1
		cycle
	endif
enddo

! Outflow
k = 0
do while (k < node_outflow)
	R = par_uni(kpar)
	if (dbug) write(nfres,'(a,f10.6)') 'out x R: ',R
	x = 1 + R*NX
	R = par_uni(kpar)
	if (dbug) write(nfres,'(a,f10.6)') 'out y R: ',R
	y = 1 + R*NY
	R = par_uni(kpar)
	if (dbug) write(nfres,'(a,f10.6)') 'out z R: ',R
	z = 1 + R*NZ        ! any z is OK to exit
	if (exit_region == EXIT_EVERYWHERE) then
		! accept it
	elseif (exit_region == EXIT_LOWERHALF) then
		if (z > z0) cycle
	endif
	indx = occupancy(x,y,z)%indx
	if (indx(1) < 0) cycle      ! OUTSIDE_TAG or DC 
	if (indx(2) > 0) then
		slot = 2
	elseif (indx(1) > 0) then
		slot = 1
	else
		cycle
	endif
	kcell = indx(slot)
!	if (SIMULATE_PERIPHERY) then
!		prob = 1/PERI_PROBFACTOR	! prob of allowing this cell to exit 
!		if (associated(cellist(kcell)%cptr)) then
!			gen = get_generation(cellist(kcell)%cptr)
!			call get_stage(cellist(kcell)%cptr,stage,region)
!			if (gen >= PERI_GENERATION .and. NDCcapable == 0) then
!				prob = 1
!			elseif (gen == 1) then	! suppress egress for undivided cognate cells.  
!				prob = 0 
!			endif
!		endif
!		R = par_uni(kpar)
!		if (R > prob) cycle
!	endif
	site = (/x,y,z/)
	call cell_exit(kcell,slot,site,left)
	if (.not.left) cycle
	k = k+1

	if (evaluate_residence_time) then
		if (cellist(kcell)%tag == RES_TAGGED_CELL) then
			noutflow_tag = noutflow_tag + 1
			restime_tot = restime_tot + tnow - cellist(kcell)%entrytime
			ihr = (tnow - cellist(kcell)%entrytime)/60. + 1
			Tres_dist(ihr) = Tres_dist(ihr) + 1
		endif
	endif
	if (track_DCvisits .and. cellist(kcell)%tag == TAGGED_CELL) then
		call recordDCvisits(kcell)
	endif
enddo
nadd_sites = nadd_sites + net_inflow
NTcells = NTcells + net_inflow
NTcellsPer = NTcellsPer + node_outflow
end subroutine

!-----------------------------------------------------------------------------------------
! To test add_Tcell()
!-----------------------------------------------------------------------------------------
subroutine add_random_cells(n, ctype, gen, stage, region)
integer :: n, ctype, gen, tag, stage, region
integer :: k, x, y, z, site(3), slots, kcell
integer :: kpar=0
logical :: cognate = .false.
logical :: ok

tag = 0
k = 0
do while (k < n)
    x = random_int(1,NX,kpar)
    y = random_int(1,NY,kpar)
    z = random_int(1,NZ,kpar)
    site = (/x,y,z/)
    slots = getslots(site)
    if (occupancy(x,y,z)%indx(1) >= 0 .and. slots < BOTH) then
        if (dbug) write(*,'(a,7i6)') 'add_random_cells: ',site,occupancy(x,y,z)%indx,slots
        call add_Tcell(site,ctype,cognate,gen,tag,stage,region,kcell,ok)
        if (dbug) write(*,'(a,7i6)') 'after add_random_cells: ',site,occupancy(x,y,z)%indx,slots
        call checkslots('add_random_cells: ',site)
        k = k+1
    endif
enddo

end subroutine

!-----------------------------------------------------------------------------------------
! Scans the cell list and builds the counts nz_sites() and nz_excess().
! nz_sites(k)  = number of available sites at the slice z = k
! nz_excess(k) = total excess of cells over sites in the paracortex zone with z >= k
! These values are used in jumper() to adjust the jump probabilities.  The probability
! of jumps in the -z direction is increased by an increment that is proportional to
! nz_excess(k)/nz_sites(k)
! Need to account for DC sites
!-----------------------------------------------------------------------------------------
subroutine scanner
integer :: x, y, z, ns, nc, nst, nct,excess, indx(2), nextra, i, k, idn, imin=0, nz1, nz2
!real, allocatable :: eratio(:), nz_sites0(:)    !, nz_totsites(:)
real :: eratio(NZ), nz_sites0(NZ)
real :: df, dfmin

if (exit_region == EXIT_EVERYWHERE) return

!allocate(eratio(NZ))
!allocate(nz_sites0(NZ))
!allocate(nz_totsites(NZ))
nz_excess = 0
excess = 0
nct = 0
nst = 0
do z = NZ,1,-1
    ns = 0
    nc = 0
    do y = 1,NY
        do x = 1,NX
            indx = occupancy(x,y,z)%indx
            if (indx(1) < 0) cycle       ! OUTSIDE_TAG or DC site (-idc)
            ns = ns+1
            if (indx(1) > 0) nc = nc+1
            if (indx(2) > 0) nc = nc+1
        enddo
    enddo
    nst = nst + ns
    nct = nct + nc
    nz_sites(z) = ns
    nz_cells(z) = nc
    excess = excess + nc - ns
    nz_excess(z) = excess
enddo

!if (mod(istep,100) == 0) then
!    write(*,*) 'nst,nct: ', nst,nct
!    write(*,*) 'nz_cells'
!    write(*,'(10i6)') nz_cells
!    write(*,*) 'nz_sites'
!    write(*,'(10i6)') nz_sites
!endif

! Note: This check fails when DC start to die, because occupancy(:)%indx is not
! updated immediately.  Therefore, we remove the check.
!if (Mnodes == 1) then
!    if (nst - nct /= NDCalive*NDCsites - nadd_sites) then
!        write(*,*) 'Site-cell imbalance: ',nst,nct,NDCalive*NDCsites-nadd_sites
!        stop
!    endif
!endif
nextra = nst - nct
! This is the imbalance between number of sites and number of T cells
! resulting from (a) DC sites and (b) sites to be added (nadd_sites)
! We need to adjust either nz_sites(:) or nz_cells(:) to bring them into balance,
! so that eratio(:) can be computed correctly to generate a drift.
! The complication arises because we want to spread the adjustment over the slices
! in a way that is proportionate to the number of sites in the slice.
! Choose to adjust nz_sites(:) (arbitrarily).

if (nextra /= 0) then
    if (nextra > 0) then
        idn = -1
    else
        idn = 1
    endif
    nz_sites0 = nz_sites    ! This conveys approx the shape of the blob - we want to maintain this
    do k = 1,abs(nextra)    ! we need to remove/add this many sites from/to nz_sites(:)
        dfmin = 1.0e10
        do i = 1,NZ
            if (nz_sites(i) == 0) cycle
            df = abs(real(nz_sites(i) + idn)/nz_sites0(i) - 1)
            if (df < dfmin) then
                dfmin = df
                imin = i
            endif
        enddo
        nz_sites(imin) = nz_sites(imin) + idn
        nst = nst + idn
    enddo
    nz_excess = 0
    excess = 0
    do z = NZ,1,-1
        excess = excess + nz_cells(z) - nz_sites(z)
        nz_excess(z) = excess
    enddo
endif
nz1 = NZ
nz2 = 1
do z = NZ,1,-1
    if (z == NZ) then
        nz_totsites(z) = nz_sites(z)
    else
        nz_totsites(z) = nz_sites(z) + nz_totsites(z+1)
    endif
    if (nz_sites(z) > 0) then
        eratio(z) = 100*nz_excess(z)/real(nz_totsites(z))
        nz1 = min(z,nz1)
        nz2 = max(z,nz2)
    else
        eratio(z) = 0
    endif
enddo
if (constant_efactor) then
    excess_factor = efactor
else
    excess_factor = efactor + (1-efactor)*exp(-etheta*(NTcells - 50000))
endif
end subroutine

!-----------------------------------------------------------------------------------------
! EGRESS_SUPPRESSION_TIME1 is the start of the ramp down to 0
! EGRESS_SUPPRESSION_TIME2 is the start of the ramp up to 1
!-----------------------------------------------------------------------------------------
real function egressFraction(tnow)
real :: tnow

if (tnow/60 < EGRESS_SUPPRESSION_TIME1) then
	egressFraction = 1
elseif (tnow/60 < EGRESS_SUPPRESSION_TIME1 + EGRESS_SUPPRESSION_RAMP) then
	egressFraction = 1 - (tnow/60 - EGRESS_SUPPRESSION_TIME1)/EGRESS_SUPPRESSION_RAMP
elseif (tnow/60 < EGRESS_SUPPRESSION_TIME2) then
	egressFraction = 0
elseif (tnow/60 < EGRESS_SUPPRESSION_TIME2 + EGRESS_SUPPRESSION_RAMP) then
	egressFraction = (tnow/60 - EGRESS_SUPPRESSION_TIME2)/EGRESS_SUPPRESSION_RAMP
else
	egressFraction = 1
endif
end function

!-----------------------------------------------------------------------------------------
! Record distributions of first DC contact for non-cognate (non-chemotactic) and
! cognate (chemotactic) cells.
!-----------------------------------------------------------------------------------------
subroutine logfirstDCcontact(cell,idc)
type(cell_type), pointer :: cell
integer :: idc
real :: tfirst, tnow
integer :: iTCchemo, iDCchemo, itfirst

tnow = istep*DELTA_T
tfirst = tnow - cell%entrytime
if (USE_DC_COGNATE) then
	if (DClist(idc)%cognate) then
		iDCchemo = 1
	else
		iDCchemo = 2
	endif
else
	iDCchemo = 1
endif
!if (DC_CHEMO_NOTRAFFIC) then	! also for use_traffic case
	if (cell%tag == TAGGED_CELL) then
!		if (cell%DCchemo < (LO_CHEMO + HI_CHEMO)/2) then
		if (cell%receptor_level(CCR1) > (LO_CHEMO + HI_CHEMO)/2) then
			iTCchemo = 1		! HI
		else
			iTCchemo = 2		! LO
		endif
	else
		return
	endif
!else
!	if (cell%tag == TAGGED_CELL) then
!		iTCchemo = 1
!	else
!		iTCchemo = 2
!	endif
!endif
firstDC_n(iTCchemo,iDCchemo) = firstDC_n(iTCchemo,iDCchemo) + 1
firstDC_tot(iTCchemo,iDCchemo) = firstDC_tot(iTCchemo,iDCchemo) + tfirst
itfirst = max(1.0,tfirst + 0.5)
firstDC_dist(iTCchemo,iDCchemo,itfirst) = firstDC_dist(iTCchemo,iDCchemo,itfirst) + 1
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine compute_stim_dist
integer :: kcell, ctype, stype, ncog, noncog
real :: s
type (cog_type), pointer :: p

ncog = 0
noncog = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    p => cellist(kcell)%cptr
    ctype = cellist(kcell)%ctype
!    stype = struct_type(ctype)
	if (associated(cellist(kcell)%cptr)) then
		stype = COG_TYPE_TAG
	else
		stype = NONCOG_TYPE_TAG
	endif
    if (stype == COG_TYPE_TAG) then
        ncog = ncog + 1
        s = p%stimulation
    elseif (stype == NONCOG_TYPE_TAG) then
        noncog = noncog + 1
    else
        write(*,*) 'ERROR: compute_stim_dist: bad stype: ',ctype,stype
        stop
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_cognate_list
integer :: k, kcell

do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell > nlist) then
        write(*,*) 'check_cognate_list: bad cognate_list entry > nlist: ',lastcogID,k,kcell,nlist
        stop
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! The aim is to display the distribution of cognate cells in the z direction, in order
! to investigate failure of cognate cell proliferation to scale with NTcells0.
! For now, just look at cognate fractions in upper and lower hemispheres.
!-----------------------------------------------------------------------------------------
subroutine get_cognate_dist(ncog1,ncog2)
integer :: z,kcell,ntot1,ntot2,ncog1,ncog2
logical :: cognate

ntot1 = 0
ntot2 = 0
ncog1 = 0
ncog2 = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    z = cellist(kcell)%site(3)
    cognate = (associated(cellist(kcell)%cptr))
    if (z < z0) then
        ntot1 = ntot1 + 1
        if (cognate) ncog1 = ncog1 + 1
    else
        ntot2 = ntot2 + 1
        if (cognate) ncog2 = ncog2 + 1
    endif
enddo
write(*,*) 'Lower: ',ntot1,ncog1,real(ncog1)/ntot1
write(*,*) 'Upper: ',ntot2,ncog2,real(ncog2)/ntot2
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine old_process_command_line
integer :: i, cnt, len, status
character :: c*(64), b*(256)
character*(64) :: progname

call get_command (b, len, status)
if (status .ne. 0) then
    write (*,*) 'get_command failed with status = ', status
    stop
end if
!write (*,*) 'command line = ', b (1:len)
call get_command_argument (0, c, len, status)
if (status .ne. 0) then
    write (*,*) 'Getting command name failed with status = ', status
    stop
end if
!write (*,*) 'command name = ', c (1:len)
progname = c(1:len)
cnt = command_argument_count ()
!write (*,*) 'number of command arguments = ', cnt
if (cnt < 1) then
    write(*,*) 'Use: ',trim(progname),' num_cpu'
    stop
!    Mnodes = 4      ! for profiling
!    write(*,*) 'Ruuning with Mnodes = ',Mnodes
endif
do i = 1, cnt
    call get_command_argument (i, c, len, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        stop
    end if
!    write (*,*) 'command arg ', i, ' = ', c (1:len)
    if (i == 1) then
!!!        read(c(1:len),'(i)') Mnodes
        read(c(1:len),*) Mnodes
        write(*,*) 'Requested threads: ',Mnodes
    elseif (i == 2) then
        inputfile = c(1:len)
        write(*,*) 'Input file: ',inputfile
    elseif (i == 3) then
        outputfile = c(1:len)
        write(*,*) 'Output file: ',outputfile
!    elseif (i == 4) then
!        resultfile = c(1:len)
!        write(*,*) 'Result file: ',resultfile 
    endif
end do
write (*,*) 'command line processed'
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_chemoactivity(cave)
real :: cave
type (cell_type), pointer :: cell
integer :: kcell

cave = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    cell => cellist(kcell)
    cave = cave + chemo_active_exit(cell)
enddo
cave = cave/NTcells
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine SaveGenDist
integer :: gendist(TC_MAX_GEN)
integer :: k, kcell, stage, region, gen, maxg
real :: fac

gendist = 0
maxg = 0
do k = 1,lastcogID
	kcell = cognate_list(k)
	if (kcell > 0) then
		if (.not.associated(cellist(kcell)%cptr)) then
			write(*,*) 'Error: SaveGenDist: cptr not associated'
			stop
		endif
		call get_stage(cellist(kcell)%cptr,stage,region)
!		if (region /= LYMPHNODE) cycle 
		gen = get_generation(cellist(kcell)%cptr)
		maxg = max(maxg,gen)
		gendist(gen) = gendist(gen) + 1
	endif
enddo
fac = 1.0/sum(gendist)
write(nflog,*) 'gendist:'
write(nflog,'(20f7.4)') fac*gendist(1:maxg)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_profile_CD69(x,y,n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_profile_cd69
use, intrinsic :: iso_c_binding
real(c_double) :: x(*)
real(c_double) :: y(*)
integer(c_int) :: n, nc
type (cog_type), pointer :: p
integer :: i, k, kcell
integer,allocatable :: cnt(:)
real :: dx

!call logger('get_profile_cd69')
n = nprofilebins
dx = 1.0/n
allocate(cnt(n))
cnt = 0
do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell == 0) cycle
    p => cellist(kcell)%cptr
    i = min(int(p%CD69/dx + 1),n)
    i = max(i,1)
    cnt(i) = cnt(i) + 1
enddo
nc = max(1,sum(cnt))
do i = 1,n
    x(i) = (i - 0.5)*dx
    y(i) = cnt(i)/real(nc)
enddo
!write(logmsg,*) 'Cognate cells: ',nc, lastcogID
!call logger(logmsg)
!write(logmsg,'(10f6.3)') x(1:n)
!call logger(logmsg)
!write(logmsg,'(10f6.3)') y(1:n)
!call logger(logmsg)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_profile_S1PR1(x,y,n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_profile_s1pr1
use, intrinsic :: iso_c_binding
real(c_double) :: x(*)
real(c_double) :: y(*)
integer(c_int) :: n, nc
type (cog_type), pointer :: p
integer :: i, k, kcell
integer,allocatable :: cnt(:)
real :: dx

n = nprofilebins
dx = 1.0/n
allocate(cnt(n))
cnt = 0
do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell == 0) cycle
    p => cellist(kcell)%cptr
    i = min(int(p%S1PR1/dx + 1),n)
    i = max(i,1)
    cnt(i) = cnt(i) + 1
enddo
nc = max(1,sum(cnt))
do i = 1,n
    x(i) = (i - 0.5)*dx
    y(i) = cnt(i)/real(nc)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_profile_CFSE(x,y,n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_profile_cfse
use, intrinsic :: iso_c_binding
real(c_double) :: x(*)
real(c_double) :: y(*)
integer(c_int) :: n, nc
type (cog_type), pointer :: p
integer :: i, k, kcell
integer,allocatable :: cnt(:)
real :: dx, cfse, logcfse
integer, parameter :: maxdiv = 20

n = nprofilebins
dx = (maxdiv+1.)/n
allocate(cnt(n))
cnt = 0
do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell == 0) cycle
    p => cellist(kcell)%cptr
    cfse = p%CFSE
    logcfse = log(cfse)/log(2.)
    i = min(int((logcfse + maxdiv + 0.5)/dx + 1),n)
!    write(nflog,'(i4,4e12.4,i3)') k,cfse,logcfse,dx,(logcfse + maxdiv + 0.5),i
    i = max(i,1)
    cnt(i) = cnt(i) + 1
enddo
!write(nflog,*) 'get_profile_CFSE:'
!write(nflog,'(20i6)') cnt
nc = max(1,sum(cnt))
do i = 1,n
    x(i) = i*dx - maxdiv - 0.5
    y(i) = cnt(i)/real(nc)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_profile_stim(x,y,n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_profile_stim
use, intrinsic :: iso_c_binding
real(c_double) :: x(*)
real(c_double) :: y(*)
integer(c_int) :: n, nc
type (cog_type), pointer :: p
integer :: i, k, kcell, region
integer,allocatable :: cnt(:)
real :: dx, stim

n = nprofilebins
dx = 1.0/n
allocate(cnt(n))
cnt = 0
do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell == 0) cycle
    p => cellist(kcell)%cptr
!    if (.not.is_activated(p)) cycle
	call get_region(p,region)
	if (region /= LYMPHNODE) cycle
    stim = p%stimulation/STIMULATION_LIMIT
    i = min(int(stim/dx + 1),n)
    i = max(i,1)
    cnt(i) = cnt(i) + 1
enddo
nc = max(1,sum(cnt))
do i = 1,n
    x(i) = (i - 0.5)*dx
    y(i) = cnt(i)/real(nc)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_profile_stimrate(x,y,n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_profile_stimrate
use, intrinsic :: iso_c_binding
real(c_double) :: x(*)
real(c_double) :: y(*)
integer(c_int) :: n, nc
type (cog_type), pointer :: p
integer :: i, k, kcell
integer,allocatable :: cnt(:)
real :: dx

n = nprofilebins
dx = 1.0/n
allocate(cnt(n))
cnt = 0
do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell == 0) cycle
    p => cellist(kcell)%cptr
    i = min(int(p%stimrate/dx + 1),n)
    i = max(i,1)
    cnt(i) = cnt(i) + 1
enddo
nc = max(1,sum(cnt))
do i = 1,n
    x(i) = (i - 0.5)*dx
    y(i) = cnt(i)/real(nc)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_profile_avidity_ln(x,y,n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_profile_avidity_ln
use, intrinsic :: iso_c_binding
real(c_double) :: x(*)
real(c_double) :: y(*)
integer(c_int) :: n, nc
type (cog_type), pointer :: p
integer :: i, k, kcell, region
integer,allocatable :: cnt(:)
real :: dx, avid

n = nprofilebins
dx = 1.0/n
allocate(cnt(n))
cnt = 0
do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell == 0) cycle
    p => cellist(kcell)%cptr
	call get_region(p,region)
	if (region /= LYMPHNODE) cycle
    avid = p%avidity/MAXIMUM_AVIDITY
    i = min(int(avid/dx + 1),n)
    i = max(i,1)
    cnt(i) = cnt(i) + 1
enddo
nc = max(1,sum(cnt))
do i = 1,n
    x(i) = (i - 0.5)*dx
    y(i) = cnt(i)/real(nc)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_profile_avidity_per(x,y,n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_profile_avidity_per
use, intrinsic :: iso_c_binding
real(c_double) :: x(*)
real(c_double) :: y(*)
integer(c_int) :: n, nc
type (cog_type), pointer :: p
integer :: i, k, kcell, region
integer,allocatable :: cnt(:)
real :: dx, avid

n = nprofilebins
dx = 1.0/n
allocate(cnt(n))
cnt = 0
do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell == 0) cycle
    p => cellist(kcell)%cptr
    if (.not.is_activated(p)) cycle
	call get_region(p,region)
	if (region /= PERIPHERY) cycle
    avid = p%avidity/MAXIMUM_AVIDITY
    i = min(int(avid/dx + 1),n)
    i = max(i,1)
    cnt(i) = cnt(i) + 1
enddo
nc = max(1,sum(cnt))
do i = 1,n
    x(i) = (i - 0.5)*dx
    y(i) = cnt(i)/real(nc)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_profile_generation_ln(x,y,n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_profile_generation_ln
use, intrinsic :: iso_c_binding
real(c_double) :: x(*)
real(c_double) :: y(*)
integer(c_int) :: n, nc
type (cog_type), pointer :: p
integer :: i, k, kcell, region, gen
integer,allocatable :: cnt(:)
real :: dx, avid

n = 20
dx = 1.0
allocate(cnt(n))
cnt = 0
do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell == 0) cycle
    p => cellist(kcell)%cptr
	call get_region(p,region)
	if (region /= LYMPHNODE) cycle
	gen = get_generation(p)
    i = min(gen,n)
    cnt(i) = cnt(i) + 1
enddo
nc = max(1,sum(cnt))
do i = 1,n
    x(i) = i*dx
    y(i) = cnt(i)/real(nc)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_profile_firstDCcontacttime(x,y,n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_profile_firstdccontacttime 
use, intrinsic :: iso_c_binding
real(c_double) :: x(*)
real(c_double) :: y(*)
integer(c_int) :: n, nc
type (cog_type), pointer :: p
integer :: i, k, kcell, gen
integer,allocatable :: cnt(:)
real :: dx, t, tnow

tnow = istep*DELTA_T
n = nprofilebins
dx = 300.0/n
allocate(cnt(n))
cnt = 0
do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell == 0) cycle
    p => cellist(kcell)%cptr
	gen = get_generation(p)
!    if (p%firstDCtime == 0) cycle
	if (gen == 1 .and. p%firstDCtime > 0 .and. (tnow-p%firstDCtime) <= 60) then
		t = (p%firstDCtime - cellist(kcell)%entrytime)
		i = min(int(t/dx + 1),n)
		i = max(i,1)
		cnt(i) = cnt(i) + 1
	endif
enddo
nc = max(1,sum(cnt))
!write(nflog,*) 'get_profile_firstDCcontacttime: ',nc
!write(nflog,'(20i4)') cnt
do i = 1,n
    x(i) = (i - 0.5)*dx
    y(i) = cnt(i)/real(nc)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_profile_DCbindtime(x,y,n) BIND(C) 
!DEC$ ATTRIBUTES DLLEXPORT :: get_profile_dcbindtime 
use, intrinsic :: iso_c_binding
real(c_double) :: x(*)
real(c_double) :: y(*)
integer(c_int) :: n
type (cog_type), pointer :: p
type (counter_type), pointer :: cp
integer :: i, nc

cp => DCbindtime_count
n = cp%nbins
do i = 1,n
	x(i) = cp%binmin + (i-0.5)*cp%binstep
enddo
y(1:n) = 0
nc = sum(cp%bincount(1:n))
if (nc > 0) then
	y(1:n) = cp%bincount(1:n)/real(nc)
else
	y(1:n) = 0
endif
cp%bincount(1:n) = 0
end subroutine

!-------------------------------------------------------------------------------- 
! Compute the clone size distribution
! For each seeding cognate cell that survives to proliferate, determine the number
! of progeny, using the %ID to identify the clones.
! Note that all cells will carry the original entrytime of the seed cell.
!
!-------------------------------------------------------------------------------- 
subroutine compute_clonesize
integer :: k, kcell, nseed, nclones, iclone, jclone, ID, gen, totalcount, i
type (cog_type), pointer :: p
type(clone_type), allocatable :: clonelist(:)
real, allocatable :: t(:)
integer, allocatable :: indx(:)

nseed = ncogseed(1) + ncogseed(2)
allocate(clonelist(nseed))

totalcount = 0
nclones = 0
do k = 1,lastcogID
	kcell = cognate_list(k)
	if (kcell > 0) then
        p => cellist(kcell)%cptr
!    	call get_region(p,region)
!    	if (region /= LYMPHNODE) cycle
	    gen = get_generation(p)
	    if (gen < 2) cycle
	    totalcount = totalcount + 1
	    ID = cellist(kcell)%ID
        iclone = 0
        do jclone = 1,nclones
            if (clonelist(jclone)%ID == ID) then
                iclone = jclone
                exit
            endif
        enddo
        if (iclone == 0) then
            nclones = nclones + 1
            clonelist(nclones)%count = 1
            clonelist(nclones)%ID = ID
            clonelist(nclones)%avidity = p%avidity
            clonelist(nclones)%entrytime = cellist(kcell)%entrytime/60
        else
            clonelist(iclone)%count = clonelist(iclone)%count + 1
        endif
    endif
enddo
write(nfout,*)
write(nfout,*) 'Clones: total: ',totalcount
write(nfout,'(a)') '    ID  entry hr   avidity  count  fraction'
do iclone = 1,nclones
    clonelist(iclone)%fraction = clonelist(iclone)%count/real(totalcount)
    write(nfout,'(i6,f10.1,f10.3,i6,f10.4)') clonelist(iclone)%ID,clonelist(iclone)%entrytime, clonelist(iclone)%avidity, &
    clonelist(iclone)%count, clonelist(iclone)%fraction
enddo
deallocate(clonelist)
return

allocate(indx(nclones))
allocate(t(nclones))
do i = 1,nclones
    indx(i) = i
    t(i) = clonelist(i)%entrytime
enddo
call qsort(t,nclones,indx)     ! sorts in increasing order 
do i = 1,nclones
    iclone = indx(i)
    write(*,*) i, iclone, clonelist(iclone)%entrytime, clonelist(iclone)%count, clonelist(iclone)%fraction
enddo

deallocate(t)
deallocate(indx)
!deallocate(clonelist)
end subroutine

!-----------------------------------------------------------------------------------------
! Runs to compute the travel time distributions are based on varying two parameters:
!	TC_COGNATE_FRACTION	! fraction of T cells that are cognate initially
!	TC_TO_DC                ! number of T cells for every DC
! The first takes N_TRAVEL_COG values, the second N_TRAVEL_DC values, while the distribution
! is evaluated at N_TRAVEL_TIME values (intervals of a minute).
! The case to be simulated corresponds to (k_travel_cog, k_travel_dc)
! NOT USED
!-----------------------------------------------------------------------------------------
subroutine set_travel_params

TC_COGNATE_FRACTION = 0
TC_TO_DC = travel_dc(k_travel_dc)
ntravel = 0
travel_dist(k_travel_cog,k_travel_dc,:) = 0
use_traffic = .false.
use_DCflux = .false.
use_cognate = .false.
write(*,*) 'TC_COGNATE_FRACTION: ',TC_COGNATE_FRACTION
write(*,*) 'TC_TO_DC: ',TC_TO_DC
end subroutine

!-----------------------------------------------------------------------------------------
! Restrictions on the distances between different kinds of sites:
!    exit portal - boundary
!    exit portal - exit portal
!
!    HEV - boundary
!    HEV - exit portal
!    HEV - HEV
!
!    DC - boundary
!    DC - exit portal
!    DC - HEV
!    DC - DC
!
! exit portals
! ------------
! Initial exit portal placement has constraints on:
!	EXIT_SITE - BDRY_SITE
!	EXIT_SITE - EXIT_SITE
! Subsequent exit portal placement (adding or moving exit portals) has
! additional constraints on:
!	EXIT_SITE - HEV_SITE
!	EXIT_SITE - DC_SITE
!
! HEV sites
! ---------
! Initial HEV placement has constraints on:
!	HEV_SITE - BDRY_SITE
!	HEV_SITE - EXIT_SITE
!	HEV_SITE - HEV_SITE
! Subsequent HEV placement (adding or moving HEV sites) has
! additional constraints on:
!	HEV_SITE - DC_SITE
!
! DC sites
! --------
! Initial DC placement has constraints on:
!	DC_SITE - BDRY_SITE
!	DC_SITE - EXIT_SITE
!	DC_SITE - HEV_SITE
!	DC_SITE - DC_SITE
! The same constraints apply if DCs are added or moved.
!
! R_limit(:,:) defines the regional constrants on placement of HEVs and DCs
! e.g. DCs must lie within the range R_limit(DC_SITE,1) < r/R < R_limit(DC_SITE,2)
! where r is the distance of the DC from the centre, R is the blob radius.
!-----------------------------------------------------------------------------------------
subroutine setup_proximity_limits

proximity_limit(BDRY_SITE,BDRY_SITE) = 0
proximity_limit(EXIT_SITE,BDRY_SITE) = 0
proximity_limit(HEV_SITE,BDRY_SITE) = 2
proximity_limit(DC_SITE,BDRY_SITE) = 3      ! was 3

proximity_limit(EXIT_SITE,EXIT_SITE) = 4    ! was 6
proximity_limit(HEV_SITE,EXIT_SITE) = 4
proximity_limit(DC_SITE,EXIT_SITE) = 2      ! was 4

proximity_limit(HEV_SITE,HEV_SITE) = 5
proximity_limit(DC_SITE,HEV_SITE) = 2

proximity_limit(DC_SITE,DC_SITE) = 2        ! was 5

proximity_limit(BDRY_SITE,EXIT_SITE) = proximity_limit(EXIT_SITE,BDRY_SITE)
proximity_limit(BDRY_SITE,HEV_SITE) = proximity_limit(HEV_SITE,BDRY_SITE)
proximity_limit(BDRY_SITE,DC_SITE) = proximity_limit(DC_SITE,BDRY_SITE)
proximity_limit(EXIT_SITE,HEV_SITE) = proximity_limit(HEV_SITE,EXIT_SITE)
proximity_limit(EXIT_SITE,DC_SITE) = proximity_limit(DC_SITE,EXIT_SITE)
proximity_limit(HEV_SITE,DC_SITE) = proximity_limit(DC_SITE,HEV_SITE)

! Testing
R_limit(HEV_SITE,1) = R_HEV_min
R_limit(HEV_SITE,2) = R_HEV_max
R_limit(DC_SITE,1)  = R_DC_min
R_limit(DC_SITE,2)  = R_DC_max

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
integer function select_CD4_CD8
integer :: kpar = 0

if (par_uni(kpar) < CTYPE_FRACTION(CD8)) then
	select_CD4_CD8 = CD8
else
	select_CD4_CD8 = CD4
endif
end function

!-----------------------------------------------------------------------------------------
! The type (cognate or non-cognate) of a T cell can be random or strictly determined.
! Is TC_COGNATE_FRACTION(CD4) the fraction of all T cells that are cognate CD4 cells,
! or is it the fraction of CD4 cells that are cognate?
! I think it should be the latter.
!-----------------------------------------------------------------------------------------
subroutine select_cell_type(ctype,cognate,kpar)
integer :: ctype, kpar
logical :: cognate
integer :: nratio

ctype = select_CD4_CD8()
!if (random_cognate) then
    if (par_uni(kpar) < TC_COGNATE_FRACTION(ctype)) then
		cognate = .true.
	else
		cognate = .false.
	endif
!else
!    ! Needs to be fixed to account for two cognate populations
!!	if (mod(istep,6*4*60) == 1) then	! update Fcognate every 6 hours - otherwise mod(k_nonrandom,nratio) messed up
!!		Fcognate = TC_COGNATE_FRACTION - (scale_factor*Ncogseed)/NTC_BODY
!!	endif
!!    nratio = 1./Fcognate
!    nratio = 1./TC_COGNATE_FRACTION(ctype)
!    k_nonrandom(ctype) = k_nonrandom(ctype) + 1
!    if (mod(k_nonrandom(ctype),nratio) == 0) then
!!        select_cell_type = COG_TYPE_TAG
!		cognate = .true.
!    else
!!        select_cell_type = NONCOG_TYPE_TAG
!		cognate = .false.
!    endif
!endif
end subroutine

!--------------------------------------------------------------------------------
! We display only T cells that are still in the lymphnode
!--------------------------------------------------------------------------------
subroutine save_cell_positions
!!!use ifport
integer :: k, kcell, site(3), j, idc, dcsite(3)
!integer :: dcstate = 1
real :: dcstate
integer :: itcstate, stype, ctype, stage, region
real :: Tcell_diam = 0.9
real :: DC_diam = 1.8
!real :: spectrum_max = 10, spectrum_freefraction = 0.9
integer :: gen, bnd(2)
logical :: ex
character*(12) :: fname = 'cell_pos.dat'
character*(9) :: removefile = 'TO_REMOVE'

if (simulation_start) then
	inquire(file=fname,exist=ex)
	if (ex) then
		call unlink(fname)
	endif
	inquire(file=removefile,exist=ex)
	if (ex) then
		call unlink(removefile)
	endif
endif
simulation_start = .false.

if (.not.clear_to_send) then
	! wait until the file called removefile exists, then remove it
	inquire(file=removefile,exist=ex)
	if (.not.ex) then
!		call logger('wait')
		do
			inquire(file=removefile,exist=ex)
			if (.not.ex) then
!				call millisleep(10) ! no good at all
			else
				exit
			endif
		enddo
	endif
	call unlink(removefile)
	clear_to_send = .true.
endif


!inquire(file=fname,exist=ex)
!if (ex) then
!	ex = .false.
!!	write(*,*) 'waiting for file removal ...'
!	call logger('waiting for file removal ...')
!	do
!		inquire(file=fname,exist=ex)
!		if (ex) then
!!!!			call sleepqq(100)
!			call logger('millisleep(10)')
!			call millisleep(10)
!		else
!			exit
!		endif
!	enddo
!!	write(*,*) 'file removed'
!	call logger('file removed')
!endif

open(nfpos,file=fname,status='new')
! DC section
if (NDC > 0) then
    do k = 1,NDC
        if (DClist(k)%alive) then
            site = DClist(k)%site
            dcstate = min(1.0,DClist(k)%density/DC_ANTIGEN_MEDIAN)
!            dcstate = 1
            ! Need dcstate to convey antigen density level (normalized to 0-1)
            write(nfpos,'(a2,i4,3i4,f4.1,f5.2)') 'D ',k-1, site, DC_diam, dcstate
!            write(logmsg,'(a2,i4,3i4,f4.1,f5.2)') 'D ',k-1, site, DC_diam, dcstate
!            call logger(logmsg)
        endif
    enddo
endif

if (.not.IV_SHOW_NONCOGNATE) then
	! T cell section
	do k = 1,lastcogID
		kcell = cognate_list(k)
		if (kcell > 0) then
			call get_stage(cellist(kcell)%cptr,stage,region)
			if (region /= LYMPHNODE) cycle
			site = cellist(kcell)%site
	!        tcstate = mod(kcell,2) + 1
			gen = get_generation(cellist(kcell)%cptr)
			bnd = cellist(kcell)%DCbound
	!        if (bnd(1) == 0 .and. bnd(2) == 0) then
	!            tcbound = 0
	!        else
	!            tcbound = 1
	!        endif
!			if (get_stage(cellist(kcell)%cptr) == NAIVE) then
			if (stage == NAIVE) then
				itcstate = 0
			else
				if (bnd(1) == 0 .and. bnd(2) == 0) then
	!				tcstate = (gen-1.0)/(TC_MAX_GEN-1.0)*spectrum_max*spectrum_freefraction
	!				tcstate = (gen/TC_MAX_GEN)*spectrum_max*spectrum_freefraction
					itcstate = gen
				else
	!				tcstate = spectrum_max
					itcstate = 99
				endif
			endif
			! Need tcstate to convey non-activated status, i.e. 0 = non-activated
			write(nfpos,'(a2,i6,3i4,f4.1,i3)') 'T ',k-1, site, Tcell_diam, itcstate
		endif
	enddo
	! Bond section
	do k = 1,lastcogID
		kcell = cognate_list(k)
		if (kcell > 0) then
			call get_stage(cellist(kcell)%cptr,stage,region)
			if (region /= LYMPHNODE) cycle
			site = cellist(kcell)%site
			do j = 1,2
				idc = cellist(kcell)%DCbound(j)
				if (idc /= 0) then
					if (DClist(idc)%capable) then
						dcsite = DClist(idc)%site
	!                    dcstate = mod(idc,2) + 1
						dcstate = 1
	!                    write(nfpos,'(a,i8,3i4,4i4)') 'Node: ',nd, site, dcsite-site, dcstate
						write(nfpos,'(a2,2i5)') 'B ',k-1,idc-1
					endif
				endif
			enddo
		endif
	enddo
else
	! T cell section
	do kcell = 1,nlist
		if (cellist(kcell)%ID == 0) cycle  ! gap
		site = cellist(kcell)%site
		bnd = cellist(kcell)%DCbound
		ctype = cellist(kcell)%ctype
!		stype = struct_type(ctype)
		if (associated(cellist(kcell)%cptr)) then
			stype = COG_TYPE_TAG
		else
			stype = NONCOG_TYPE_TAG
		endif
		if (stype == NONCOG_TYPE_TAG) then
			itcstate = -1
		else
			call get_stage(cellist(kcell)%cptr,stage,region)
			if (region /= LYMPHNODE) cycle
			gen = get_generation(cellist(kcell)%cptr)
!			if (get_stage(cellist(kcell)%cptr) == NAIVE) then
			if (stage == NAIVE) then
				itcstate = 0
			else
				if (bnd(1) == 0 .and. bnd(2) == 0) then
					itcstate = gen
				else
					itcstate = 99
				endif
			endif
		endif
		! Need tcstate to convey non-activated status, i.e. 0 = non-activated
		write(nfpos,'(a2,i6,3i4,f4.1,i3)') 'T ',kcell-1, site, Tcell_diam, itcstate
	enddo
	! Bond section
	do kcell = 1,nlist
		if (cellist(kcell)%ID == 0) cycle  ! gap
		if (associated(cellist(kcell)%cptr)) then
			call get_stage(cellist(kcell)%cptr,stage,region)
			if (region /= LYMPHNODE) cycle
		endif
		site = cellist(kcell)%site
		do j = 1,2
			idc = cellist(kcell)%DCbound(j)
			if (idc /= 0) then
				if (DClist(idc)%capable) then
					dcsite = DClist(idc)%site
!                    dcstate = mod(idc,2) + 1
					dcstate = 1
!                    write(nfpos,'(a,i8,3i4,4i4)') 'Node: ',nd, site, dcsite-site, dcstate
					write(nfpos,'(a2,2i5)') 'B ',kcell-1,idc-1
				endif
			endif
		enddo
	enddo
endif
write(nfpos,'(a2,i6)') 'E ',istep
close(nfpos)

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine save_exnode
integer :: k, kcell, site(3), nd, j, idc, dcsite(3)
integer :: dcstate = 1
real :: tcstate
real :: Tcell_diam = 0.9
real :: DC_diam = 1.8
real :: spectrum_max = 10, spectrum_freefraction = 0.6
integer :: gen, stage, region, bnd(2)
character*(64) :: fname = '\CMGUI\DCU\dcu_00000.exnode'

write(fname(16:20),'(i5.5)') istep
write(*,*) fname
open(nfcmgui,file=fname,status='replace')

nd = 0

! DC section
if (NDC > 0) then
    write(nfcmgui,*) 'Group name : DC'
    write(nfcmgui,*) '  #Fields=3'
    write(nfcmgui,*) '  1) coordinates, coordinate, rectangular cartesian, #Components=3'
    write(nfcmgui,*) '    x.  Value index= 1, #Derivatives= 0'
    write(nfcmgui,*) '    y.  Value index= 2, #Derivatives= 0'
    write(nfcmgui,*) '    z.  Value index= 3, #Derivatives= 0'
    write(nfcmgui,*) '  2) diameter, field, real, #Components=1'
    write(nfcmgui,*) '    value.'
    write(nfcmgui,*) '  3) dcstate, field, real, #Components=1'
    write(nfcmgui,*) '    value.'

    do k = 1,NDC
        if (DClist(k)%alive) then
            nd = nd+1
            site = DClist(k)%site
            dcstate = 1
            write(nfcmgui,'(a,i8,3i4,f4.1,i3)') 'Node: ',nd, site, DC_diam, dcstate
        endif
    enddo
endif

! T cell section
write(nfcmgui,*) 'Group name : Tcell'
write(nfcmgui,*) '  #Fields=3'
write(nfcmgui,*) '  1) coordinates, coordinate, rectangular cartesian, #Components=3'
write(nfcmgui,*) '    x.  Value index= 1, #Derivatives= 0'
write(nfcmgui,*) '    y.  Value index= 2, #Derivatives= 0'
write(nfcmgui,*) '    z.  Value index= 3, #Derivatives= 0'
write(nfcmgui,*) '  2) diameter, field, real, #Components=1'
write(nfcmgui,*) '    value.'
write(nfcmgui,*) '  3) tcstate, field, real, #Components=1'
write(nfcmgui,*) '    value.'
!write(nfcmgui,*) '  3) tcgen, field, integer, #Components=1'
!write(nfcmgui,*) '    value.'
!write(nfcmgui,*) '  4) tcbound, field, integer, #Components=1'
!write(nfcmgui,*) '    value.'

do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell > 0) then
		call get_stage(cellist(kcell)%cptr,stage,region)
		if (region /= LYMPHNODE) cycle
        nd = nd+1
        site = cellist(kcell)%site
!        tcstate = mod(kcell,2) + 1
        gen = get_generation(cellist(kcell)%cptr)
        bnd = cellist(kcell)%DCbound
!        if (bnd(1) == 0 .and. bnd(2) == 0) then
!            tcbound = 0
!        else
!            tcbound = 1
!        endif
        if (bnd(1) == 0 .and. bnd(2) == 0) then
            tcstate = (gen-1.0)/(TC_MAX_GEN-1.0)*spectrum_max*spectrum_freefraction
        else
            tcstate = spectrum_max
        endif
!        write(nfcmgui,'(a,i8,3i4,f4.1,i3,i2)') 'Node: ',nd, site, Tcell_diam, gen, tcbound
        write(nfcmgui,'(a,i8,3i4,f4.1,f6.2)') 'Node: ',nd, site, Tcell_diam, tcstate
    endif
enddo

! Bond section
write(nfcmgui,*) 'Group name : Bond'
write(nfcmgui,*) '  #Fields=3'
write(nfcmgui,*) '  1) coordinates, coordinate, rectangular cartesian, #Components=3'
write(nfcmgui,*) '    x.  Value index= 1, #Derivatives= 0'
write(nfcmgui,*) '    y.  Value index= 2, #Derivatives= 0'
write(nfcmgui,*) '    z.  Value index= 3, #Derivatives= 0'
write(nfcmgui,*) '  2) vector, coordinate, rectangular cartesian, #Components=3'
write(nfcmgui,*) '    x.  Value index= 1, #Derivatives= 0'
write(nfcmgui,*) '    y.  Value index= 2, #Derivatives= 0'
write(nfcmgui,*) '    z.  Value index= 3, #Derivatives= 0'
write(nfcmgui,*) '  3) dcstate, field, real, #Components=1'
write(nfcmgui,*) '    value.'

do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell > 0) then
        site = cellist(kcell)%site
        do j = 1,2
            idc = cellist(kcell)%DCbound(j)
            if (idc /= 0) then
                if (DClist(idc)%capable) then
                    nd = nd+1
                    dcsite = DClist(idc)%site
!                    dcstate = mod(idc,2) + 1
                    dcstate = 1
                    write(nfcmgui,'(a,i8,3i4,4i4)') 'Node: ',nd, site, dcsite-site, dcstate
                endif
            endif
        enddo
    endif
enddo

close(nfcmgui)
end subroutine

!--------------------------------------------------------------------------------
! The ability of a DC to provide TCR stimulation is conveyed by the current
! value of antigen density, %density.  The level decays over time, and possibly
! is reduced by every episode of TCR stimulation, at the rate DC_DENS_BY_STIM.
! When t > %dietime - DC_ACTIV_TAPER the density reduces each time step by a factor
! such that when %dietime is reached the product of the factors = 0.1
! tfactor**n = 0.1 where n = DC_ACTIV_TAPER/DELTA_T
! tfactor = (0.1)**(1/n)
! A DC loses its TCR stimulating capability when %density falls below the
! limiting value STIM_HILL_THRESHOLD.
! A DC is scheduled to die when t > %dietime.  At this point %capability is set
! to .false., but the DC waits until %nbound = 0 before it dies.
!--------------------------------------------------------------------------------
subroutine update_DCstate(ok)
logical :: ok
real :: tnow, tfactor, decay_factor
integer :: idc, nalive, nbound, ncbound, ncapable

if (NDCalive == 0) then
	NDCcapable = 0
	ok = .true.
	return
endif
ncapable = 0
tfactor = 0.1**(DELTA_T/(60*DC_ACTIV_TAPER))
tnow = istep*DELTA_T
decay_factor = 1 - DCdecayrate*DELTA_T
nalive = 0
do idc = 1,NDC
    if (DClist(idc)%alive) then
        if (.not.IN_VITRO .and. DC_outside(idc)) then	! A kludge to remove DCs stranded outside the blob
            DClist(idc)%dietime = min(DClist(idc)%dietime,tnow)
		endif
        nbound = DClist(idc)%nbound
        ncbound = DClist(idc)%ncogbound
        if (save_DCbinding) then
            dcbind(ncbound) = dcbind(ncbound) + 1
        endif
        if (tnow > DClist(idc)%dietime) then
            DClist(idc)%capable = .false.
            if (nbound /= 0) then
            elseif (nbound == 0) then
!				write(logmsg,'(a,i4,f8.2,i6)') 'DC dies: ',idc,DClist(idc)%dietime,NDCalive
!				call logger(logmsg)
                DClist(idc)%alive = .false.
                ndeadDC = ndeadDC + 1
                DCdeadlist(ndeadDC) = idc
                ! later we need to adjust occupancy to reflect the death of DC idc
            else
                write(logmsg,'(a,2i6)') 'Error: update_DCstatus: nbound < 0: ',idc,nbound
                call logger(logmsg)
				ok = .false.
				return
            endif
        endif
        if (DClist(idc)%capable) then
            DClist(idc)%density = DClist(idc)%density*decay_factor - DC_DENS_BY_STIM*DClist(idc)%stimulation
            DClist(idc)%stimulation = 0     ! reset the tally of stimulation delivered to zero (OMP)
            if (tnow > DClist(idc)%dietime - 60*DC_ACTIV_TAPER) then	! antigen density decays over DC_ACTIV_TAPER
                DClist(idc)%density = DClist(idc)%density*tfactor
            endif
            ! Do not make DCs die when the antigen level drops to zero (Philippe)
!            if (DClist(idc)%density < STIM_HILL_THRESHOLD) then
!                DClist(idc)%capable = .false.
!                if (incapable_DC_dies) then
!                    ! A non-capable DC might as well die - but this would eliminate low-antigen DCs
!                    DClist(idc)%dietime = tnow + 60     ! give it 1 hour to die
!                endif
!!                write(*,*) 'DC incapable: ',idc
!            endif
        endif
    endif
    if (DClist(idc)%alive) nalive = nalive + 1
    if (DClist(idc)%capable) ncapable = ncapable + 1
enddo
NDCalive = nalive
NDCcapable = ncapable
!write(nflog,*) 'update_DCstate: live DCs: ',nalive
if (nalive == 0) then
    write(logmsg,*) 'No live DC'
    call logger(logmsg)
endif
!if (ncapable == 0) then
!	call logger("No capable DC")
!	return
!endif
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! If a DC is left outside the blob as it contracts, the simplest thing to do is to make
! it die.  Better if the DC moves.
!-----------------------------------------------------------------------------------------
logical function DC_outside(idc)
integer :: idc
real :: d

d = cdistance(DClist(idc)%site)
if (d - Radius > 1) then
	DC_outside = .true.
else
	DC_outside = .false.
endif
end function

!-----------------------------------------------------------------------------------------
! When a DC dies it is labelled as dead in the DC list (%alive = .false.) and the count
! of DCs that have died is incremented.
! Later, after gather_data(), the DCs that have died are cleared out of occupancy()%indx().
! The reassignment of occupancy()%DC values is done all at once in
! reassign_DC(), called from balancer().
!-----------------------------------------------------------------------------------------
subroutine clearDC(idc)
integer :: idc
integer :: k, site0(3), site(3)

site0 = DClist(idc)%site
occupancy(site0(1),site0(2),site0(3))%indx = 0
occupancy(site0(1),site0(2),site0(3))%DC = 0
if (NDCsites > 1) then  ! need to clear non-central sites occupied by the DC
    do k = 2,NDCsites
        site = site0 + DCoffset(:,k)
        if (occupancy(site(1),site(2),site(3))%indx(1) == -idc) then
            occupancy(site(1),site(2),site(3))%indx = 0
            occupancy(site(1),site(2),site(3))%DC = 0
        endif
    enddo
endif
DClist(idc)%ID = 0
DClist(idc)%density = 0
DClist(idc)%stimulation = 0
DClist(idc)%dietime = 0
DClist(idc)%site = 0
DClist(idc)%nbound = 0
DClist(idc)%ncogbound = 0
deallocate(DClist(idc)%cogbound)
if (IN_VITRO) then	! Need to adjust zrange2D()
	call resetDCzrange(site0(1),site0(2))
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Set the zrange2D() values for sites near a DC at (xdc,ydc).
! For now, just use DCstack().
! Note that zrange2D was initialized to (1,1) for all interior sites before DC placement.
! Since the ranges of influence of DCs may overlap, both lower and upper z limits can
! only increase those values set for previously encountered DCs.
!-----------------------------------------------------------------------------------------
subroutine setDCzrange(xdc,ydc)
integer :: xdc, ydc
integer :: x,y,z,dx,dy,zval(3)

do dx = -4,4
	do dy = -4,4
		x = xdc + dx
		y = ydc + dy
		zval = DCstack(dx,dy,:)
		if (zval(1) == 0) cycle
		do z = 1,3
			if (zval(z) > 0) then
				zrange2D(x,y,1) = max(z,zrange2D(x,y,1))
				exit
			endif
		enddo
		do z = 3,1,-1
			if (zval(z) > 0) then
				zrange2D(x,y,2) = max(z,zrange2D(x,y,2))
				exit
			endif
		enddo
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! When a DC dies the zrange2D() values in its vicinity must be adjusted.
!-----------------------------------------------------------------------------------------
subroutine resetDCzrange(xdc,ydc)
integer :: xdc, ydc
integer :: x,y,z,dx,dy,zval(3)

do dx = -4,4
	do dy = -4,4
		x = xdc + dx
		y = ydc + dy
		zval = DCstack(dx,dy,:)
		if (zval(1) == 0) cycle
		if (zrange2D(x,y,2) == 1) cycle
		do z = 3,1,-1
			if (zval(z) > 1) then
				if (zval(z) == zrange2D(x,y,2)) then
					zrange2D(x,y,:) = 1
					write(nflog,*) 'resetDCzrange: ',xdc,ydc,dx,dy,zval(z)
					exit
				endif
			endif
		enddo
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Determine the number of DC that enter the paracortex between times t1 and t2.
! For now just use a rate of arrivals/min., declining to zero between times
! tdc1, tdc2 (hours).
! This rate must be scaled to be proportional to the steady-state number of
! T cells, i.e. NTcells0
! A rate of 0.1 is not bad for NTcells0 = 100k
! It also must be scaled by the specified rate of occurrence of DCs, i.e. to
! take account of NT_TO_DC.  The value NT_TO_DC = 200 can be used as a baseline.
! THIS IS HARD_WIRED, NEEDS TO DEPEND ON STATE OF INFECTION IN THE TISSUE
! In the case of DCs injected into experimental animals, the schedule of DCs
! estimated to enter the paracortex in each hour is read from an input file.
! The number of DCs entering is precomputed for each hour, scaled by the initial
! T cell population.
!--------------------------------------------------------------------------------
integer function DCinflux(t1,t2,kpar)
integer :: kpar
integer :: ih1, ih2
real :: t1, t2
real :: rate, temp, dn
real :: DCrate0
real(DP) :: R
!real, save :: dn_last = 0

!write(logmsg,*) 'DCinflux: dn_last ',dn_last
!call logger(logmsg)
if (DC_INJECTION) then
	ih1 = t1/60
	ih2 = t2/60
	if (ih1 < ih2) then
		DCinflux = DCinjected(ih2)
	else
		DCinflux = 0
	endif
else
	DCrate0 = DC_FACTOR*DCrate_100k*(NTcells0/1.0e5)     ! DC influx is scaled by the initial T cell population
	!DCrate0 = DCrate0*(200./TC_TO_DC)                              ! and scaled by 1/TC_TO_DC (TRY CANCELLING THIS)
	if (t1 > T_DC2*60) then
		DCinflux = 0
		return
	elseif (t1 < T_DC1*60) then
		rate = DCrate0      ! rate /min
	else
		rate = DCrate0*(T_DC2 - t1/60)/(T_DC2 - T_DC1)
	endif
	if (RANDOM_DCFLUX) then
		temp = rate*(t2-t1)
		DCinflux = temp
		dn = temp - DCinflux
		R = par_uni(kpar)
		if (R < dn) then
			DCinflux = DCinflux + 1
		endif
	else
		! Try making the DC influx deterministic
		temp = rate*(t2-t1) + DCinflux_dn_last
		DCinflux = temp
		DCinflux_dn_last = temp - DCinflux
	endif
endif
!write(logmsg,*) 'DC rate: ',rate,DCinflux,dn_last
!call logger(logmsg)
end function

!--------------------------------------------------------------------------------
! Determines the degree to which a cell is subject to chemotaxis.
! (Determined by CD69 level, or S1PR1 level.)
! If use_exit_chemotaxis is true, i.e. chemotaxis is used to control cell exit, any cell
! that gets close enough to the exit will leave the paracortex.  Is this
! acceptable?
! For a noncognate cell, should rise from 0 to 1 in an hour or so.
!--------------------------------------------------------------------------------
real function chemo_active_exit(cell)
type(cell_type), pointer :: cell
real :: tnow, t

if (turn_off_chemotaxis) then
    chemo_active_exit = 0
    return
endif

if (RELAX_INLET_EXIT_PROXIMITY) then
    tnow = istep*DELTA_T
    t = tnow - cell%entrytime
    if (t > CHEMO_K_RISETIME) then
		chemo_active_exit = 1
	else
		chemo_active_exit = t/CHEMO_K_RISETIME
	endif
	return
endif

if (TAGGED_EXIT_CHEMOTAXIS) then
    tnow = istep*DELTA_T
    t = tnow - cell%entrytime
    if (cell%tag == RES_TAGGED_CELL) then     ! testing effect of S1PR1
        chemo_active_exit = (1 - exp(-K1_S1PR1*t))*TAGGED_CHEMO_ACTIVITY
    else
		chemo_active_exit = 0
	endif
	return
endif

if (TAGGED_LOG_PATHS) then
	if (cell%tag == CHEMO_TAGGED_CELL) then
		write(*,*) 'chemo_active_exit: CHEMO_K_EXIT is no longer used'
		stop
		chemo_active_exit = CHEMO_K_EXIT
	else
		chemo_active_exit = 0
	endif
	return
endif

if (associated(cell%cptr)) then     ! cognate cell
    chemo_active_exit = cell%cptr%S1PR1
!    if (CD69 < CD69_threshold) then
!        chemo_active = 1 - CD69/CD69_threshold
!    else
!        chemo_active = 0
!    endif
else
    tnow = istep*DELTA_T
    t = tnow - cell%entrytime
!    if (cell%tag == RES_TAGGED_CELL) then     ! testing effect of S1PR1
!        chemo_active = 1 - exp(-K1_S1PR1*0.1*t)
!        chemo_active = 0
!    else
        chemo_active_exit = 1 - exp(-K1_S1PR1*t)
!    endif
endif
end function

!--------------------------------------------------------------------------------
! Returns the level of CCR7 ligand (i.e. CCL19/21) at a distance r sites from an
! exit site, within the exit SOI.  The value must be in the range (0,1).
! The parameters CCR7_R1 and CCR7_R2 specify the range of r over which the
! ligand level ranges linearly from 0 to 1.  This is a simple way to program a
! decreased level of CCR7 near an exit, following Cyster's ideas.
! Note: currently this always returns 1
!--------------------------------------------------------------------------------
real function CCR7_ligand(r)
real :: r
real, parameter :: CCR7_R1 = 0, CCR7_R2 = 0

if (r < CCR7_R1) then
    CCR7_ligand = 0
elseif (r < CCR7_R2) then
    CCR7_ligand = (r - CCR7_R1)/(CCR7_R2 - CCR7_R1)
else
    CCR7_ligand = 1
endif
end function

!--------------------------------------------------------------------------------
! Returns cell's susceptibility to DC chemotaxis
!--------------------------------------------------------------------------------
real function chemo_active_DC(cell)
type(cell_type), pointer :: cell

if (turn_off_chemotaxis) then
    chemo_active_DC = 0
    return
endif
!if (associated(cell%cptr)) then     ! cognate cell
!    chemo_active_DC = cell%cptr%DCchemo
!else
!	chemo_active_DC = BASE_DCchemo
!endif
!chemo_active_DC = cell%DCchemo
chemo_active_DC = cell%receptor_level(CCR1)
end function

!--------------------------------------------------------------------------------
! Determines e(:), the location of the nearest exit, if there is a close one.
! Otherwise returns e(:) = 0.
! Note: currently only one exit number is stored in occupancy()
!--------------------------------------------------------------------------------
subroutine nearest_exit(site,in_exit_SOI,e)
integer :: site(3), e(3)
logical :: in_exit_SOI
integer :: iexit

iexit = occupancy(site(1),site(2),site(3))%exitnum
if (iexit == 0) then
    in_exit_SOI = .false.
    e = (/0,0,0/)
    return
endif
in_exit_SOI = .true.
e = exitlist(abs(iexit))%site
!write(*,*) 'exit site: ',e
end subroutine

!--------------------------------------------------------------------------------
! The criterion for a near exit is based on chemo_radius
!--------------------------------------------------------------------------------
subroutine near_exits(site,ne,ee)
integer :: site(3), ne, ee(3,*)
integer :: iexit

iexit = occupancy(site(1),site(2),site(3))%exitnum
if (iexit == 0) then
	ne = 0
    return
endif
ne = 1
ee(:,1) = exitlist(abs(iexit))%site
!write(*,*) 'exit site: ',e
end subroutine

!--------------------------------------------------------------------------------
! This determines which DCs (if any) the site is within chemotactic range of.
! The criterion for a near DC is based on first sites within DC_RADIUS (%DC)
! then sites not in the %DC list but within chemo_radius.
!--------------------------------------------------------------------------------
subroutine near_DCs(site,nd,near)
integer :: site(3), nd, near(*)
integer :: i
integer(2) :: DC(0:max(DCDIM,cDCDIM)-1)     ! DC(0) = number of near DCs, DC(k) = ID of kth near DC

nd = 0
DC(0:DCDIM-1) = occupancy(site(1),site(2),site(3))%DC
if (DC(0) /= 0) then
	do i = 1,DC(0)
		if (DClist(DC(i))%alive) then
			nd = nd + 1
			near(nd) = DC(i)
			if (nd == MAX_DC_CHEMO) return
		endif
	enddo
endif
DC(0:cDCDIM-1) = occupancy(site(1),site(2),site(3))%cDC
if (DC(0) /= 0) then
	do i = 1,DC(0)
		if (DClist(DC(i))%alive) then
			nd = nd + 1
			near(nd) = DC(i)
			if (nd == MAX_DC_CHEMO) return
		endif
	enddo
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_big_occ
integer :: x,y,z,k,kcell,indx(2)
logical :: OK

OK = .true.
do x = 1,NX
    do y = 1,NY
        do z = 1,NZ
            indx = occupancy(x,y,z)%indx
            do k = 1,2
                kcell = indx(k)
                if (kcell < 0) then
                    if (kcell /= OUTSIDE_TAG .and. -kcell > NDC) then
                        OK = .false.
                        write(*,*) x,y,z,k,kcell
                    endif
                endif
            enddo
        enddo
    enddo
enddo
if (.not.OK) stop
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checkcell(str)
character*(*) :: str
integer :: ictest = 12

write(*,*) 'checkcell: ',str,'  ',ictest
!write(*,*) cellist(ictest)%ID,cellist(ictest)%PTR1%entrytime
end subroutine

!-----------------------------------------------------------------------------------------
! Checks global cell locations for tag validity.
!-----------------------------------------------------------------------------------------
subroutine check_tagged
integer :: k,d2,n,site(3)
type(cell_type) :: cell

do k = 1,nlist
    cell = cellist(k)
    if (cell%tag == TAGGED_CELL) then
        n = n+1
        site = cell%site
        if (.not.taggable(site)) then
            write(*,*) 'Bad tagged cell: ',site,d2,Radius*Radius
        endif
    endif
enddo
write(*,*) 'did check_tagged'
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
real function get_IL2store(p)
type(cog_type), pointer :: p
integer :: IL2store_index

if (.not.use_cytokines) then
    get_IL2store = 0
    return
endif
IL2store_index = NP_offset(cyt_seq(IL2_TAG)) + IL2_Store
get_IL2store = p%IL_state(IL2store_index)
end function

!-----------------------------------------------------------------------------------------
! Convert one of the sites near a DC into a DC peripheral site for DC idc.
! The index k indicates which peripheral site of 2:NDCsites to convert.
! A site is a candidate for a DC site provided:
! (1) It is not outside the blob, or already a DC site
! (2) It does not hold two T cells
! (3) It is not too close the the blob boundary
! If a candidate site is free, it is used.  If it holds a single T cell, the cell
! is moved if possible to a neighbouring site.
! Note that when a T cell is bumped it retains its binding (if any) to a DC, even though
! it may have been moved to a site that is - strictly speaking - not within the SOI of
! the DC.
!-----------------------------------------------------------------------------------------
subroutine addDCsite(idc,site0,k,err)
integer :: idc, site0(3), k, err
integer :: indx(2), jslot, kcell, freeslot, site1(3), site2(3)
real :: prox
logical :: OK

site1 = site0 + DCoffset(:,k)
indx = occupancy(site1(1),site1(2),site1(3))%indx
if (indx(1) < 0) then                       ! OUTSIDE_TAG or DC
    err = 1
    return
endif
if (indx(1) /= 0 .and. indx(2) /= 0) then   ! two T cells at this site
    err = 2
    return
endif
!if (tooNearDC(site1,NDC)) then    ! too close to another DC
!    err = 3
!    return
!endif
!prox = 0.5*DC_DCprox*DC_RADIUS
prox = bdry_DCprox*DC_RADIUS
if (tooNearBdry(site1,prox)) then                ! too close to the blob boundary
    err = 4
    return
endif

OK = .false.
jslot = 0
if (indx(1) /= 0) then
    jslot = 1
elseif (indx(2) /= 0) then
    jslot = 2
endif
if (jslot == 0) then ! free site, use it
    OK = .true.
else  ! one T cell here in jslot, it must be moved
    ! Can the T cell be bumped to a neighbour site?
    call get_free_slot(occupancy,NX,site1,site2,freeslot)
    if (freeslot == 0) then    ! cannot be bumped
        err = 5
        return
    endif
   ! Move the cell in site1/jslot to site2/freeslot
    kcell = indx(jslot)
    occupancy(site1(1),site1(2),site1(3))%indx = 0
    occupancy(site2(1),site2(2),site2(3))%indx(freeslot) = kcell
    cellist(kcell)%site = site2
    ! Now site1 is free to use
    OK = .true.
endif
!if (OK) then
    err = 0
    occupancy(site1(1),site1(2),site1(3))%indx = -idc
!else
!    err = 6
!endif
end subroutine

!-----------------------------------------------------------------------------------------
! A DC is allowed to move only if it has grown to its full extent,
! i.e. if DClist()%nsites = NDCsites
! For now allow only moves in the 6 principal directions (Neumann).
! Currently valid only for NDCsites = 7.
! Cells in the 5 sites in the path of the DC step are moved.
! For 4 sites, the shift is by one site, for the one site in line with the DC centre the
! shift is 3 sites.
!-----------------------------------------------------------------------------------------
subroutine moveDC(idc,dir)
integer :: idc, dir
integer :: k, i, kcell, indx(2), site0(3), site1(3), site2(3), step(3)

if (DClist(idc)%nsites /= NDCsites) return
step = neumann(:,dir)
site0 = DClist(idc)%site
!write(*,*) 'moveDC: ',idc,dir,'  ',site0,'  ',step
do k = 2,NDCsites
    if (all(DCoffset(:,k) == step)) then       ! this is in the direction of step
        ! move site contents by 3 sites in opposite direction to step
        site1 = site0 - step       ! old DC site
        site2 = site0 + 2*step     ! new DC site
!        write(*,'(a,7i4)') 'step dir: ',k,site1,site2
        indx = occupancy(site2(1),site2(2),site2(3))%indx
        if (any(indx == OUTSIDE_TAG)) then
            write(*,*) 'moveDC: site2 outside!'
            stop
        endif
        occupancy(site1(1),site1(2),site1(3))%DC = occupancy(site2(1),site2(2),site2(3))%DC
        occupancy(site1(1),site1(2),site1(3))%indx = indx
        do i = 1,2
            kcell = indx(i)
            if (kcell /= 0) then
                cellist(kcell)%site = site1
            endif
        enddo
        occupancy(site2(1),site2(2),site2(3))%indx = -idc
        occupancy(site0(1),site0(2),site0(3))%indx = -idc
        site2 = site0  + step
        occupancy(site2(1),site2(2),site2(3))%indx = -idc
    elseif (all(DCoffset(:,k) == -step)) then   ! this is the direction opposite to step
        ! do nothing
    else                                        ! one of the other 4 directions
        ! move site contents by 1 site in opposite direction to step
        site1 = site0 + DCoffset(:,k)       ! old DC site, new T cell site
        site2 = site1 + step                ! new DC site, old T cell site
!        write(*,'(a,7i4)') 'other dir: ',k,site1,site2
        indx = occupancy(site2(1),site2(2),site2(3))%indx
        if (any(indx == OUTSIDE_TAG)) then
            write(*,*) 'moveDC: site2 outside!'
            stop
        endif
        occupancy(site1(1),site1(2),site1(3))%DC = occupancy(site2(1),site2(2),site2(3))%DC
        occupancy(site1(1),site1(2),site1(3))%indx = indx
        do i = 1,2
            kcell = indx(i)
            if (kcell /= 0) then
                cellist(kcell)%site = site1
            endif
        enddo
        occupancy(site2(1),site2(2),site2(3))%indx = -idc
    endif
enddo
DClist(idc)%site = site0 + step
end subroutine

!-----------------------------------------------------------------------------------------
! All associations of sites with DC are recomputed.
! That is for every site (x,y,z), create the list of DC that are near this site:
! occupancy(x,y,z)%DC(1:3), where occupancy(x,y,z)%DC(0) = number of nearby DC (if >= 0)
! To avoid having to explicitly select the closest DC (when there are more than DCDIM-1 near
! a site), the order of scanning the DClist is randomized.
!-----------------------------------------------------------------------------------------
subroutine reassign_DC(kpar,ok)
integer :: kpar
logical :: ok
integer, allocatable :: perm(:)
integer :: xdc, ydc, zdc, xmin, xmax, ymin, ymax, zmin, zmax
integer :: idc, kdc, k, x, y, z, y2, z2, d2, site(3), nassigned
logical :: added

!write(*,*) 'reassign_DC'

occupancy(:,:,:)%DC(0) = 0
occupancy(:,:,:)%cDC(0) = 0
allocate(perm(MAX_DC))
do k = 1, NDC
    perm(k) = k
enddo

call permute(perm,NDC,kpar)
NDCalive = 0
do k = 1,NDC
    idc = perm(k)
    if (.not.DClist(idc)%alive) cycle
!    write(*,*) 'idc: ',idc
!	write(*,*) '(4) cell 16243: ',cellist(16243)%site,occupancy(70,63,88)%indx
	if (DClist(idc)%nsites < NDCsites) then
		call assignDCsites(idc,nassigned,ok)
		if (.not.ok) return
!		write(*,*) 'assigned DC sites for: ',idc,nassigned
		DClist(idc)%nsites = nassigned
	endif
    NDCalive = NDCalive + 1
    site = DClist(idc)%site
    xdc = site(1)
    ydc = site(2)
    zdc = site(3)
    xmin = xdc - DC_RADIUS
    xmax = xdc + DC_RADIUS
    ymin = ydc - DC_RADIUS
    ymax = ydc + DC_RADIUS
    zmin = zdc - DC_RADIUS
    zmax = zdc + DC_RADIUS
    xmin = max(1,xmin)
    xmax = min(NX,xmax)
    ymin = max(1,ymin)
    ymax = min(NY,ymax)
    zmin = max(1,zmin)
    zmax = min(NZ,zmax)
    do z = zmin,zmax
        z2 = (z-zdc)*(z-zdc)
        do y = ymin,ymax
            y2 = (y-ydc)*(y-ydc)
            do x = xmin,xmax
                d2 = (x-xdc)*(x-xdc) + y2 + z2
                added = .false.
                ! The following procedure is correct only if DC_RADIUS <= chemo_radius
	            if (d2 <= DC_RADIUS*DC_RADIUS) then
	                kdc = occupancy(x,y,z)%DC(0)
	                if (kdc < 0) cycle			! this is a DC site
		            if (kdc < DCDIM-1) then		! not all possible DC assigned
						kdc = kdc + 1
						occupancy(x,y,z)%DC(kdc) = idc
						occupancy(x,y,z)%DC(0) = kdc
						added = .true.
					endif
				endif
	            if (.not.added .and. (d2 <= chemo_radius*chemo_radius)) then
	                kdc = occupancy(x,y,z)%cDC(0)
		            if (kdc == cDCDIM-1) cycle   ! all possible DC assigned
		            kdc = kdc + 1
		            occupancy(x,y,z)%cDC(kdc) = idc
	                occupancy(x,y,z)%cDC(0) = kdc
	            endif
            enddo
        enddo
    enddo
enddo
deallocate(perm)
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! Try to expand any DCs that are not yet occupying their full allocation of sites.
!-----------------------------------------------------------------------------------------
subroutine growDC
integer :: k, idc, err, site0(3), site(3)

do idc = 1,NDC
    if (.not.DClist(idc)%alive) cycle
    if (DClist(idc)%nsites < NDCsites) then
        site0 = DClist(idc)%site
        do k = 2,NDCsites
            site = DClist(idc)%site + DCoffset(:,k)
            if (occupancy(site(1),site(2),site(3))%indx(1) /= -idc) then
                call addDCsite(idc,site0,k,err)
                if (err == 0) then
                    DClist(idc)%nsites = DClist(idc)%nsites + 1
                endif
            endif
        enddo
!        write(*,*) 'growDC: ',idc,DClist(idc)%nsites
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Test moveDC() with a randomly selected DC, random step direction.
!-----------------------------------------------------------------------------------------
subroutine test_moveDC
integer :: idc, dir
integer :: kpar=0

write(*,*) 'test_moveDC'
if (NDCsites /= 7) then
    write(*,*) 'Valid only for NDCsites = 7'
    stop
endif
do
    idc = random_int(1,NDC,kpar)
    if (DClist(idc)%alive) exit
enddo
dir = random_int(1,6,kpar)
call moveDC(idc,dir)
end subroutine

!-----------------------------------------------------------------------------------------
! Depends on use_HEV_portals (currently hard-coded)
! The number of HEV portals should increase as Ncells increases, or with inflowtotal
!-----------------------------------------------------------------------------------------
subroutine PlaceHEVPortals(ok)
logical :: ok
real(DP) :: R
integer :: kpar = 0
integer :: ihev, k, x, y, z, site(3)
integer :: nhev_required = 200	! testing - what criterion to use?

!write(nfout,*) 'HEV sites'
allocate(HEVlist(nhev_required))
do ihev = 1,nhev_required
	k = 0
	do
		k = k+1
		if (k > 10000) then
			call logger('PlaceHEVPortals: too many iterations')
			ok = .false.
			return
		endif
        R = par_uni(kpar)
        x = 1 + R*NX
        R = par_uni(kpar)
        y = 1 + R*NY
        R = par_uni(kpar)
        z = 1 + R*NZ
        site = (/x,y,z/)
        if (.not.portalOK(site)) cycle
		call CheckSite(HEV_SITE,site,ok)
		if (ok) exit
	enddo
	HEVlist(ihev)%site = site
	occupancy(site(1),site(2),site(3))%hevnum = ihev
!	write(nfout,'(i6,2x,3i4)') ihev,site
enddo
NHEV = nhev_required
end subroutine

!-----------------------------------------------------------------------------------------
! Checks to see if a site can be used as a portal site, either exit portal or HEV (if used)
!-----------------------------------------------------------------------------------------
logical function portalOK(site)
integer :: site(3)

portalOK = .false.
if (occupancy(site(1),site(2),site(3))%indx(1) < 0) return
if (occupancy(site(1),site(2),site(3))%exitnum < 0) return
if (occupancy(site(1),site(2),site(3))%hevnum /= 0) return
portalOK = .true.
end function


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine assignDCsites(idc,nassigned,ok)
integer :: idc, nassigned
logical :: ok
integer :: k, site(3), site1(3)

ok = .false.
site = DClist(idc)%site
occupancy(site(1),site(2),site(3))%indx = -idc     ! This site holds a DC (centre)
occupancy(site(1),site(2),site(3))%DC(0) = -idc    ! This site holds a DC
nassigned = 1
!write(*,*) 'assignDCsites: ',site
do k = 2,NDCsites
    site1 = site + DCoffset(:,k)
    if (occupancy(site1(1),site1(2),site1(3))%indx(1) == -idc) then
		nassigned = nassigned + 1
		cycle
	endif
    ! Reduce available blob sites only if the site is within the blob.
    ! A DC site should never fall outside, I think.  Check this.
    if (.not.inside_xyz(site1)) then
		write(*,*) 'Error: assignDCsites: site outside grid: '
		write(*,*) 'DC site: ',site
		write(*,*) 'DCoffset: ',k,DCoffset(:,k)
		write(*,*) 'site1: ',site1
		call logger("Error: assignDCsites: site outside grid")
		stop
	endif
    if (occupancy(site1(1),site1(2),site1(3))%indx(1) == OUTSIDE_TAG) then
		if (.not.IN_VITRO .or. (IN_VITRO .and. site1(3) == 1)) then
			call logger('assignDCsites: DC hits boundary')
			return
		endif
	endif
    if (IN_VITRO) then
		if (site1(3) == 1 .and. occupancy(site1(1),site1(2),site1(3))%indx(1) /= OUTSIDE_TAG) then
			nassigned = nassigned + 1
		endif
	else
		if (occupancy(site1(1),site1(2),site1(3))%indx(1) /= 0) cycle
		if (occupancy(site1(1),site1(2),site1(3))%indx(2) /= 0) cycle
		nassigned = nassigned + 1
	endif
    occupancy(site1(1),site1(2),site1(3))%indx = -idc     ! This site holds a DC (soma)
    occupancy(site1(1),site1(2),site1(3))%DC(0) = -idc    ! This site holds a DC
enddo
ok = .true.
end subroutine

!-----------------------------------------------------------------------------
! Number of exit portals required for the current cell population, ncells. 
! For SURFACE_PORTALS, linearly proportional to surface area, i.e. to
! ncells^(2/3), factor = exit_fraction
! Modified to use an adjustment that was computed with:
!	Tres = 12
!	chemotaxis factor = 1
!	chemotaxis radius = 100 um
! exit_fraction = mN + c
! where
!	N = ncells
!	m = 1.607E-5
!	c = 0.00602
! 
! Better is a quadratic fit (from Excel, steady_chemo_0_1.xls):
! exit_fraction = a.x^2 + b.x + c
! where x = Ncells/1000
!
! The number of exit portals is based on the residence time for CD4 cells
!-----------------------------------------------------------------------------
integer function requiredExitPortals(ncells)
integer :: ncells
real :: a, b, c, x, Fe
real, parameter :: pow = 2./3.
!real, parameter :: m = 1.607E-8, c = 0.00602
! parameters for chemotaxis
!real, parameter :: a_chemo_12h = -2.228E-08
!real, parameter :: b_chemo_12h = 2.624E-05
!real, parameter :: c_chemo_12h = 5.510E-03

! parameters for no chemotaxis, INLET_R_FRACTION = 1.0
! Best fit:
! y = -3.873E-08x2 + 5.150E-05x + 1.029E-02
!real, parameter :: a_nochemo_24h = -3.873E-08
!real, parameter :: b_nochemo_24h = 5.150E-05
!real, parameter :: c_nochemo_24h = 1.029E-02

! parameters for no chemotaxis, INLET_R_FRACTION = 0.7
! Best fit:
!y = -4.058E-08x2 + 5.590E-05x + 1.006E-02
real, parameter :: a_nochemo_24h = -4.058E-08
real, parameter :: b_nochemo_24h = 5.590E-05
real, parameter :: c_nochemo_24h = 1.006E-02

! This is to adjust Ne when there is general exit chemotaxis, to generate steady-state
real, parameter :: K_Ke = 1.0		! 0.68 for Ke = 1.0, 1.0 for Ke = 0.0 (Pe = 0.02)

! TESTING
!if (use_exit_chemotaxis) then 
	a = a_nochemo_24h
	b = b_nochemo_24h
	c = c_nochemo_24h
! try making it constant
!	a = 0
!	b = 0
! Calibration...
!	c = 12.5E-03*24/residence_time	! for R=23 (50k)
!	c = 15.5E-03*24/residence_time	! for R=29 (100k)
!	c = 19.3E-03*24/residence_time	! for R=36.3 (200k)
!	c = 22.0E-03*24/residence_time	! for R=41.5 (300k)
!	c = 24.5E-03*24/residence_time	! for R=45.7 (400k)
!	c = 26.5E-03*24/residence_time	! for R=49.2 (500k)
!else
!	a = a_nochemo_12h
!	b = b_nochemo_12h
!	c = c_nochemo_12h
!endif
if (TAGGED_LOG_PATHS) then
	requiredExitPortals = 1
elseif (FIXED_NEXITS) then
	x = NTcells0/1000
	exit_fraction = (a*x**2 + b*x + c)*24/residence_time(CD4)
	requiredExitPortals = exit_fraction*NTcells0**pow + 0.5
elseif (RELAX_INLET_EXIT_PROXIMITY) then	! just for ncells = 100k
	x = ncells/1000
	Fe = 2.50E-03*x**3.595E-01		! power law fit (8 points, 51k - 1.1m cells)
	exit_fraction = Fe*24.0/residence_time(CD4)
	requiredExitPortals = K_Ke*exit_fraction*ncells**pow + 0.5 
!	write(logmsg,'(a,f7.4)') 'exit_fraction: ',exit_fraction 
!	call logger(logmsg)
else
	x = ncells/1000
!	Fe = (a*x**2 + b*x + c)
!	Fe = 3.028E-03*x**3.571E-01		! power law fit (7 points)
	Fe = 2.990E-03*x**3.595E-01		! power law fit (8 points, 51k - 1.1m cells)
	exit_fraction = Fe*24.0/residence_time(CD4)
	requiredExitPortals = K_Ke*exit_fraction*ncells**pow + 0.5 
!	write(logmsg,'(a,f7.4)') 'exit_fraction: ',exit_fraction 
!	call logger(logmsg)
endif
end function

!---------------------------------------------------------------------
! Currently, exit locations are distributed randomly through the blob. 
! Sites within the SOI of an exit are labelled with %exitnum, unless
! the site is also within the SOI of another exit, in which case the
! site is labelled with the exit index of the nearest exit, or in
! the case of a tie the choice is made randomly.
! When a site is an exit portal location (centre), %exitnum = -(the exit index)
! Note:
! When the blob grows, more exits will be added, and when the blob
! shrinks again exits must be removed.  When an exit is removed from
! the blob we set exitlist(iexit)%ID = 0, and all sites within the SOI
! must have %exitnum either set to 0 or, if within the SOI of another
! exit, set to the ID of the nearest exit.
! This is analogous to the treatment of DCs
!
! PlaceExitPortals > AddExitPortal > ChoosePortalSite > PlaceExitPortal
!---------------------------------------------------------------------
subroutine PlaceExitPortals(ok)
integer :: Nex, iexit, site(3)
logical :: ok
logical :: testing = .false.

!if (exit_region /= EXIT_CHEMOTAXIS) then
!    write(*,*) 'Error: placeExits: not EXIT_CHEMOTAXIS'
!    stop
!endif
call logger('PlaceExitPortals')
ok = .true.
if (TAGGED_LOG_PATHS) then
	Nex = 1
    max_exits = 10*Nex
    allocate(exitlist(max_exits))       ! Set the array size to 10* the initial number of exits
    Nexits = 1
    site = Centre
	call PlaceExitPortal(1,site)
	return
endif
if (testing) then
    lastexit = 1
    Nexits = lastexit
    allocate(exitlist(lastexit))
    exitlist(1)%ID = 1
    exitlist(1)%site = (/NX/2,NY/2,NZ/2/)
    return
endif
lastexit = 0
Nexits = 0
if (use_traffic) then
	Nex = requiredExitPortals(NTcells0)
    max_exits = 10*Nex
    write(logmsg,*) 'NTcells0, Nex, max_exits: ',NTcells0,Nex,max_exits
    call logger(logmsg)
    allocate(exitlist(max_exits))       ! Set the array size to 10* the initial number of exits
else
    Nex = 0
    write(logmsg,*) 'No exits'
    call logger(logmsg)
    return
endif
do iexit = 1,Nex
	call AddExitPortal
!	write(logmsg,*) 'Placed exit portal: ',iexit,Nex
!	call logger(logmsg)
enddo
write(logmsg,*) 'Number of exit portals: ',Nex
call logger(logmsg)
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine getExitNum(iexit)
integer :: iexit
integer :: i

iexit = 0
! First look for an unused index
do i = 1,lastexit
	if (exitlist(i)%ID == 0) then
		iexit = i
		exit
	endif
enddo
if (iexit == 0) then
	lastexit = lastexit + 1
	iexit = lastexit
endif
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine AddExitPortal
integer :: site(3)
integer :: iexit

if (dbug) write(nflog,*) 'ChoosePortalSite'
call ChoosePortalSite(site)
if (dbug) write(nflog,*) 'got site: ',site
if (dbug) write(nflog,*) 'getExitNum'
call getExitNum(iexit)
if (dbug) write(nflog,*) 'got iexit: ',iexit
Nexits = Nexits + 1
if (lastexit > max_exits) then
	write(logmsg,*) 'Error: AddExitPortal: too many exits: need to increase max_exits: ',max_exits
	call logger(logmsg)
	stop
endif
if (dbug) write(nflog,*) 'PlaceExitPortal: ',iexit,site
call PlaceExitPortal(iexit,site)
if (dbug) write(nflog,*) 'placed'
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
integer function blobNeighbours(x,y,z)
integer :: x,y,z
integer :: nb, k, xx, yy, zz

nb = 0
blobNeighbours = 0
if (x<1 .or. x>NX) return
if (y<1 .or. y>NY) return
if (z<1 .or. z>NZ) return
if (occupancy(x,y,z)%indx(1) < 0) return
do k = 1,27
	if (k == 14) cycle
	xx = x + jumpvec(1,k)
	yy = y + jumpvec(2,k)
	zz = z + jumpvec(3,k)
	if (occupancy(xx,yy,zz)%indx(1) >= 0) nb = nb + 1
enddo
blobNeighbours = nb
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine getBoundarySite(u,site,ok)
real :: u(3)
integer :: site(3)
logical :: ok
integer :: k, x, y, z, nb
real :: r0, r

ok = .false.		
r0 = radius
do k = -10,10
	r = r0 + k*0.5
	x = x0 + u(1)*r
	y = y0 + u(2)*r
	z = z0 + u(3)*r
	nb = blobNeighbours(x,y,z)
	if (nb == 0) then
		exit
	elseif (nb <= 17) then
		ok = .true.		
		exit
	endif
enddo
site = (/x,y,z/)
end subroutine

!-----------------------------------------------------------------------------
! Exit sites (portals) are either on the blob surface, or within the blob.
! If this is the initial placement (istep = 0) then only blob boundary and
! other exit portals need be considered.
! Otherwise HEV and DC locations need to be avoided as well.
!-----------------------------------------------------------------------------
subroutine ChoosePortalSite(site)
integer :: site(3)
integer :: xex, yex, zex
real(DP) :: R
real :: prox, u(3), theta, rx
integer :: kpar = 0
logical :: ok

!call logger('ChoosePortalSite')
if (SURFACE_PORTALS) then
	! randomly choose direction in 3D, then locate sites near this vector on the blob boundary
	! Need to restrict the sinus interface by removing the follicle interface.
	! Assume that the follicle interface is a cap from (normalized) x = XFOLLICLE (< 1)
	do
		R = par_uni(kpar)
		u(1) = 2*R - 1
		if (u(1) > XFOLLICLE) cycle
		rx = sqrt(1-u(1)*u(1))
		R = par_uni(kpar)
		theta = 2*PI*R
		u(2) = rx*cos(theta)
		u(3) = rx*sin(theta)
		call getBoundarySite(u,site,ok)
		if (.not.ok) then
		    if (dbug) write(nflog,*) 'getBoundarySite not OK'
		    cycle
		endif
		ok = portalOK(site)
		if (.not.ok) then
		    if (dbug) write(nflog,*) 'portalOK not OK'
		    cycle
		endif
		call CheckSite(EXIT_SITE,site,ok)
		if (.not.ok) then
		    if (dbug) write(nflog,*) 'CheckSite not OK'
		    cycle
		endif
!		if (use_DC) then
!			if (tooNearDC(site,exit_DCprox)) then		! exit_DCprox is min distance in sites
!				call logger('tooNearDC')
!				cycle
!			endif
!		endif
!		prox = exit_prox*chemo_N				! chemo_N is chemo_radius in units of sites
!		prox = 0.5*prox
!!		if (USE_PORTAL_EGRESS) then
!!			prox = 0.3*prox
!!		endif
!		if (tooNearExit(site,prox)) then	! too near another exit
!!			call logger('tooNearExit')
!			cycle
!		endif	
		exit
	enddo
else	! blob portals
	do
		R = par_uni(kpar)
		xex = 1 + R*NX
		R = par_uni(kpar)
		yex = 1 + R*NY
		R = par_uni(kpar)
		zex = 1 + R*NZ
		site = (/xex,yex,zex/)
		if (.not.portalOK(site)) cycle
		if (use_DC) then
	!        write(*,*) 'use_DC?'
	!        stop
			if (tooNearDC(site,exit_DCprox)) cycle    ! exit_DCprox is min distance in sites
		endif
		prox = exit_prox*chemo_N				! chemo_N is chemo_radius in units of sites
		if (tooNearExit(site,prox)) cycle       ! too near another exit
		prox = 0.5*exit_prox*chemo_N
		if (tooNearBdry(site,prox)) cycle       ! too near the blob boundary
		exit
	enddo
endif
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine PlaceExitPortal(iexit,site)
integer :: iexit, site(3)
integer :: kexit, x, y, z, x2, y2, z2, ek(3), vk(3), ns
integer :: xex, yex, zex, xmin, xmax, ymin, ymax, zmin, zmax
real(DP) :: R
real :: d2, d2k
integer :: kpar = 0

if (iexit == lastexit + 1) then
	lastexit = iexit
endif
ns = 0
exitlist(iexit)%ID = iexit      ! when an exit site is lost (because the blob retracted) set %ID = 0
exitlist(iexit)%site = site
xex = site(1)
yex = site(2)
zex = site(3)
occupancy(xex,yex,zex)%exitnum = -iexit     ! This site holds an exit 
!write(nfout,'(a17,i6,2x,3i4)') 'Place portal: ',iexit,site

xmin = xex - chemo_N
xmax = xex + chemo_N
ymin = yex - chemo_N
ymax = yex + chemo_N
zmin = zex - chemo_N
zmax = zex + chemo_N
xmin = max(1,xmin)
xmax = min(NX,xmax)
ymin = max(1,ymin)
ymax = min(NY,ymax)
zmin = max(1,zmin)
zmax = min(NZ,zmax)
do x = xmin,xmax
    x2 = (x-xex)*(x-xex)
    do y = ymin,ymax
        y2 = (y-yex)*(y-yex)
        do z = zmin,zmax
	        z2 = (z-zex)*(z-zex)
	        d2 = x2 + y2 + z2
	        if (d2 <= chemo_N*chemo_N) then
		        if (d2 > 0) then   ! don't touch exit site
     ! NOTE!!!  If this site is already marked as within another exit's SOI, we need to
     ! determine which exit is closer, and set %exitnum to this exit's ID
                    kexit = occupancy(x,y,z)%exitnum
                    if (kexit > 0) then
                        ek = exitlist(kexit)%site
                        vk = (/x,y,z/) - ek
                        d2k = dot_product(vk,vk)    ! distance from the other exit (kexit)
!                            write(*,*) 'other exit: ',iexit,kexit,d2,d2k
                        if (d2k < d2) then
!	                            write(*,*) 'closer exit: ',iexit,kexit,d2,d2k
                            cycle
                        endif
                        if (d2k == d2) then     ! toss a coin
                            R = par_uni(kpar)
                            if (R < 0.5) cycle
                        endif
                        occupancy(x,y,z)%exitnum = iexit  ! the site is closer to iexit than to kexit
                        ns = ns + 1
                    elseif (kexit == 0) then
                        occupancy(x,y,z)%exitnum = iexit    ! this site is closest to exit iexit
                        ns = ns + 1
                    endif
                endif
	        endif
        enddo
    enddo
enddo
!write(*,*) 'near sites: ',iexit,ns
end subroutine


!---------------------------------------------------------------------
! Remove exits that are close to others.
!---------------------------------------------------------------------
subroutine RemoveExitPortals(nremex)
integer :: nremex
integer :: Nex, k, iexit, kexit, site1(3), site2(3)
real :: r(3), sum, d2, cmin
real, allocatable :: closeness(:)

Nex = lastexit
allocate(closeness(Nex))
closeness = 0
do iexit = 1,Nex
	sum = 0
    if (exitlist(iexit)%ID == 0) cycle
    site1 = exitlist(iexit)%site
    do kexit = 1,Nex
		if (exitlist(kexit)%ID == 0) cycle
		if (iexit == kexit) cycle
	    site2 = exitlist(kexit)%site
	    r = site1 - site2
	    d2 = norm2(r)
	    sum = sum + 1/d2
	enddo
	closeness(iexit) = sum
enddo

do k = 1,nremex
	cmin = 1.0e10
	do kexit = 1,Nex
		if (closeness(kexit) == 0) cycle
		if (closeness(kexit) < cmin) then
			cmin = closeness(kexit)
			iexit = kexit
		endif
	enddo
	closeness(iexit) = 0
    site1 = exitlist(iexit)%site
    call RemoveExitPortal(site1)
enddo
deallocate(closeness)
!call checkExits("after RemoveExitPortals")
end subroutine

!---------------------------------------------------------------------
! Remove the exit portal at site.
!---------------------------------------------------------------------
subroutine RemoveExitPortal(site)
integer :: site(3)
integer :: iexit,xex,yex,zex,xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,ek(3),vk(3),k,kmin,i
real :: d2k, d2min

xex = site(1)
yex = site(2)
zex = site(3)
iexit = -occupancy(xex,yex,zex)%exitnum
if (iexit <= 0) then
	write(logmsg,*) 'Error: RemoveExitPortal: no exit at: ',site,iexit
	call logger(logmsg)
	do i = 1,lastexit
		site = exitlist(i)%site
		write(logmsg,*) 'exit: ',i,site,occupancy(site(1),site(2),site(3))%exitnum
		call logger(logmsg)
	enddo
	stop
endif
exitlist(iexit)%ID = 0
occupancy(xex,yex,zex)%exitnum = iexit      ! to ensure that the site is processed in the next section
											! effectively it is treated like a site within the portal's SOI
xmin = xex - chemo_N
xmax = xex + chemo_N
ymin = yex - chemo_N
ymax = yex + chemo_N
zmin = zex - chemo_N
zmax = zex + chemo_N
xmin = max(1,xmin)
xmax = min(NX,xmax)
ymin = max(1,ymin)
ymax = min(NY,ymax)
zmin = max(1,zmin)
zmax = min(NZ,zmax)
do x = xmin,xmax
    do y = ymin,ymax
        do z = zmin,zmax
            if (occupancy(x,y,z)%exitnum == iexit) then
                occupancy(x,y,z)%exitnum = 0
                ! at this point we need to check to see if (x,y,z) is within the SOI of another exit
                d2min = 9999
                kmin = 0
                do k = 1,lastexit
                    if (exitlist(k)%ID == 0) cycle
                    ek = exitlist(k)%site
                    vk = (/x,y,z/) - ek
                    d2k = dot_product(vk,vk)    ! distance from the other exit (kexit)
		            if (d2k <= chemo_N*chemo_N) then
		                if (d2k < d2min) then
		                    d2min = d2k
		                    kmin = k
		                endif
		            endif
		        enddo
		        if (kmin > 0) then
		            occupancy(x,y,z)%exitnum = kmin
		        endif
            endif
        enddo
    enddo
enddo
if (iexit == lastexit) then
	lastexit = lastexit - 1
endif
Nexits = Nexits - 1
end subroutine

!---------------------------------------------------------------------
! Move the exit portal to a site closer to the blob centre.
!---------------------------------------------------------------------
subroutine MoveExitPortalInwards(site0)
integer :: site0(3)
integer :: site(3), site1(3), iexit0, iexit1, dx, dy, dz
real :: d0, dc, d2, d2min, r(3)
real, parameter :: dlim = 1.95

iexit0 = -occupancy(site0(1),site0(2),site0(3))%exitnum
! First, remove the exit
call RemoveExitPortal(site0)
!call checkexits("in moveExitPortalInwards")
! Now add the exit back at a closer site.  Choose a new site at
! least a distance of dlim closer to blob centre.  Choose the
! site closest to site0 from among the neighbours within the
! specified distance from the centre.
! NOTE: need to check that there is not a DC or another exit portal nearby.  
! If there is then a completely new portal site will need to be chosen, 
! using ChoosePortalSite()
d0 = cdistance(site0)
d2min = 1.0e10
site1 = 0
do dx = -3,3
	do dy = -3,3
		do dz = -3,3
			site = site0 + (/dx,dy,dz/)
			if (.not.portalOK(site)) cycle
			dc = cdistance(site)
			if (d0-dc >= dlim) then
				r = site - site0
				d2 = norm2(r)
				if (d2 < d2min) then
					d2min = d2
					site1 = site
				endif
			endif
		enddo
	enddo
enddo
if (site1(1) == 0) then
	write(logmsg,*) 'Error: moveExitPortalInwards: no site to move to: ',iexit0,site0
	call logger(logmsg)
	stop
endif
call getExitNum(iexit1)
Nexits = Nexits + 1
call PlaceExitPortal(iexit1,site1)
if (exitlist(iexit1)%site(1) /= site1(1) .or. &
	exitlist(iexit1)%site(2) /= site1(2) .or. &
	exitlist(iexit1)%site(3) /= site1(3)) then
	write(logmsg,*) 'Error: moveExitPortalInwards: bad site: ',iexit1,exitlist(iexit1)%site,site1
	call logger(logmsg)
	stop
endif
end subroutine

!---------------------------------------------------------------------
! The exit portal iexit needs to be moved by one site in the direction 
! closest to that given by the unit vector v(:).
!---------------------------------------------------------------------
subroutine MoveExitPortal(iexit0,v)
integer :: iexit0
real :: v(3)
integer :: iexit1, site0(3), site1(3), site(3), k, kmax
real :: proj, pmax, jump(3)

site0 = exitlist(iexit0)%site
pmax = 0
do k = 1,27
	if (k == 14) cycle
	jump = jumpvec(:,k)
	site = site0 + jump
	if (.not.portalOK(site)) cycle
	proj = dot_product(jump,v)/norm(jump)
	if (proj > pmax) then
		pmax = proj
		kmax = k
	endif
enddo
site1 = site0 + jumpvec(:,kmax)
!write(logmsg,'(a,i4,a,3i4,a,3i4)') 'moveExitPortal: ',iexit0,' from: ',site0,' to: ',site1
!call logger(logmsg)
!write(nfout,'(a17,i6,2x,3i4)') 'Move portal out: ',iexit0,site1
call RemoveExitPortal(site0)
call getExitNum(iexit1)
Nexits = Nexits + 1
call PlaceExitPortal(iexit1,site1)
end subroutine

!---------------------------------------------------------------------
! Check the suitability of site(:) for locating sitetype:
! EXIT_SITE, HEV_SITE, or DC_SITE
!---------------------------------------------------------------------
subroutine CheckSite(sitetype,site,ok)
integer :: sitetype, site(3)
logical :: ok
real :: proxlimit, r

select case(sitetype)
case(EXIT_SITE)
	proxlimit = proximity_limit(EXIT_SITE,EXIT_SITE)
	if (tooNearExit(site,proxlimit)) then	! too near another exit
	    if (dbug) write(nflog,*) 'tooNearExit'
		ok = .false.
		return
	endif
	proxlimit = proximity_limit(EXIT_SITE,DC_SITE)
	if (tooNearDC(site,proxlimit)) then	! too near a DC
	    if (dbug) write(nflog,*) 'tooNearDC'
		ok = .false.
		return
	endif
	proxlimit = proximity_limit(EXIT_SITE,HEV_SITE)
	if (tooNearHEV(site,proxlimit)) then	! too near an HEV
	    if (dbug) write(nflog,*) 'tooNearHEV'
		ok = .false.
		return
	endif
case(HEV_SITE)
	r = cdistance(site)
	if (r/radius < R_limit(HEV_SITE,1) .or. r/radius > R_limit(HEV_SITE,2)) then
		ok = .false.
		return
	endif
   proxlimit = proximity_limit(HEV_SITE,HEV_SITE)
    if (tooNearHEV(site,proxlimit)) then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(HEV_SITE,EXIT_SITE)
    if (tooNearExit(site,proxlimit))  then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(HEV_SITE,BDRY_SITE)
    if (tooNearBdry(site,proxlimit))  then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(HEV_SITE,DC_SITE)
	if (tooNearDC(site,proxlimit)) then	! too near an HEV
		ok = .false.
		return
	endif
case(DC_SITE)
	r = cdistance(site)
	if (r/radius < R_limit(DC_SITE,1) .or. r/radius > R_limit(DC_SITE,2)) then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(DC_SITE,DC_SITE)
    if (tooNearDC(site,proxlimit)) then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(DC_SITE,EXIT_SITE)
    if (tooNearExit(site,proxlimit))  then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(DC_SITE,BDRY_SITE)
    if (tooNearBdry(site,proxlimit))  then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(DC_SITE,HEV_SITE)
	if (tooNearHEV(site,proxlimit)) then	! too near an HEV
		ok = .false.
		return
	endif
end select
ok = .true.
end subroutine

!---------------------------------------------------------------------
! Display exit portal location info
!---------------------------------------------------------------------
subroutine displayExits
integer :: Nex, k, iexit, kexit, site1(3), site2(3), count
real :: r(3), sum, d2

Nex = lastexit
do iexit = 1,Nex
	sum = 0
    if (exitlist(iexit)%ID == 0) cycle
    site1 = exitlist(iexit)%site
    do kexit = 1,Nex
		if (exitlist(kexit)%ID == 0) cycle
		if (iexit == kexit) cycle
	    site2 = exitlist(kexit)%site
	    r = site1 - site2
	    d2 = norm2(r)
	    sum = sum + 1/d2
	enddo
	count = neighbourhoodCount(site1)
	write(logmsg,'(5i4,3f10.3)') iexit,site1,count,cdistance(site1),(PI/180)*atan2(site1(2)-Centre(2),site1(3)-Centre(3)),1000*sum
	call logger(logmsg)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! need to ensure that exits do not get too close together
!-----------------------------------------------------------------------------------------
subroutine CheckExitSpacing(msg)
character*(*) :: msg
integer :: iexit1, iexit2, site1(3), site2(3), im1=0, im2=0
real :: r(3), d, dmin, f(3), df(3)
logical :: ok = .true.
real :: dlim

dlim = exit_prox*chemo_radius
dmin = 1.0e10
do iexit1 = 1,lastexit
	if (exitlist(iexit1)%ID == 0) cycle
	site1 = exitlist(iexit1)%site
	f = 0
	do iexit2 = 1,lastexit
		if (iexit1 == iexit2) cycle
		if (exitlist(iexit2)%ID == 0) cycle
		site2 = exitlist(iexit2)%site
		r = site1 - site2
		d = norm(r)
		if (d < dmin) then
			dmin = d
			im1 = iexit1
			im2 = iexit2
		endif
		if (d < dlim) then
			write(logmsg,'(2i4,f6.2)') iexit1,iexit2,d
			call logger(logmsg)
			df = r/d		! unit vector in direction of force
			df = df/d		! scale by inverse distance
			f = f + df
		endif
	enddo
	if (f(1) /= 0 .or. f(2) /= 0 .or. f(3) /= 0) then
		f = f/norm(f)
		call MoveExitPortal(iexit1,f)
	endif
enddo
write(logmsg,*) msg,': Min exit spacing: ',im1,im2,dmin
call logger(logmsg)
end subroutine


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine init_cytokine
integer :: x, y, z, icyt, iseq, tag, offset
integer :: np(MAX_CYT)
real :: dt, dx, mols_pM

Ncytokines = N_CYT
cyt_seq = 0
do iseq = 1,Ncytokines
    tag = cyt_tag(iseq)
    cyt_seq(tag) = iseq
enddo
do iseq = 1,Ncytokines
    tag = cyt_tag(iseq)
    select case(tag)
    case (IL2_TAG)
        np(iseq) = IL2_NP
    case (IL4_TAG)
        np(iseq) = IL4_NP
    case (IL7_TAG)
        np(iseq) = IL7_NP
    case (IL9_TAG)
        np(iseq) = IL9_NP
    case (IL15_TAG)
        np(iseq) = IL15_NP
    case (IL21_TAG)
        np(iseq) = IL21_NP
    end select
enddo
offset = 0
do iseq = 1,Ncytokines
    NP_offset(iseq) = offset
    offset = offset + np(iseq)
    NP_offset(iseq+1) = offset
enddo

allocate(cyt_constit(Ncytokines))
if (calibrate_diffusion) then
    cyt_constit = 0
else
    cyt_constit = 0
    iseq = cyt_seq(IL2_TAG)
    cyt_constit(iseq) = IL2_CONSTIT_RATE    ! units molecules/um^3/min
endif

! For testing
K_diff = 2.0E-12                ! m^2/s
dt = DELTA_T*60.0/NDIFFSTEPS    ! s
dx = DELTA_X*1.0E-6             ! m
do icyt = 1,Ncytokines
    delta_diff(icyt) = K_diff(icyt)*dt/(dx*dx)
    if (6*delta_diff(icyt) >= 1.0) then     ! was 6*
        write(*,*) 'Bad delta_diff: need to increase NDIFFSTEPS: ',delta_diff(icyt)
!        stop
    endif
    if (calibrate_diffusion) then
        cyt_init(icyt) = 1.0
    else
        cyt_init(icyt) = cyt0(icyt)
    endif
enddo

allocate(cyt(NX,NY,NZ,Ncytokines))
allocate(cyt_mols(Ncytokines))
allocate(dcyt_mols(Ncytokines))
dcyt_mols = 0
!write(*,*) 'Initial cytokine conc., mols'
do icyt = 1,Ncytokines
    if (calibrate_diffusion) then
        do x = 1,NX/2
            cyt(x,:,:,icyt) = cyt_init(icyt)
        enddo
        do x = NX/2,NX
            cyt(x,:,:,icyt) = 0
        enddo
    else
        do z = 1,NZ
            do y = 1,NY
                do x = 1,NX
                    if (occupancy(x,y,z)%indx(1) >= 0) then
                        cyt(x,y,z,icyt) = cyt_init(icyt)
                    else
                        cyt(x,y,z,icyt) = 0
                    endif
                enddo
            enddo
        enddo
        ! need to make # of mols correspond to cyt_init
        mols_pM = L_um3*M_pM/(NTcells*Vc*Navo)
        cyt_mols(icyt) = cyt_init(icyt)/mols_pM   ! -> total number of molecules
        write(*,'(a,f8.3,f12.0)') cyt_name(cyt_tag(icyt)),cyt_init(icyt),cyt_mols(icyt)
    endif
enddo

end subroutine

!-----------------------------------------------------------------------------------------
! If "mean field" concentrations are used for cytokines (i.e. not use_diffusion)
! the net number of molecules released into the blob in one time step (in updater()) is
! accumulated in dcyt_mols(:)
! The overall total count of molecules in cyt_mols(:) is updated every time step and
! broadcast to the nodes.
!-----------------------------------------------------------------------------------------
subroutine molsynch
integer :: k, kIL2  !, source, dest, ierr, status(MPI_STATUS_SIZE)
real :: mols_pM
!real, allocatable :: dmols(:)

!allocate(dmols(Ncytokines))
!if (me == 0) then
    kIL2 = cyt_seq(IL2_TAG)
    cyt_mols(kIL2) = cyt_mols(kIL2)*(1 - IL2_DECAY_RATE*DELTA_T)                ! decay of IL2
    mols_pM = L_um3*M_pM/(NTcells*Vc*Navo)
!    write(*,'(a,4f8.1)') 'molsynch: ',cyt_mols(kIL2),dcyt_mols(kIL2), &
!        cyt_constit(kIL2)*DELTA_T*NTcells*Vc,mols_pM*cyt_mols(kIL2)

    cyt_mols = cyt_mols + dcyt_mols + cyt_constit*DELTA_T*NTcells*Vc  ! constitutive production of cytokines
    dcyt_mols = 0
!    if (Mnodes > 1) then
!        do source = 1,Mnodes-1
!            call MPI_RECV(dmols, Ncytokines, MPI_REAL, source, MOL1_TAG, MPI_COMM_WORLD, status, ierr)
!            cyt_mols = cyt_mols + dmols
!        enddo
!        do dest = 1,Mnodes-1
!            call MPI_SEND(cyt_mols, Ncytokines, MPI_REAL, dest, MOL2_TAG, MPI_COMM_WORLD, ierr)
!        enddo
!    endif
!elseif (Mnodes > 1) then
!    dest = 0
!    call MPI_SEND(dcyt_mols, Ncytokines, MPI_REAL, dest, MOL1_TAG, MPI_COMM_WORLD, ierr)
!    dcyt_mols = 0
!    source = 0
!    call MPI_RECV(cyt_mols, Ncytokines, MPI_REAL, source, MOL2_TAG, MPI_COMM_WORLD, status, ierr)
!endif
!deallocate(dmols)
do k = 1,Ncytokines
    if (cyt_mols(k) < 0) then
        write(*,*) 'molsynch: cyt_mols < 0: ',k,cyt_mols(k)
        stop
    endif
enddo
end subroutine

if (use_DC) then
    allocate(DClist(MAX_DC))
    allocate(DCdeadlist(MAX_DC))
    allocate(DCvisits(0:MAX_DC+1,2,2,2))
    allocate(DCtotvisits(0:1000,2,2,2))		! a guess
    z1 = -DC_RADIUS
    z2 = DC_RADIUS
    if (IN_VITRO) then
		z1 = 0
		z2 = NZ - 1
	endif
    k = 0
    do x = -DC_RADIUS,DC_RADIUS
	    do y = -DC_RADIUS,DC_RADIUS
		    do z = z1, z2
		        rr = (/x,y,z/)
		        d = norm(rr)
			    if (d <= DC_RADIUS) then
				    k = k+1
			    endif
		    enddo
	    enddo
    enddo
    k = k - NDCsites
    nbindmax = min(k,MAX_TC_BIND)
    nbind1 = ABIND1*nbindmax
    nbind2 = ABIND2*nbindmax
!    write(*,*) 'Space for T cells/DC: ',k
!    write(*,*) 'nbind1,nbind2,nbindmax: ',nbind1,nbind2,nbindmax
    if (.not.IN_VITRO) then
	    DCoffset(:,1) = (/ 0,0,0 /)
	    if (NDCsites > 1 .and. NDCsites <= 7) then
	        DCoffset(:,2:NDCsites) = neumann(:,1:NDCsites-1)
	    endif
	else
		! need to create DCoffset() for IN_VITRO case
		if (NDCsites == 9) then
			k = 0
			do x = -1,1
				do y = -1,1
					k = k+1
					DCoffset(:,k) = (/ x, y, 0 /)
				enddo
			enddo
		elseif (NDCsites >= 5 .and. NDCsites <=8) then
			DCoffset(:,1) = (/ 0, 0, 0 /)
			DCoffset(:,2) = (/-1, 0, 0 /)
			DCoffset(:,3) = (/ 1, 0, 0 /)
			DCoffset(:,4) = (/ 0,-1, 0 /)
			DCoffset(:,5) = (/ 0, 1, 0 /)
			if (NDCsites > 5) DCoffset(:,6) = (/ 0, 0, 1 /)
			if (NDCsites > 6) DCoffset(:,7) = (/-1, 0, 1 /)
			if (NDCsites > 7) DCoffset(:,8) = (/ 1, 0, 1 /)
		endif
	endif
    ndeadDC = 0
	DClist(:)%nbound = 0		! moved from initial_binding
	DClist(:)%ncogbound = 0
endif

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine save_pos
!integer :: k, kcell, site(3)
character*(8) :: msg
integer :: error

!write(nfout,*) 'STEP: ',istep,lastcogID
!do k = 1,lastcogID
!    kcell = cognate_list(k)
!    if (kcell > 0) then
!        site = cellist(kcell)%site
!        gen = get_generation(cellist(kcell)%cptr)
!    else
!        site = 0
!        gen = 0
!    endif
!    write(nfout,'(i8,3i4)') kcell,site,gen
!enddo
if (save_pos_cmgui) then
	call save_exnode
elseif (use_TCP) then
	call save_cell_positions
	msg = 'VTK'
	clear_to_send = .false.
    call winsock_send(awp_1,msg,len_trim(msg),error)
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine initial_binding
integer :: kcell, site(3), dc(0:3), k, idc
real(DP) :: R
real :: bindtime
type(cell_type) :: tcell
integer :: kpar=0
logical :: cognate

do kcell = 1,nlist
    tcell = cellist(kcell)
    cognate = (associated(tcell%cptr))
    site = tcell%site
!    write(*,*) 'initial_binding: ',kcell,nlist,' ',cognate,' ',site
    dc = occupancy(site(1),site(2),site(3))%DC
    if (dc(0) > 0) then
        do k = 1,dc(0)
            idc = dc(k)
            if (.not.DClist(idc)%capable) cycle
            if (cognate .and. DClist(idc)%ncogbound == MAX_COG_BIND) cycle
            if (bindDC(idc,kpar)) then
!                call random_number(R)
                R = par_uni(kpar)
                bindtime = 1 + R*5.     ! mean = 3.5 min
!                call random_number(R)
                R = par_uni(kpar)
                cellist(kcell)%DCbound(1) = idc
                cellist(kcell)%unbindtime(1) = bindtime*R
                DClist(idc)%nbound = DClist(idc)%nbound + 1
                if (cognate) then
                    DClist(idc)%ncogbound = DClist(idc)%ncogbound + 1
                    call AddCogbound(idc,kcell)
                endif
                exit
            endif
        enddo
    endif
enddo
!write(*,*) 'DC occupancy:'
!write(*,'(8i6)') DClist(1:NDC)%nbound
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_DCproximity
integer :: x,y,z,k,cnt(0:DCDIM-1)
integer :: dc(DCDIM-1)

cnt = 0
do x = 1,NX
    do y = 1,NY
        do z = 1,NZ
            k = occupancy(x,y,z)%DC(0)
            if (k >= 0) then
                cnt(k) = cnt(k) + 1
                dc = occupancy(x,y,z)%DC(1:DCDIM-1)
                if (k >= 2) then
                    call randomizeDC(dc,k)
                    occupancy(x,y,z)%DC(1:DCDIM-1) = dc
                endif
            endif
        enddo
    enddo
enddo
!write(*,*) 'DC proximity counts: ',cnt

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine randomizeDC(dc,k)
integer :: k, dc(:)
integer :: i,p(DCDIM-1),tmp(DCDIM-1)
integer :: kpar=0

do i = 1,k
    p(i) = i
enddo
call permute(p,k,kpar)
tmp(1:k) = dc(1:k)
do i = 1,k
    dc(i) = tmp(p(i))
enddo
!write(*,*) tmp(1:k),dc(1:k)
end subroutine

!-----------------------------------------------------------------------------------------
! Is site near a DC?
! The criterion for a DC site might be different from an exit site.
! prox = DC_DCprox*DC_RADIUS for DC - DC
!-----------------------------------------------------------------------------------------
logical function tooNearDC(site,prox)
integer :: site(3)
real :: prox
integer :: idc
real :: r(3), d

if (NDCalive == 0) then
    tooNearDC = .false.
    return
endif
do idc = 1,NDC
    if (.not.DClist(idc)%alive) cycle
    r = site - DClist(idc)%site
    d = norm(r)     ! units sites
    if (d < prox) then
        tooNearDC = .true.
        return
    endif
enddo
tooNearDC = .false.
end function

!-----------------------------------------------------------------------------------------
! prox is the minimum distance from the boundary (sites)
!-----------------------------------------------------------------------------------------
logical function tooNearBdry(site,prox)
integer :: site(3)
real :: prox
real :: d, dmin

!dmin = proxfac*DC_RADIUS
if (use_blob) then
    d = cdistance(site)
    if (Radius - d < prox) then
        tooNearBdry = .true.
        return
    endif
!    write(*,*) 'tooNearBdry: ',Radius,d,dmin
else
    if (site(1) < dmin .or. (NX - site(1)) < dmin) then
        tooNearBdry = .true.
        return
    endif
    if (site(2) < dmin .or. (NY - site(2)) < dmin) then
        tooNearBdry = .true.
        return
    endif
    if (site(3) < dmin .or. (NZ - site(3)) < dmin) then
        tooNearBdry = .true.
        return
    endif
endif
tooNearBdry = .false.
end function

!-----------------------------------------------------------------------------------------
! Is site near an exit?  prox is the minimum separation (sites)
!-----------------------------------------------------------------------------------------
logical function tooNearExit(site,prox)
integer :: site(3)
real :: prox
integer :: iexit
real :: r(3), d

if (Nexits == 0) then
    tooNearExit = .false.
    return
endif
if (.not.allocated(exitlist)) then
	call logger('exitlist not allocated')
	stop
endif
do iexit = 1,lastexit
    if (exitlist(iexit)%ID == 0) cycle  ! exit no longer exists
    r = site - exitlist(iexit)%site
    d = norm(r)
    if (d < prox) then
        tooNearExit = .true.
        return
    endif
enddo
tooNearExit = .false.
end function

!-----------------------------------------------------------------------------------------
! Is site near an HEV?  prox is the minimum separation (sites)
!-----------------------------------------------------------------------------------------
logical function tooNearHEV(site,prox)
integer :: site(3)
real :: prox
integer :: ihev
real :: r(3), d

if (NHEV == 0) then
    tooNearHEV = .false.
    return
endif
do ihev = 1,NHEV
    r = site - HEVlist(ihev)%site
    d = norm(r)
    if (d < prox) then
        tooNearHEV = .true.
        return
    endif
enddo
tooNearHEV = .false.
end function

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
logical function bindDC(idc,kpar)
integer :: idc,kpar
integer :: n
real(DP) :: p, R

n = DClist(idc)%nbound
if (n >= nbind2) then
    bindDC = .false.
elseif (n <= nbind1) then
    bindDC = .true.
else
    p = real(nbind2 - n)/(nbind2 - nbind1)
!    call random_number(R)
    R = par_uni(kpar)
    if (R < p) then
        bindDC = .true.
    else
        bindDC = .false.
    endif
endif
end function

!-----------------------------------------------------------------------------------------
! Add a cognate T cell index to the list of bound cognate cells for DC idc
!-----------------------------------------------------------------------------------------
subroutine AddCogbound(idc,kcell)
integer :: idc, kcell
integer :: i

do i = 1,MAX_COG_BIND
	if (DClist(idc)%cogbound(i) == 0) then
		DClist(idc)%cogbound(i) = kcell
		return
	endif
enddo
write(logmsg,*) 'Error: AddCogbound: ',idc,kcell
call logger(logmsg)
stop
end subroutine

!-----------------------------------------------------------------------------------------
! Remove a cognate T cell index from the list of bound cognate cells for DC idc
!-----------------------------------------------------------------------------------------
subroutine RemoveCogbound(idc,kcell)
integer :: idc, kcell
integer :: i

do i = 1,MAX_COG_BIND
	if (DClist(idc)%cogbound(i) == kcell) then
		DClist(idc)%cogbound(i) = 0
		return
	endif
enddo
write(logmsg,*) 'Error: RemoveCogbound: ',idc,kcell
call logger(logmsg)
stop
end subroutine

!-----------------------------------------------------------------------------------------
! Pre-programmed inflammation signal is flat at Inflammation_100k (for NTcells0 = 100k)
! for Inflammation_days1 then ends at Inflammation_days2
!-----------------------------------------------------------------------------------------
real function get_inflammation()
real :: tnow, plateau

tnow = istep*DELTA_T    ! mins
tnow = tnow/(24*60)     ! days
plateau = Inflammation_level*NTcells0
if (tnow < Inflammation_days1) then
    get_inflammation = plateau
elseif (tnow > Inflammation_days2) then
    get_inflammation = 0
else
    get_inflammation = plateau*(Inflammation_days2 - tnow)/(Inflammation_days2 - Inflammation_days1)
endif
end function

!-----------------------------------------------------------------------------------------
! Computes the total DC activity summed over all DC that are capable.
! For this to work there must be a model for:
! (a) influx of cognate DC during the response
! (b) the decline of DC stimulating ability (%density) either as a function of time,
!     or in proportion to TCR stimulation delivered.  This is provided by update_DCstate()
!-----------------------------------------------------------------------------------------
real function get_DCactivity()
integer :: idc, ncap
real :: a, tnow, td

if (test_vascular) then
    ! A simple variation in DC activity to test vascular development
    tnow = istep*DELTA_T    ! mins
    td = tnow/(60*24)       ! days
    if (td < 2) then
        a = 1
    elseif (td < 5) then
        a = 1 - (td-2)/3
    else
        a = 0
    endif
    get_DCactivity = a*NDC
    return
endif

ncap = 0
a = 0
do idc = 1,NDC
    if (DClist(idc)%capable) then
        ncap = ncap + 1
        a = a + DClist(idc)%density
    endif
enddo
get_DCactivity = a
end function

subroutine generate_traffic(inflow0)
real :: inflow0
real :: act, expansion, actfactor, tnow
real :: inflow, outflow

if (suppress_egress) then
    call blocked_egress_inflow(inflow)
    InflowTotal = inflow
    OutflowTotal = 0
    return
endif
    
if (traffic_mode == TRAFFIC_MODE_1) then    ! naive
    write(*,*) 'generate_traffic: do not use TRAFFIC_MODE_1'
    stop
    act = get_DCactivity()
    act = act*100/NTcells ! to make act into a concentration, and Bflow approx 1.0
    expansion = real(NTcells)/NTcells0 - 1
    actfactor = 1 + (Kflow3-1)*act**Nflow/(act**Nflow + Bflow**Nflow)
    write(*,'(a,4f6.2)') 'act,actfactor,expansion,radius: ',act,actfactor,expansion,Radius
    inflow = inflow0*actfactor*(1 + Kflow1*expansion)
    inflow = max(inflow,inflow0)
    outflow = inflow0*(1 + Kflow2*expansion)
    outflow = max(outflow,inflow0)
else        !traffic_mode == TRAFFIC_MODE_2 or TRAFFIC_MODE_3
	! Note: if inflammation signal = 0 the vascularity (and inflow) should be constant
    tnow = istep*DELTA_T
    inflow = inflow0*Vascularity   ! level of vascularity (1 = steady-state)
    outflow = NTcells*DELTA_T/(ave_RESIDENCE_TIME*60)
endif

if (L_selectin) then
	if (mod(istep,1000) == 0) then
		call logger('Using L_selectin!')
	endif
    OutflowTotal = NTcells*DELTA_T/(ave_RESIDENCE_TIME*60)
    InflowTotal = 0
    return
endif

!DCactivity = act
InflowTotal = inflow
! This is a kludge to induce a return to steady-state maintenance when NTcells drops
! back to "close enough" to the steady-state value.
! I think this is no longer needed.
! Try removing it
!if (use_exit_chemotaxis .and. NTcells < 0.99*NTcells0) then
!    OutflowTotal = InflowTotal      ! => steady-state with chemotaxis
!    steadystate = .true.
!else
    OutflowTotal = outflow
!endif
!if (mod(istep,240) == 0) then
!	write(logmsg,*) 'generate_traffic: inflow: ',inflow0,Vascularity,InflowTotal
!	call logger(logmsg)
!endif
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine check_bindings(str)
character*(*) :: str
integer :: kcell, k, bnd(2), idc
integer :: n26, cells(1000)
integer, allocatable :: nbnd(:)

!call MPI_BARRIER ( MPI_COMM_WORLD, ierr )
write(*,*) 'check_bindings: ',str,' ',istep
allocate(nbnd(NDC))
nbnd = 0
n26 = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    bnd = cellist(kcell)%DCbound
    if (bnd(1) == 0 .and. bnd(2) /= 0) then
        write(*,*) 'Bad DCbound order: ',kcell,bnd
        stop
    endif
    do k = 1,MAX_DC_BIND
        idc = bnd(k)
        if (idc > 0) then
            nbnd(idc) = nbnd(idc) + 1
            if (idc == 26) then
                n26 = n26 + 1
                cells(n26) = kcell
            endif
        endif
    enddo
enddo
deallocate(nbnd)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine big_check_bindings
integer :: kcell, k, bnd(2), idc
integer, allocatable :: nbnd(:)

allocate(nbnd(NDC))
nbnd = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    bnd = cellist(kcell)%DCbound
    do k = 1,MAX_DC_BIND
        idc = bnd(k)
        if (idc > 0) then
            nbnd(idc) = nbnd(idc) + 1
        endif
    enddo
enddo
do idc = 1,NDC
    if (.not.DClist(idc)%alive) cycle
    if (nbnd(idc) /= DClist(idc)%nbound) then
        write(*,*) 'big_check_bindings: inconsistent nbnd: ',idc,nbnd(idc),DClist(idc)%nbound
        stop
    endif
enddo
deallocate(nbnd)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine read_fixed_params(ok)
logical :: ok
logical :: ext
logical :: use_DCvisits_params
character*(128) :: DCparamfile

inquire(file=fixedfile,exist=ext)
if (.not.ext) then
	call logger("Fixed parameter input file not found")
	ok = .false.
	return
endif
open(nfcell,file=fixedfile,status='old',err=99)
ok = .true.
! T cell parameters
!read(nfcell,*) use_traffic			!= .false.
read(nfcell,*) TCR_splitting		!= .false.  ! enable sharing of integrated TCR signal between progeny cells
read(nfcell,*) transient_stagetime	!= 1.0
read(nfcell,*) clusters_stagetime	!= 13.0
read(nfcell,*) transient_bindtime	!= 10.0
read(nfcell,*) clusters_bindtime	!= 180.0
read(nfcell,*) swarms_bindtime		!= 20.0
read(nfcell,*) TC_life_median1		!= 196				! median lifetime of naive T cells
read(nfcell,*) TC_life_median2		!= 196				! median lifetime of activated T cells
read(nfcell,*) TC_life_shape		!= 1.4				! shape parameter for lifetime of T cells
read(nfcell,*) NTC_LN				!= 3.0e07			! number of T cells in a LN
read(nfcell,*) NTC_BODY				!= 1.6e09			! number of circulating T cells in the whole body
read(nfcell,*) NLN_RESPONSE								! number of LNs in the response
!write(*,*) 'NTC_LN, NTC_BODY: ',NTC_LN, NTC_BODY


!CD25/IL2
read(nfcell,*) optionA                      ! 1,2
read(nfcell,*) optionB                      ! 1,2,3
read(nfcell,*) optionC                      ! 1,2
read(nfcell,*) IL2_PRODUCTION_TIME			! duration of IL2/CD25 production
read(nfcell,*) CD25_DIVISION_THRESHOLD	    ! CD25 store level needed for division of activated cell (optionA = 1)
read(nfcell,*) CD25_SURVIVAL_THRESHOLD		! CD25 store level needed for survival of activated cell (optionC = 1)

read(nfcell,*) TC_RADIUS					! radius of T cell (um)
read(nfcell,*) TC_STIM_WEIGHT				! contribution of stimulation to act level
read(nfcell,*) TC_MAX_GEN				    ! maximum T cell generation
read(nfcell,*) CD8_DIVTIME_FACTOR			! CD8 divide time as multiple of CD4 divide time
read(nfcell,*) DC_MULTIBIND_PROB			! reducing factor to bind prob for each current DC binding
read(nfcell,*) MAX_DC_BIND					! number of DC that a cell can bind to simultaneously
read(nfcell,*) divide_dist1%class
read(nfcell,*) divide_dist2%class
read(nfcell,*) GAMMA						! controls crowding

read(nfcell,*) DC_ACTIV_TAPER				! time (hours) over which DC activity decays to zero
read(nfcell,*) DC_DENS_BY_STIM              ! rate of reduction of density by TCR stimulation
read(nfcell,*) DC_FACTOR                    ! multiplying factor for DC number (initial and influx)

read(nfcell,*) efactor                      ! If constant_efactor = true, this is the factor for the p correction
read(nfcell,*) VEGF_MODEL                   ! 1 = VEGF signal from inflammation, 2 = VEGF signal from DCactivity

read(nfcell,*) fix_avidity                  ! true if avidity takes discrete values, false if a distribution is used
read(nfcell,*) avidity_logscale             ! true if actual avidity = 10^(avidity_min + i*avidity_step)
read(nfcell,*) avidity_nlevels              ! If fix_avidity, number of discrete avidity values
read(nfcell,*) avidity_min                  ! minimum value
read(nfcell,*) avidity_step                 ! step between equi-spaced values

! Parameters added to control HEV and DC placement
read(nfcell,*) R_HEV_min
read(nfcell,*) R_HEV_max
read(nfcell,*) R_DC_min
read(nfcell,*) R_DC_max
! Parameters added for DCvisits simulations, now in special case input file
!read(nfcell,'(L)') use_DCvisits_params
!if (use_DCvisits_params) then
!	track_DCvisits = .true.
!	read(nfcell,'(a)') DCparamfile
!	close(nfcell)
!	open(nfcell,file=DCparamfile,status='old')
!	read(nfcell,'(L)') TAGGED_DC_CHEMOTAXIS		! DC visits are logged for tagged T cells only
!	read(nfcell,*) istep_DCvisits	! 24*60*4	! delay before counting visits (when use_traffic) (steps)
!	read(nfcell,*) t_log_DCvisits	! 2*24*60	! time to log DC visits of retained cells (mins)
!	read(nfcell,'(L)') USE_DC_COGNATE 		! if only a fraction (0.5) of DCs bear cognate antigen, secrete chemokine
!	read(nfcell,*) DC_COGNATE_FRACTION		! fraction of DCs that bear cognate antigen
!	read(nfcell,*) RETAIN_TAGLIMIT_FRACTION ! 0.5 (when use_traffic)
!	read(nfcell,'(L)') DC_CHEMO_NOTRAFFIC   ! cells are tagged initially, no use_traffic
!	read(nfcell,*) DC_CHEMO_FRACTION 		! fraction of cells tagged when DC_CHEMO_NOTRAFFIC
!	read(nfcell,*) HI_CHEMO_FRACTION 		! fraction of tagged cells with full chemotactic susceptibility
!	read(nfcell,*) HI_CHEMO 				! chemotactic susceptibility of the HI_CHEMO_FRACTION cells
!	read(nfcell,*) LO_CHEMO 				! chemotactic susceptibility of the other tagged cells
!endif
close(nfcell)
return
99	continue
ok = .false.
end subroutine

!----------------------------------------------------------------------------------------
! Special cases:
! 1 
! To quantify DC visits for an artificial situation of no trafficking.
! The T cell population does not change, a small fraction of T cells are tagged.  These 
! tagged cells are assigned either HI or LO chemotactic susceptibility.  There is no TCR
! stimulation, no T cell activation.  DC visits are counted.
! 2
! To quantify DC visits with trafficking.
! 3
! DC chemokine secretion is linked to contact with cognate CD4  cells.
! 4
! Discrete levels of TCR avidity and/or discrete DC antigens
!----------------------------------------------------------------------------------------
subroutine read_specialcase(icase,casefile,ok)
integer :: icase
character*(*) :: casefile
logical :: ok

!if (use_DCvisits_params) then
ok = .true.
write(logmsg,*) 'read_specialcase: icase, casefile: ',icase,'  ',casefile
call logger(logmsg)
if (icase == 1) then
	if (receptor(CCR1)%strength == 0) then
		call logger('Error: need to specify CCR1 receptor strength for special case: TAGGED_DC_CHEMOTAXIS')
		ok = .false.
		return
	endif
	use_traffic = .false.
	open(nfcell,file=casefile,status='old')
	read(nfcell,'(L)') USE_DC_CHEMOTAXIS		! DC chemotaxis
	read(nfcell,'(L)') TAGGED_DC_CHEMOTAXIS		! DC visits are logged for tagged T cells only
	read(nfcell,*) istep_DCvisits	! 24*60*4	! delay before counting visits (when use_traffic) (steps)
	read(nfcell,*) t_log_DCvisits	! 2*24*60	! time to log DC visits of retained cells (mins)
	read(nfcell,'(L)') USE_DC_COGNATE 			! if only a fraction (0.5) of DCs bear cognate antigen, secrete chemokine
	read(nfcell,*) DC_COGNATE_FRACTION			! fraction of DCs that bear cognate antigen
	read(nfcell,'(L)') RETAIN_TAGGED_CELLS		! applies with use_traffic
	read(nfcell,'(L)') DC_CHEMO_NOTRAFFIC		! cells are tagged initially, no use_traffic
	read(nfcell,*) DC_CHEMO_FRACTION 			! fraction of cells tagged when DC_CHEMO_NOTRAFFIC
	read(nfcell,*) HI_CHEMO_FRACTION 			! fraction of tagged cells with full chemotactic susceptibility
	read(nfcell,*) HI_CHEMO 					! chemotactic susceptibility of the HI_CHEMO_FRACTION cells
	read(nfcell,*) LO_CHEMO 					! chemotactic susceptibility of the other tagged cells
	close(nfcell)
elseif (icase == 2) then
	track_DCvisits = .true.
	use_traffic = .true.
	open(nfcell,file=casefile,status='old')
	read(nfcell,'(L)') USE_DC_CHEMOTAXIS		! DC chemotaxis
	read(nfcell,'(L)') TAGGED_DC_CHEMOTAXIS		! DC visits are logged for tagged T cells only
	read(nfcell,*) istep_DCvisits	! 24*60*4	! delay before counting visits (when use_traffic) (steps)
	read(nfcell,*) t_log_DCvisits	! 2*24*60	! time to log DC visits of retained cells (mins)
	read(nfcell,'(L)') USE_DC_COGNATE 		! if only a fraction (0.5) of DCs bear cognate antigen, secrete chemokine
	read(nfcell,*) DC_COGNATE_FRACTION		! fraction of DCs that bear cognate antigen
	read(nfcell,*) RETAIN_TAGGED_CELLS		! applies with use_traffic
	read(nfcell,*) RETAIN_TAGLIMIT_FRACTION ! 0.5 (when use_traffic)
	read(nfcell,'(L)') DC_CHEMO_NOTRAFFIC   ! cells are tagged initially, no use_traffic
	read(nfcell,*) TAGGED_CHEMO_FRACTION 		! fraction of cells tagged when DC_CHEMO_NOTRAFFIC
	read(nfcell,*) HI_CHEMO_FRACTION 		! fraction of tagged cells with full chemotactic susceptibility
	read(nfcell,*) HI_CHEMO 				! chemotactic susceptibility of the HI_CHEMO_FRACTION cells
	read(nfcell,*) LO_CHEMO 				! chemotactic susceptibility of the other tagged cells
	close(nfcell)

endif
end subroutine

!-----------------------------------------------------------------------------------------
! Read the DC injection schedule file, determine the total number of DCs entering the LN
! in the experiment.  The total for the simulation run is determined by scaling the
! experiment total by the ratio of the initial T cell population numbers for the simulation
! and the experiment.  The total NDCrun must then allocated to the hourly intervals in
! the appropriate proportions.
!-----------------------------------------------------------------------------------------
subroutine setup_DCinjected
integer :: nthours, NTCexpt, NDCexpt, ntimes, i, ihour, n, NDCrun, nsum, ndeficit
real :: ratio
integer, allocatable :: DCexpt(:)

open(nfDCinfo,file=DCinjectionfile,status='old')
nthours = days*24 + 1
allocate(DCinjected(nthours))
allocate(DCexpt(nthours))
DCinjected = 0
DCexpt = 0
read(nfDCinfo,*) NTCexpt
read(nfDCinfo,*) T_DC_INJECTION
read(nfDCinfo,*) ntimes
NDCexpt = 0
do i = 1,ntimes
	read(nfDCinfo,*) ihour, n
	if (ihour > nthours) exit
	DCexpt(ihour) = n
	NDCexpt = NDCexpt + n
enddo
close(nfDCinfo)
ratio = real(NTcells0)/NTCexpt
NDCrun = ratio*NDCexpt + 0.5
write(logmsg,*) 'NTcells0, NTCexpt, ratio, NDCrun: ',NTcells0,NTCexpt,ratio,NDCrun
call logger(logmsg)
nsum = 0
do ihour = 1,nthours
	DCinjected(ihour) = ratio*DCexpt(ihour) + 0.5
	nsum = nsum + DCinjected(ihour)
enddo
ndeficit = NDCrun - nsum
do ihour = 1,nthours
	if (ndeficit == 0) exit
	if (DCinjected(ihour) > 0) then
		if (ndeficit > 0) then
			DCinjected(ihour) = DCinjected(ihour) + 1
			ndeficit = ndeficit - 1
		elseif (ndeficit < 0) then
			DCinjected(ihour) = DCinjected(ihour) - 1
			ndeficit = ndeficit + 1
		endif
	endif
enddo
deallocate(DCexpt)
do ihour = 1,nthours
	write(logmsg,*) ihour,DCinjected(ihour)
	call logger(logmsg)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine summary
integer :: idc

if (use_DC) then
    do idc = 1,NDC
        if (DClist(idc)%alive) then
            write(*,*) idc,DClist(idc)%nbound,real(DClist(idc)%nbound)/nbindmax
        else
            write(*,*) idc,' dead'
        endif
    enddo
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Using the complete list of cells, cellist(), extract info about the current state of the
! paracortex.  This info must be supplemented by counts of cells that have died and cells that
! have returned to the circulation.
!-----------------------------------------------------------------------------------------
subroutine show_snapshot(ok)
logical :: ok
integer :: kcell, ctype, stype, ncog, noncog, ntot, stage, region, i, iseq, error
integer :: gen, ngens, neffgens, teffgen, dNdead, Ndead
real :: stim(FINISHED), IL2sig(FINISHED), tgen, tnow, fac, act, cyt_conc, mols_pM
type (cog_type), pointer :: p
integer :: nst(FINISHED)
integer, allocatable :: gendist(:)
integer, allocatable :: div_gendist(:)  ! cells that are capable of dividing
character*(6) :: numstr
character*(256) :: msg

ok = .true.
allocate(gendist(TC_MAX_GEN))
allocate(div_gendist(TC_MAX_GEN))
tnow = istep*DELTA_T
noncog = 0
ncog = 0
ntot = 0
nst = 0
stim = 0
IL2sig = 0
gendist = 0
div_gendist = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    p => cellist(kcell)%cptr
    ntot = ntot + 1
    ctype = cellist(kcell)%ctype
!    stype = struct_type(ctype)
	if (associated(cellist(kcell)%cptr)) then
		stype = COG_TYPE_TAG
	else
		stype = NONCOG_TYPE_TAG
	endif
    if (stype == COG_TYPE_TAG) then
        ncog = ncog + 1
!        stage = get_stage(p)
		call get_stage(p,stage,region)
        nst(stage) = nst(stage) + 1
        stim(stage) = stim(stage) + p%stimulation
!        IL2sig(stage) = IL2sig(stage) + get_IL2store(p)
        gen = get_generation(p)
        if (gen < 0 .or. gen > TC_MAX_GEN) then
            write(logmsg,'(a,2i6)') 'show_snapshot: bad gen: ',kcell,gen
            call logger(logmsg)
            ok = .false.
            return
        endif
        gendist(gen) = gendist(gen) + 1
        if ((gen == 1 .and. p%stimulation > FIRST_DIVISION_THRESHOLD(1)) .or. &
			(gen > 1 .and. p%stimulation > DIVISION_THRESHOLD(1))) then
			div_gendist(gen) = div_gendist(gen) + 1
        endif
        max_TCR = max(p%stimulation,max_TCR)
    elseif (stype == NONCOG_TYPE_TAG) then
        noncog = noncog + 1
    else
        write(*,*) 'ERROR: show_snapshot: bad stype: ',ctype,stype
        stop
    endif
enddo
do i = 1,FINISHED
    if (nst(i) > 0) then
        stim(i) = stim(i)/nst(i)
        IL2sig(i) = IL2sig(i)/nst(i)
    else
        stim(i) = 0
        IL2sig(i) = 0
    endif
enddo
tgen = sum(gendist)
do i = TC_MAX_GEN,1,-1
    if (gendist(i) /= 0) exit
enddo
ngens = i

teffgen = sum(totalres%N_EffCogTCGen(1:TC_MAX_GEN))
do i = TC_MAX_GEN,1,-1
    if (totalres%N_EffCogTCGen(i) /= 0) exit
enddo
neffgens = i

dNdead = totalres%dN_Dead
Ndead = totalres%N_Dead
mols_pM = L_um3*M_pM/(NTcells*Vc*Navo)

if (teffgen > 0) then
    fac = 1/real(teffgen)
else
    fac = 0
endif
if (.not.use_TCP .and. use_cognate) then
write(*,'(a)') '----------------------------------------------------------------------'
write(*,*) 'use_cognate: ',use_cognate
write(*,'(a,i6,4i8,a,2i8)') 'snapshot: ',istep,ntot,ncogseed,ncog,'     dead: ',dNdead,Ndead
write(*,'(a,7i7)')   '# in stage:  ',nst
write(*,'(a,7f7.0)') 'stimulation: ',stim
!write(*,'(a,7f7.0)') 'IL2 signal:  ',IL2sig
write(*,'(a,2i8,4x,i8)') 'Recent efferent: ',totalres%dN_EffCogTC(2:3),sum(totalres%dN_EffCogTC)
write(*,'(a,2i8,4x,i8)') 'Total efferent:  ',totalres%N_EffCogTC(2:3),teffgen
write(*,'(a,10i6)')   'gen dist: ',(i,i=1,10)
write(*,'(a)')        'In node:  '
write(*,'(10x,10f6.3)') gendist(1:ngens)/tgen
write(*,'(a)')        'Efferent: '
write(*,'(10x,10f6.3)') fac*totalres%N_EffCogTCGen(1:neffgens)
write(*,'(a,i6,a,f6.0)') 'Live DC: ',NDCalive,'  max_TCR: ',max_TCR
!if (use_cytokines) then
!    do iseq = 1,Ncytokines
!        if (use_diffusion) then
!            cyt_conc = cyt_mean(iseq)
!        else
!            cyt_conc = cyt_mols(iseq)*mols_pM   ! -> conc in pM
!        endif
!        write(*,'(3a,f8.4)') 'Mean cytokine conc: ',cyt_name(cyt_tag(iseq)),'  ',cyt_conc
!    enddo
!endif
write(*,'(a)') '----------------------------------------------------------------------'

!call check_cognate_list
!kcog = 1
!kcell = cognate_list(kcog)
!if (kcell > 0) then
!    write(*,'(2i6,f8.4,f8.1)') kcog,kcell,cellist(kcell)%cptr%stimrate,cellist(kcell)%cptr%CD69
!endif
write(*,*) '========= Average time to IL2 threshold: ',nIL2thresh,tIL2thresh/max(1,nIL2thresh)
write(*,'(a)') '----------------------------------------------------------------------'
endif

!call get_cognate_dist(ncog1,ncog2)

write(nfout,'(2(i8,f8.0),7i8,25f7.4)') istep,tnow/60,NDCalive,act,ntot,ncogseed,ncog,dNdead,Ndead,teffgen, &
    fac*totalres%N_EffCogTCGen(1:TC_MAX_GEN)
if (use_tcp) then
!    if (.not.awp_1%is_open) then
!        call logger("in show_snapshot: awp_1 is not open")
!    endif
!    write(msg,'(2(i6,f8.0),5i8)') istep,tnow,NDCalive,act,ntot,ncogseed,ncog,Ndead,teffgen
!    call winsock_send(awp_1,msg,len_trim(msg),error)
    msg = ''
    do i = 1,ngens
		write(numstr,'(i6)') gendist(i)
		msg = trim(msg)//trim(adjustl(numstr))
		msg = trim(msg)//'('
		write(numstr,'(i6)') div_gendist(i)
		msg = trim(msg)//trim(adjustl(numstr))
		msg = trim(msg)//')-'
	enddo
    call logger(msg)
endif

! To plot outflow variation with time
!write(nfout,'(2f8.2)') tnow/60,OutflowTotal

!write(nfout,'(a,7f7.0)') 'stimulation: ',stim
!write(nfout,'(a,7f7.0)') 'IL2 signal:  ',IL2sig
deallocate(gendist)
totalres%dN_EffCogTC = 0
totalres%dN_EffCogTCGen = 0
totalres%dN_Dead = 0

if (save_DCbinding) then
    ncog = sum(dcbind(0:MAX_COG_BIND))
    write(*,'(21i6)') dcbind(0:MAX_COG_BIND)
    write(*,*) 'ncog: ',ncog
    write(nfdcbind,'(f8.2,i6,30f7.4)') tnow,ncog,dcbind(0:MAX_COG_BIND)/real(ncog)
    dcbind = 0
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Distance from the iexit exit (units = grids)
!-----------------------------------------------------------------------------------------
real function ExitDistance(site,iexit)
integer :: site(3),iexit
real :: r(3)

r = site - exitlist(iexit)%site
ExitDistance = norm(r)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checkExits(msg)
character*(*) :: msg
integer :: iexit,site(3)
logical :: ok = .true.

write(logmsg,'(a,i8,a,a)') 'checkExits: ',istep,'  ',msg
call logger(logmsg)
do iexit = 1,lastexit
    if (exitlist(iexit)%ID == 0) cycle
	site = exitlist(iexit)%site
	if (occupancy(site(1),site(2),site(3))%exitnum /= -exitlist(iexit)%ID) then
		write(logmsg,*) 'checkExits: ',iexit,site,occupancy(site(1),site(2),site(3))%exitnum
		call logger(logmsg)
		ok = .false.
	endif
enddo
if (.not.ok) stop
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine check_exit(iexit)
integer :: iexit
integer :: x,y,z,a,nc,nctot(0:10),site(3)
integer :: n = 5

site = exitlist(iexit)%site
nctot(0) = sitecells(site,0,0,0)
do a = 1,n
    nc = 0
    do z = -a,a,2*a
        do x = -a,a
            do y = -a,a
                nc = nc + sitecells(site,x,y,z)
            enddo
        enddo
    enddo
    do z = -a+1,a-1
        y = -a
        do x = -a+1,a
            nc = nc + sitecells(site,x,y,z)
        enddo
        x = a
        do y = -a+1,a
            nc = nc + sitecells(site,x,y,z)
        enddo
        y = a
        do x = -a,a-1
            nc = nc + sitecells(site,x,y,z)
        enddo
        x = -a
        do y = -a,a-1
            nc = nc + sitecells(site,x,y,z)
        enddo
    enddo
    nctot(a) = nc
enddo
write(nfout,'(10i6)') nctot(0:n),sum(nctot(0:n))
end subroutine

!-----------------------------------------------------------------------------------------
! The values returned, for setting the size of arrays TC_list[] etc. in the GUI, are not
! the precise numbers, rather they are upper bounds.
! This code assumes that only cognate cells will be displayed.
!-----------------------------------------------------------------------------------------
subroutine get_scene_dimensions(nTCsize, nDCsize, nbondsize) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: nTCsize, nDCsize, nbondsize
integer :: k, kcell, site(3), jb, idc

nTCsize = lastcogID
nDCsize = NDC
	k = 0
	do kcell = 1,nlist
		if (cellist(kcell)%ID == 0) cycle  ! gap
		site = cellist(kcell)%site
!		do jb = 1,2
!			idc = cellist(kcell)%DCbound(jb)
!			if (idc /= 0) then
!				if (DClist(idc)%capable) then
!					k = k+1
!				endif
!			endif
!		enddo
	enddo
	nbondsize = k
end subroutine

!-----------------------------------------------------------------------------------------
! The global variable values are computed from the global data:
!    NTcells
!    NDC
!    Radius
!    LNactivation
!    Inflow
!    Outflow
!    FlowFraction(:)
!-----------------------------------------------------------------------------------------
subroutine set_globalvar
real :: inflow0

!if (TAGGED_LOG_PATHS) then
!	InflowTotal = LOG_PATH_FACTOR*NTcells0*DELTA_T/(ave_residence_time*60)
!	OutflowTotal = InflowTotal
!	return
!endif
if (IN_VITRO) then
	InflowTotal = 0
	OutflowTotal = 0
	return
endif
!if (use_traffic) then
!    inflow0 = NTcells0*DELTA_T/(ave_residence_time*60)
!else
!    inflow0 = 0
!endif

!if (.not.steadystate) then     ! surrogate for modeling an immune response
!    call generate_traffic(inflow0)
!else
!    InflowTotal = inflow0
!    OutflowTotal = inflow0
!endif
if (istep == 1 .and. .not.use_TCP) then
	write(*,'(a,i8,4f8.2)') 'NTcells,Inflow,Outflow: ',NTcells0,InflowTotal,OutflowTotal
endif
end subroutine

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function f_blocked_egress_inflow(t,a)
real :: t, a
real :: inflow0

inflow0 = (DELTA_T/60)*NTcells0/ave_residence_time
f_blocked_egress_inflow = inflow0*(1 + a*(1-exp(-eblock%k1*t)))*exp(-eblock%k2*t)
!write(*,'(i6,4f8.1)') NTcells0,a,ave_residence_time,inflow0,f_blocked_egress_inflow
end function

!--------------------------------------------------------------------------------------
! When egress is suppressed, inflow initially increases, then decreases to 0.
! The inflow is non-zero for a specified duration of T_inflow (h).
! The steady-state inflow is NTcells0/ave_residence_time per hour.
! The total inflow in T_inflow is a specified multiple of NTcells0: expansion*NTcells0
! The functional variation for t: 0 -> T_inflow is given by f(t)
! We want the inflow quantities per timestep to add up to expansion*NTcells0
!--------------------------------------------------------------------------------------
subroutine setup_blocked_egress_inflow
real :: t, fsum, a
integer :: i, n

eblock%k1 = 0.8
eblock%k2 = 0.15
eblock%expansion = 1.3
eblock%duration = 24*60     ! min
do i = 1,20
    a = 1 + i*0.2
    n = 0
    fsum = 0
    t = 0
    do while (t < eblock%duration)
        fsum = fsum + f_blocked_egress_inflow(t,a)
        t = t + DELTA_T
        n = n+1
    enddo
    write(*,'(2i6,4e12.3)') i,n,a,t,fsum,fsum/((eblock%expansion-1)*NTcells0)
enddo
end subroutine

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine blocked_egress_inflow(inflow)
real :: inflow
real :: scale, t

t = istep*DELTA_T/60
inflow = f_blocked_egress_inflow(t,eblock%amp)

end subroutine

!--------------------------------------------------------------------------------------
! Interim
!--------------------------------------------------------------------------------------
subroutine set_stage(p,stage)
type(cog_type), pointer :: p
integer :: stage
integer :: status, oldstage, region
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

call get_stage(p,oldstage,region)
status = p%status
!statusbyte(STAGE_BYTE) = stage + (region-1)*STAGELIMIT
p%status = status
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine set_stage_region(p,stage,region)
type(cog_type), pointer :: p
integer :: stage, region
integer :: status
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

status = p%status
!statusbyte(STAGE_BYTE) = stage + (region-1)*STAGELIMIT
p%status = status
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine get_stage(p,stage,region)
type(cog_type), pointer :: p
integer :: stage, region
integer :: status
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

status = p%status
stage = statusbyte(STAGE_BYTE)
!region = LYMPHNODE
!if (stage > STAGELIMIT) then
!	stage = stage - STAGELIMIT
!	region = PERIPHERY
!endif
!write(logmsg,*) 'get_stage: ',stage,region
!call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine get_region(p,region)
type(cog_type), pointer :: p
integer :: region
integer :: stage, status
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

status = p%status
stage = statusbyte(STAGE_BYTE)
!region = LYMPHNODE
!if (stage > STAGELIMIT) then
!	region = PERIPHERY
!endif
end subroutine

!--------------------------------------------------------------------------------------
! Interim
!--------------------------------------------------------------------------------------
subroutine set_generation(p,gen)
type(cog_type), pointer :: p
integer :: gen
integer :: status
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

status = p%status
statusbyte(GENERATION_BYTE) = gen
p%status = status
end subroutine

!----------------------------------------------------------------------------------------
! Interim
!----------------------------------------------------------------------------------------
integer function get_generation(p)
type(cog_type), pointer :: p
integer :: status
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

status = p%status
get_generation = statusbyte(GENERATION_BYTE)
end function

!--------------------------------------------------------------------------------------
! Interim
!--------------------------------------------------------------------------------------
subroutine set_activation(p,active)
type(cog_type), pointer :: p
integer :: active
integer :: status
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

status = p%status
statusbyte(ACTIVATION_BYTE) = active
p%status = status
end subroutine

!----------------------------------------------------------------------------------------
! Interim
!----------------------------------------------------------------------------------------
integer function get_activation(p)
type(cog_type), pointer :: p
integer :: status
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

status = p%status
get_activation = statusbyte(ACTIVATION_BYTE)
end function

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
logical function is_activated(p)
type(cog_type), pointer :: p

is_activated = (get_activation(p) == 1)
end function

!-----------------------------------------------------------------------------------------
! Display the properties of a cognate cell.  The calling program must ascertain that kcell
! is cognate.
!-----------------------------------------------------------------------------------------
subroutine show_cognate_cell(kcell)
integer :: kcell
type (cog_type), pointer :: p
!integer :: cogID
type(cell_type) :: tcell
integer :: gen, stage, region

tcell = cellist(kcell)
if (.not.associated(tcell%cptr)) then
    write(*,*) 'ERROR: show_cognate_cell: cptr not associated: ',kcell
    stop
endif
p => tcell%cptr
write(*,*) 'Cognate cell: ',p%cogID,kcell,cellist(kcell)%ID
write(*,'(a,i10,a,3i4,a,i2)') '  ID: ',tcell%ID,' site: ',tcell%site,' ctype: ',tcell%ctype
gen = get_generation(p)
!stage = get_stage(p)
call get_stage(p,stage,region)
write(*,'(a,i8,a,i2,a,i2)') '   cogID: ',p%cogID,' gen: ',gen,' stage: ', stage
write(*,'(a,4f10.2)') '   times: entry,die,div,stage: ',tcell%entrytime,p%dietime,p%dividetime,p%stagetime

end subroutine

!-----------------------------------------------------------------------------------------
! Called whenever balancer carries out add_sites or removeSites.
! cognate_list(k) = 0 when the kth cognate cell has gone (left or died).
! If Mnodes = 1 this is called only once, after PlaceCells.  After that the list
! is maintained directly when cognate cells arrive or leave.
!-----------------------------------------------------------------------------------------
subroutine make_cognate_list(ok)
logical :: ok
integer :: kcell, ctype, stype, cogID
type (cog_type), pointer :: p

!write(*,*) 'make_cognate_list: ', lastcogID
ok = .true.
cognate_list(1:lastcogID) = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    ctype = cellist(kcell)%ctype
!    stype = struct_type(ctype)
	if (associated(cellist(kcell)%cptr)) then
		stype = COG_TYPE_TAG
	else
		stype = NONCOG_TYPE_TAG
	endif
    if (stype == COG_TYPE_TAG) then
        p => cellist(kcell)%cptr
        cogID = p%cogID
        if (cogID == 0) then
            lastcogID = lastcogID + 1
            if (lastcogID > MAX_COG) then
                write(logmsg,'(a,i6)') 'Error: make_cognate_list: cognate_list dimension exceeded: ',MAX_COG
                call logger(logmsg)
                ok = .false.
                return
            endif
            cogID = lastcogID
            p%cogID = cogID
        endif
        cognate_list(cogID) = kcell
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_log_paths
integer :: ic, k, kcell, nc(2), site(3)

nc = 0
do ic = 1,2
	do k = 1,n_log_path(ic)
		if (log_path(k,ic)%in .and. log_path(k,ic)%kstep < MAX_PATH_STEPS) then
			nc(ic) = nc(ic) + 1
			kcell = log_path(k,ic)%kcell
			site = cellist(kcell)%site
			log_path(k,ic)%kstep = log_path(k,ic)%kstep + 1
			log_path(k,ic)%pos(:,log_path(k,ic)%kstep) = site
		endif
	enddo
enddo
if (nc(2) == 0) then
	write(*,*) 'update_log_paths: all chemo cells have exited'
	call terminate_run(0)
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine end_log_path(kcell,ic)
integer :: kcell, ic
integer :: k
logical :: found

found = .false.
do k = 1,n_log_path(ic)
	if (log_path(k,ic)%kcell == kcell) then
		log_path(k,ic)%in = .false.
		found = .true.
		exit
	endif
enddo
if (.not.found) then
	write(*,*) 'Error: end_log_path: cell not found: ',kcell,ic
	stop
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine write_log_paths
character*(32) :: fname
real :: d(MAX_LOG_PATHS), r(3)
integer :: kstep, ic, k
logical :: done

write(*,*) 'Enter filename for paths'
read(*,'(a)') fname
open(nfpath,file=fname,status='replace')
do kstep = 1,MAX_PATH_STEPS
	ic = 2
	d = 0
	done = .true.
	do k = 1,MAX_LOG_PATHS
		if (kstep > log_path(k,ic)%kstep) cycle
		done = .false.
!		d(k) = max(0.5,ExitDistance(log_path(k,ic)%pos(:,kstep),1))
	enddo
	write(nfpath,'(i6,100f5.1)') kstep,d
	if (done) exit
enddo
close(nfpath)	
end subroutine

!--------------------------------------------------------------------------------
! Makes an approximate count of the number of sites of the spherical blob that
! are in the xth slice.  Uses the area of the slice.
! The blob centre is at (x0,y0,z0), and the blob radius is R = Radius
!--------------------------------------------------------------------------------
integer function slice_count(x)
integer :: x
real :: r2

slice_count = 1
!r2 = Radius**2 - (x-x0)**2
!if (r2 < 0) then
!    slice_count = 0
!else
!    slice_count = PI*r2
!endif
end function


!read(nfcell,*) TC_AVIDITY_MEAN              ! mean of avidity distribution (only if fix_avidity = false)
!read(nfcell,*) TC_AVIDITY_SHAPE			    ! shape -> 1 gives normal dist with small variance
!read(nfcell,*) TC_CD8_FRACTION				! fraction of all T cells that are CD8
!read(nfcell,*) TC_COGNATE_FRACTION(1)		! fraction of CD4 T cells that are cognate
!read(nfcell,*) TC_COGNATE_FRACTION(2)		! fraction of CD8 T cells that are cognate
!!read(nfcell,*) TC_CONVERSION_TIME			! time (in hours) over which T cells become cognate
!read(nfcell,*) TC_STIM_RATE_CONSTANT		! rate const for TCR stimulation (-> molecules/min)
!!read(nfcell,*) TC_STIM_WEIGHT				! contribution of stimulation to act level
!read(nfcell,*) TC_STIM_HALFLIFE				! halflife of T cell stimulation (hours)
!!read(nfcell,*) TC_MAX_GEN				    ! maximum T cell generation
!!read(nfcell,*) divide_dist1%p1
!!read(nfcell,*) divide_dist1%p2
!!read(nfcell,*) divide_dist2%p1
!!read(nfcell,*) divide_dist2%p2
!read(nfcell,*) divide_mean1
!read(nfcell,*) divide_shape1
!read(nfcell,*) divide_mean2
!read(nfcell,*) divide_shape2
!read(nfcell,*) BETA							! speed: 0 < beta < 1		(0.65)
!read(nfcell,*) RHO							! persistence: 0 < rho < 1	(0.95)
!
!read(nfcell,*) DC_ANTIGEN_MEAN				! mean DC antigen density
!read(nfcell,*) DC_ANTIGEN_SHAPE				! DC antigen density shape param
!read(nfcell,*) DC_LIFETIME_MEAN				! days
!read(nfcell,*) DC_LIFETIME_SHAPE 			! days
!!read(nfcell,*) DC_ACTIV_TAPER				! time (hours) over which DC activity decays to zero
!read(nfcell,*) DC_BIND_DELAY				! delay after unbinding before next binding (min)
!!read(nfcell,*) DC_BIND_ALFA					! binding prob parameter
!!read(nfcell,*) DC_MULTIBIND_PROB			! reducing factor to bind prob for each current DC binding
!!read(nfcell,*) DC_DENS_BY_STIM              ! rate of reduction of density by TCR stimulation
!read(nfcell,*) DC_DENS_HALFLIFE             ! base half-life of DC activity (hours)
!read(nfcell,*) MAX_TC_BIND						! number of T cells that can bind to a DC
!read(nfcell,*) MAX_COG_BIND					! number of cognate T cells that can bind to a DC simultaneously
!
!!read(nfcell,*) optionA                      ! 1,2
!!read(nfcell,*) optionB                      ! 1,2,3
!!read(nfcell,*) optionC                      ! 1,2
!!read(nfcell,*) IL2_PRODUCTION_TIME			! duration of IL2/CD25 production
!
!read(nfcell,*) IL2_THRESHOLD			    ! stimulation needed to initiate IL-2/CD25 production
!read(nfcell,*) ACTIVATION_THRESHOLD			! stimulation needed for activation
!read(nfcell,*) FIRST_DIVISION_THRESHOLD(1)	! activation level needed for first division
!read(nfcell,*) DIVISION_THRESHOLD(1)		! activation level needed for subsequent division
!read(nfcell,*) EXIT_THRESHOLD(1)			! activation level below which exit is permitted
!read(nfcell,*) STIMULATION_LIMIT			! maximum activation level
!read(nfcell,*) THRESHOLD_FACTOR             ! scales all threshold values
!read(nfcell,*) STAGED_CONTACT_RULE			! 2 = CT_HENRICKSON
!read(nfcell,*) STIM_HILL_THRESHOLD	        ! Normalized stim rate threshold for TCR signalling
!read(nfcell,*) STIM_HILL_N                  ! Parameters of Hill function for stimulation rate
!read(nfcell,*) STIM_HILL_C                  ! as function of x = (avidity/max avidity)*(pMHC/max pMHC)
!read(nfcell,*) ACTIVATION_MODE              ! indicates selection of STAGED_MODE or UNSTAGED_MODE
!read(nfcell,*) BINDTIME_HILL_THRESHOLD      ! potential normalized stimulation rate required for a cognate DC interaction
!read(nfcell,*) BINDTIME_HILL_N              ! N parameter for Hill function that determines bind duration
!read(nfcell,*) BINDTIME_HILL_C              ! C parameter for Hill function that determines bind duration
!read(nfcell,*) BINDTIME_MIN                 ! minimum cognate bind duration, i.e. kinapse (mins)
!read(nfcell,*) BINDTIME_MAX                 ! maximum cognate bind duration, i.e. synapse (mins converted from input hrs)
!read(nfcell,*) UNSTAGED_MIN_DIVIDE_T        ! minimum time elapsed before start of 1st division (mins or hrs?)
!read(nfcell,*) MAXIMUM_AVIDITY              ! maximum TCR avidity, used to normalize T cell avidity levels
!read(nfcell,*) MAXIMUM_ANTIGEN              ! maximum DC antigen density, used to normalize DC antigen density levels
!read(nfcell,*) K1_CD69
!read(nfcell,*) K2_CD69
!read(nfcell,*) K1_S1PR1
!read(nfcell,*) K2_S1PR1
!read(nfcell,*) S1PR1_EXIT_THRESHOLD
!
!!read(nfcell,*) CD25_DIVISION_THRESHOLD	    ! CD25 store level needed for division of activated cell (optionA = 1)
!!read(nfcell,*) CD25_SURVIVAL_THRESHOLD		! CD25 store level needed for survival of activated cell (optionC = 1)
!
!!read(nfcell,*) divide_dist1%dclass
!!read(nfcell,*) divide_dist2%dclass
!!read(nfcell,*) CD8_DIVTIME_FACTOR			! CD8 divide time as multiple of CD4 divide time
!
!read(nfcell,*) NX							! rule of thumb: about 4*BLOB_RADIUS
!read(nfcell,*) BLOB_RADIUS					! initial T cell blob size (sites)
!read(nfcell,*) TC_FRACTION					! T cell fraction of paracortex volume
!!read(nfcell,*) TC_RADIUS						! radius of T cell (um)
!read(nfcell,*) FLUID_FRACTION				! fraction of paracortex that is fluid
!read(nfcell,*) real_DCradius						! radius of DC sphere of influence (um)
!
!!read(nfcell,*) GAMMA						! controls crowding
!
!read(nfcell,*) TC_TO_DC						! number of T cells for every DC
!!read(nfcell,*) DC_FACTOR                    ! multiplying factor for DC number (initial and influx)
!!read(nfcell,*) MAX_DC_BIND					! number of DC that a cell can bind to simultaneously
!read(nfcell,*) DCrate_100k                  ! DC influx rate corresponding to NTcells0 = 100k
!read(nfcell,*) T_DC1                        ! Duration of constant DC influx (hours)
!read(nfcell,*) T_DC2                        ! Time of cessation of DC influx (hours)
!read(nfcell,*) DCsinjected					! select DC injection into experimental animals
!read(nfcell,*) DCinjectionfile				! file with schedule of DC injection
!!read(nfcell,*) T_DC_INJECTION				! Time of DC injection
!!read(nfcell,*) DC_FRACTION					! Fraction of DCs that are bearing antigen
!
!read(nfcell,*) usetraffic					! use T cell trafficking
!read(nfcell,*) usehev					    ! use HEV portals
!read(nfcell,*) useexitchemo                 ! use exit chemotaxis
!read(nfcell,*) useDCchemo					! use DC chemotaxis
!read(nfcell,*) cognateonly				    ! simulate only cognate cells
!read(nfcell,*) halveCD69				    ! halve CD69 on cell division
!read(nfcell,*) CD8_effector_prob            ! probability of CD8 effector switch on division
!read(nfcell,*) EXIT_RULE					! rule controlling egress of cognate cells
!read(nfcell,*) RESIDENCE_TIME(CD4)          ! CD4 T cell residence time in hours -> inflow rate
!read(nfcell,*) RESIDENCE_TIME(CD8)          ! CD8 T cell residence time in hours -> inflow rate
!! Vascularity parameters
!read(nfcell,*) Inflammation_days1	        ! Days of plateau level - parameters for VEGF_MODEL = 1
!read(nfcell,*) Inflammation_days2	        ! End of inflammation
!read(nfcell,*) Inflammation_level	        ! This is the level of inflammation
!
!	read(nfcell,*) receptor(CCR1)%strength      ! relative strength of CCL3-CCR1 chemotactic influence (chemo_K_DC)
!    read(nfcell,*) receptor(CCR1)%level(1)
!    read(nfcell,*) receptor(CCR1)%level(2)
!    read(nfcell,*) receptor(CCR1)%level(3)
!    read(nfcell,*) receptor(CCR1)%saturation_threshold	! CCL3 signal level to saturate CCR1 receptor
!    read(nfcell,*) receptor(CCR1)%refractory_time	! Time for CCR1 receptor to recover sensitivity after desensitization
!    read(nfcell,*) chemo(CCL3)%diff_coef
!    read(nfcell,*) chemo(CCL3)%halflife
!    read(nfcell,*) useCCL3_0
!    read(nfcell,*) useCCL3_1
!    read(nfcell,*) chemo(CCL3)%bdry_rate
!    read(nfcell,*) chemo(CCL3)%bdry_conc
!    read(nfcell,*) chemo(CCL3)%radius		! Sites within this radius of the DC receive CCL3 concentration (+ all DC sites)
!
!!read(nfcell,*) fix_avidity                  ! true if avidity takes discrete values, false if a distribution is used
!!read(nfcell,*) avidity_logscale             ! true if actual avidity = 10^(avidity_min + i*avidity_step)
!!read(nfcell,*) avidity_nlevels              ! If fix_avidity, number of discrete avidity values
!!read(nfcell,*) avidity_min                  ! minimum value
!!read(nfcell,*) avidity_step                 ! step between equi-spaced values
!
!read(nfcell,*) days							! number of days to simulate
!read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
!read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
!read(nfcell,*) ncpu_dummy					! just a placeholder for ncpu, not used currently
!read(nfcell,*) NT_GUI_OUT					! interval between GUI outputs (timesteps)
!read(nfcell,*) facs_h						! interval between FACS outputs (h)
!read(nfcell,*) SPECIES						! animal species source of T cells
!read(nfcell,*) invitro 						! select in vivo or in vitro simulation
!read(nfcell,*) IV_WELL_DIAMETER				! diameter of in vitro well (mm)
!read(nfcell,*) IV_NTCELLS					! initial T cell population in vitro
!read(nfcell,*) IV_COGNATE_FRACTION			! fraction of in vitro cells that are cognate for DC antigen
!read(nfcell,*) shownoncog			        ! display non-cognate T cells
!read(nfcell,*) ispecial						! special case
!read(nfcell,*) specialfile					! special case input data file
!read(nfcell,*) fixedfile					! file with "fixed" parameter values


! T cell parameters
logical :: TCR_splitting = .false.      ! enable sharing of integrated TCR signal between progeny cells
real(REAL_KIND) :: transient_stagetime				! stagetime for TRANSIENT (Stage 1)
real(REAL_KIND) :: clusters_stagetime				! stagetime for CLUSTERS (Stage 2)
real(REAL_KIND) :: transient_bindtime				! bindtime for TRANSIENT (Stage 1)
real(REAL_KIND) :: clusters_bindtime				! bindtime for CLUSTERS (Stage 2)
real(REAL_KIND) :: swarms_bindtime					! bindtime for SWARMS (Stage 3)
real(REAL_KIND) :: TC_life_median1					! median lifetime of naive T cells
real(REAL_KIND) :: TC_life_median2					! median lifetime of activated T cells
real(REAL_KIND) :: TC_life_shape					! shape parameter for lifetime of T cells
integer :: NTC_LN = 3.0e07				! number of T cells in a LN
integer :: NTC_BODY = 1.6e09			! number of circulating T cells in the whole body
integer :: NLN_RESPONSE					! number of LNs in the response
real(REAL_KIND) :: K1_S1PR1 != 0.005				! S1PR1/CD69 system parameters
real(REAL_KIND) :: K2_S1PR1 != 0.05
real(REAL_KIND) :: K1_CD69 != 0.04
real(REAL_KIND) :: K2_CD69 != 0.01
real(REAL_KIND) :: S1PR1_EXIT_THRESHOLD

!--------------------------------------------------------------------------------
! Generates the arrays wz(), zoffset() and zdomain().
! The domains (slices) are numbered 0,...,2*Mnodes-1
! wz(k) = width of the slice for kth domain
! zoffset(k) = offset of kth domain occupancy array in the occupancy array.
! zdomain(x) = domain that global z lies in.
! The kth domain (slice) extends from z = zoffset(k)+1 to z = zoffset(k+1)
! The idea is to set the domain boundaries such that each domain has roughly the
! same number of available sites.
! This is the initial split, which will continue to be OK if:
! not using a blob, or Mnodes <= 2
! If IN_VITRO, cells lie in the x-y plane (z=1), and the split is based on x value,
! => xoffset(k), xdomain(k)
!--------------------------------------------------------------------------------
subroutine make_split(force)
logical :: force
integer :: k, wsum, kdomain, nsum, Ntot, N, last, x, y, z
integer, allocatable :: scount(:)
integer, allocatable :: wz(:), ztotal(:)
integer :: Mslices
real(REAL_KIND) :: dNT, diff1, diff2
logical :: show = .false.

if (IN_VITRO) then
	xdomain = -1
	if (Mnodes == 1) then
		xdomain = 0
		Mslices = 1
		return
	endif
	Mslices = 2*Mnodes
	! Divide up a disk into Mslices equal pieces
	! First find total number of sites inside
	Ntot = 0
	z = 1
	do x = 1,NX
		do y = 1,NY
			if (occupancy(x,y,z)%indx(1) /= OUTSIDE_TAG) Ntot = Ntot + 1
		enddo
	enddo
	k = 0
	nsum = 0
	do x = 1,NX
		do y = 1,NY
			if (occupancy(x,y,z)%indx(1) /= OUTSIDE_TAG) then
				if (k == 0) then
					xoffset(0) = x - 1
					k = 1
				endif
				nsum = nsum + 1
			endif
		enddo
		if (k > 0 .and. nsum >= (k*Ntot)/Mslices) then
			xoffset(k) = x
			k = k+1
			if (k > Mslices) exit
		endif
	enddo
	do k = 0,Mslices-1
		do x = xoffset(k)+1,xoffset(k+1)
			xdomain(x) = k
		enddo
	enddo
	write(logmsg,'(a,12i4)') 'xoffset: ',xoffset
	call logger(logmsg)
	write(nflog,'(a,12i4)') 'xoffset: ',xoffset
	write(nflog,*) 'xdomain:'
	write(nflog,'(10i5)') xdomain(1:NX)
	return
endif

if (Mnodes == 1) then
    Mslices = 1
    zdomain = 0
else
	Mslices = 2*Mnodes
endif
dNT = abs(NTcells - lastNTcells)/real(lastNTcells+1)
if (.not.force .and. dNT < 0.03) then
    return
endif
lastNTcells = NTcells
if (Mslices > 1) then
	allocate(wz(0:Mslices))
	allocate(ztotal(0:Mslices))
	allocate(scount(NX))
endif
blobrange(:,1) = 99999
blobrange(:,2) = 0
nsum = 0
do z = 1,NZ
    k = 0
    do y = 1,NY
        do x = 1,NX
            if (occupancy(x,y,z)%indx(1) /= OUTSIDE_TAG) then
                k = k + 1
                blobrange(1,1) = min(blobrange(1,1),x)
                blobrange(1,2) = max(blobrange(1,2),x)
                blobrange(2,1) = min(blobrange(2,1),y)
                blobrange(2,2) = max(blobrange(2,2),y)
                blobrange(3,1) = min(blobrange(3,1),z)
                blobrange(3,2) = max(blobrange(3,2),z)
            endif
        enddo
    enddo
    if (Mslices > 1) then
	    scount(z) = k
	    nsum = nsum + scount(z)
	endif
enddo
if (Mslices == 1) return

Ntot = nsum
N = Ntot/Mslices
nsum = 0
last = 0
k = 0
do z = 1,NZ
    nsum = nsum + scount(z)
    if (nsum >= (k+1)*N) then
        diff1 = nsum - (k+1)*N
        diff2 = diff1 - scount(z)
        if (abs(diff1) < abs(diff2)) then
            wz(k) = z - last
            last = z
        else
            wz(k) = z - last - 1
            last = z - 1
        endif
        k = k+1
        if (k == Mslices-1) exit
    endif
enddo
wz(Mslices-1) = NZ - last
if (show) then
    write(*,*) 'Ntot, N: ',Ntot,N
    write(*,'(10i6)') scount
endif

zoffset(0) = 0
do k = 1,Mslices-1
    zoffset(k) = zoffset(k-1) + wz(k-1)
enddo
zoffset(Mslices) = NZ
z = 0
do kdomain = 0,Mslices-1
    do k = 1,wz(kdomain)
        z = z+1
        zdomain(z) = kdomain      ! = kpar with two sweeps
    enddo
enddo
if (show) then
    write(*,*) 'zoffset: ',zoffset
    write(*,*) 'wz:      ',wz
    write(*,*) 'zdomain: '
    write(*,'(10i4)') zdomain
endif
ztotal = 0
do k = 0,2*Mnodes-1
    do z = zoffset(k)+1,zoffset(k+1)
        ztotal(k) = ztotal(k) + scount(z)
    enddo
    if (show) write(*,*) k,ztotal(k)
enddo
deallocate(wz)
deallocate(ztotal)
deallocate(scount)
end subroutine

!-----------------------------------------------------------------------------------------
! In this version the parallel section has been made into a subroutine.
! Try varying the sweep order, i.e. 0,1 and 1,0
!-----------------------------------------------------------------------------------------
subroutine mover(ok)
logical :: ok
integer :: kpar=0, sweep, nsweeps, sweep1, sweep2, dsweep

if (Mnodes > 1) then
!DEC$ IF .NOT. DEFINED (_OPENMP)
!    stop
!DEC$ ENDIF
endif

if (Mnodes == 1) then
    nsweeps = 1
else
    nsweeps = 2
endif

if (mod(istep,2) == 0) then
	sweep1 = 0
	sweep2 = nsweeps-1
	dsweep = 1
else
	sweep2 = 0
	sweep1 = nsweeps-1
	dsweep = -1
endif
do sweep = sweep1,sweep2,dsweep
	if (IN_VITRO) then
		if (Mnodes > 1) then
		!$omp parallel do
			do kpar = 0,Mnodes-1
				call par_mover2D(sweep,kpar,ok)
			enddo
		!$omp end parallel do
		else
			call par_mover2D(sweep,kpar,ok)
		endif
		if (.not.ok) return
	else
		if (Mnodes > 1) then
		!$omp parallel do
			do kpar = 0,Mnodes-1
				call par_mover(sweep,kpar)
			enddo
		!$omp end parallel do
		else
			call par_mover(sweep,kpar)
		endif
	endif
enddo
ok = .true.
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine par_mover(sweep,kpar)
integer :: sweep
integer :: kpar
integer :: site1(3), kcell, indx(2), slot, z, slice, site(3), stage, region
logical :: go
type(cell_type), pointer :: cell
integer :: z_lo,z_hi

slice = sweep + 2*kpar
do kcell = 1,nlist
	if (Mnodes == 1) then
	    z_lo = 1
	    z_hi = NZ
	else
	    z_lo = zoffset(slice) + 1
	    z_hi = zoffset(slice+1)
	endif
    cell => cellist(kcell)
    if (cell%ID == 0) cycle             ! skip gaps in the list
    if (cell%step == istep) cycle
    site1 = cell%site
    z = site1(3)
    if (zdomain(z) /= slice) cycle      ! not in the slice for this processor
    indx = occupancy(site1(1),site1(2),site1(3))%indx
    if (indx(1) < 0) then
        write(*,*) 'Error: par_mover: OUTSIDE_TAG or DC: ',kcell,site1,indx
        stop
    endif
    if (kcell == indx(1)) then
        slot = 1
    elseif (kcell == indx(2)) then
        slot = 2
    else
        write(*,'(a,6i8)') 'Error: par_mover: bad indx: ',kcell,site1,indx
        stop
    endif
	call jumper(kcell,indx,slot,go,kpar)
    cell%step = istep
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! For (mainly) 2D motion, i.e. IN_VITRO simulation.
!-----------------------------------------------------------------------------------------
subroutine par_mover2D(sweep,kpar,ok)
integer :: sweep
integer :: kpar
logical :: ok
integer :: site1(3), kcell, indx(2), slot, x, slice, site(3)
logical :: go
type(cell_type), pointer :: cell
integer :: x_lo,x_hi

slice = sweep + 2*kpar
do kcell = 1,nlist
	if (Mnodes == 1) then
	    x_lo = 1
	    x_hi = NX
	else
	    x_lo = xoffset(slice) + 1
	    x_hi = xoffset(slice+1)
	endif
    cell => cellist(kcell)

    if (cell%ID == 0) cycle             ! skip gaps in the list
    if (cell%step == istep) cycle
    site1 = cell%site
    x = site1(1)
    if (xdomain(x) /= slice) cycle      ! not in the slice for this processor
    indx = occupancy(site1(1),site1(2),site1(3))%indx
    if (indx(1) < 0) then
        write(logmsg,'(a,6i6)') 'Error: par_mover2D: OUTSIDE_TAG or DC: ',kcell,site1,indx
        call logger(logmsg)
        ok = .false.
        return
    endif
    if (kcell == indx(1)) then
        slot = 1
    elseif (kcell == indx(2)) then
        slot = 2
    else
        write(logmsg,'(a,6i8)') 'Error: par_mover2D: bad indx: ',kcell,site1,indx
        call logger(logmsg)
        ok = .false.
        return
    endif
    call jumper2D(kcell,indx,slot,go,kpar)
    cell%step = istep
enddo
ok = .true.
end subroutine

!--------------------------------------------------------------------------------
! Computes the jump probabilities (absolute directions) accounting for chemotaxis
! towards exit (or DC).  On input:
!   p(:) holds the jump probabilities not accounting for  chemotaxis
!   v(:) is the site offset relative to the exit (or DC)
! CHANGED
!   v(:) is offset the exit (or DC) offset relative to the site
! OR
!   v(:) is the concentration gradient vector
!   f is the amount of chemotactic influence
!   c is the amount of CCR7 ligand influence (Cyster).
! Note that f incorporates both the distance from the exit (or DC) and the cell's
! susceptibility to chemotaxis, which may be the S1PR1 level of the T cell.
! On return p(:) holds the modified jump probabilities.
! Note: should have p remaining unchanged as chemo_K -> 0
! Note: when njumpdirs = 27, jump 14 corresponds to (0,0,0) - unused.
! Note: code modifications now base v and f on net chemotactic attraction of 
! multiple exits and DCs.
!--------------------------------------------------------------------------------
subroutine chemo_probs(p,v,f,c)
real(REAL_KIND) :: p(:)
integer :: v(:)
real(REAL_KIND) :: f, c
integer :: k
real(REAL_KIND) :: pc(MAXRELDIR+1)

if (f == 0 .and. c == 1) then
    return
endif
if (f > 1) then
	write(logmsg,*) 'ERROR: chemo_probs: f > 0: ',f
	call logger(logmsg)
	return
endif
p = p/sum(p)
if (f > 0) then
    pc(1:njumpdirs) = chemo_p(v(1),v(2),v(3),:)
else
    pc = 0
endif
do k = 1,njumpdirs
    if (p(k) > 0) then      ! prevents jumps in disallowed directions
        p(k) = (1-f)*c*p(k) + f*pc(k)
    endif
enddo
end subroutine



