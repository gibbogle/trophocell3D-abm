! Code that is currently not used

!----------------------------------------------------------------------------------------- 
! This is the original version, using files to communicate signals and data with the Qt main.
! NOT USED
!-----------------------------------------------------------------------------------------
subroutine simulate(ok)
!!!use ifport
logical :: ok
real :: cave,tres
real(8) :: t1,t2,td
integer :: save_step1, save_step2, save_count=0
integer :: k, nd, kpar = 0
real(DP) :: R
!type(cog_type), pointer :: p

t1 = secnds(tres)
!write(logmsg,'(a,f10.2)') 't1: ',t1
!call logger(logmsg)

!t1 = mpi_wtime()
!t1 = timef()
ok = .true.
write(logmsg,'(a,i6)') 'Start simulation: Nsteps: ',Nsteps
call logger(logmsg)
istep = 0
if (use_cognate .and. use_DC) then
    call show_snapshot(ok)
    if (.not.ok) return
endif
clear_to_send = .true.
simulation_start = .true.
dbug = .false.
do istep = 1,Nsteps
!    call check_xyz(0)
    write(nfres,'(i6,f10.6)') istep,par_uni(kpar)

	inquire(file=stopfile,exist=stopped)
	if (stopped) then
		write(logmsg,'(a)') 'Stop order received'
		call logger(logmsg)
		ok = .true.
		return
	endif
	call check_pause

    if (idbug > 0) write(*,'(4i6,f10.6)') istep,cellist(idbug)%site    !,par_uni(kpar)
!        if (evaluate_residence_time .and. istep == istep_res1) then
!            call tag_cells
!        endif
	! Transfer cell location data to GUI
	if (mod(istep,NT_GUI_OUT) == 0 .and. use_TCP) then
!		write(logmsg,'(a,i6)') 'call save_pos: ',istep
!		call logger(logmsg)
		call save_pos
	endif
	if (mod(istep,240) == 0) then
        globalvar%Radius = (globalvar%NTcells*3/(4*PI))**0.33333
	    if (use_cognate .and. use_DC) then
	        call show_snapshot(ok)
	        if (.not.ok) return
	    endif
        if (log_traffic) then
            call check_chemoactivity(cave)
            write(nftraffic,'(5i8,3f8.3)') istep, globalvar%NTcells, globalvar%Nexits, total_in, total_out, &
                globalvar%InflowTotal, globalvar%vascularity,cave
        endif
        total_in = 0
	    total_out = 0
!            call arrayview(NX)
	endif
	if (log_results .and. mod(istep,ntres) == 0) then
 	    call write_results
	endif
	
	! Start step simulation
	
    if (use_cytokines) then
        if (use_diffusion) then
            call diffuser
        else
            call molsynch
        endif
    endif

	if (use_traffic .and. mod(istep,SCANNER_INTERVAL) == 0) then
		if (dbug) write(*,*) 'call scanner'
		call scanner
		if (dbug) write(*,*) 'did scanner'
	endif
!	if (dbug) write(*,*) 'call mover'
	call mover(ok)
    if (dbug) write(nfres,'(a,i6,f10.6)') 'mover: ',istep,par_uni(kpar)
	if (.not.ok) return
!	if (dbug) write(*,*) 'did mover'
    if (dbug) call check_xyz(1)

    if (use_DC .and. globalvar%NDCalive > 0) then
		if (dbug) write(*,*) 'call binder'
        call binder(ok)
        if (.not.ok) return
		if (dbug) write(*,*) 'did binder'
        if (.not.track_DCvisits .and. .not.evaluate_residence_time) then
			if (dbug) write(*,*) 'call update_DCstate'
            call update_DCstate(ok)
            if (.not.ok) return
			if (dbug) write(*,*) 'did update_DCstate'
        endif
    endif
    if (.not.track_DCvisits .and. .not.evaluate_residence_time) then
!		if (dbug) write(*,*) 'call updater'
        call updater(ok)
        if (dbug) write(nfres,'(a,i6,f10.6)') 'updater: ',istep,par_uni(kpar)
        if (.not.ok) return
!		if (dbug) write(*,*) 'did updater'
    endif
    if (dbug) call check_xyz(2)

    if (use_traffic) then
        if (vary_vascularity) then	! There is a problem with this system
            call vascular
        endif
        if (use_portal_egress) then
			if (dbug) write(*,*) 'call portal_traffic'
            call portal_traffic(ok)
            if (.not.ok) return
			if (dbug) write(*,*) 'did portal_traffic'
        else
            call traffic(ok)
            if (dbug) write(nfres,'(a,i6,f10.6)') 'traffic: ',istep,par_uni(kpar)
            if (.not.ok) return
        endif
    endif
    if (dbug) call check_xyz(3)

	if (dbug) write(*,*) 'call balancer'
    call balancer(ok)
    if (.not.ok) return
	if (dbug) write(*,*) 'did balancer'
	if (dbug) write(*,*) 'call set_globalvar'
    call set_globalvar
	if (dbug) write(*,*) 'did set_globalvar'
	
	! End step simulation
	
    if (evaluate_residence_time) then
        if (noutflow_tag > 0) then
            if (mod(istep,10) == 0) then
                write(nfout,'(f8.1,i8,f8.1)') istep*DELTA_T, noutflow_tag, restime_tot/noutflow_tag
                if (ninflow_tag > 0 .and. noutflow_tag == ninflow_tag) then
                    write(logmsg,*) 'All tagged cells have exited'
                    call logger(logmsg)
                    call write_Tres_dist
                    stop
                endif
            endif
        endif
    endif
!        if (Mnodes == 1 .and. save_pos_cmgui) then
    if (save_pos_cmgui) then
        if (save_count == 0) then
            save_step1 = 1
            save_step2 = save_length_hours*60*nsteps_per_min
        else
            save_step2 = save_count*save_interval_hours*60*nsteps_per_min
            save_step1 = save_step2 - save_length_hours*60*nsteps_per_min + 1
        endif
        if (istep >= save_step1 .and. istep <= save_step2) then
            call save_pos
            if (istep == save_step2) then
                save_count = save_count + 1
            endif
        endif
    endif
enddo

!write(logmsg,'(a,f10.2)') 't1: ',t1
!call logger(logmsg)

t2 = secnds(tres)
!t2 = mpi_wtime()
!t2 = timef()
write(logmsg,'(a,3f10.2)') 'Execution completed, time (seconds): ',t2-t1
call logger(logmsg)

if (log_results) then
    write(nfout,*) 'Total efferent avidity distribution:'
    write(nfout,'(10f7.4)') real(avid_count_total%ndist)/sum(avid_count_total%ndist)
endif

if (.not.use_TCP) then
	nd = 0
	td = 0
	write(*,*) 'Average times between divisions'
	do k = 1,MMAX_GEN
		if (ndivided(k) > 0) then
			write(*,'(i3,i6,f6.1)') k,ndivided(k),tdivided(k)/(60*ndivided(k))
			if (k > 1) then
				nd = nd + ndivided(k)
				td = td + tdivided(k)
			endif
		else
			exit
		endif
	enddo
	write(*,'(a,f6.1)') 'Average for all later divisions: ',td/max(1,nd)
endif
ok = .true.
end subroutine

!--------------------------------------------------------------------------------
! Parallel
! This doesn't work well for Mnodes=4 because the DC sphere of influence may span
! more than a stripe width, therefore there can be simultaneous writing to the
! data for a DC from two threads.
! This problem is less likely to occur as the number of T cells increases, in
! fact what is relevant is the slice width, which is proportional to the
! (cube root of NTcells)/Mnodes, i.e. to Rblob/Mnodes.
!--------------------------------------------------------------------------------
subroutine binder2
integer :: nsweeps, sweep

if (Mnodes == 1) then
    nsweeps = 1
else
    nsweeps = 2
endif

do sweep = 0,nsweeps-1
!$omp parallel
!    call par_binder(sweep)     ! NOT USED
!$omp end parallel
enddo
end subroutine

!--------------------------------------------------------------------------------
! Not used
!--------------------------------------------------------------------------------
subroutine par_binder(sweep)
integer :: sweep
integer :: kcell, nbnd, k, site(3), ctype, stage, region, neardc(0:DCDIM-1), idc, DCbound(2)
integer :: nadd, nsub
real :: tnow, bindtime, unbindtime(2), pMHC
logical :: unbound, cognate
type(cell_type), pointer :: cell
type(cog_type), pointer :: cog_ptr
integer :: kpar = 0
integer :: z_lo,z_hi

kpar = omp_get_thread_num()
nadd = 0
nsub = 0
tnow = istep*DELTA_T
do kcell = 1,nlist
	if (Mnodes == 1) then
	    z_lo = 1
	    z_hi = NX
	else
	    z_lo = zoffset(sweep+2*kpar) + 1
	    z_hi = zoffset(sweep+2*kpar+1)
	endif
    cell => cellist(kcell)
!    if (cellist(kcell)%ID == 0) cycle
    if (cell%ID == 0) cycle
    site = cell%site
    if (site(3) < z_lo .or. site(3) > z_hi) cycle      ! not in the slice for this processor

    ! Handle unbinding first
    nbnd = 0
!    DCbound = cellist(kcell)%DCbound
!    unbindtime = cellist(kcell)%unbindtime
    DCbound = cell%DCbound
    unbindtime = cell%unbindtime
    unbound = .false.
    do k = 1,MAX_DC_BIND
        idc = DCbound(k)
        if (idc == 0) cycle
        nbnd = nbnd + 1
        if (tnow >= unbindtime(k) .or. .not.DClist(idc)%capable) then   ! unbind
            DClist(idc)%nbound = DClist(idc)%nbound - 1
            if (DClist(idc)%nbound < 0) then
                write(logmsg,*) 'binder: DClist(idc)%nbound < 0: ',kcell,nbnd,idc,DCbound,unbindtime,tnow
				call logger(logmsg)
                stop
            endif
            DCbound(k) = 0
            nbnd = nbnd - 1
            unbound = .true.
           nsub = nsub + 1
        endif
     enddo
     if (unbound) then  ! ensure that a remaining binding is stored in slot 1
        if (DCbound(1) == 0 .and. DCbound(2) /= 0) then
            DCbound(1) = DCbound(2)
            DCbound(2) = 0
            unbindtime(1) = unbindtime(2)
        endif
        unbindtime(2) = tnow
    endif

    ! Handle new binding
    if (nbnd < MAX_DC_BIND .and. tnow >= unbindtime(2) + DC_BIND_DELAY) then
        if (evaluate_residence_time .or. track_DCvisits) then
            cognate = .false.
        else
!            cog_ptr => cellist(kcell)%cptr
            cog_ptr => cell%cptr
            cognate = (associated(cog_ptr))
            if (cognate) then
                if (.not.bindable(cog_ptr)) then
!                    cellist(kcell)%DCbound = DCbound
!                    cellist(kcell)%unbindtime = unbindtime
                    cell%DCbound = DCbound
                    cell%unbindtime = unbindtime
                    cycle
                endif
!                stage = get_stage(cog_ptr)
				call get_stage(cog_ptr,stage,region)
                if (stage == NAIVE) then    ! on first binding a cognate cell is ready for next stage
                    cog_ptr%stagetime = tnow
                endif
            endif
        endif
!        site = cellist(kcell)%site
!        site = cell%site
        neardc = occupancy(site(1),site(2),site(3))%DC
        if (neardc(0) /= 0) then
            do k = 1,neardc(0)
                idc = neardc(k)
                if (.not.DClist(idc)%capable) cycle
                if (idc == DCbound(1)) cycle      ! valid for MAX_DC_BIND <= 2
                if (bindDC(idc,kpar)) then
                    nbnd = nbnd + 1
                    if (DCbound(nbnd) /= 0) then
                        write(logmsg,*) 'binder: DCbound(nbnd) /= 0: ',DCbound(nbnd)
						call logger(logmsg)
                        stop
                    endif
                    DCbound(nbnd) = idc
!                    ctype = cellist(kcell)%ctype
                    ctype = cell%ctype
                    pMHC = DClist(idc)%density
                    bindtime = get_bindtime(cog_ptr,cognate,ctype,pMHC,kpar)
                    unbindtime(nbnd) = tnow + bindtime
                    nadd = nadd + 1
                    DClist(idc)%nbound = DClist(idc)%nbound + 1
                    if (DCbound(1) == 0 .and. DCbound(2) /= 0) then
                        write(logmsg,*) 'DCbound order: ',kcell,DCbound
						call logger(logmsg)
                        stop
                    endif
					if (track_DCvisits .and. ctype == TAGGED_CELL) then
						if (revisit(kcell,idc)) then
							Nrevisits = Nrevisits + 1
!							cellist(kcell)%revisits = cellist(kcell)%revisits + 1
							cell%revisits = cell%revisits + 1
						else
							Nvisits = Nvisits + 1
!							cellist(kcell)%visits = cellist(kcell)%visits + 1
!							cellist(kcell)%ndclist = cellist(kcell)%ndclist + 1
!							cellist(kcell)%dclist(cellist(kcell)%ndclist) = idc
							cell%visits = cell%visits + 1
							cell%ndclist = cell%ndclist + 1
							cell%dclist(cell%ndclist) = idc
						endif
					endif
                    if (nbnd == MAX_DC_BIND) exit
                endif
            enddo
        endif
    endif
!    cellist(kcell)%DCbound = DCbound
!    cellist(kcell)%unbindtime = unbindtime
    cell%DCbound = DCbound
    cell%unbindtime = unbindtime
enddo
!write(*,*) 'sweep, kpar, nadd, nsub: ',sweep,kpar,nadd,nsub
end subroutine

!-----------------------------------------------------------------------------------------
! In the event-based formulation, only cognate T cells are considered, and there is no
! spatial aspect.
! NOT USED
!-----------------------------------------------------------------------------------------
subroutine event_binder
integer :: kcell, nbnd, k, ctype, stage, region, idc=0
integer :: nadd, nsub
!integer :: neardc(0:DCDIM-1)
real :: tnow, bindtime, pMHC
logical :: unbound, cognate
type(cell_type), pointer :: cell
type(cog_type), pointer :: cog_ptr
integer :: kpar = 0

nadd = 0
nsub = 0
tnow = istep*DELTA_T
do kcell = 1,nlist
    cell => cellist(kcell)
    if (cell%ID == 0) cycle

    cog_ptr => cell%cptr
!    cognate = (associated(cog_ptr))

    ! Handle unbinding first
    nbnd = 0
    unbound = .false.
!    do k = 1,MAX_DC_BIND
	k = 1
!        idc = cell%DCbound(k)
!        if (idc == 0) cycle
        nbnd = nbnd + 1
        if (tnow >= cell%unbindtime(k)) then   ! unbind
!            DClist(idc)%nbound = DClist(idc)%nbound - 1
!            if (DClist(idc)%nbound < 0) then
!                write(*,*) 'binder: DClist(idc)%nbound < 0: ',kcell,nbnd,idc,cell%DCbound,cell%unbindtime,tnow
!                stop
!            endif
!            if (cognate) then
!                DClist(idc)%ncogbound = DClist(idc)%ncogbound - 1
!                if (DClist(idc)%ncogbound < 0) then
!                    write(*,*) 'binder: DClist(idc)%ncogbound < 0: ',kcell,nbnd,idc,cell%DCbound,cell%unbindtime,tnow
!                    stop
!                endif
!            endif
            cell%DCbound(k) = 0
            nbnd = nbnd - 1
            unbound = .true.
!           nsub = nsub + 1
        endif
!     enddo
     if (unbound) then  ! ensure that a remaining binding is stored in slot 1
!        if (cell%DCbound(1) == 0 .and. cell%DCbound(2) /= 0) then
!            cell%DCbound(1) = cell%DCbound(2)
!            cell%DCbound(2) = 0
!            cell%unbindtime(1) = cell%unbindtime(2)
!        endif
        cell%unbindtime(2) = tnow
!!        DCencounter(kcell) = tnow + get_traveltime()
    endif

    ! Handle new binding
!    if (nbnd < MAX_DC_BIND .and. tnow >= cell%unbindtime(2) + DC_BIND_DELAY) then
    if (nbnd < 1) then	! .and. tnow >= DCencounter(kcell)) then
!        if (cognate) then

!            if (.not.bindable(cog_ptr)) cycle
!            stage = get_stage(cog_ptr)
			call get_stage(cog_ptr,stage,region)
            if (stage == NAIVE) then    ! on first binding a cognate cell is ready for next stage
                cog_ptr%stagetime = tnow
            endif
!        endif
!        site = cell%site
!        neardc = occupancy(site(1),site(2),site(3))%DC      ! list of DC near this site
!        if (neardc(0) /= 0) then
!            do k = 1,neardc(0)
!                idc = neardc(k)
!                if (.not.DClist(idc)%capable) cycle
!                if (idc == cell%DCbound(1)) cycle      ! valid for MAX_DC_BIND <= 2
!                if (DClist(idc)%ncogbound > MAX_COG_BIND) then      ! ERROR CHECKING
!                    write(*,*) 'binder: ncogbound: ',idc,DClist(idc)%ncogbound
!                    stop
!                endif
!                if (cognate .and. DClist(idc)%ncogbound == MAX_COG_BIND) cycle
!                if (bindDC(idc,kpar)) then
!					if (compute_travel_time .and. cognate .and. nbnd == 1) then
!					if (tnow > 2*60 .and. compute_travel_time .and. nbnd == 1) then
!						t_travel = tnow - cell%unbindtime(2)
!						ntravel = ntravel + 1
!						i = t_travel + 0.5
!						i = max(i,1)
!						i = min(i,N_TRAVEL_DIST)
!						travel_dist(k_travel_cog,k_travel_dc,i) = travel_dist(k_travel_cog,k_travel_dc,i) + 1
!					endif
                    nbnd = nbnd + 1
!                    if (cell%DCbound(nbnd) /= 0) then
!                        write(*,*) 'binder: DCbound(nbnd) /= 0: ',cell%DCbound(nbnd)
!                        stop
!                    endif
!                    cell%DCbound(nbnd) = idc
                    ctype = cell%ctype
                    pMHC = DClist(idc)%density
                    bindtime = get_bindtime(cog_ptr,cognate,ctype,pMHC,kpar)
                    cell%unbindtime(nbnd) = tnow + bindtime
                    nadd = nadd + 1
                    DClist(idc)%nbound = DClist(idc)%nbound + 1
!                   if (cognate) then
                        DClist(idc)%ncogbound = DClist(idc)%ncogbound + 1
!                   endif
!                   if (cell%DCbound(1) == 0 .and. cell%DCbound(2) /= 0) then
!                       write(*,*) 'DCbound order: ',kcell,cell%DCbound
!                       stop
!                   endif
					if (track_DCvisits .and. ctype == TAGGED_CELL) then
						if (revisit(kcell,idc)) then
							Nrevisits = Nrevisits + 1
							cell%revisits = cell%revisits + 1
						else
							Nvisits = Nvisits + 1
							cell%visits = cell%visits + 1
							cell%ndclist = cell%ndclist + 1
							cell%dclist(cell%ndclist) = idc
						endif
					endif
                    if (nbnd == MAX_DC_BIND) exit
                endif
            enddo
!        endif
!    endif
!enddo
!write(*,*) 'me, nadd, nsub: ',me,nadd,nsub,nadd-nsub
end subroutine