! Implementing a list of cells:
! In this version the cells in the domain are stored in a list, while the
! occupancy array holds the indices of cells in the list.  When a cell
! leaves the domain or dies a gap is created in the list.
! The locations of such gaps are stored in the gaplist, the total number
! of such gaps is ngaps.  A cell entering the domain is allocated an index
! from the tail of this list, if ngaps > 0, or else it is added to the end of the cell list.

module tropho_mod
use global
use behaviour
use packer
use nbr
use Mesh_Generate
use sparse     ! removed for use of the SuperLU_MT solver
use m_unista
use fmotion
use winsock

implicit none

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine rng_initialisation
integer, allocatable :: zig_seed(:)
integer :: i, n, R
integer :: kpar = 0
integer :: npar, grainsize = 32

!do i = 1,8
!    my_seed(i) = i
!enddo
!call random_seed(size = m)
!write(*,*) 'random_number seed size: ',m
!my_seed(1:2) = seed(1:2)
!call random_seed(put=my_seed(1:m))

npar = Mnodes
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
par_zig_init = .true.

end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine omp_initialisation(ok)
logical :: ok
integer :: npr, nth

ok = .true.
if (Mnodes == 1) return
!!DEC$ IF ( DEFINED (_OPENMP) .OR. DEFINED (IBM))
#if defined(OPENMP) || defined(_OPENMP)
write(logmsg,'(a,i2)') 'Requested Mnodes: ',Mnodes
call logger(logmsg)
npr = omp_get_num_procs()
write(logmsg,'(a,i2)') 'Machine processors: ',npr
call logger(logmsg)

nth = omp_get_max_threads()
write(logmsg,'(a,i2)') 'Max threads available: ',nth
call logger(logmsg)
if (nth < Mnodes) then
    Mnodes = nth
    write(logmsg,'(a,i2)') 'Setting Mnodes = max thread count: ',nth
	call logger(logmsg)
endif

call omp_set_num_threads(Mnodes)
!$omp parallel
nth = omp_get_num_threads()
write(logmsg,*) 'Threads, max: ',nth,omp_get_max_threads()
call logger(logmsg)
!$omp end parallel
#endif
call logger('did omp_initialisation')
end subroutine

!-----------------------------------------------------------------------------------------
! Geometry:
! NX = tube length
! NY = NZ = 2*(tube radius) in units of DELTA_X (= cell diameter = 20 um)
!-----------------------------------------------------------------------------------------
subroutine array_initialisation(ok)
logical :: ok
integer :: x,y,z,k
integer :: MAXX

ok = .false.
call rng_initialisation

nsteps_per_min = 1.0/DELTA_T
ngaps = 0
nlist = 0

if (allocated(perm_index)) deallocate(perm_index)

allocate(cell_list(MAX_NLIST))
allocate(gaplist(MAX_NGAPS))
allocate(perm_index(MAX_NLIST))

lastID = 0
k_nonrandom = 0

ok = .true.

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function hill(x,b,n)
real(REAL_KIND) :: x, b
integer :: n
hill = x**n/(x**n + b**n)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine testqsort(n)
integer :: n
integer :: i
integer :: kpar = 0
real(REAL_KIND), allocatable :: a(:)
integer, allocatable :: t(:)

call rng_initialisation

allocate(a(n))
allocate(t(n))
do i = 1,n
    a(i) = par_uni(kpar)
    t(i) = i
enddo
call qsort(a,n,t)
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine testrnor
integer :: n = 100000
integer :: k, j
integer :: kpar = 0
real(REAL_KIND) :: r, rmin=1.0e10, rmax = -1.0e10

do k = 1,n
    do j = 1,n
        r = par_rnor(kpar)
        rmin = min(r,rmin)
        rmax = max(r,rmax)
    enddo
    write(*,'(i12,2e12.4)') k,rmin,rmax
enddo
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine save_positions
integer :: kcell, site(3)
real(REAL_KIND) :: centre(3)

write(nflog,*) 'n_cell_positions: ',n_cell_positions
write(nflog,'(i6,$)') istep
do kcell = 1,n_cell_positions
    centre = cell_list(kcell)%centre(:,1)
    write(nflog,'(3f10.3,$)') centre(1:3)
enddo
write(nflog,*)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_cum_prob
real(REAL_KIND) :: m = 30, s = 2.0
real(REAL_KIND) :: p1, p2, a
integer :: i

p1 = log(m)
p2 = log(s)
do i = 1,30
    a = 2*i
    write(*,'(i4,2f8.3)') i,a,1-cum_prob_lognormal(a,p1,p2)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Using the complete list of cells, cell_list(), extract info about the current state of the
! paracortex.  This info must be supplemented by counts of cells that have died and cells that
! have returned to the circulation.
! We now store stim() and IL2sig() for cells in the periphery.
!-----------------------------------------------------------------------------------------
subroutine get_summary(summaryData) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*)
real(REAL_KIND) :: tnow
logical :: ok

if (.not.use_TCP) then
endif

tnow = istep*DELTA_T
summaryData(1:3) = [ int(tnow/60), istep, Ncells ]
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine write_header
write(nfout,'(a)') '========================================================================================='
write(nfout,*)
write(nfout,'(a10,$)') '      Hour'
write(nfout,'(a10,$)') '  Timestep'
write(nfout,'(a10,$)') '  N_Tcells'
write(nfout,*)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine copyname(name,name_array,n)
character*(*) :: name
character :: name_array(*)
integer :: n
integer :: k

do k = 1,n
	name_array(k) = name(k:k)
enddo
end subroutine


!--------------------------------------------------------------------------------
! Pass a list of cell positions and associated data
! Position is passes as an integer (um)
!--------------------------------------------------------------------------------
subroutine get_scene(nTC_list,TC_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene
use, intrinsic :: iso_c_binding
integer(c_int) :: nTC_list, TC_list(*)
integer :: k, kcell, j, site(3)
integer :: itcstate, ctype

! Cell section
k = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state /= ALIVE) cycle
	k = k+1
	j = 5*(k-1)
	site = cell_list(kcell)%centre(:,1) + 0.5
	itcstate = 0
	TC_list(j+1) = kcell-1
	TC_list(j+2:j+4) = site
	TC_list(j+5) = itcstate
enddo
nTC_list = k
end subroutine

!-----------------------------------------------------------------------------------------
! Now this is used only to set use_TCP = .false.
! The lines
!    call get_command (b, len, status)
!    call get_command_argument (0, c, len, status)
! were failing with gfortran (don't know why), but in any case there was no need to
! get the command line arguments in this way.
!-----------------------------------------------------------------------------------------
subroutine process_command_line(ncpu,infile,outfile)
!DEC$ ATTRIBUTES DLLEXPORT :: process_command_line
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"PROCESS_COMMAND_LINE" :: process_command_line
integer :: i, cnt, len, status
integer :: ncpu
character :: c*(64), b*(256)
character*(64) :: infile,outfile
character*(64) :: progname

!write(*,*) 'process_command_line'
use_TCP = .false.   ! because this is called from para_main()							! --> use_TCP

return

ncpu = 3
!infile = 'omp_para.inp'
!outfile = 'omp_para.out'

call get_command (b, len, status)
if (status .ne. 0) then
    write (logmsg,'(a,i4)') 'get_command failed with status = ', status
    call logger(logmsg)
    stop
end if
call logger('command: ')
call logger(b)
c = ''
call get_command_argument (0, c, len, status)
if (status .ne. 0) then
    write (*,*) 'Getting command name failed with status = ', status
    write(*,*) c
    stop
end if
progname = c(1:len)
cnt = command_argument_count ()
if (cnt < 1) then
    write(*,*) 'Use: ',trim(progname),' num_cpu'
    stop
endif

do i = 1, cnt
    call get_command_argument (i, c, len, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        stop
    end if
    if (i == 1) then
!        read(c(1:len),'(i)') ncpu
        read(c(1:len),*) ncpu															! --> ncpu
        write(*,*) 'Requested threads: ',ncpu
    elseif (i == 2) then
        infile = c(1:len)																! --> infile
        write(*,*) 'Input file: ',infile
    elseif (i == 3) then
        outfile = c(1:len)																! --> outfile
        write(*,*) 'Output file: ',outfile
    endif
end do

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step
use, intrinsic :: iso_c_binding
integer(c_int) :: res
integer :: hour, nit, nt_hour, kpar=0
real(REAL_KIND) :: tnow, dt
logical :: ok, done, changed

res = 0
if (Ncells == 0) then
	call logger('Cells all gone')
	res = 1
	return
endif
if (nocells) then
	write(*,*) 'SIMULATING NOCELLS CASE TO TEST VELOCITY FIELD GENERATION!!!'
	write(nflog,*) 'SIMULATING NOCELLS CASE TO TEST VELOCITY FIELD GENERATION!!!'
	nlist = 0
	ncells = 0
endif

dbug = .false.
nt_hour = 3600/DELTA_T
ok = .true.
istep = istep + 1

if (n_cell_positions > 0) then
	call save_positions
endif
tnow = istep*DELTA_T
use_settling = (settle_hrs > 0)
if (use_settling) then
	settling = (tnow < settle_hrs*3600)
endif
if (mod(istep,nt_hour) == 0) then
	write(logmsg,'(a,i6,a,i6,a,f8.1)') 'simulate_step: ',istep, '   Ncells: ',Ncells, '  hour: ',tnow/3600
	call logger(logmsg)
    if (TAGGED_LOG_PATHS) then
		call add_log_paths
	endif
endif
if (TAGGED_LOG_PATHS .and. mod(istep,1) == 0) then
	call update_log_paths
endif

call make_perm_index(ok)
if (.not.ok) then
	call logger('make_perm_index error')
	res = 4
	return
endif

if (.not.calibration_run) then
	call FEsolve
	call logger('Done FEsolve')
!	call test_getPointVel
!	call makeVelDist(500.d0)
	call GetVel
endif
call update_all_nbrlists

t_fmover = 0
nit = 0
done = .false.
do while (.not.done)
	nit = nit + 1
	call fmover(dt,done,ok)
	if (.not.ok) then
		call logger('mover error')
		res = 1
		return
	endif
	t_fmover = t_fmover + dt
	changed = .false.
	if (changed) then
		call make_perm_index(ok)
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine connection(awp,port,error)
TYPE(winsockport) :: awp
integer :: port, error
integer :: address = 0
!!!character*(64) :: ip_address = "127.0.0.1"C      ! need a portable way to make a null-terminated C string
character*(64) :: host_name = "localhost"

if (.not.winsock_init(1)) then
    call logger("winsock_init failed")
    stop
endif

awp%handle = 0
awp%host_name = host_name
awp%ip_port = port
awp%protocol = IPPROTO_TCP
call Set_Winsock_Port (awp,error)

if (.not.awp%is_open) then
    write(nflog,*) 'Error: connection: awp not open: ',port
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
call connection(awp_0,TCP_PORT_0,error)
if (awp_0%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_0 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_0%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_0'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_0  '
call logger(logmsg)

if (use_CPORT1) then
	call connection(awp_1,TCP_PORT_1,error)
	if (awp_1%handle < 0 .or. error /= 0) then
		write(logmsg,'(a)') 'TCP connection to TCP_PORT_1 failed'
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (.not.awp_1%is_open) then
		write(logmsg,'(a)') 'No connection to TCP_PORT_1'
		call logger(logmsg)
		ok = .false.
		return
	endif
	write(logmsg,'(a)') 'Connected to TCP_PORT_1  '
	call logger(logmsg)
endif
! Allow time for completion of the connection
call sleeper(2)
end subroutine

!-----------------------------------------------------------------------------------------
! This subroutine is called to initialize a simulation run.
! ncpu = the number of processors to use
! infile = file with the input data
! outfile = file to hold the output
! runfile = file to pass info to the master program (e.g. Python) as the program executes.
!-----------------------------------------------------------------------------------------
subroutine setup(ncpu,infile,outfile,ok)
integer :: ncpu
!character*(*) :: infile, outfile, filename
character*(*) :: infile, outfile
logical :: ok
character*(64) :: msg
integer :: error


ok = .true.
initialized = .false.
par_zig_init = .false.
Mnodes = ncpu
!Mnodes = 1
inputfile = infile
outputfile = outfile
write(logmsg,*) 'ncpu: ',Mnodes
call logger(logmsg)

#if defined(OPENMP) || defined(_OPENMP)
    call logger("OPENMP defined")
    call omp_initialisation(ok)
    if (.not.ok) return
#else
    call logger("OPENMP NOT defined")
    if (Mnodes > 1) then
        write(logmsg,'(a)') 'No OpenMP, using one thread only'
        call logger(logmsg)
        Mnodes = 1
    endif
#endif

call logger("read_cell_params")
call read_cell_params(ok)
if (.not.ok) return
call logger("did read_cell_params")

call setup_force_parameters
call logger("did setup_force_parameters")

call array_initialisation(ok)
if (.not.ok) return
call logger('did array_initialisation')

call PlaceCells(ok)
if (ok) then
	call logger('did PlaceCells: OK')
else
	call logger('did PlaceCells: not OK')
	stop
endif
if (.not.ok) return

call FEsetup
call logger('did FEsetup')
!call test_getPointElement

nwallcells = 0
!call make_wall

d_nbr_limit = 1.5*2*Raverage	!*Rdivide0	! 1.5 is an arbitrary choice - was 1.2
call setup_nbrlists

return

if (save_input) then
    call save_inputfile(inputfile)
!	call save_inputfile(fixedfile)
    call save_parameters
    call write_header
endif
!call chemokine_setup
firstSummary = .true.
initialized = .true.
write(logmsg,'(a,i6)') 'Startup procedures have been executed: initial T cell count: ',Ncells
call logger(logmsg)

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine wrapup
integer :: ic, ierr
logical :: isopen

call logger('doing wrapup ...')
ierr = 0
!if (allocated(Tres_dist)) deallocate(Tres_dist)
if (allocated(cell_list)) deallocate(cell_list,stat=ierr)
if (ierr /= 0) then
    write(*,*) 'cellist deallocate error: ',ierr
    stop
endif
ierr = 0
if (allocated(gaplist)) deallocate(gaplist,stat=ierr)
if (allocated(life_dist)) deallocate(life_dist)
!if (allocated(divide_dist)) deallocate(divide_dist)
if (allocated(chemo_p)) deallocate(chemo_p)

! Close all open files
inquire(unit=nfout,OPENED=isopen)
if (isopen) then
	close(nfout)
	call logger('closed nfout')
endif
inquire(nfres,OPENED=isopen)
if (isopen) close(nfres)
inquire(nftraffic,OPENED=isopen)
if (isopen) close(nftraffic)
inquire(nfchemo,OPENED=isopen)
if (isopen) close(nfchemo)

if (par_zig_init) then
	call par_zigfree
endif

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run
use, intrinsic :: iso_c_binding
integer(c_int) :: res
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

!call SaveGenDist
!if (evaluate_residence_time) then
!	call write_Tres_dist
!endif
if (TAGGED_LOG_PATHS) then
	call write_log_paths
endif
!call write_FACS

call wrapup

if (res >= 0) then
	call logger(' Execution successful!')
else
	call logger('  === Execution failed ===')
	call sleeper(1)
endif

if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    if (use_CPORT1) call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
!	    call logger("closed PORT_0")
		if (use_CPORT1) then
			call winsock_send(awp_1,quit,8,error)
			call winsock_close(awp_1)
!			call logger("closed PORT_1")
		endif
	endif
endif

end subroutine

!-----------------------------------------------------------------------------------------
! This is the DLL procedure that can be called from an external non-Fortran program to
! make a simulation run.
! Called from Python with:
!     mydll.EXECUTE(byref(ncpu),infile,n1,outfile,n2,resfile,n3,runfile,n4)
! Note that the arguments n1,n2,n3,n4, the lengths of the filename strings, are
! hidden arguments (they are not explicitly in the Fortran subroutine argument list).
! Every Python string needs to be followed by a hidden length parameter, and unless
! the declared length of a Fortran string matches the actual length of that passed from
! Python, the form character*(*) must be used.
!-----------------------------------------------------------------------------------------
subroutine execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128), outfile_array(128)
integer(c_int) :: ncpu, inbuflen, outbuflen
character*(128) :: infile, outfile
logical :: ok, success, isopen
integer :: i, res

use_CPORT1 = .false.	! DIRECT CALLING FROM Fortran, C++
infile = ''
do i = 1,inbuflen
	infile(i:i) = infile_array(i)
enddo
outfile = ''
do i = 1,outbuflen
	outfile(i:i) = outfile_array(i)
enddo

inquire(unit=nflog,OPENED=isopen)
inquire(unit=nflog2,OPENED=isopen)
if (.not.isopen) then
    open(nflog,file='tropho.log',status='replace')
    open(nflog2,file='flow.log',status='replace')
endif
awp_0%is_open = .false.
awp_1%is_open = .false.

#ifdef GFORTRAN
    write(logmsg,'(a)') 'Built with GFORTRAN'
	call logger(logmsg)
#endif

logmsg = 'OS??'
#ifdef LINUX
    write(logmsg,'(a)') 'OS is Linux'
#endif
#ifdef OSX
    write(logmsg,'(a)') 'OS is OS-X'
#endif
#ifdef _WIN32
    write(logmsg,'(a)') 'OS is Windows'
#endif
#ifdef WINDOWS
    write(logmsg,'(a)') 'OS is Windows'
#endif
call logger(logmsg)

!#ifdef OPENMP
#if defined(OPENMP) || defined(_OPENMP)
    write(logmsg,'(a)') 'Executing with OpenMP'
	call logger(logmsg)
#endif

write(logmsg,*) 'inputfile:  ', infile
call logger(logmsg)
write(logmsg,*) 'outputfile: ', outfile
call logger(logmsg)
if (use_tcp) then
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
endif
call setup(ncpu,infile,outfile,ok)
if (ok) then
	clear_to_send = .true.
	simulation_start = .true.
	istep = 0
	res = 0
else
	call logger('=== Setup failed ===')
	res = 1
	stop
endif
return

end subroutine

!------------------------------------------------------
!------------------------------------------------------
subroutine FEsetup

integer :: i, j, k, m, e, count1,row_size,sample_size,kcell
integer :: count_empty, ND, Nvert,II(6), signum(6)
integer, parameter :: out_unit=20
real (REAL_KIND), allocatable :: vertices(:,:)
integer, allocatable :: InternalFaces(:,:), INLET(:,:)
integer, allocatable :: OUTLET(:,:), WALLs(:,:)
integer, allocatable :: ELEMENT(:,:), ELEMENT2(:,:)
integer, allocatable :: n2f(:,:), e2f(:,:),f2e(:,:)
real (REAL_KIND), allocatable :: NODES(:,:)
integer, allocatable :: All_Surfaces(:,:)

character * ( 255 ) :: n2f_filename, e2f_filename, f2e_filename
character * ( 255 ) :: CylinCube_filename, Centroid_filename
character * ( 255 ) :: fileplace

integer, allocatable :: element2face(:,:)
integer, allocatable :: face2element(:,:)
integer, allocatable :: nodes2face(:,:)
!integer, allocatable :: CylinCub_list(:,:)

!real (REAL_KIND), allocatable :: Centroids(:,:)
real(REAL_KIND), allocatable :: BB(:)
real(REAL_KIND), allocatable :: CC(:)
integer, allocatable :: IB(:), JB(:)
integer, allocatable :: IC(:), JC(:)
!real (REAL_KIND), allocatable :: B_Matrix(:,:)

integer :: nel,nel_z !number of elements

logical :: answer1, answer2
real (REAL_KIND) :: vertix_vec(8)

!---------------------------------
! We generate mesh nodes here
! and set boundary conditions
!---------------------------------
!************** here we read in mesh properties form .msh file and .k files ************
	call readmshfile(vertices,InternalFaces,INLET,OUTLET,WALLs)
	call readKfile(NODES,ELEMENT)
	nel = size(ELEMENT,1)

	allocate(ELEMENT2(size(ELEMENT,1),8))
	do e=1,size(ELEMENT,1)
		if (e<=1800 .OR. e>2700) then
			ELEMENT2(e,1)=ELEMENT(e,1)
			ELEMENT2(e,2)=ELEMENT(e,5)
			ELEMENT2(e,3)=ELEMENT(e,6)
			ELEMENT2(e,4)=ELEMENT(e,2)
			ELEMENT2(e,5)=ELEMENT(e,4)
			ELEMENT2(e,6)=ELEMENT(e,8)
			ELEMENT2(e,7)=ELEMENT(e,7)
			ELEMENT2(e,8)=ELEMENT(e,3)
		else
			ELEMENT2(e,1)=ELEMENT(e,1)
			ELEMENT2(e,2)=ELEMENT(e,2)
			ELEMENT2(e,3)=ELEMENT(e,3)
			ELEMENT2(e,4)=ELEMENT(e,4)
			ELEMENT2(e,5)=ELEMENT(e,5)
			ELEMENT2(e,6)=ELEMENT(e,6)
			ELEMENT2(e,7)=ELEMENT(e,7)
			ELEMENT2(e,8)=ELEMENT(e,8)
		endif
	enddo
	ELEMENT=ELEMENT2
	!print *, "number of elements=", nel
	row_size=size(InternalFaces,1)+size(INLET,1)+&
					size(OUTLET,1)+size(WALLs,1)
	allocate(All_Surfaces(row_size,6))
	do i=1,size(InternalFaces,1)
		All_Surfaces(i,:)= InternalFaces(i,:)
	enddo
	i=size(InternalFaces,1)
	do j=1,size(INLET,1)
		All_Surfaces(i+j,:)= INLET(j,:)
	enddo
	j=size(INLET,1)
	do k=1,size(OUTLET,1)
		All_Surfaces(i+j+k,:)= OUTLET(k,:)
	enddo
	k=size(OUTLET,1)
	do m=1,size(WALLs,1)
		All_Surfaces(i+j+k+m,:)= WALLs(m,:)
	enddo
	nofaces = size(All_Surfaces,1)

	!do i=1,row_size
     !   write(nflog, '(a,6i6)') 'all_faces ', All_Surfaces(i,:)
    !enddo

	allocate(nodes2face(6*nel,5))
	allocate(element2face(nel,6))
	allocate(face2element(nofaces,6))
	!call face(NODES,ELEMENT,All_Surfaces,,nodes2face,element2face,face2element)
										   !this part takes a lot of time to run so we ran it
										   !in MATLAB and saved the outputs in .txt files and
										   !read them here to use the data later in the code:
	n2f_filename='nodes2face.txt'
	!print *, "nel*6", nel*6
	j=5 !number of columns should be read from the n2f file
	call Read_face(n2f_filename,nel*6,j,n2f)
	nodes2face=n2f
	!do i=1,6*nel
	!	print *, "n2f", nodes2face(i,:)
	!enddo

	e2f_filename='element2face.txt'
	j=6 !number of columns should be read from the n2f file
	call Read_face(e2f_filename,nel,j,e2f)
	element2face=e2f
	!print *, "e2f", element2face(nel,:)

	f2e_filename='face2element.txt'
	j=6 !number of columns should be read from the n2f file
	call Read_face(f2e_filename,nofaces,j,f2e)
	face2element=f2e
	CommonFaceNo=0
	do j = 1,size(face2element,1)
		if (face2element(j,6)/= 0) then
			CommonFaceNo=CommonFaceNo+1
		endif
	enddo
	!print *, "CommonFaceNo", CommonFaceNo
	!print *, "f2e", face2element(1,:)

	if (allocated(DCInlet)) deallocate(DCInlet)
	if (allocated(DCOutlet)) deallocate(DCOutlet)
	if (allocated(NCWall)) deallocate(NCWall)
	if (allocated(NCInlet)) deallocate(NCInlet)
!set boundary condition
	!1)set up dirichlet surfaces
	allocate(DCInlet(size(INLET,1)))
	allocate(DCOutlet(size(OUTLET,1)))
	!2)set up Neumann surfaces - normal flux in/out = neumann wall
	allocate(NCWall(size(WALLs,1)))
	allocate(NCInlet(size(INLET,1)))
	call set_boundary_condition(nel, INLET,OUTLET,WALLs,All_Surfaces)	!, &
!	DCInlet,DCOutlet, NCWall, NCInlet)
	!print *, "inlet face=", DCInlet
	!print *, "outlet face=", DCOutlet
	!write (*,*) "wall face=", NCWall


	CylinCube_filename='CylinCub_list.txt'
	call Read_face(CylinCube_filename,nel,2,CylinCub_list)

	Centroid_filename='Centroid.txt'
	fileplace= ""
	call Read_real(fileplace,Centroid_filename,nel,3,Centroids)
	!print *, "Centroids=", Centroids(1,:)

	call unique(REAL(NINT(vertices(:,3)), REAL_KIND),e_Z)
	!print *, "e_Z=", e_Z(:)
	nel_z=size(e_Z)-1
	if (allocated(layered_el)) deallocate(layered_el)
	allocate(layered_el(nel_z,INT(nel/nel_z)))
	do k=2,nel_z+1
		j=0
		do e=1,size(ELEMENT,1)
			do i=1,8
				vertix_vec(i)=REAL(NINT(vertices(ELEMENT(e,i),3)),REAL_KIND)
			enddo
		answer1=isempty_find(vertix_vec,e_Z(k))
		answer2=isempty_find(vertix_vec,e_Z(k-1))
			if ( .not.answer1 .AND. .not.answer2 ) then
				j=j+1
				layered_el(k-1,j)=e
			endif
		enddo
	enddo
	!print *, "el=", el(1,:)

if (allocated(B_matrix)) deallocate(B_matrix)
if (allocated(ElToFace)) deallocate(ElToFace)
if (allocated(AllSurfs)) deallocate(AllSurfs)
if (allocated(Cylinder_element)) deallocate(Cylinder_element)
if (allocated(Cylinder_vertices)) deallocate(Cylinder_vertices)
if (allocated(u_cell_x)) deallocate(u_cell_x)
if (allocated(u_cell_y)) deallocate(u_cell_y)
if (allocated(u_cell_z)) deallocate(u_cell_z)

allocate(B_MATRIX(nel))
call READ_Bmatrix(nel)

allocate(ElToFace(nel,6))
ElToFace=element2face
!allocate(FaToEl(nofaces,6))
!FaToEl=face2element
allocate(AllSurfs(row_size,6))
AllSurfs=All_Surfaces
NofElements=nel
NofElements_z=nel_z
allocate(Cylinder_ELEMENT(nel,8))
Cylinder_ELEMENT=ELEMENT
allocate(Cylinder_vertices(size(vertices,1),3))
Cylinder_vertices=vertices

allocate(u_cell_x(ncells))
allocate(u_cell_y(ncells))
allocate(u_cell_z(ncells))
u_cell_x(:)=0.0D0
u_cell_y(:)=0.0D0
u_cell_z(:)=0.0D0
!write(*,*) 'initial cell velocityz:  ', u_cell_z(1)
end subroutine

end module

