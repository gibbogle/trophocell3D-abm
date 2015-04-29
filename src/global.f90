!==========================================================================================
!==========================================================================================

module global

use omp_lib
!use CD69
!use IL2
!use IL7
!use IL_dummy
use par_zig_mod
use winsock

implicit none

!INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND( 12, 60 )
integer, parameter :: SP = kind(1.0), DP = kind(1.0d0)
integer, parameter :: REAL_KIND = SP

! Files
integer, parameter :: nfcell = 10, nfout = 11, nfvec = 12, nfpath = 13, nfres = 14, nftraffic = 16, nfrun = 17, &
					  nftravel=18, nfcmgui=19, nfpos=20, nflog=21, nfchemo=22, nffacs = 24

! General parameters
integer, parameter :: BIG_INT = 2**30
integer, parameter :: NEUMANN_MODEL = 1
integer, parameter :: MOORE18_MODEL = 2
integer, parameter :: MOORE26_MODEL = 3
integer, parameter :: MODEL = MOORE18_MODEL
integer, parameter :: MAXRELDIR = 26
integer, parameter :: MAXRELDIR2D = 8

integer, parameter :: OUTSIDE_TAG = -99
integer, parameter :: TAGGED_CELL   = 99
!integer, parameter :: FINISHED = 10

integer, parameter :: TCP_PORT_0 = 5000		! main communication port (logging)
integer, parameter :: TCP_PORT_1 = 5001		! data transfer port (plotting)
integer, parameter :: TCP_PORT_2 = 5002
integer, parameter :: TCP_PORT_3 = 5003

!integer, parameter :: STAGE_BYTE = 1
!integer, parameter :: GENERATION_BYTE = 2
!integer, parameter :: ACTIVATION_BYTE = 3
integer, parameter :: SLOT_NUM1 = 1
integer, parameter :: SLOT_NUM2 = 2
integer, parameter :: BOTH = 3

integer, parameter :: CHEMO_1 = 1
integer, parameter :: CHEMO_2 = 2

integer, parameter :: NCTYPES = 1
integer, parameter :: TROPHO_CELL = 1

integer, parameter :: neumann(3,6) = reshape((/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /), (/3,6/))
real(REAL_KIND), parameter :: jumpvec2D(3,8) = reshape((/ 1,0,0, 1,1,0, 0,1,0, -1,1,0, -1,0,0, -1,-1,0, 0,-1,0, 1,-1,0 /), (/3,8/))

real(REAL_KIND), parameter :: BIG_TIME = 100000
logical, parameter :: use_add_count = .true.    ! keep count of sites to add/remove, do the adjustment at regular intervals
logical, parameter :: save_input = .true.

integer, parameter :: NZ_2D = 1

! Diffusion parameters
logical, parameter :: use_cytokines = .false.       ! to use IL-2
logical, parameter :: use_diffusion = .false.		! For cytokines
integer, parameter :: NDIFFSTEPS = 6    ! divisions of DELTA_T for diffusion computation

! Chemokine/receptor parameters
integer, parameter :: MAX_CHEMO = 2
integer, parameter :: MAX_RECEPTOR = 2
logical, parameter :: USE_GENERAL_CODE = .true.	! this code covers dynamic and SS solvers, conc and secretion b.c.s
logical, parameter :: use_ODE_diffusion = .false.	! for general case, use ODE system (dynamic solver)
logical, parameter :: USE_ORIGINAL_CODE = .not.USE_GENERAL_CODE	! simple case, fixed secretion into DC neighborhood

real(REAL_KIND), parameter :: CFSE_std = 0.05


! Data above this line almost never change
!==============================================================================================================

! Run parameters

! Parameters and switches for calibration etc.
logical, parameter :: calibrate_motility = .false.
logical, parameter :: motility_param_range = .false.
logical, parameter :: motility_save_paths = .false.
logical, parameter :: calibrate_diffusion = .false.
integer, parameter :: n_multiple_runs = 1

! To investigate how chemotaxis influences T cell paths.
! Set up a single exit at the centre of the blob.  Tag cells starting at a specified distance from the exit,
! and log their paths until they exit.  Cells are either chemotactic or not.
logical, parameter :: TAGGED_LOG_PATHS = .false.
integer, parameter :: MAX_LOG_PATHS = 100		! number of paths to track for non-chemo (1) and chemo (2) cells
integer, parameter :: MAX_LOG_PATH_SITES = 1000		! max number of tracking start locations
integer :: log_path_site(3,MAX_LOG_PATH_SITES)		! list of tracking start locations
integer :: n_log_path_sites						! number of tracking start locations
integer :: n_log_path(2)						! current numbers of cells being tracked
real(REAL_KIND), parameter :: R_LOG_PATH = 7				! distance of start location from exit
real(REAL_KIND), parameter :: LOG_PATH_FACTOR = 0.11		! inflow factor for a single exit (depends on blob size, 18.1 -> 0.11)
integer, parameter :: MAX_PATH_STEPS = 10000	! maximum number of steps to log on a path
type path_type
	integer :: kcell					! cell number
	integer :: id						! cell ID
	integer :: kstep					! number of last location logged (kstep = 1 is the start location)
	logical :: in						! flag for cell in the blob
	integer :: pos(3,MAX_PATH_STEPS)	! sequence of cell locations
end type
type (path_type) :: log_path(MAX_LOG_PATHS,2)

! Parameters for controlling data capture for graphical purposes etc.
integer, parameter :: save_interval_hours = 48
real(REAL_KIND), parameter :: save_length_hours = 0.5      ! 30 minutes
integer, parameter :: ntaglimit_base = 200000
integer, parameter :: ntres = 60    ! 15 min
logical, parameter :: log_results = .false.

! Type definitions
integer :: NTcells0
integer :: NTcells
integer :: Nsites
real (REAL_KIND) :: InflowTotal

type cell_type
    integer :: ID
    integer :: site(3)
    integer :: step
	real(REAL_KIND) :: receptor_level(MAX_RECEPTOR)
	real(REAL_KIND) :: receptor_saturation_time(MAX_RECEPTOR)
    integer :: tag
    real(REAL_KIND) :: entrytime       ! time that the cell entered the paracortex (by HEV or cell division)
!    type(cog_type),    pointer :: cptr    ! because NULL is used by winsock (from ifwinty).  NULLIFY() instead.
    integer(2) :: ctype
	integer(2) :: lastdir
	integer :: dtotal(3)
end type

type occupancy_type
    integer :: indx(2)
!    real(REAL_KIND) :: chemo_conc		! chemokine concentration
!    real(REAL_KIND) :: chemo_grad(3)	! chemokine gradient vector
    integer :: isrc						! index into source list
end type

type dist_type
	integer :: dclass
	real(REAL_KIND) :: p1, p2, p3
end type

type counter_type
    logical :: logscale
    integer :: nsamples
    integer :: nbins
    real(REAL_KIND) :: binmin, binstep
    integer, allocatable :: bincount(:)
    real(REAL_KIND) :: total
end type

type source_type
	integer :: bc
	integer :: site(3)
	real(REAL_KIND) :: level(MAX_CHEMO)
end type

real(REAL_KIND) :: TC_RADIUS

real(REAL_KIND) :: GAMMA                          ! controls crowding
real(REAL_KIND) :: BETA                            ! speed: 0 < beta < 1
real(REAL_KIND) :: RHO                             ! persistence: 0 < rho < 1

logical :: use_traffic = .false.

real(REAL_KIND) :: days                            ! number of days to simulate
integer :: seed(2)                      ! seed vector for the RNGs
integer :: FACS_INTERVAL				! interval between FACs plot outputs (h)
!integer :: SPECIES						! animal species source of T cells
!real(REAL_KIND) :: IV_WELL_DIAMETER				! diameter of in vitro well (mm)
!integer :: IV_NTCELLS					! initial T cell population in vitro
!real(REAL_KIND) :: IV_COGNATE_FRACTION				! fraction of in vitro cells that are cognate for DC antigen
!logical :: IV_SHOW_NONCOGNATE			! display non-cognate T cells
!character*(128) :: fixedfile

!---------------------------------------------------
! end of parameters to read from input file
!---------------------------------------------------

!-------------------------------------------------------
! More input parameters to be read from fixed input file
!-------------------------------------------------------

!---------------------------------------------------------
! end of more parameters to be read from fixed input file
!---------------------------------------------------------

! Geometry data
logical, parameter :: SIMULATE_2D = .false.
integer :: NX, NY, NZ, Ncells
integer :: nlength, nradius, nplug
!integer, allocatable :: xdomain(:),xoffset(:),zdomain(:),zoffset(:)
!integer :: blobrange(3,2)
real(REAL_KIND) :: Radius0
real(REAL_KIND) :: DELTA_X, PI
real(REAL_KIND) :: TagRadius
real(REAL_KIND) :: x0,y0,z0   ! centre in global coordinates (units = grids)
real(REAL_KIND) :: Centre(3)
real(REAL_KIND) :: Vc, Ve

! Motility data
real(REAL_KIND) :: DELTA_T       ! minutes
integer :: nreldir, njumpdirs
integer :: jumpvec(3,27)    ! 14 is no-jump case (0,0,0)
real(REAL_KIND) :: jumpdist(27)
!integer :: jumpvec2D(3,8)
integer :: reldir(6,MAXRELDIR)
real(REAL_KIND) :: dirprob(0:MAXRELDIR)

integer :: nreldir2D, njumpdirs2D
integer :: reldir2D(8,8)
real(REAL_KIND) :: dirprob2D(0:8)
logical :: diagonal_jumps
logical :: use_wrapping
integer :: n_cell_positions

real(REAL_KIND) :: Vmax
real(REAL_KIND) :: chemo_attraction
real(REAL_KIND) :: Kdrag, Kadhesion, Kstay

! Chemotaxis data
integer, parameter :: chemo_N = 2
real(REAL_KIND), allocatable :: chemo_p(:,:,:,:)
real(REAL_KIND) :: grad_amp(2), grad_dir(2,3)
! Background flow
real(REAL_KIND) :: BG_flow_amp

! Cell data
type(occupancy_type), allocatable :: occupancy(:,:,:)
type(cell_type), allocatable, target :: cell_list(:)
type(source_type), allocatable :: sourcelist(:)
integer :: nsources
integer, allocatable :: gaplist(:)
!integer, allocatable :: zrange2D(:,:,:)	! for 2D case
integer :: lastID, nlist, n2Dsites, ngaps, ntaglimit, ntagged=0, ntagged_left
integer :: lastNTcells, k_nonrandom(2)
integer :: max_nlist, max_ngaps
integer :: nleft
integer :: nadd_sites, ndivisions
real(REAL_KIND) :: lastbalancetime

integer :: noutflow_tag, ninflow_tag
logical :: firstSummary

! Miscellaneous data
logical :: initialized, steadystate
integer :: total_in, total_out
integer :: Nsteps, nsteps_per_min, istep
integer :: Mnodes
integer :: IDtest

character*(128) :: inputfile
character*(128) :: outputfile
character*(2048) :: logmsg
TYPE(winsockport) :: awp_0, awp_1, awp_2, awp_3
logical :: use_TCP = .true.         ! turned off in para_main()
logical :: use_CPORT1
logical :: stopped, clear_to_send, simulation_start, par_zig_init

logical :: dbug = .false.

!!DEC$ ATTRIBUTES DLLEXPORT :: ntravel, N_TRAVEL_COG, N_TRAVEL_DC, N_TRAVEL_DIST, k_travel_cog, k_travel_dc
!!DEC$ ATTRIBUTES DLLEXPORT :: travel_dc, travel_cog, travel_dist
!DEC$ ATTRIBUTES DLLEXPORT :: nsteps	!istep
contains

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
integer :: n1,n2,kpar
integer :: k,R

!k = irand()     ! intrinsic
if (n1 == n2) then
    random_int = n1
elseif (n1 > n2) then
    write(logmsg,*) 'ERROR: random_int: n1 > n2: ',n1,n2
    call logger(logmsg)
    stop
endif
R = par_shr3(kpar)
if (R < -2147483647) R = par_shr3(kpar)
k = abs(R)
random_int = n1 + mod(k,(n2-n1+1))
end function

!--------------------------------------------------------------------------------
! Returns a permutation of the elements of a()
!--------------------------------------------------------------------------------
subroutine permute(a,n,kpar)
integer :: a(*),n,kpar
integer :: i,k,tmp

do i = 1,n
    k = random_int(1,n,kpar)
	tmp = a(i)
	a(i) = a(k)
	a(k) = tmp
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine waste_time(n,dummy)
integer :: k, n
real(REAL_KIND) :: dummy
real(REAL_KIND) :: rsum,R
integer :: kpar=0

rsum = 0
do k = 1,n
!    call random_number(R)
    R = par_uni(kpar)
    rsum = rsum + R
enddo
dummy = rsum
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function norm(r)
real(REAL_KIND) :: r(3)

norm = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function norm2(r)
real(REAL_KIND) :: r(3)

norm2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine normalize(r)
real(REAL_KIND) :: r(3)

r = r/norm(r)
end subroutine

!-----------------------------------------------------------------------------------------
! Distance from the blob centre (units = grids)
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function cdistance(site)
integer :: site(3)
real(REAL_KIND) :: r(3)

r = site - Centre
cdistance = norm(r)
end function


!--------------------------------------------------------------------------------
! A site is taggable if it is less than a specified distance TagRadius from
! the centre.
!--------------------------------------------------------------------------------
logical function taggable(site)
integer :: site(3)
real(REAL_KIND) :: d, r(3)

r = site - Centre
d = norm(r)
if (d <= TagRadius) then
    taggable = .true.
else
    taggable = .false.
endif
end function

!----------------------------------------------------------------------------------------
! Check to see if (x,y,z) is outside the grid
!----------------------------------------------------------------------------------------
logical function outside_xyz(x,y,z)
integer :: x, y, z
outside_xyz = .true.
if (x < 1 .or. x > NX) return
if (y < 1 .or. y > NY) return
if (z < 1 .or. z > NZ) return
outside_xyz = .false.
end function

!-----------------------------------------------------------------------------------------
! Squeeze gaps out of cellist array, adjusting occupancy array.
!-----------------------------------------------------------------------------------------
subroutine squeezer(force)
logical :: force
integer :: last, kcell, site(3), indx(2), i, j, idc, n, region

!write(*,*) 'squeezer'
if (ngaps == 0) return
if (.not.force .and. (ngaps < max_ngaps/2)) return
if (dbug) write(nflog,*) 'squeezer: ',ngaps,max_ngaps,nlist

last = nlist
kcell = 0
n = 0
do
    kcell = kcell+1
    if (cell_list(kcell)%ID == 0) then    ! a gap
        if (kcell == last) exit
        do
            if (last == 0) then
                write(nflog,*) 'last = 0: kcell: ',kcell
                stop
            endif
            if (cell_list(last)%ID == 0) then
                last = last-1
                n = n+1
                if (n == ngaps) exit
            else
                exit
            endif
        enddo
        if (n == ngaps) exit
        call copycell2cell(cell_list(last),cell_list(kcell),kcell)
!        cell_list(kcell) = cell_list(last)
! Note that the DCs that were in cell_list(last)%DCbound(:), which are now in cell_list(kcell)%DCbound(:),
! have in the DClist(:)%cogbound(:) the old cell index (last), but should have the new index (kcell)
! Ths is only an issue for a cognate cell
        last = last-1
        n = n+1
    endif
    if (n == ngaps) exit
enddo
nlist = nlist - ngaps
ngaps = 0
if (dbug) write(nflog,*) 'squeezed: ',n,nlist

end subroutine

!-----------------------------------------------------------------------------------------
! Copy the contents of cell_list(kfrom) to the entry cell_list(kto)
! Need to fix up the corresponding cognate_list() entry as well.
! Note that:
!   cell = cell_list(kcell)
!   kcog = cell%cptr%cogID
!   k = cognate_list(kcog)
! =>
!   k = kcell
!-----------------------------------------------------------------------------------------
subroutine copycell2cell(cell_from,cell_to,kcell)
integer :: kcell
type(cell_type) :: cell_from, cell_to
integer :: ctype, stype, kcog

!write(*,*) 'copycell2cell: ',kcell
ctype = cell_from%ctype
!if (associated(cell_from%cptr)) then
!	stype = COG_TYPE_TAG
!else
!	stype = NONCOG_TYPE_TAG
!endif
!if (stype == NONCOG_TYPE_TAG .and. associated(cell_to%cptr)) then
!    deallocate(cell_to%cptr)
!endif
!if (stype == COG_TYPE_TAG) then
!    if (.not.associated(cell_to%cptr)) then
!        allocate(cell_to%cptr)
!    endif
!    cell_to%cptr = cell_from%cptr
!    kcog = cell_to%cptr%cogID
!    cognate_list(kcog) = kcell
!elseif (stype /= NONCOG_TYPE_TAG) then
!    write(*,*) 'ERROR: copycell2cell: istep, ID, ctype, stype: ',istep,cell_from%ID,ctype,stype
!    stop
!endif
cell_to%ID = cell_from%ID
cell_to%site = cell_from%site
cell_to%receptor_level = cell_from%receptor_level
cell_to%receptor_saturation_time = cell_from%receptor_saturation_time
cell_to%tag = cell_from%tag
cell_to%ctype = cell_from%ctype
cell_to%lastdir = cell_from%lastdir
cell_to%entrytime = cell_from%entrytime
if (cell_from%ctype == 0) then
    write(*,*) 'ERROR: copycell2cell: ctype = 0'
    stop
endif

end subroutine

!--------------------------------------------------------------------------------
! Returns:
! 0 if no slots are occupied
! 1 if slot 1 is occupied
! 2 if slot 2 is occupied
! 3 if both slots are occupied
!--------------------------------------------------------------------------------
integer function getslots(site)
integer :: site(3)
integer :: k

getslots = 0
do k = 1,2
    if (occupancy(site(1),site(2),site(3))%indx(k) > 0) then
        getslots = getslots + k
    endif
enddo
end function

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
integer function neighbourhoodCount(site)
integer :: site(3)
integer :: site2(3), k, count

count = 0
do k = 1,27
	site2 = site + jumpvec(:,k)
	if (site2(1) < 1 .or. site2(1) > NX) cycle
	if (site2(2) < 1 .or. site2(2) > NY) cycle
	if (site2(3) < 1 .or. site2(3) > NZ) cycle
	if (occupancy(site2(1),site2(2),site2(3))%indx(1) >= 0) then
		count = count + 1
	endif
enddo
neighbourhoodCount = count
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
integer function cell_count()
integer :: kcell, ntot
ntot = 0
do kcell = 1,nlist
    if (cell_list(kcell)%ID == 0) cycle             ! skip gaps in the list
    ntot = ntot + 1
enddo
cell_count = ntot
end function

!--------------------------------------------------------------------------------
! Make a random choice of an integer from 1 - N on the basis of probabilities in
! the array p(:) (assumed to be normalized).
!--------------------------------------------------------------------------------
integer function random_choice(p,N,kpar)
integer :: N,kpar
real(REAL_KIND) :: p(:)
integer :: k
real(REAL_KIND) :: R, psum

!call random_number(R)
R = par_uni(kpar)
psum = 0
do k = 1,N
    psum = psum + p(k)
    if (R <= psum) then
        random_choice = k
        return
    endif
enddo
write(logmsg,*) 'ERROR: random_choice: ',N,p
call logger(logmsg)
stop
end function

!-----------------------------------------------------------------------------------------
! Generate a random value for CFSE from a distribution with mean = average
! In the simplest case we can allow a uniform distribution about the average.
! Multiplying factor in the range (1-a, 1+a)
! Better to make it a Gaussian distribution:
!  = average*(1+s*R)
! where R = N(0,1), s = std deviation
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function generate_CFSE(average)
real(REAL_KIND) :: average, std
integer :: kpar = 0
real(REAL_KIND) :: R

! Uniform distribution
!R = par_uni(kpar)
!generate_CFSE = (1 - a + 2*a*R)*average
! Gaussian distribution
R = par_rnor(kpar)	! N(0,1)
generate_CFSE = (1 + CFSE_std*R)*average
end function

!-----------------------------------------------------------------------------------------
! Given c (0 < c < 0.5), the fraction of the sphere volume that is cut off, compute
! the distance of the cut from the centre, d, as x = d/R, where R = sphere radius.
!-----------------------------------------------------------------------------------------
subroutine get_slice(c,x)
real(REAL_KIND) :: c,x,xp
real(REAL_KIND), parameter :: epsilon = 0.0001
integer :: k

x = 0.5
xp = x
k = 0
do
    k = k+1
    if (k > 100) then
        write(logmsg,*) 'ERROR: get_slice: failed to converge'
		call logger(logmsg)
        stop
    endif
    x = 1 - 4*c/(2-(1+x)*x**3)
    if (abs(x-xp) < epsilon) exit
    xp = x
enddo

end subroutine



!----------------------------------------------------------------------------------------
! Convert a half life in hours to a decay rate /min
!----------------------------------------------------------------------------------------
real(REAL_KIND) function DecayRate(halflife)
real(REAL_KIND) :: halflife

DecayRate = log(2.0)/(halflife*60)    ! rate/min
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine logger(msg)
character*(*) :: msg
integer :: error
logical :: isopen
character*(1) :: LF = char(94)

error = 0
if (use_TCP) then
    if (awp_0%is_open) then
        call winsock_send(awp_0,trim(msg)//LF,len_trim(msg)+1,error)
    endif
else
	write(*,*) trim(msg)
endif
inquire(unit=nflog,OPENED=isopen)
if (isopen) then
	write(nflog,*) 'msg: ',trim(msg)
	if (error /= 0) then
	    write(nflog,'(a,i4)') 'winsock_send error: ',error
	    close(nflog)
	endif
endif
!if (error /= 0) stop
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_counter(counter, nbins, binmin, binstep, logscale)
type(counter_type) :: counter
integer :: nbins
real(REAL_KIND) :: binmin, binstep
logical :: logscale

counter%nbins = nbins
counter%binmin = binmin
counter%binstep = binstep
counter%logscale = logscale
if (allocated(counter%bincount)) then
	deallocate(counter%bincount)
endif
allocate(counter%bincount(counter%nbins))
counter%bincount = 0
counter%nsamples = 0
counter%total = 0
end subroutine

!-----------------------------------------------------------------------------------------
! The ndist(i) count is incremented if binmin + (i-1)*binstep <= val < binmin + i*binstep
!-----------------------------------------------------------------------------------------
subroutine log_count(counter, sample)
type(counter_type) :: counter
real(REAL_KIND) :: sample
real(REAL_KIND) :: val
integer :: i

if (counter%logscale) then
    val = log10(sample)
else
	val = sample
endif
if (counter%nbins == 1) then
    i = 1
else
    i = (val - counter%binmin)/counter%binstep + 1
    i = max(i,1)
    i = min(i,counter%nbins)
endif
counter%bincount(i) = counter%bincount(i) + 1
counter%nsamples = counter%nsamples + 1
counter%total = counter%total + sample
end subroutine


!--------------------------------------------------------------------------------
!     NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
!     BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
!     SINGLE PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
!--------------------------------------------------------------------------------
SUBROUTINE qsort(a, n, t)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL(REAL_KIND), INTENT(INOUT) :: a(n)
INTEGER, INTENT(INOUT) :: t(n)

!     Local Variables

INTEGER :: i, j, k, l, r, s, stackl(15), stackr(15), ww
REAL(REAL_KIND) :: w, x

s = 1
stackl(1) = 1
stackr(1) = n

!     KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

10 CONTINUE
l = stackl(s)
r = stackr(s)
s = s - 1

!     KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

20 CONTINUE
i = l
j = r
k = (l+r) / 2
x = a(k)

!     REPEAT UNTIL I > J.

DO
  DO
    IF (a(i).LT.x) THEN                ! Search from lower end
      i = i + 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  DO
    IF (x.LT.a(j)) THEN                ! Search from upper end
      j = j - 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  IF (i.LE.j) THEN                     ! Swap positions i & j
    w = a(i)
    ww = t(i)
    a(i) = a(j)
    t(i) = t(j)
    a(j) = w
    t(j) = ww
    i = i + 1
    j = j - 1
    IF (i.GT.j) EXIT
  ELSE
    EXIT
  END IF
END DO

IF (j-l.GE.r-i) THEN
  IF (l.LT.j) THEN
    s = s + 1
    stackl(s) = l
    stackr(s) = j
  END IF
  l = i
ELSE
  IF (i.LT.r) THEN
    s = s + 1
    stackl(s) = i
    stackr(s) = r
  END IF
  r = j
END IF

IF (l.LT.r) GO TO 20
IF (s.NE.0) GO TO 10

RETURN
END SUBROUTINE qsort

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_xyz(k)
integer :: k
integer :: kcell, xyzsum(3)

xyzsum = 0
do kcell = 1,nlist
    if (cell_list(kcell)%ID == 0) cycle
    xyzsum = xyzsum + cell_list(kcell)%site
enddo
write(nfres,'(2i6,3i12)') istep,k,xyzsum
end subroutine


!--------------------------------------------------------------------------------
! Returns the number of cells on the site offset by (x,y,z) from site(:)
!--------------------------------------------------------------------------------
integer function sitecells(site,x,y,z)
integer :: site(3),x,y,z
integer :: k,indx(2)

indx = occupancy(site(1)+x,site(2)+y,site(3)+z)%indx
sitecells = 0
do k = 1,2
    if (indx(k) > 0) sitecells = sitecells + 1
enddo
end function

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine check_zdistribution
real(REAL_KIND) :: tot(NZ)
integer :: kcell, site(3)

tot = 0
do kcell = 1,nlist
    if (cell_list(kcell)%ID == 0) cycle
    site = cell_list(kcell)%site
    tot(site(3)) = tot(site(3)) + 1
enddo
tot = tot/nlist
!write(*,'(a,i8,2f6.3)') 'z distribution: ',nlist,tot(NZ/2+10),tot(NZ/2-10)
!write(*,'(10f6.3)') tot
end subroutine


!-----------------------------------------------------------------------------------------
! For a location xyz, check that the occupancy() info is consistent with the info in
! the cell_list() entries corresponding to occupancy()%indx (if any).  There could be
! 0, 1, or 2 cell_list() entries.
!-----------------------------------------------------------------------------------------
subroutine checkslots(msg,xyz)
character*(*) :: msg
integer :: xyz(3)
integer :: indx(2), k, cells, slots, site(3)
logical :: occ(2)

occ = .false.
cells = 0
slots = getslots(xyz)
indx = occupancy(xyz(1),xyz(2),xyz(3))%indx
do k = 1,2
    if (indx(k) > 0) then
        cells = cells + k
        occ(k) = .true.
        site = cell_list(indx(k))%site
        if (xyz(1) /= site(1) .or. xyz(2) /= site(2) .or. xyz(3) /= site(3)) then
            write(*,'(a,a,8i6)') msg,' checkslots: site error: ',k,xyz,site
            stop
        endif
    elseif (indx(1) < 0) then
        write(*,*) msg,' checkslots: indx: ',xyz,indx
        stop
    endif
enddo
if (slots /= cells) then
    write(*,'(a,a,6i4,2L2)') msg,' checkslots: mismatch: ',xyz,slots,cells,occ
    stop
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check(msg,x,y,z)
character*(*) :: msg
integer :: x,y,z,slots,k,indx(2),site(3)
logical :: occ(2)

site = (/x,y,z/)
slots = getslots(site)
occ = .false.
indx = occupancy(x,y,z)%indx
do k = 1,2
    if (cell_list(indx(k))%ctype > 0) occ(k) = .true.
enddo
if (slots == 1 .and. (.not.occ(1) .or. occ(2))) then
    write(*,*) 'Bad check: ',msg,x,y,z,slots,occ
    stop
endif
if (slots == 2 .and. (.not.occ(2) .or. occ(1))) then
    write(*,*) 'Bad check: ',msg,x,y,z,slots,occ
    stop
endif
if (slots == 3 .and. .not.(occ(1) .and. occ(2))) then
    write(*,*) 'Bad check: ',msg,x,y,z,slots,occ
    stop
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checker
integer :: ic,x,y,z,slots,slot,cells,k,indx(2),site(3)
integer, allocatable :: tot(:)
logical :: occ(2)

allocate(tot(NX))
do x = 1,NX
    tot(x) = 0
    do y = 1,NY
        do z = 1,NZ
            site = (/x,y,z/)
            slots = getslots(site)
            cells = 0
            occ = .false.
            indx = occupancy(x,y,z)%indx
            do k = 1,2
                if (indx(k) > 0) then
                    tot(x) = tot(x) + 1
                    cells = cells + k
                    occ(k) = .true.
                    site = cell_list(indx(k))%site
                    if (x /= site(1) .or. y /= site(2) .or. z /= site(3)) then
                        write(*,'(a,2i2,2i7,6i4)') 'checker: site error: ',k,indx(k),cell_list(indx(k))%ID,x,y,z,site
                        stop
                    endif
                endif
            enddo
            if (slots /= cells) then
                write(*,'(a,6i4,2L2)') 'checker: mismatch: ',x,y,z,slots,cells,occ
                stop
            endif
        enddo
    enddo
enddo

do ic = 1,nlist
    if (cell_list(ic)%ID == 0) cycle  ! gap
    site = cell_list(ic)%site
    indx = occupancy(site(1),site(2),site(3))%indx
    if (ic == indx(1)) then
        slot = 1
    elseif (ic == indx(2)) then
        slot = 2
    else
        write(*,'(a,7i6)') 'ERROR: checker: bad indx: ',ic,site,indx
        stop
    endif
enddo
deallocate(tot)
write(*,*) 'checked OK: ',' nlist: ',nlist
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checkcellsite(kcell)
integer :: kcell
integer :: id,site(3)

site = cell_list(kcell)%site
id = cell_list(kcell)%ID
write(*,'(a,2i8,3i4,4i8)') 'cell: site,indx: ',kcell,id,site,occupancy(site(1),site(2),site(3))%indx
site = cell_list(kcell)%site
id = cell_list(kcell)%ID
write(*,'(a,2i8,3i4,4i8)') 'big_: site,indx: ',kcell,id,site,occupancy(site(1),site(2),site(3))%indx
end subroutine


!-----------------------------------------------------------------------------------------
! Checks the two sites (before and after jump) to see if the one free slot is always #2
! HAPPENS ALL THE TIME
!-----------------------------------------------------------------------------------------
subroutine check_site_indx(site1,site2)
integer :: site1(3),site2(3)
integer :: indx(2)

indx = occupancy(site1(1),site1(2),site1(3))%indx
if (indx(1) == 0 .and. indx(2) /= 0) write(*,*) 'check_site_indx: ',site1,indx
indx = occupancy(site2(1),site2(2),site2(3))%indx
if (indx(1) == 0 .and. indx(2) /= 0) write(*,*) 'check_site_indx: ',site1,indx
end subroutine


end module


