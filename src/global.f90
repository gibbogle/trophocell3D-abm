!==========================================================================================
!==========================================================================================

module global

use real_kind_mod
use omp_lib
use par_zig_mod
use winsock

implicit none

! Files
integer, parameter :: nfinput=1, nfcell = 10, nfout = 11, nfvec = 12, nfpath = 13, nfres = 14, nftraffic = 16, nfrun = 17, &
					  nftravel=18, nfcmgui=19, nfpos=20, nflog=21, nfchemo=22, nffacs = 24, nflog2=26, nfvel=27

! General parameters
integer, parameter :: BIG_INT = 2**30

integer, parameter :: WALL = 0
integer, parameter :: ALIVE = 1
integer, parameter :: DEAD = 2
integer, parameter :: GONE_BACK = 3
integer, parameter :: GONE_THROUGH = 4
integer, parameter :: TAGGED_CELL   = 99

integer, parameter :: TCP_PORT_0 = 5000		! main communication port (logging)
integer, parameter :: TCP_PORT_1 = 5001		! data transfer port (plotting)
integer, parameter :: TCP_PORT_2 = 5002
integer, parameter :: TCP_PORT_3 = 5003

integer, parameter :: CHEMO_1 = 1
integer, parameter :: CHEMO_2 = 2

integer, parameter :: NCTYPES = 1
integer, parameter :: TROPHO_CELL = 1
integer, parameter :: MAX_CELLTYPES = 1
integer, parameter :: MAX_NBRS = 100

integer, parameter :: TERMINAL_MITOSIS   = 1
integer, parameter :: CONTINUOUS_MITOSIS = 2
integer, parameter :: MITOSIS_MODE = TERMINAL_MITOSIS

real(REAL_KIND), parameter :: BIG_TIME = 100000
!real(REAL_KIND), parameter :: small_d = 0.1e-4		! 0.1 um -> cm
real(REAL_KIND), parameter :: small_d = 0.1		! 0.1 um
logical, parameter :: save_input = .true.

integer, parameter :: ndt_max = 30

integer, parameter :: MAX_NLIST = 100000
integer, parameter :: MAX_NGAPS = 1000

! Chemokine/receptor parameters
integer, parameter :: MAX_CHEMO = 2
integer, parameter :: MAX_RECEPTOR = 2

!-------flow parameters:--------------------------
real (REAL_KIND), parameter :: outletPressure = 10/(7.5)*10**(3)!*133.322
real (REAL_KIND), parameter :: blood_viscosity= 0.003

! Data above this line almost never change
!==============================================================================================================

! Run parameters

! Parameters and switches for calibration etc.
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

!Node data
type Node_type
    integer :: ID
    real (REAL_KIND) :: site(3)
end type

type element_type
    integer :: nel, nel_x ,nel_y ,nel_z, nnp
    real (REAL_KIND) :: x_min, x_max, y_min, y_max, z_min, z_max
    real (REAL_KIND) :: x_width, y_width, z_width, volume
    real (REAL_KIND), allocatable :: vol_fraction(:)
    integer :: nsd, nen ,ndof, neq
end type


type neighbour_type
	integer :: indx
!	logical :: incontact
	logical*1 :: contact(2,2)
end type


type cell_type
    integer :: ID
	integer :: state
    integer :: site(3)
	logical :: Iphase
    integer :: nspheres             ! =1 for Iphase, =2 for Mphase
	real(REAL_KIND) :: V			! actual volume (um3)
	real(REAL_KIND) :: dVdt
	real(REAL_KIND) :: radius(2)	! sphere radii (um)
	real(REAL_KIND) :: centre(3,2)  ! sphere centre positions
	real(REAL_KIND) :: d			! centre separation distance (um)
	real(REAL_KIND) :: birthtime
	real(REAL_KIND) :: t_start_mitosis
	real(REAL_KIND) :: mitosis		! level of mitosis (0 - 1)
	real(REAL_KIND) :: V_divide
	real(REAL_KIND) :: d_divide		! centre separation distance at the end of mitosis
    integer :: step
	real(REAL_KIND) :: receptor_level(MAX_RECEPTOR)
	real(REAL_KIND) :: receptor_saturation_time(MAX_RECEPTOR)
    integer :: tag
    integer(2) :: ctype
	integer(2) :: lastdir
	integer :: dtotal(3)
	integer :: nbrs
	type(neighbour_type) :: nbrlist(100)
end type

type occupancy_type
    integer :: indx(2)
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

real(REAL_KIND) :: days                            ! number of days to simulate
integer :: seed(2)                      ! seed vector for the RNGs

!---------------------------------------------------
! end of parameters to read from input file
!---------------------------------------------------

! Geometry data
integer :: Ncells, nwallcells, nlist, ndt
real(REAL_KIND) :: DELTA_X, PI
real(REAL_KIND) :: TagRadius
real(REAL_KIND) :: Vc, Ve

real(REAL_KIND) :: Raverage
real(REAL_KIND) :: tube_radius, tube_length, plug_zmin, plug_zmax, plug_hmax

!mesh data
integer :: nofaces, NofElements,NofElements_z
integer :: CommonFaceNo
integer, allocatable :: Cylinder_ELEMENT(:,:)
real (REAL_KIND), allocatable :: Cylinder_vertices(:,:)
integer, allocatable :: ElToFace(:,:)
integer, allocatable :: AllSurfs(:,:)
real (REAL_KIND), allocatable :: Centroids(:,:)
integer, allocatable :: CylinCub_list(:,:)

real (REAL_KIND), allocatable :: e_Z(:)
integer, allocatable :: layered_el(:,:)

!Q-P Matrix
real (REAL_KIND), allocatable :: Q_P(:)
!Body Force Matrix
real (REAL_KIND), allocatable :: BodyForce(:)
integer, allocatable :: FreeFaces(:)

type Bmatrix_cell
    real (REAL_KIND), ALLOCATABLE :: Inside_B(:,:)
end type
type(Bmatrix_cell), allocatable:: B_MATRIX(:)

!Boundary Condition
integer, allocatable :: DCInlet(:)
integer, allocatable ::DCOutlet(:)
integer, allocatable ::NCWall(:)
integer, allocatable ::NCInlet(:)

!cell velocity
real (REAL_KIND), allocatable :: u_cell_x(:), u_cell_y(:), u_cell_z(:)
real (REAL_KIND), allocatable :: F_shear(:,:)

! Motility data
real(REAL_KIND) :: DELTA_T       ! minutes
real(REAL_KIND) :: BETA, RHO
integer :: n_cell_positions

real(REAL_KIND) :: Vmax
real(REAL_KIND) :: chemo_attraction
real(REAL_KIND) :: Kadhesion, Kstay	!,Kdrag

! Force parameters
real(REAL_KIND) :: a_separation, kdrag, frandom
real(REAL_KIND) :: a_force, c_force_cell, c_force_wall, x0_force, x1_force, xcross1_force, xcross2_force
real(REAL_KIND) :: t_fmover, delta_tmove, dt_min, delta_min, delta_max

! Chemotaxis data
integer, parameter :: chemo_N = 2
real(REAL_KIND), allocatable :: chemo_p(:,:,:,:)
real(REAL_KIND) :: grad_amp(2), grad_dir(2,3)
! Background flow
real(REAL_KIND) :: BG_flow_amp
real(REAL_KIND):: inletPressure != 65/(7.5)*10**(3)
integer :: chemo_coef1, chemo_coef2

! Cell data
type(cell_type), allocatable, target :: cell_list(:)
type(source_type), allocatable :: sourcelist(:)
integer :: nsources
integer, allocatable :: gaplist(:)
integer, allocatable :: perm_index(:)
integer :: lastID, n2Dsites, ngaps, ntaglimit, ntagged=0, ntagged_left
integer :: k_nonrandom(2)
integer :: nleft
integer :: nadd_sites, ndivisions
real(REAL_KIND) :: lastbalancetime
real(REAL_KIND) :: settle_hrs

real(REAL_KIND) :: alpha_v, k_detach
real(REAL_KIND) :: dr_mitosis, mitosis_hours, mitosis_duration
real(REAL_KIND) :: Vdivide0, dVdivide, Rdivide0, MM_THRESHOLD, medium_volume0, total_volume
real(REAL_KIND) :: divide_time_median(MAX_CELLTYPES), divide_time_shape(MAX_CELLTYPES), divide_time_mean(MAX_CELLTYPES), celltype_fraction(MAX_CELLTYPES)
type(dist_type) :: divide_dist(MAX_CELLTYPES)
real(REAL_KIND) :: d_nbr_limit


logical :: firstSummary

integer :: kcell_debug

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
logical :: use_hysteresis = .false.
logical :: use_permute = .true.
logical :: use_packing = .false.
logical :: use_loosepack = .false.
logical :: use_makeRing = .true.
logical :: use_settling = .true.
logical :: settling
logical :: calibration_run = .false.
logical :: no_cells = .true.

logical :: dbug = .false.

!DEC$ ATTRIBUTES DLLEXPORT :: nsteps	!istep
contains

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
integer :: n1,n2,kpar
integer :: k,R

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

!-----------------------------------------------------------------------------------------
! This needs to be given radial symmetry
!-----------------------------------------------------------------------------------------
subroutine get_random_dr(dr)
real(REAL_KIND) :: dr(3)
integer :: kpar=0

dr(1) = 2*(par_uni(kpar) - 0.5)
dr(2) = 2*(par_uni(kpar) - 0.5)
dr(3) = 2*(par_uni(kpar) - 0.5)
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
subroutine normalize(r)
real(REAL_KIND) :: r(3)

r = r/norm(r)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!real(REAL_KIND) function norm2(r)
!real(REAL_KIND) :: r(3)
!
!norm2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
!end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
integer function cell_count()
integer :: kcell, ntot
ntot = 0
do kcell = 1,nlist
    if (cell_list(kcell)%state /= ALIVE) cycle             ! skip gaps in the list
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
! The perm_index captures only ALIVE cells, i.e. excludes WALL cells
!-----------------------------------------------------------------------------------------
subroutine make_perm_index(ok)
logical :: ok
integer :: np, kcell, kpar=0

np = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == ALIVE .or. cell_list(kcell)%state == GONE_BACK .or. cell_list(kcell)%state == GONE_THROUGH) then!cycle
        np = np + 1
        perm_index(np) = kcell
    else
        cycle
    endif
enddo
if (np /= ncells-nwallcells) then
	write(logmsg,*) 'Error: make_perm_index: np /= Ncells: ',np,ncells-nwallcells,nlist
	call logger(logmsg)
	ok = .false.
	return
endif
if (use_permute) then
	call permute(perm_index,np,kpar)
endif
ok = .true.
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
    if (cell_list(kcell)%state /= ALIVE) cycle
    xyzsum = xyzsum + cell_list(kcell)%site
enddo
write(nfres,'(2i6,3i12)') istep,k,xyzsum

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
subroutine ax_st ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! AX_ST computes A*x for a matrix stored in sparset triplet form.
!
!  Discussion:
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) W(N), the value of A*X.
!
  implicit none

  integer :: n
  integer :: nz_num

  real (REAL_KIND) :: a(nz_num)
  integer :: i
  integer :: ia(nz_num)
  integer :: j
  integer :: ja(nz_num)
  integer :: k
  real (REAL_KIND) :: w(n)
  real (REAL_KIND) :: x(n)

  w(1:n) = 0.0D+00

  do k = 1, nz_num
    i = ia(k)
    j = ja(k)
    w(i) = w(i) + a(k) * x(j)
  end do

  return
end subroutine

!------------------------------------------------
!------------------------------------------------
subroutine get_Nsteps(N)
integer :: N
N=Nsteps
end subroutine

end module


