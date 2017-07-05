!=======================================================================
!== umf4hb =============================================================
!=======================================================================

module umf_solver

implicit none

integer, parameter :: nfHB = 10, nfrhs = 11

integer, parameter :: SOLVE_NO_ITERATION = 1
integer, parameter :: SOLVE_ITERATION    = 2
integer, parameter :: SOLVE_STAGEWISE    = 3

logical :: verbose

contains

!-----------------------------------------------------------------------
! Read the Harwell/Boeing matrix
! (I think the HB file is always 1-based)
!-----------------------------------------------------------------------
subroutine readHB(HBfilename, Aval, Ap, Ai, n, nz)
character*(*) :: HBfilename
integer :: n, nz
integer, pointer :: Ap(:), Ai(:)
real(8), pointer :: Aval(:)

integer :: totcrd, ptrcrd, indcrd, valcrd, rhscrd, ncol, nrow, nrhs, nzrhs, nel, i
character*(72) :: title
character*(30) :: key
character*(3) :: atype, rhstyp
character*(16) :: ptrfmt, indfmt
character*(20) :: valfmt, rhsfmt

open(nfHB,file=HBfilename,status='old')
read (nfHB, 10, err = 998) &
	   title, key, &
	   totcrd, ptrcrd, indcrd, valcrd, rhscrd, &
	   atype, nrow, ncol, nz, nel, &
	   ptrfmt, indfmt, valfmt, rhsfmt
if (rhscrd > 0) then
!  new Harwell/Boeing format:
   read (nfHB, 20, err = 998) rhstyp, nrhs, nzrhs
endif
10      format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
20      format (a3, 11x, 2i14)

write(*,*) 'Matrix key: ', key

n = nrow
if (atype /= 'RUA' .or. nrow /= ncol) then
   write(*,*) 'Error: can only handle square RUA matrices'
   stop
endif
allocate(Ap(n+1))
allocate(Ai(nz))
allocate(Aval(nz))
! read the matrix (1-based)
read (nfHB, ptrfmt, err = 998) (Ap(i), i = 1, ncol+1)
read (nfHB, indfmt, err = 998) (Ai(i), i = 1, nz)
read (nfHB, valfmt, err = 998) (Aval(i), i = 1, nz)
close(nfHB)
return
998 write(*,*) 'Read error: Harwell/Boeing matrix'
stop
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine solve_umf(isolve_mode, n, Ap, Ai, Ax, b, x)
integer :: isolve_mode, n
real(8) :: Ax(:), b(:), x(:)
real(8), allocatable :: r(:)
integer :: Ai(:), Ap(:)

real(8) :: control(20), info(90)
integer :: numeric, symbolic, status, sys, filenum

! set default parameters
write(*,*) 'call umf4def'
call umf4def (control)

! print control parameters.  set control (1) to 1 to print
! error messages only
if (verbose) then
	control(1) = 2
else
	control(1) = 1
endif
write(*,*) 'call umf4pcon'
call umf4pcon (control)

! pre-order and symbolic analysis
call umf4sym (n, n, Ap, Ai, Ax, symbolic, control, info)

if (verbose) then
!   print statistics computed so far
!   call umf4pinf (control, info) could also be done.
	write(*,80) info(1), info(16), &
	   (info(21) * info(4)) / 2**20, &
	   (info(22) * info(4)) / 2**20, &
	   info(23), info(24), info(25)
80  format ('symbolic analysis:',/, &
	   '   status:  ', f5.0, /, &
	   '   time:    ', e10.2, ' (sec)'/, &
	   '   estimates (upper bound) for numeric LU:', /, &
	   '   size of LU:    ', f10.2, ' (MB)', /, &
	   '   memory needed: ', f10.2, ' (MB)', /, &
	   '   flop count:    ', e10.2, / &
	   '   nnz (L):       ', f10.0, / &
	   '   nnz (U):       ', f10.0)
endif

! check umf4sym error condition
if (info(1) < 0) then
	write(*,*) 'Error occurred in umf4sym: ', info(1)
	stop
endif

! numeric factorization
write(*,*) 'call umf4num'
call umf4num (Ap, Ai, Ax, symbolic, numeric, control, info)

if (verbose) then
!   print statistics for the numeric factorization
!   call umf4pinf (control, info) could also be done.
	write(*,90) info(1), info(66), &
	   (info (41) * info(4)) / 2**20, &
	   (info (42) * info(4)) / 2**20, &
	   info (43), info(44), info(45)
90   format ('numeric factorization:',/, &
	   '   status:  ', f5.0, /, &
	   '   time:    ', e10.2, /, &
	   '   actual numeric LU statistics:', /, &
	   '   size of LU:    ', f10.2, ' (MB)', /, &
	   '   memory needed: ', f10.2, ' (MB)', /, &
	   '   flop count:    ', e10.2, / &
	   '   nnz (L):       ', f10.0, / &
	   '   nnz (U):       ', f10.0)
endif

! check umf4num error condition
if (info(1) < 0) then
	write(*,*) 'Error occurred in umf4num: ', info(1)
	stop
endif
! free the symbolic analysis
call umf4fsym (symbolic)

if (isolve_mode == SOLVE_NO_ITERATION) then
	! solve Ax=b, without iterative refinement
	sys = 0
	write(*,*) 'call umf4sol'
	call umf4sol (sys, x, b, numeric, control, info)
	if (info(1) < 0) then
		write(*,*) 'Error occurred in umf4sol: ', info(1)
		stop
	endif
elseif (isolve_mode == SOLVE_ITERATION) then
	! solve Ax=b, with iterative refinement
	sys = 0
	write(*,*) 'call umf4solr'
        call umf4solr (sys, Ap, Ai, Ax, x, b, numeric, control, info)
	if (info(1) < 0) then
		write(*,*) 'Error occurred in umf4solr: ', info(1)
		stop
	endif
elseif (isolve_mode == SOLVE_STAGEWISE) then
! the factorization is PAQ=LU, PRAQ=LU, or P(R\A)Q=LU.
!   x = R*b (or x=R\b, or x=b, as appropriate)
	call umf4scal (x, b, numeric, status)
	if (status < 0) then
		write(*,*) 'Error occurred in umf4scal: ', status
		stop
	endif

!   solve P'Lr=x for r (using r as workspace)
	allocate(r(n))
	sys = 3
	call umf4sol (sys, r, x, numeric, control, info)
	if (info(1) < 0) then
		write(*,*) 'Error occurred in umf4sol: ', info(1)
		stop
	endif

!   solve UQ'x=r for x
	sys = 9
	call umf4sol (sys, x, r, numeric, control, info)
	if (info(1) < 0) then
		write(*,*) 'Error occurred in umf4sol: ', info(1)
		stop
	endif
	deallocate(r)
endif

! free the numeric factorization
write(*,*) 'call umf4fnum'
call umf4fnum (numeric)

! No LU factors (symbolic or numeric) are in memory at this point.

! print final statistics
if (verbose) then
	call umf4pinf (control, info)
endif

end subroutine

!=======================================================================
!== resid ==============================================================
! Compute the residual, r = Ax-b, its max-norm, and print the max-norm
! Note that A is zero-based.
!=======================================================================
subroutine resid (n, Ap, Ai, Ax, x, b, r)
integer :: n, Ap(:), Ai(:), j, i, p
double precision :: Ax(:), x(:), b(:), r(:), rmax, aij

r(1:n) = -b(1:n)

do j = 1,n
	do p = Ap(j) + 1, Ap(j+1)
		i = Ai(p) + 1
		aij = Ax(p)
		r(i) = r(i) + aij * x(j)
	enddo
enddo
rmax = 0
do i = 1, n
	rmax = max(rmax, r(i))
enddo

write(*,*) 'norm (A*x-b): ', rmax
end subroutine

end module


!-----------------------------------------------------------------------
! UMFPACK Copyright (c) 2005-2012 by Timothy A. Davis,
! http://www.suitesparse.com.  All Rights Reserved.  See ../Doc/License
! for License.
!-----------------------------------------------------------------------

! umf4hb:
! read a sparse matrix in the Harwell/Boeing format, factorizes
! it, and solves Ax=b.  Also saves and loads the factors to/from a
! file.  Saving to a file is not required, it's just here to
! demonstrate how to use this feature of UMFPACK.  This program
! only works on square RUA-type matrices.
!
! This is HIGHLY non-portable.  It may not work with your C and
! FORTRAN compilers.  See umf4_f77wrapper.c for more details.

!Program main
!use umf_solver
!implicit none

!integer, pointer :: Ap(:), Ai(:)
!integer :: n, nz, i, j, p
!real(8) :: aij, xj
!real(8), pointer :: Ax(:)
!real(8), allocatable :: x(:), b(:), r(:)
!character*(64) :: HBfilename, rhsfilename
!character :: ans
!integer :: isolve_mode
!real(4) :: t1, t2

!verbose = .false.

!----------------------------------------------------------------
! read the Harwell/Boeing matrix
!----------------------------------------------------------------
!write(*,*) 'Enter HB file name'
!read(*,'(a)') HBfilename
!call readHB(HBfilename, Ax, Ap, Ai, n, nz)
!allocate(x(n),b(n),r(n))

!write(*,*) 'Read the rhs (Y/N)?'
!read(*,'(a)') ans
!if (ans == 'Y' .or. ans == 'y') then
!	write(*,*) 'Enter rhs file'
!	read(*,'(a)') rhsfilename
!	open(nfrhs,file=rhsfilename,status='old')
!	read(nfrhs,*) b(1:n)
!	close(nfrhs)
!else
	!----------------------------------------------------------------
	! create the right-hand-side, assume x (i) = 1 + i/n
	!----------------------------------------------------------------
!	b = 0
	! b = A*x
!	do j = 1,n
!		xj = j
!		xj = 1 + xj / n
!		do p = Ap(j), Ap(j+1)-1
!			i = Ai(p)
!			aij = Ax(p)
!			b(i) = b(i) + aij * xj
!		enddo
!	enddo
!endif

!----------------------------------------------------------------
! convert from 1-based to 0-based
!----------------------------------------------------------------
!do j = 1, n+1
!	Ap(j) = Ap(j) - 1
!enddo
!do p = 1, nz
!	Ai(p) = Ai(p) - 1
!enddo

!isolve_mode = SOLVE_NO_ITERATION
!call cpu_time(t1)
!call solve(isolve_mode, n, Ap, Ai, Ax, b, x)
!call cpu_time(t2)
!write(*,'(a,f8.1)') 'Solve time (sec): ',t2-t1

! print the residual.
!call resid (n, Ap, Ai, Ax, x, b, r)

! write(*,'(10f7.0)') (x(i),i=1,n)

!end program

