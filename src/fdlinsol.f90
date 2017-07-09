module slu_solver

implicit none

integer, parameter :: nfHB = 10, nfrhs = 11

logical :: verbose

contains

!-----------------------------------------------------------------------
subroutine readHBsize(HBfilename, n, nz)
character*(*) :: HBfilename
integer :: n, nz

integer :: totcrd, ptrcrd, indcrd, valcrd, rhscrd, ncol, nrow, nrhs, nzrhs, nel, i
character*(72) :: title
character*(30) :: key
character*(3) :: atype, rhstyp
character*(16) :: ptrfmt, indfmt
character*(20) :: valfmt, rhsfmt

write(*,*) 'HBfilename: ',trim(HBfilename)
write(*,*)
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
close(nfHB)
return
998 write(*,*) 'Read error: Harwell/Boeing matrix'
stop
end subroutine

!-----------------------------------------------------------------------
! Read the Harwell/Boeing matrix
! (I think the HB file is always 1-based)
!-----------------------------------------------------------------------
subroutine readHB(HBfilename, Ap, Ai, Aval, n, nz)
character*(*) :: HBfilename
integer :: n, nz
integer :: Ap(:), Ai(:)
real(8) :: Aval(:)

integer :: totcrd, ptrcrd, indcrd, valcrd, rhscrd, ncol, nrow, nrhs, nzrhs, nel, i
character*(72) :: title
character*(30) :: key
character*(3) :: atype, rhstyp
character*(16) :: ptrfmt, indfmt
character*(20) :: valfmt, rhsfmt

write(*,*) 'HBfilename: ',trim(HBfilename)
write(*,*)
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
!allocate(Ap(n+1))
!allocate(Ai(nz))
!allocate(Aval(nz))
! read the matrix (1-based)
read (nfHB, ptrfmt, err = 998) (Ap(i), i = 1, ncol+1)
read (nfHB, indfmt, err = 998) (Ai(i), i = 1, nz)
read (nfHB, valfmt, err = 998) (Aval(i), i = 1, nz)
close(nfHB)
return
998 write(*,*) 'Read error: Harwell/Boeing matrix'
stop
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

