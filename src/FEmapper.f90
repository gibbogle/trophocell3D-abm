module FEmapper

use nleq1_mod

implicit none

integer :: NDim
integer :: Nvertices
real(8) :: vertex(8,3)
real(8) :: xglobal(3)

contains

!---------------------------------------------------------------------------
! Determine the local (reference) coords x(:) for a point xg(:) in a hex
! element with vertices vert(:,:)
!---------------------------------------------------------------------------
subroutine mapper(ND,Nvert,vert,xg,x)
integer :: ND, Nvert
real(8) :: vert(8,3)
real(8) :: xg(3), x(3)
real(8) :: rtol = 1.0d-5
integer :: ierr

Ndim = ND
Nvertices = Nvert
vertex(1:Nvert,1:Ndim) = vert(1:Nvert,1:Ndim)
xglobal(1:Ndim) = xg(1:Ndim)

rtol = 1.0d-5
x(1:Ndim) = 0.0	! initial guess
CALL NLEQ1E(FCN,dummyjac,Ndim,x,rtol,ierr)
end subroutine

!---------------------------------------------------------------------------
! This function is used by nleq1e to solve the set of nonlinear equations
! given by eqtn 21 in Matringe(2007).  In the case of a hex element there are
! 8 simultaneous equations, and F(1:8) holds the deviation from 0 for each
! equation.
!---------------------------------------------------------------------------
subroutine FCN(N,X,F,IFAIL)
integer :: N, IFAIL
real(8) :: X(N),F(N)
integer :: i
real(8) :: fsum(3)

fsum = 0
do i = 1,Nvertices
	fsum(1:N) = fsum(1:N) + vertex(i,1:N)*fshape(Nvertices,i,X)
enddo
F(1:N) = fsum(1:N) - xglobal(1:N)
end subroutine
	
!-----------------------------------------------------------------
! Note the error in Matringe 2007, Eq. 23
! In N3, N4, N7 and N8 the signs of the x^ terms all need to be 
! switched to be consistent with the vertex numbering in Fig. 2
!-----------------------------------------------------------------
real(8) function fshape(nvert,i,X)
integer :: nvert, i
real(8) :: X(3)

if (Ndim == 2 .and. nvert == 4) then	! 2D quad element
	if     (i == 1) then
		fshape = 0.25*(1 - X(1))*(1 - X(2))
	elseif (i == 2) then
		fshape = 0.25*(1 + X(1))*(1 - X(2))
	elseif (i == 3) then
		fshape = 0.25*(1 + X(1))*(1 + X(2))
	elseif (i == 4) then
		fshape = 0.25*(1 - X(1))*(1 + X(2))
	endif
elseif (nvert == 6) then				! 3D wedge element
	if     (i == 1) then
		fshape = 0.5*(1 - X(1) - X(2))*(1 - X(3))
	elseif (i == 2) then
		fshape = 0.5*X(1)*(1 - X(3))
	elseif (i == 3) then
		fshape = 0.5*X(2)*(1 - X(3))
	elseif (i == 4) then
		fshape = 0.5*(1 - X(1) - X(2))*(1 + X(3))
	elseif (i == 5) then
		fshape = 0.5*X(1)*(1 + X(3))
	elseif (i == 6) then
		fshape = 0.5*X(2)*(1 + X(3))
	endif
elseif (nvert == 8) then				! 3D hexahedral element
	if     (i == 1) then
		fshape = 0.125*(1 - X(1))*(1 - X(2))*(1 - X(3))
	elseif (i == 2) then
		fshape = 0.125*(1 + X(1))*(1 - X(2))*(1 - X(3))
	elseif (i == 3) then
		fshape = 0.125*(1 + X(1))*(1 + X(2))*(1 - X(3))
	elseif (i == 4) then
		fshape = 0.125*(1 - X(1))*(1 + X(2))*(1 - X(3))
	elseif (i == 5) then
		fshape = 0.125*(1 - X(1))*(1 - X(2))*(1 + X(3))
	elseif (i == 6) then
		fshape = 0.125*(1 + X(1))*(1 - X(2))*(1 + X(3))
	elseif (i == 7) then
		fshape = 0.125*(1 + X(1))*(1 + X(2))*(1 + X(3))
	elseif (i == 8) then
		fshape = 0.125*(1 - X(1))*(1 + X(2))*(1 + X(3))
	endif
endif
end function
	
subroutine dummyjac(N,LDJAC,X,DFDX,IFAIL)
integer :: N, LDJAC, IFAIL
real(8) :: X(N), DFDX(LDJAC,N)
end subroutine

end module