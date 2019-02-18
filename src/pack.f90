module packer

use real_kind_mod

implicit none

type cloc_type
    real(REAL_KIND) :: centre(3)
    real(REAL_KIND) :: d2
end type

type xyz
	integer :: x, y, z
end type

type(cloc_type), allocatable :: cloc(:,:,:)
real(REAL_KIND), allocatable :: cdist(:)
integer, allocatable :: t(:)
type(xyz), allocatable :: xyz_lookup(:)

contains

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine addHCPSheet(j, grid_x, grid_z, y, radius, global_x, global_z)
integer :: j, grid_x, grid_z;
real(REAL_KIND) :: y, radius, global_x, global_z
  !grid_x,		//number of particles in x direction
  !grid_z,		//number of particles in z direction
  !radius,		//radius of spheres
  !global_x,	//global offset of sheet in x
  !global_z		//global offset of sheet in z
real(REAL_KIND) :: offset, x, z;
integer :: i, k

do i = 0,grid_x-1
    do k = 0, grid_z-1
        !need to offset alternate rows by radius
        !offset = (k % 2 != 0) ? radius : 0;
        if (mod(k,2) == 0) then
            offset = 0
        else
            offset = radius;
        endif
        !x position, shifted to center
        x = i * 2 * radius + offset  - grid_x * 2 * radius / 2.0 + global_x
        !z position shifted to center
        z = k * (sqrt(3.0) * radius)  - grid_z * sqrt(3.0) * radius / 2.0 + global_z
        ! x, y, z contain coordinates for sphere position
        cloc(i+1,j+1,k+1)%centre = [x, y, z]
    enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine addHCPCube(grid_x, grid_y, grid_z, radius, global_x, global_y, global_z)
integer :: grid_x, grid_y, grid_z
real(REAL_KIND) :: radius, global_x, global_y, global_z
  !grid_x,      //number of particles in x direction
  !grid_y,      //number of particles in y direction
  !grid_z,      //number of particles in z direction
  !radius,		//radius of sphere
  !global_x,	//global offset in x
  !global_y,	//global offset in y
  !global_z)	//global offset in z
!{
real(REAL_KIND) :: offset_x, offset_z, y
integer :: j
do j = 0, grid_y-1
    y = j * (sqrt(3.0) * radius)
    !need to offset each alternate layer by radius in both x and z direction
    if (mod(j,2) == 0) then
        offset_x = 0
        offset_z = 0
    else
        offset_x = radius
        offset_z = radius
    endif
    call addHCPSheet(j, grid_x, grid_z, y + global_y, radius, offset_x+global_x, offset_z+global_z)
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_cellcentre(i, centre)
integer :: i
real(REAL_KIND) :: centre(3)
integer :: ix, iy, iz

ix = xyz_lookup(t(i))%x
iy = xyz_lookup(t(i))%y
iz = xyz_lookup(t(i))%z
centre = cloc(ix,iy,iz)%centre
end subroutine

!--------------------------------------------------------------------------------
!     NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
!     BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
!     SINGLE PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
!--------------------------------------------------------------------------------
SUBROUTINE qqsort(a, n, t)
IMPLICIT NONE

INTEGER, INTENT(IN)    :: n
REAL(REAL_KIND), INTENT(INOUT)    :: a(n)
INTEGER, INTENT(INOUT) :: t(n)

!     Local Variables

INTEGER         :: i, j, k, l, r, s, stackl(15), stackr(15), ww
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
END SUBROUTINE qqsort

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine SelectCellLocations(nxx, nyy, nzz, cellradius, block_centre)
real(REAL_KIND) :: xrng, zrng, cellradius, block_centre(3)
integer :: npts
integer :: nxx, nyy, nzz
integer :: i, ix, iy, iz
real(REAL_KIND) :: ave(3), global_x=0, global_y=0, global_z=0
real(REAL_KIND) :: xyzmin(3), xyzmax(3)

npts = nxx*nyy*nzz
allocate(cloc(nxx,nyy,nzz))
call addHCPCube(nxx,nyy,nzz,cellradius,global_x,global_y,global_z)

allocate(xyz_lookup(npts))
allocate(cdist(npts))

xyzmin = 1.0e10
xyzmax = -xyzmin
ave = 0
do ix = 1,nxx
	do iy = 1,nyy
		do iz = 1,nzz
			ave = ave + cloc(ix,iy,iz)%centre
			do i = 1,3
				xyzmin(i) = min(cloc(ix,iy,iz)%centre(i),xyzmin(i))
				xyzmax(i) = max(cloc(ix,iy,iz)%centre(i),xyzmax(i))
			enddo
		enddo
	enddo
enddo
ave = ave/npts
write(*,'(a,3f8.1)') 'xyzmin: ',xyzmin
write(*,'(a,3f8.1)') 'xyzmax: ',xyzmax

i = 0
do ix = 1,nxx
	do iy = 1,nyy
		do iz = 1,nzz
			i = i+1
			xyz_lookup(i)%x = ix
			xyz_lookup(i)%y = iy
			xyz_lookup(i)%z = iz
!			cloc(x,y,z)%d2 = (cloc(x,y,z)%centre(1) - ave(1))**2 &
!							+(cloc(x,y,z)%centre(2) - ave(2))**2 &
!							+(cloc(x,y,z)%centre(3) - ave(3))**2 
!			cdist2(i) = cloc(x,y,z)%d2
			cdist(i) = sqrt((cloc(ix,iy,iz)%centre(1) - ave(1))**2 + (cloc(ix,iy,iz)%centre(2) - ave(2))**2)
		enddo
	enddo
enddo
block_centre = ave

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine FreeCellLocations
deallocate(cloc)
deallocate(xyz_lookup)
deallocate(cdist)
!deallocate(t)
end subroutine

end module

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!program main
!use packer
!implicit none
!integer :: i, nblob=1000
!real(REAL_KIND) :: blobcentre(3), centre(3), ave(3)
!real(REAL_KIND) :: radius = 1.0
!
!blobcentre = [0, 0, 0]
!call SelectCellLocations(nblob, radius)
!ave = 0
!write(*,*) 'Cell centres:'
!do i = 1,nblob
!	call get_cellcentre(i,centre)
!!	write(*,'(3f8.2)') centre
!	ave = ave + centre
!enddo
!ave = ave/nblob
!write(*,*)
!write(*,'(a,3f8.3)') 'Average cell centre: ',ave
!do i = 1,nblob
!	call get_cellcentre(i,centre)
!	centre = centre - ave + blobcentre
!enddo
!call FreeCellLocations
!
!end program

    