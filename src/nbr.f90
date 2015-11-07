module nbr

use global

implicit none

contains

!-----------------------------------------------------------------------------------------
! Dumb version to start
!-----------------------------------------------------------------------------------------
subroutine setup_nbrlists
type(cell_type), pointer :: cp1, cp2
real(REAL_KIND) :: R1, c1(3), R2, c2(3), v(3), d2, d, Rsum, dfactor
integer :: kcell, k2, k, isphere1, nspheres1, isphere2, nspheres2, nbrs
type(neighbour_type) :: nbrlist(100)
logical :: near, incontact, contact(2,2)
logical :: dbug = .false.

do kcell = 1,nlist
	cp1 => cell_list(kcell)
	if (cp1%state == DEAD) cycle
	!if (cp1%Iphase) then
	!	nspheres1 = 1
	!else
	!	nspheres1 = 2
	!endif
    nspheres1 = cp1%nspheres
	nbrs = 0
	do k2 = 1,nlist
		if (k2 == kcell) cycle
		cp2 => cell_list(k2)
		if (cp2%state == DEAD) cycle
		!if (cp2%Iphase) then
		!	nspheres2 = 1
		!else
		!	nspheres2 = 2
		!endif
        nspheres2 = cp2%nspheres
		near = .false.
		incontact = .false.
		contact = .false.
		do isphere1 = 1,nspheres1
			R1 = cp1%radius(isphere1)
			c1 = cp1%centre(:,isphere1)
			do isphere2 = 1,nspheres2
				R2 = cp2%radius(isphere2)
				c2 = cp2%centre(:,isphere2)
				Rsum = R1 + R2
				v = c2 - c1
				d2 = dot_product(v,v)
				d = sqrt(d2)
				if (use_hysteresis) then
					if (d < d_nbr_limit) then
						if (dbug) write(*,'(5f8.2)') R1,R2,Rsum,d,d_nbr_limit
						dfactor = 1
						do k = 1,cp1%nbrs
							if (cp1%nbrlist(k)%indx == k2) then
								if (cp1%nbrlist(k)%contact(isphere1,isphere2)) dfactor = k_detach
								exit
							endif
						enddo
						near = .true.
						if (d < dfactor*Rsum) then
							incontact = .true.
							contact(isphere1,isphere2) = .true.
						endif
					endif
				else
					if (d < d_nbr_limit) then
						near = .true.
						exit
					endif
				endif
			enddo
		enddo
		if (near) then
			nbrs = nbrs + 1
			if (nbrs > MAX_NBRS) then
				write(logmsg,'(a,5i4)') 'Size of nbrlist exceeded: istep, kcell, ncells: ',istep,kcell,ncells,nbrs,MAX_NBRS
				call logger(logmsg)
				write(*,'(10i6)') nbrlist(:)%indx
				stop
			endif
			nbrlist(nbrs)%indx = k2
			nbrlist(nbrs)%contact = contact
		endif
	enddo
	cp1%nbrs = nbrs
	cp1%nbrlist(1:nbrs) = nbrlist(1:nbrs)
    !if (kcell == 10 .or. kcell == 30) then  !!!!!!!!!! #10 does not see #30, and #30 does not see #10 !!!!!! this is the problem
    !	write(*,'(a,3i4)') 'kcell: ',kcell,nspheres1,nbrs
	   ! write(*,'(10i6)') nbrlist(1:nbrs)%indx 
    !endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_all_nbrlists
integer :: kcell, max_nbrs

!write(nflog,*) 'update_all_nbrlists'
max_nbrs = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	call update_nbrlist(kcell)
	max_nbrs = max(max_nbrs,cell_list(kcell)%nbrs)
enddo
!write(nflog,*) 'max_nbrs: ',max_nbrs
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_nbrlist(kcell)
integer :: kcell
type (cell_type), pointer :: cp1, cp2, cp3
integer :: nbrlist(1000), t(1000)
real(REAL_KIND) :: nbr_d(1000)
integer :: nbrs, ns1, k, k2, k3, knbr, i
real(REAL_KIND) :: c1(3,2), d
logical :: inlist

cp1 => cell_list(kcell)
if (cp1%Iphase) then
    ns1 = 1
    c1(:,1) = cp1%centre(:,1)
else
    ns1 = 2
    c1(:,1) = cp1%centre(:,1)
    c1(:,2) = cp1%centre(:,2)
endif
nbrs = 0
do k2 = 1,cp1%nbrs
    knbr = cp1%nbrlist(k2)%indx
    if (knbr == kcell) cycle
    cp2 => cell_list(knbr)
	if (cp2%state == DEAD) cycle
    call get_dist(c1,ns1,cp2,d)
    if (d < d_nbr_limit) then
        call add_to_list(knbr,d,nbrlist,nbr_d,nbrs)
    endif
    do k3 = 1,cp2%nbrs
        knbr = cp2%nbrlist(k3)%indx
        if (knbr == kcell) cycle
        cp3 => cell_list(knbr)
		if (cp3%state == DEAD) cycle
        call get_dist(c1,ns1,cp3,d)
        if (d < d_nbr_limit) then
            call add_to_list(knbr,d,nbrlist,nbr_d,nbrs)
        endif
    enddo
enddo
if (nbrs == 0) return
! sort nbrlist(:) by nbr_d
do i = 1,nbrs
    t(i) = i
enddo
call qqsort(nbr_d,nbrs,t)     ! sort in increasing order
! Now the ith nearest is nbrlist(t(i))
! set cp1%nbrlist(:) as the closest #
cp1%nbrs = min(nbrs,100)
do i = 1,cp1%nbrs
	cp1%nbrlist(i)%indx = nbrlist(t(i))
enddo
! Need to ensure that if knbr is in nbrlist then kcell is in the list for knbr
do i = 1,cp1%nbrs
	knbr = cp1%nbrlist(i)%indx
    cp2 => cell_list(knbr)
	if (cp2%state == DEAD) cycle
	inlist = .false.
    do k2 = 1,cp2%nbrs
		if (kcell == cp2%nbrlist(k2)%indx) then
			inlist = .true.
			exit
		endif
	enddo
	if (.not.inlist) then
		cp2%nbrs = cp2%nbrs + 1
		cp2%nbrlist(cp2%nbrs)%indx = kcell
	endif
enddo
	
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_dist(c1,ns1,cp2,d)
integer :: ns1
type (cell_type), pointer :: cp2
real(REAL_KIND) :: c1(3,2), d
real(REAL_KIND) :: r2min, r(3)
integer :: is1, is2, ns2

if (cp2%Iphase) then
    ns2 = 1
else
    ns2 = 2
endif
r2min = 1.0e10
do is1 = 1,ns1
    do is2 = 1,ns2
        r = c1(:,is1) - cp2%centre(:,is2)
        r2min = min(r2min,dot_product(r,r))
    enddo
enddo
d = sqrt(r2min)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine add_to_list(knbr,d,nbrlist,nbr_d,nbrs)
integer :: knbr, nbrs, nbrlist(:)
real(REAL_KIND) :: d,nbr_d(:)
integer :: k

do k = 1,nbrs
    if (nbrlist(k) == knbr) return
enddo
nbrs = nbrs + 1
nbrlist(nbrs) = knbr
nbr_d(nbrs) = d
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

end module
