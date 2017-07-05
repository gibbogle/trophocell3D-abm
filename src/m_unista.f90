Module m_unista

use global
!use ISO_FORTRAN_ENV, only : REAL32, REAL64, REAL128

contains

!----------------------------------------------
subroutine find_cell_in_CylElem(e_Z, nel_z,N_Cell,Element,vertices,el,cell_elem)
integer :: nel_z, N_Cell,i, m, k,e
integer :: Element(:,:), el(:,:)
integer :: cell_elem(:)
integer :: base(6,4)
real (REAL_KIND):: vertices(:,:),e_Z(:)
real (REAL_KIND):: normal(6,3),dP(6)

do i=1,N_Cell
	!print *, "!!", cell_list(i)%centre(:,1)
    do k=2,nel_z+1
        if ( cell_list(i)%centre(3,1)<=e_Z(k) .AND. &
			cell_list(i)%centre(3,1)>=e_Z(k-1) ) then
			!print *, "!!", el(k-1,:)
            do m=1, size(el,2)
                e=el(k-1,m)
                base(1,1) = Element(e,1)
				base(1,2) = Element(e,4)
				base(1,3) = Element(e,8)
				base(1,4) = Element(e,5)
                normal(1,:) = cross(&
				(vertices(base(1,2),:)-vertices(base(1,1),:)), &
				(vertices(base(1,3),:)-vertices(base(1,1),:)))
                dP(1) = dot_product(normal(1,:),(cell_list(i)%centre(:,1) &
				-vertices(base(1,1),:)))/sqrt(NORM2(normal(1,:)))
				!*******************************************************
				!NOTE that norm2 here is not using the built-in fortran funmction
				!instead it is using norm2 function defined in global
				!because of that I needed to take sqrt of the norm2
				!************************************************************

                base(2,1) = Element(e,1)
				base(2,2) = Element(e,5)
				base(2,3) = Element(e,6)
				base(2,4) = Element(e,2)
                normal(2,:) = cross(&
				(vertices(base(2,2),:)-vertices(base(2,1),:)), &
				(vertices(base(2,3),:)-vertices(base(2,1),:)))
                dP(2) = dot_product(normal(2,:),(cell_list(i)%centre(:,1) &
				-vertices(base(2,1),:)))/sqrt(NORM2(normal(2,:)))

                base(3,1) = Element(e,4)
				base(3,2) = Element(e,3)
				base(3,3) = Element(e,7)
				base(3,4) = Element(e,8)
                normal(3,:) = cross(&
				(vertices(base(3,2),:)-vertices(base(3,1),:)), &
				(vertices(base(3,3),:)-vertices(base(3,1),:)))
                dP(3) = dot_product(normal(3,:),(cell_list(i)%centre(:,1) &
				-vertices(base(3,1),:)))/sqrt(NORM2(normal(3,:)))

                base(4,1) = Element(e,1)
				base(4,2) = Element(e,2)
				base(4,3) = Element(e,3)
				base(4,4) = Element(e,4)
                normal(4,:) = cross(&
				(vertices(base(4,2),:)-vertices(base(4,1),:)), &
				(vertices(base(4,3),:)-vertices(base(4,1),:)))
                dP(4) = dot_product(normal(4,:),(cell_list(i)%centre(:,1) &
				-vertices(base(4,1),:)))/sqrt(NORM2(normal(4,:)))

                base(5,1) = Element(e,5)
				base(5,2) = Element(e,8)
				base(5,3) = Element(e,7)
				base(5,4) = Element(e,6)
                normal(5,:) = cross(&
				(vertices(base(5,2),:)-vertices(base(5,1),:)), &
				(vertices(base(5,3),:)-vertices(base(5,1),:)))
                dP(5) = dot_product(normal(5,:),(cell_list(i)%centre(:,1) &
				-vertices(base(5,1),:)))/sqrt(NORM2(normal(5,:)))

                base(6,1) = Element(e,2)
				base(6,2) = Element(e,6)
				base(6,3) = Element(e,7)
				base(6,4) = Element(e,3)
                normal(6,:) = cross(&
				(vertices(base(6,2),:)-vertices(base(6,1),:)), &
				(vertices(base(6,3),:)-vertices(base(6,1),:)))
                dP(6) = dot_product(normal(6,:),(cell_list(i)%centre(:,1) &
				-vertices(base(6,1),:)))/sqrt(NORM2(normal(6,:)))

				!if ((i==1 .or. i==2) .and. (e==101)) then
					!print *, "Element", Element(e,:)
					!print *, "base", base(1,:)
					!print *, "normal", normal(1,:)
					!print *, "!!", dp(1),dp(2),dp(3),dp(4),dp(5),dp(6)
				!elseif (i==3 .and. e==111) then
				!	print *, "!!", dp(1),dp(2),dp(3),dp(4),dp(5),dp(6)
				!endif

              if (NINT(dP(1))>=0 .AND. NINT(dP(2))>=0 .AND. NINT(dP(3))>=0 &
				.AND. NINT(dP(4))>=0 .AND. NINT(dP(5))>=0 .AND.  NINT(dP(6))>=0) then
                  cell_elem(i)=e
              endif
            enddo
        endif
    enddo
enddo


end subroutine
!----------------------------------------------
subroutine unique(vec,vec_unique)
! Return only the unique values from vec.

implicit none

real (REAL_KIND),intent(in) :: vec(:)
real (REAL_KIND),allocatable,intent(out) :: vec_unique(:)

integer :: i,num
logical,dimension(size(vec)) :: mask

mask = .false.

do i=1,size(vec)

    !count the number of occurrences of this element:
    num = count( vec(i)==vec )

    if (num==1) then
        !there is only one, flag it:
        mask(i) = .true.
    else
        !flag this value only if it hasn't already been flagged:
        if (.not. any(vec(i)==vec .and. mask) ) mask(i) = .true.
    end if

end do

!return only flagged elements:
IF (ALLOCATED (vec_unique)) DEALLOCATE (vec_unique)
allocate( vec_unique(count(mask)) )
vec_unique = pack( vec, mask )

!if you also need it sorted, then do so.
! For example, with slatec routine:
!call ISORT (vec_unique, [0], size(vec_unique), 1)

end subroutine unique
!----------------------------------------------
subroutine grouping(N_Cell,uniquex,uniquez,Cellin,Cellonline,Cellonedge, &
	Celloncorner, size_Cellin,size_Cellonline,size_Cellonedge,size_Celloncornet)

real (REAL_KIND) :: uniquex(:), uniquez(:)
integer :: N_Cell,i, jj, kk, ll, mm, nn
integer :: size_Cellin, size_Cellonline
integer :: size_Cellonedge, size_Celloncornet
real (REAL_KIND), allocatable :: Cellin(:,:)
real (REAL_KIND), allocatable :: Cellboundary(:,:)
real (REAL_KIND), allocatable :: Cellonline(:,:)
real (REAL_KIND), allocatable :: Cellonedge(:,:)
real (REAL_KIND), allocatable :: Celloncorner(:,:)
logical :: answer1, answer2, answer3

jj=0
kk=0
ll=0
mm=0
nn=0
IF (ALLOCATED (Cellin)) DEALLOCATE (Cellin)
allocate(Cellin(N_Cell,3))
IF (ALLOCATED (Cellboundary)) DEALLOCATE (Cellboundary)
allocate(Cellboundary(N_Cell,3))
IF (ALLOCATED (Cellonline)) DEALLOCATE (Cellonline)
allocate(Cellonline(N_Cell,3))
IF (ALLOCATED (Cellonedge)) DEALLOCATE (Cellonedge)
allocate(Cellonedge(N_Cell,3))
IF (ALLOCATED (Celloncorner)) DEALLOCATE (Celloncorner)
allocate(Celloncorner(N_Cell,3))
do i=1,N_Cell
	answer1=isempty_find(uniquex,cell_list(i)%centre(1,1))
	answer2=isempty_find(uniquex,cell_list(i)%centre(2,1))
	answer3=isempty_find(uniquez,cell_list(i)%centre(3,1))
	if (answer1 .AND. answer2 .AND. answer3) then
		jj=jj+1
		!print *, "!!", i
		Cellin(jj,:) = cell_list(i)%centre(:,1)
	!elseif ( abs(cell_list(i)%centre(1,1))==200 .OR. &
	!		abs(cell_list(i)%centre(2,1))==200 ) then
	!	mm=mm+1
	!	Cellboundary(mm,:) =cell_list(i)%centre(:,1)
	elseif (.not.answer1 .AND. answer2 .AND. answer3) then
		kk=kk+1
		Cellonline(kk,:) = cell_list(i)%centre(:,1)
	elseif (answer1 .AND. .not.answer2 .AND. answer3) then
		kk=kk+1
		Cellonline(kk,:) = cell_list(i)%centre(:,1)
	elseif (.not.answer1 .AND. .not.answer2 .AND. .not.answer3) then
		ll=ll+1
		Celloncorner(ll,:)=cell_list(i)%centre(:,1)
	elseif (.not.answer1 .AND. .not.answer2 .AND. answer3) then
		nn=nn+1
		Cellonedge(nn,:) = cell_list(i)%centre(:,1)
	else
		if (answer1 .AND. answer2) then
            kk=kk+1
			Cellonline(kk,:) = cell_list(i)%centre(:,1)
        else
			nn=nn+1
			Cellonedge(nn,:) = cell_list(i)%centre(:,1)
		endif
	endif
end do
size_Cellin=jj
size_Cellonline=kk
size_Cellonedge=nn
size_Celloncornet=ll
!print *, "size_Cellin", size_Cellin
!print *, "size_Cellonline", size_Cellonline
!print *, "size_Cellonedge", size_Cellonedge
!print *, "size_Celloncornet", size_Celloncornet
end subroutine
!----------------------------------------------
function isempty_find(vec,num) result(answer)

real (REAL_KIND) :: vec(:), num
integer :: i,j
logical :: answer
!print *, "size(vec)", size(vec)
j=0
do i=1,size(vec)
	if (vec(i)==num) then
		j=j+1
	endif
enddo

if (j==0) then
	answer = .true.
else
	answer = .false.
end if
end function
!----------------------------------------------
subroutine cellcount(EPm,EpmCoor,EPmIEN,Cell,sizeCell,countedCells)

type(element_type) :: EPm
real (REAL_KIND) :: EpmCoor(:,:), Cell(:,:)
real (REAL_KIND) :: sum_PyVol,Pyramids_volume
integer :: sizeCell,i,j,k,e
integer :: EPmIEN(:,:)
integer, allocatable :: countedCells(:)
integer, allocatable :: base(:,:)
real (REAL_KIND), allocatable :: normal(:,:)
real (REAL_KIND), allocatable :: dP(:)
real (REAL_KIND), allocatable :: face_area(:)

!integer, parameter :: out_unit=20

!open (unit=out_unit,file="CellZCoor.txt",action="write",status="replace")
IF (ALLOCATED (countedCells)) DEALLOCATE (countedCells)
allocate(countedCells(EPm%nel))
countedCells(:) = 0
IF (ALLOCATED (base)) DEALLOCATE (base)
allocate(base(6,4))
IF (ALLOCATED (normal)) DEALLOCATE (normal)
allocate(normal(6,3))
IF (ALLOCATED (dP)) DEALLOCATE (dP)
allocate(dP(6))
IF (ALLOCATED (face_area)) DEALLOCATE (face_area)
allocate(face_area(6))


!print *, "!!", EPm%z_width
do i=1,sizeCell
!i=1
	do k=1,EPm%nel_z
		if (Cell(i,3)<=k*EPm%z_width .AND. Cell(i,3)>=(k-1)*EPm%z_width) then
			!write (*,*) "cell z coor!!", i,Cell(i,3)
			do e=(k-1)*EPm%nel_x*Epm%nel_y+1 , (k)*EPm%nel_x*EPm%nel_y
				sum_PyVol=0.0D+00

				base(1,1) = EPmIEN(e,1)
				base(1,2) = EPmIEN(e,5)
				base(1,3) = EPmIEN(e,6)
				base(1,4) = EPmIEN(e,2)
                normal(1,:) = cross( &
				EpmCoor(base(1,2),:)- EpmCoor(base(1,1),:), &
				EpmCoor(base(1,3),:)-EpmCoor(base(1,1),:) )
				!print *, "normal", normal(1,:)
                dP(1) = dot_product(normal(1,:),(Cell(i,:)- &
									EpmCoor(base(1,1),:) ) )/sqrt(NORM2(normal(1,:)))
                !*******************************************************
				!NOTE that norm2 here is not using the built-in fortran funmction
				!instead it is using norm2 function defined in global
				!because of that I needed to take sqrt of the norm2
				!************************************************************
				!print *, "dp", dP(1)
                face_area(1)=EPm%x_width*EPm%z_width

				base(2,1) = EPmIEN(e,1)
				base(2,2) = EPmIEN(e,3)
				base(2,3) = EPmIEN(e,7)
				base(2,4) = EPmIEN(e,5)
                normal(2,:) = cross( &
				EpmCoor(base(2,2),:)-EpmCoor(base(2,1),:), &
				EpmCoor(base(2,3),:)-EpmCoor(base(2,1),:) )
                dP(2) = dot_product(normal(2,:),(Cell(i,:)- &
									EpmCoor(base(2,1),:) ) )/sqrt(NORM2(normal(2,:)))
                face_area(2)=EPm%x_width*EPm%z_width

				base(3,1) = EPmIEN(e,2)
				base(3,2) = EPmIEN(e,6)
				base(3,3) = EPmIEN(e,8)
				base(3,4) = EPmIEN(e,4)
                normal(3,:) = cross( &
				EpmCoor(base(3,2),:)-EpmCoor(base(3,1),:), &
				EpmCoor(base(3,3),:)-EpmCoor(base(3,1),:) )
                dP(3) = dot_product(normal(3,:),(Cell(i,:)- &
									EpmCoor(base(3,1),:) ) )/sqrt(NORM2(normal(3,:)))
                face_area(3)=EPm%x_width*EPm%z_width

				base(4,1) = EPmIEN(e,1)
				base(4,2) = EPmIEN(e,2)
				base(4,3) = EPmIEN(e,4)
				base(4,4) = EPmIEN(e,3)
                normal(4,:) = cross( &
				EpmCoor(base(4,2),:)-EpmCoor(base(4,1),:), &
				EpmCoor(base(4,3),:)-EpmCoor(base(4,1),:) )
                dP(4) = dot_product(normal(4,:),(Cell(i,:)- &
									EpmCoor(base(4,1),:) ) )/sqrt(NORM2(normal(4,:)))
                face_area(4)=EPm%x_width*EPm%y_width

				base(5,1) = EPmIEN(e,5)
				base(5,2) = EPmIEN(e,7)
				base(5,3) = EPmIEN(e,8)
				base(5,4) = EPmIEN(e,6)
                normal(5,:) = cross( &
				EpmCoor(base(5,2),:)-EpmCoor(base(5,1),:), &
				EpmCoor(base(5,3),:)-EpmCoor(base(5,1),:) )
                dP(5) = dot_product(normal(5,:),(Cell(i,:)- &
									EpmCoor(base(5,1),:) ) )/sqrt(NORM2(normal(5,:)))
                face_area(5)=EPm%x_width*EPm%y_width

				base(6,1) = EPmIEN(e,3)
				base(6,2) = EPmIEN(e,4)
				base(6,3) = EPmIEN(e,8)
				base(6,4) = EPmIEN(e,7)
                normal(6,:) = cross( &
				EpmCoor(base(6,2),:)-EpmCoor(base(6,1),:), &
				EpmCoor(base(6,3),:)-EpmCoor(base(6,1),:) )
                dP(6) = dot_product(normal(6,:),(Cell(i,:)- &
									EpmCoor(base(6,1),:) ) )/sqrt(NORM2(normal(6,:)))
                face_area(6)=EPm%x_width*EPm%z_width

				do j=1,6
					Pyramids_volume = abs(dP(j))*face_area(j)/3.
					!print *, "!!", Pyramids_volume
                    sum_PyVol = sum_PyVol+ Pyramids_volume
				enddo

				!if (i==21 .AND. e==17) then
						!print *, "dp", dP(:)
						!print *, "normal 1!", normal(1,:)
						!print *, "normal 2!", normal(2,:)
						!print *, "normal 3!", normal(3,:)
						!print *, "normal 4!", normal(4,:)
						!print *, "normal 5!", normal(5,:)
						!print *, "normal 6!", normal(6,:)
						!print *, "face_area", face_area(:)
						!print *, "!!", NINT(sum_PyVol)
						!print *, "Cell!", i, Cell(i,:)
				!endif

				if (NINT(sum_PyVol)==EPm%volume) then
                   countedCells(e)=countedCells(e)+1;
					!write (out_unit,*) "e, i, cell coor!!", e, i,Cell(i,:)
				!else
				 !  print *, "sum_PyVol!!", e, sum_PyVol
                endif
			enddo
		endif
	enddo
enddo
!close (out_unit)

!print *, "EPm%volume", EPm%volume
end subroutine
!-----------------------------------------
FUNCTION cross(a, b)
  real (REAL_KIND), DIMENSION(3) :: cross
  real (REAL_KIND), DIMENSION(3), INTENT(IN) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross
!---------------------------------
subroutine refined_k_grid(coordinates,fitted_k,EP,EPCoordinates, &
							EPm_Elemlist, k_new)
real (REAL_KIND) :: coordinates(:,:), EPCoordinates(:,:)
integer :: EPm_Elemlist(:,:)
real (REAL_KIND) :: fitted_k(:,:), basis(8), mat(8)
type(element_type) :: EP
real (REAL_KIND) :: eta1, eta2, eta3
real (REAL_KIND) :: k_new(size(coordinates,1))
integer :: e, i

!integer, parameter :: out_unit=20

!open (unit=out_unit,file="K_ref.txt",action="write",status="replace")

do i = 1,size(coordinates,1)

    call FIND_ELEM(EP%x_min,EP%y_min,EP%z_min,EP%x_width,EP%y_width, &
	EP%z_width,EP%nel_x,EP%nel_y,EP%nel_z,coordinates(i,:), e)
!	use new coordinates to find the element number of the original grid

    eta1 = (coordinates(i,1) - EPCoordinates(EPm_Elemlist(e,1),1))/EP%x_width
    eta2 = (coordinates(i,2) - EPCoordinates(EPm_Elemlist(e,1),2))/EP%y_width
    eta3 = (coordinates(i,3) - EPCoordinates(EPm_Elemlist(e,1),3))/EP%z_width

	!print *, " eta1=", eta1
	!print *, " eta2=", eta2
	!print *, " eta3=", eta3

    if (eta1 < 10**(-10)) then
        eta1 = 0
    endif

    if (eta2 < 10**(-10)) then
        eta2 = 0
    endif

    if (eta3 < 10**(-10)) then
        eta3 = 0
    endif

    basis(1) = (1-eta1)*(1-eta2)*(1-eta3)
	basis(2) = eta1*(1-eta2)*(1-eta3)
	basis(3) = (1-eta1)*eta2*(1-eta3)
	basis(4) = eta1*eta2*(1-eta3)
	basis(5) = (1-eta1)*(1-eta2)*eta3
	basis(6) = eta1*(1-eta2)*eta3
	basis(7) = (1-eta1)*eta2*eta3
	basis(8) = eta1*eta2*eta3

	mat(1) = fitted_k(EPm_Elemlist(e,1),1)
	mat(2) = fitted_k(EPm_Elemlist(e,2),1)
	mat(3) = fitted_k(EPm_Elemlist(e,3),1)
	mat(4) = fitted_k(EPm_Elemlist(e,4),1)
	mat(5) = fitted_k(EPm_Elemlist(e,5),1)
	mat(6) = fitted_k(EPm_Elemlist(e,6),1)
	mat(7) = fitted_k(EPm_Elemlist(e,7),1)
	mat(8) = fitted_k(EPm_Elemlist(e,8),1)

    k_new(i) = dot_product(basis,mat)
	!write (out_unit,*) " k_new=", k_new(i)

enddo
!close (out_unit)

end subroutine
!---------------------------------
subroutine FIND_K_Cylinder(k_elem,ep, EP_Elemlist, EP_coordinates, Cent,CylinCub,k_Cylinder)
real (REAL_KIND) :: k_elem(:), EP_coordinates(:,:)
real (REAL_KIND) ::  Cent(:,:), k_Cylinder(:), basis(8), mat(8) !CyCuK_list(:,:)
integer :: CylinCub(:,:), EP_Elemlist(:,:)
type(element_type) :: EP
real (REAL_KIND) :: coef, x_width, y_width, z_width, eta1,eta2,eta3
integer :: i, e

!integer, parameter :: out_unit=20

!open (unit=out_unit,file="k_Cylinder.txt",action="write",status="replace")
do i =1,size(Cent,1)
    k_Cylinder(i)= k_elem(CylinCub(i,2))
    if (k_Cylinder(i)<0.285) then
        k_Cylinder(i)=0.285
    endif
enddo
!coef=1/(ep%x_width*ep%y_width*ep%z_width);
!x_width=ep%x_width
!y_width=ep%y_width
!z_width=ep%z_width

!do i =1,size(Centroids,1)
 !   e=CylinCub_list(i,2)
 !       eta1 = Centroids(i,1)-min(EP_coordinates(EP_Elemlist(e,1),1), &
!		EP_coordinates(EP_Elemlist(e,2),1), EP_coordinates(EP_Elemlist(e,3),1), &
!		EP_coordinates(EP_Elemlist(e,4),1), EP_coordinates(EP_Elemlist(e,5),1), &
!		EP_coordinates(EP_Elemlist(e,6),1), EP_coordinates(EP_Elemlist(e,7),1), &
!		EP_coordinates(EP_Elemlist(e,8),1))
 !       eta2 = Centroids(i,2)-min(EP_coordinates(EP_Elemlist(e,1),2), &
!		EP_coordinates(EP_Elemlist(e,2),2), EP_coordinates(EP_Elemlist(e,3),2), &
!		EP_coordinates(EP_Elemlist(e,4),2), EP_coordinates(EP_Elemlist(e,5),2), &
!		EP_coordinates(EP_Elemlist(e,6),2), EP_coordinates(EP_Elemlist(e,7),2), &
!		EP_coordinates(EP_Elemlist(e,8),2))
 !       eta3 = Centroids(i,3)-min(EP_coordinates(EP_Elemlist(e,1),3), &
!		EP_coordinates(EP_Elemlist(e,2),3), EP_coordinates(EP_Elemlist(e,3),3), &
!		EP_coordinates(EP_Elemlist(e,4),3), EP_coordinates(EP_Elemlist(e,5),3), &
!		EP_coordinates(EP_Elemlist(e,6),3), EP_coordinates(EP_Elemlist(e,7),3), &
!		EP_coordinates(EP_Elemlist(e,8),3))

!        basis(1) =  coef*(x_width-eta1)*(y_width-eta2)*(z_width-eta3)
!		basis(2) =  coef* eta1*(y_width-eta2)*(z_width-eta3)
 !       basis(3) =  coef*(x_width-eta1)*eta2*(z_width-eta3)
!		basis(4) =  coef*eta1*eta2*(z_width-eta3)
!		basis(5) =  coef*(x_width-eta1)*(y_width-eta2)*eta3
!        basis(6) =  coef* eta1*(y_width-eta2)*eta3
!		basis(7) =  coef* (x_width-eta1)*eta2*eta3
!		basis(8) =  coef* eta1*eta2*eta3

!		mat(1) = k_ref(EP_Elemlist(e,1))
!		mat(2) = k_ref(EP_Elemlist(e,2))
!		mat(3) = k_ref(EP_Elemlist(e,3))
!		mat(4) = k_ref(EP_Elemlist(e,4))
!		mat(5) = k_ref(EP_Elemlist(e,5))
!		mat(6) = k_ref(EP_Elemlist(e,6))
!		mat(7) = k_ref(EP_Elemlist(e,7))
!		mat(8) = k_ref(EP_Elemlist(e,8))

!       k_Cylinder(i)= dot_product(basis,mat)
!       !CyCuK_list(i,:)=[i,e,k_Cylinder(i)]
!	   !write (out_unit,*) " k_Cylinder=", k_Cylinder(i)
!enddo
!close (out_unit)
end subroutine
!----------------------------------
subroutine find_k_mid(k_ref,EP,EP_Elemlist, k_mid)
real (REAL_KIND) :: k_ref(:), k_mid(:)
real (REAL_KIND) :: basis(8), mat(8), eta1,eta2,eta3
integer :: EP_Elemlist(:,:)
type(element_type) :: EP
real (REAL_KIND) :: coef
integer :: e

eta1 = 0.5* EP%x_width
eta2 = 0.5* EP%y_width
eta3 = 0.5* EP%z_width
coef = 1/(EP%x_width*EP%y_width*EP%z_width)
 basis(1) =  coef*(EP%x_width-eta1)*(EP%y_width-eta2)*(EP%z_width-eta3)
 basis(2) =  coef* eta1*(EP%y_width-eta2)*(EP%z_width-eta3)
 basis(3) =  coef*(EP%x_width-eta1)*eta2*(EP%z_width-eta3)
 basis(4) =  coef*eta1*eta2*(EP%z_width-eta3)
 basis(5) =  coef*(EP%x_width-eta1)*(EP%y_width-eta2)*eta3
 basis(6) =  coef* eta1*(EP%y_width-eta2)*eta3
 basis(7) =  coef* (EP%x_width-eta1)*eta2*eta3
 basis(8) =  coef* eta1*eta2*eta3

 do e=1, EP%nel_x*EP%nel_y*EP%nel_z
    mat(1) = k_ref(EP_Elemlist(e,1))
	mat(2) = k_ref(EP_Elemlist(e,2))
	mat(3) = k_ref(EP_Elemlist(e,3))
	mat(4) = k_ref(EP_Elemlist(e,4))
	mat(5) = k_ref(EP_Elemlist(e,5))
	mat(6) = k_ref(EP_Elemlist(e,6))
	mat(7) = k_ref(EP_Elemlist(e,7))
	mat(8) = k_ref(EP_Elemlist(e,8))
    k_mid(e)= dot_product(basis,mat)
enddo

end subroutine
!----------------------------------
subroutine FIND_ELEM(x_min,y_min,z_min,x_width,y_width, &
	z_width,nel_x,nel_y,nel_z,datapoints, e)

real (REAL_KIND) :: datapoints(:)
real (REAL_KIND) :: x_min,y_min,z_min,x_width,y_width,z_width
integer :: nel_x,nel_y,nel_z, e
integer :: num_x,num_y,num_z

	num_x = FLOOR((datapoints(1)-x_min)/x_width)+1;
	num_y = FLOOR((datapoints(2)-y_min)/y_width)+1;
	num_z = FLOOR((datapoints(3)-z_min)/z_width)+1;

if (num_x > nel_x) then
    num_x = nel_x
endif

if (num_y >= nel_y) then
    num_y = nel_y
endif

if (num_z >= nel_z) then
    num_z = nel_z
endif

e = (num_z-1)*nel_x*nel_y + (num_y-1)*nel_x + num_x;

end subroutine
!----------------------------------
end module m_unista
