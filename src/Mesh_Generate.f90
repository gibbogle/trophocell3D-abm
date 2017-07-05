module Mesh_Generate

use global
!use ISO_FORTRAN_ENV, only : REAL32, REAL64, REAL128

contains
!-----------------------------------------------
subroutine face(NODES,ELEMENT,All_Surfaces,&
				nodes2face,element2face,face2element)
integer :: ELEMENT(:,:)
real (REAL_KIND) :: NODES(:,:)
integer :: All_Surfaces(:,:), nodes2face(:,:)
integer :: element2face(:,:),face2element(:,:)
integer :: facecount, rowcount, i, j, k, e, f, row
logical :: answer1, answer2,  answer3, answer4

! nodes2face:
facecount = 0
rowcount = 0
!allocate(nodes2face(6*size(ELEMENT,1),5))

do e = 1,size(ELEMENT,1)

    nodes2face(1+(e-1)*6, 1)=ELEMENT(e,1)
	nodes2face(1+(e-1)*6, 2)=ELEMENT(e,2)
	nodes2face(1+(e-1)*6, 3)=ELEMENT(e,3)
	nodes2face(1+(e-1)*6, 4)=ELEMENT(e,4)
    do i=1,size(All_Surfaces,1)
		answer1=notempty_find(All_Surfaces(i,1:4),nodes2face(1+(e-1)*6,1))
		answer2=notempty_find(All_Surfaces(i,1:4),nodes2face(1+(e-1)*6,2))
		answer3=notempty_find(All_Surfaces(i,1:4),nodes2face(1+(e-1)*6,3))
		answer4=notempty_find(All_Surfaces(i,1:4),nodes2face(1+(e-1)*6,4))
		if (answer1 .AND. answer2 .AND. answer3 .AND. answer4) then
		    nodes2face(1+(e-1)*6, 5) = i
		endif
    enddo

	nodes2face(2+(e-1)*6,1)=ELEMENT(e,1)
	nodes2face(2+(e-1)*6,2)=ELEMENT(e,4)
	nodes2face(2+(e-1)*6,3)=ELEMENT(e,8)
	nodes2face(2+(e-1)*6,4)=ELEMENT(e,5)
	do i=1,size(All_Surfaces,1)
		answer1=notempty_find(All_Surfaces(i,1:4),nodes2face(2+(e-1)*6,1))
		answer2=notempty_find(All_Surfaces(i,1:4),nodes2face(2+(e-1)*6,2))
		answer3=notempty_find(All_Surfaces(i,1:4),nodes2face(2+(e-1)*6,3))
		answer4=notempty_find(All_Surfaces(i,1:4),nodes2face(2+(e-1)*6,4))
		if (answer1 .AND. answer2 .AND. answer3 .AND. answer4) then
		    nodes2face(2+(e-1)*6, 5) = i
		endif
    enddo


    nodes2face(3+(e-1)*6,1)=ELEMENT(e,2)
	nodes2face(3+(e-1)*6,2)=ELEMENT(e,6)
	nodes2face(3+(e-1)*6,3)=ELEMENT(e,7)
	nodes2face(3+(e-1)*6,4)=ELEMENT(e,3)
    do i=1,size(All_Surfaces,1)
		answer1=notempty_find(All_Surfaces(i,1:4),nodes2face(3+(e-1)*6,1))
		answer2=notempty_find(All_Surfaces(i,1:4),nodes2face(3+(e-1)*6,2))
		answer3=notempty_find(All_Surfaces(i,1:4),nodes2face(3+(e-1)*6,3))
		answer4=notempty_find(All_Surfaces(i,1:4),nodes2face(3+(e-1)*6,4))
		if (answer1 .AND. answer2 .AND. answer3 .AND. answer4) then
		    nodes2face(3+(e-1)*6, 5) = i
		endif
    enddo


	nodes2face(4+(e-1)*6,1)=ELEMENT(e,1)
	nodes2face(4+(e-1)*6,2)=ELEMENT(e,5)
	nodes2face(4+(e-1)*6,3)=ELEMENT(e,6)
	nodes2face(4+(e-1)*6,4)=ELEMENT(e,2)
	do i=1,size(All_Surfaces,1)
		answer1=notempty_find(All_Surfaces(i,1:4),nodes2face(4+(e-1)*6,1))
		answer2=notempty_find(All_Surfaces(i,1:4),nodes2face(4+(e-1)*6,2))
		answer3=notempty_find(All_Surfaces(i,1:4),nodes2face(4+(e-1)*6,3))
		answer4=notempty_find(All_Surfaces(i,1:4),nodes2face(4+(e-1)*6,4))
		if (answer1 .AND. answer2 .AND. answer3 .AND. answer4) then
		    nodes2face(4+(e-1)*6, 5) = i
		endif
    enddo

    nodes2face(5+(e-1)*6,1)=ELEMENT(e,4)
	nodes2face(5+(e-1)*6,2)=ELEMENT(e,3)
	nodes2face(5+(e-1)*6,3)=ELEMENT(e,7)
	nodes2face(5+(e-1)*6,4)=ELEMENT(e,8)
    do i=1,size(All_Surfaces,1)
		answer1=notempty_find(All_Surfaces(i,1:4),nodes2face(5+(e-1)*6,1))
		answer2=notempty_find(All_Surfaces(i,1:4),nodes2face(5+(e-1)*6,2))
		answer3=notempty_find(All_Surfaces(i,1:4),nodes2face(5+(e-1)*6,3))
		answer4=notempty_find(All_Surfaces(i,1:4),nodes2face(5+(e-1)*6,4))
		if (answer1 .AND. answer2 .AND. answer3 .AND. answer4) then
		    nodes2face(5+(e-1)*6, 5) = i
		endif
    enddo

    nodes2face(6+(e-1)*6,1)=ELEMENT(e,5)
	nodes2face(6+(e-1)*6,2)=ELEMENT(e,8)
	nodes2face(6+(e-1)*6,3)=ELEMENT(e,7)
	nodes2face(6+(e-1)*6,4)=ELEMENT(e,6)
	do i=1,size(All_Surfaces,1)
		answer1=notempty_find(All_Surfaces(i,1:4),nodes2face(6+(e-1)*6,1))
		answer2=notempty_find(All_Surfaces(i,1:4),nodes2face(6+(e-1)*6,2))
		answer3=notempty_find(All_Surfaces(i,1:4),nodes2face(6+(e-1)*6,3))
		answer4=notempty_find(All_Surfaces(i,1:4),nodes2face(6+(e-1)*6,4))
		if (answer1 .AND. answer2 .AND. answer3 .AND. answer4) then
		    nodes2face(6+(e-1)*6, 5) = i
		endif
    enddo
enddo

nofaces = size(All_Surfaces,1)
!print *, "number of total faces", nofaces

!element2face:
!allocate(element2face(size(ELEMENT,1),6))

do e = 1,size(ELEMENT,1)
    element2face(e,1) = nodes2face(3+(e-1)*6,5)
    element2face(e,2) = nodes2face(2+(e-1)*6,5)
    element2face(e,3) = nodes2face(6+(e-1)*6,5)
    element2face(e,4) = nodes2face(1+(e-1)*6,5)
    element2face(e,5) = nodes2face(5+(e-1)*6,5)
    element2face(e,6) = nodes2face(4+(e-1)*6,5)
enddo


!face2element:
!allocate(face2element(nofaces,6))
do e = 1,nofaces
	face2element(e,:)=(/0,0,0,0,0,0/)
	!print *, "face2element", face2element(e,:)
enddo

do e = 1,size(ELEMENT,1)
    do f = 1,6
		row=nodes2face(f+(e-1)*6,5)
        if (face2element(row,1)==0) then
            face2element(row,1:4) = nodes2face(f+(e-1)*6,1:4)
            face2element(row,5) = e
        else
            face2element(row,6) = e
        endif
    enddo
enddo

CommonFaceNo=0
do j = 1,size(face2element,1)
	if (face2element(j,6)/= 0) then
		CommonFaceNo=CommonFaceNo+1
	endif
enddo
end subroutine
!-----------------------------------------------
function notempty_find(vec,num) result(answer)
integer :: vec(4), num, i
logical :: answer
!print *, "vec", vec
do i=1,4
	if (vec(1)==num .OR. vec(2)==num .OR. &
		vec(3)==num .OR. vec(4)==num) then
		answer = .true.
	else
		answer = .false.
	endif
enddo
end function
!----------------------------------------------
!       ***for generating our own mesh**
! Here we make mesh nodes having the number of
! elements in each direction (x,y,z) and having
! the max and min lengths in mm.
! Node%ID = Node number
! Node%site = (x,y,z) coordinates
!----------------------------------------------
subroutine MESH_Node_3D (Ep,coordinates)

type(element_type) :: Ep
type(Node_type), allocatable, target :: Nodelist(:)
real (REAL_KIND), allocatable :: coordinates(:,:)
type(Node_type), pointer :: Node
integer :: i, j, k, counter
integer :: x, y, z, site(3)
integer :: nel_x ,nel_y , nel_z
integer :: nnp
real (REAL_KIND) :: x_min, y_min, z_min, x_max, y_max, z_max


counter=0
allocate(Nodelist(Ep%nnp))
allocate(coordinates(Ep%nnp,3))
do k= 0,Ep%nel_z
	z=(Ep%z_max - Ep%z_min)/Ep%nel_z*k + Ep%z_min
	do j= 0,Ep%nel_y
		y=(Ep%y_max - Ep%y_min)/Ep%nel_y*j + Ep%y_min
		do i= 0,Ep%nel_x
			x=(Ep%x_max - Ep%x_min)/Ep%nel_x*i + Ep%x_min
			counter=counter+1
			site=(/x,y,z/)
			Node => Nodelist(counter)
			!Node%ID = counter
			Node%site=site
			coordinates(counter,:)=site
			!print *, "Coordinates!=", Node%site
		enddo
	enddo
enddo

!nlist = counter
end subroutine
!--------------------------------------------------------
!             ***for generating our own mesh**
! Here we make elements with node numbers:
! element1 = 8 node numbers arranged in a specific order
!--------------------------------------------------------
subroutine MESH_Element_3D(Ep,Elementlist)

!type(Element_type), pointer  :: Element
type(element_type) :: Ep
integer, allocatable :: Elementlist(:,:)
!integer :: Elist
integer :: nel_x ,nel_y , nel_z, nel
integer :: i, j, k, EN1, EN2, EN3, EN4, EN5, EN6, EN7, EN8
integer :: NodeNo(8), ElCount




allocate(Elementlist(Ep%nel,8))
ElCount=0
do k =	1,Ep%nel_z
	do j = 1,Ep%nel_y
		do i = 1,Ep%nel_x
		    EN1 = i+(Ep%nel_x+1)*(j-1)+(Ep%nel_x+1)*(Ep%nel_y+1)*(k-1)
            EN2 = EN1+1
            EN3 = EN1+Ep%nel_x+1
            EN4 = EN3+1
            EN5 = EN1+(Ep%nel_x+1)*(Ep%nel_y+1)
            EN6 = EN2+(Ep%nel_x+1)*(Ep%nel_y+1)
            EN7 = EN3+(Ep%nel_x+1)*(Ep%nel_y+1)
            EN8 = EN4+(Ep%nel_x+1)*(Ep%nel_y+1)

            ElCount = ElCount+1
			NodeNo=(/EN1, EN2, EN3, EN4, EN5, EN6, EN7, EN8/)
			Elementlist(ElCount,:)= NodeNo
			!print *, "elements Node Numbers!=", Element%NodeNumbers
		enddo
	enddo
enddo
!Elist=ElCount

!do j=1,Elist
!	print *, "elements Node Numbers!=", Elementlist(j,:)
!enddo

end subroutine
!-------------------------------------------------------------------------
!                     ***for generating our own mesh**
!Here we define boundaries of the geometry e.g. in a cube we have
!6 sides: Top, Bottom, Back, Front, Right hand side, and Left hand side
!--------------------------------------------------------------------------
subroutine boundaries(nel, nel_x,nel_y,nel_z,nodes2face,element2face, &
						face2element, Elementlist)

integer, allocatable :: TopFace(:,:)
integer, allocatable :: RHSFace(:,:)
integer, allocatable :: LHSFace(:,:)
integer, allocatable :: BackFace(:,:)
integer, allocatable :: BottomFace(:,:)
integer, allocatable :: FrontFace(:,:)
integer, allocatable :: Elementlist(:,:)
integer :: nodes2face(:,:), element2face(:,:),face2element(:,:)
integer :: count1, count2, count3, count4, count5, count6
integer :: top1,top2, top3, top4, topNodes(4)
integer :: RHS1,RHS2, RHS3, RHS4, RHSNodes(4)
integer :: LHS1,LHS2, LHS3, LHS4, LHSNodes(4)
integer :: Back1,Back2, Back3, Back4, BackNodes(4)
integer :: Bottom1,Bottom2, Bottom3, Bottom4, BottomNodes(4)
integer :: Front1,Front2, Front3, Front4, FrontNodes(4)
integer :: i, j, k
integer :: nel, nel_x ,nel_y , nel_z

integer :: nfT  ! number of nodal nodes at Top surface
integer :: nfLHS  ! number of nodal nodes at LHS surface
integer :: nfRHS  ! number of nodal nodes at RHS surface
integer :: nfBa  ! number of nodal nodes at Back surface
integer :: nfBo  ! number of nodal nodes at Bottom surface
integer :: nfF  ! number of nodal nodes at Front surface

nfT = (nel_x)*(nel_z)
nfLHS = (nel_y)*(nel_z)
nfRHS = (nel_y)*(nel_z)
nfBa = (nel_x)*(nel_y)
nfBo = (nel_x)*(nel_z)
nfF = (nel_x)*(nel_y)

	!Find top faces
	count1=0
	allocate(TopFace(nfT,4))
	do k = 1,nel_z
		do i = 1,nel_x
			count1 = count1+1
			top1 = (nel_x+1)*(nel_y+1)*(k-1)+(nel_x+1)*nel_y+i
			top2 = (nel_x+1)*(nel_y+1)*(k-1)+(nel_x+1)*nel_y+i+1
			top3 = (nel_x+1)*(nel_y+1)*(k)+(nel_x+1)*nel_y+i+1
			top4 = (nel_x+1)*(nel_y+1)*(k)+(nel_x+1)*nel_y+i
			topNodes=(/top1,top2,top3,top4/)
			TopFace(count1,1:4)=topNodes
		enddo
	enddo

	!Find RHS faces
	count2=0
	allocate(RHSFace(nfRHS,4))
	do k = 1,nel_z
		do j = 1,nel_y
			count2 = count2+1
			RHS1 = (k-1)*(nel_x+1)*(nel_y+1) + (nel_x+1)*j
			RHS2 =  k*(nel_x+1)*(nel_y+1) + (nel_x+1)*j
			RHS3 =  k*(nel_x+1)*(nel_y+1) + (nel_x+1)*(j+1)
			RHS4 = (k-1)*(nel_x+1)*(nel_y+1) + (nel_x+1)*(j+1)
			RHSNodes=(/RHS1,RHS2,RHS3,RHS4/)
			RHSFace(count2,1:4) = RHSNodes
		enddo
	enddo

	!Find LHS faces
	count3=0
	allocate(LHSFace(nfLHS,4))
	do k = 1,nel_z
		do j = 1,nel_y
			count3 = count3+1
			LHS1 = (k-1)*(nel_x+1)*(nel_y+1)+(j-1)*(nel_x+1)+1
			LHS2 =  (k-1)*(nel_x+1)*(nel_y+1)+(j)*(nel_x+1)+1
			LHS3 =   k*(nel_x+1)*(nel_y+1)+(j)*(nel_x+1)+1
			LHS4 =  k*(nel_x+1)*(nel_y+1)+(j-1)*(nel_x+1)+1
			LHSNodes=(/LHS1,LHS2,LHS3,LHS4/)
			LHSFace(count3,1:4) = LHSNodes
		enddo
	enddo

	!Find Back faces
	count4=0
	allocate(BackFace(nfBa,4))
	do j = 1,nel_y
		do i = 1,nel_x
			count4 = count4+1
			Back1 = i + (nel_x+1)*(j-1)
			Back2 = i + (nel_x+1)*(j-1) + 1
			Back3 = i + (nel_x+1)*(j) + 1
			Back4 = i + (nel_x+1)*(j)
			BackNodes=(/Back1,Back2,Back3,Back4/)
			BackFace(count4,1:4) = BackNodes
		enddo
	enddo

	!Find Bottom faces
	count5=0
	allocate(BottomFace(nfBo,4))
	do k = 1,nel_z
		do i = 1,nel_x
			count5 = count5+1
			Bottom1 = i+(k-1)*(nel_x+1)*(nel_y+1)
			Bottom2 = i+(nel_x+1)*(nel_y+1)+(k-1)*(nel_x+1)*(nel_y+1)
			Bottom3 = i+1+(nel_x+1)*(nel_y+1)+(k-1)*(nel_x+1)*(nel_y+1)
			Bottom4 = i+1+(k-1)*(nel_x+1)*(nel_y+1)
			BottomNodes=(/Bottom1,Bottom2,Bottom3,Bottom4/)
			BottomFace (count5,1:4)= BottomNodes
		enddo
	enddo

	!Find Front faces
	count6=0
	allocate(FrontFace(nfF,4))
	do j = 1,nel_y
		do i = 1,nel_x
			count6 = count6+1
			Front1 = (nel_x+1)*(nel_y+1)*nel_z+(nel_x+1)*(j-1)+i
			Front2 = (nel_x+1)*(nel_y+1)*nel_z+(nel_x+1)*(j)+i
			Front3 = (nel_x+1)*(nel_y+1)*nel_z+(nel_x+1)*(j)+i+1
			Front4 = (nel_x+1)*(nel_y+1)*nel_z+(nel_x+1)*(j-1)+i+1
			FrontNodes=(/Front1,Front2,Front3,Front4/)
			FrontFace (count6,1:4)= FrontNodes
		enddo
	enddo

	!do j=1,nfT
	!	print *, "boundary faces: TOP", TopFace(:,:)
	!enddo
	call MESH_Face_3D(nel,nel_x,nel_y,nel_z,nodes2face,element2face, &
						face2element, Elementlist)
	!do j=1,nofaces
		!N2Srf => nodes2face(j)
	!	print *, "face2element", face2element(j,:)
	!enddo
	call set_BC(TopFace,RHSFace,LHSFace,BackFace,BottomFace,FrontFace,count1, &
	count2, count3, count4, count5, count6,nodes2face,element2face,face2element, nel)

end subroutine
!-------------------------------------------------------------------------------
!                     ***for generating our own mesh**
! nodes2face: here we list the surfaces of each element with their node numbers
!            (in total we would have ElementNo*6 rows where 6 is the 6 surfaces
!             of a hexahedral element):
!              element1: - - - -, -  => First 4 numbers represent node numbers
!                        - - - -, -     of the surface (RH rule) and 5th one shows
!                        - - - -, -     the face number. Some elements share the
!                        .              same face number but node numbering is
!                        .              reverse for them.
!                        .
! element2face: Each element has 6 faces. Here we list elements' face numbers
! face2element: Each surface has 4 nodes. Here we list faces' node number with
!               neighbouring elements:
!               Face1: - - - -,- - First four numbers represent node numbers
!                      The last 2 numbers show the neighbouring elements (e+,e-)
!--------------------------------------------------------------------------------
subroutine MESH_Face_3D(nel, nel_x,nel_y,nel_z,nodes2face,element2face, &
						face2element,Elementlist)

integer :: facecount, rowcount, i, j, k, e, f, row
integer :: nodes2face(:,:), element2face(:,:),face2element(:,:)
integer :: nel_x ,nel_y , nel_z, nel
integer, allocatable :: Elementlist(:,:)

! nodes2face:
facecount = 0
rowcount = 0
!allocate(nodes2face(6*nel,5))

do k = 1,nel_z
    do j = 1,nel_y

        do i = 1,nel_x
            e = i + nel_x*(j-1) + nel_x*nel_y*(k-1)
            nodes2face(1+(e-1)*6, 1)=Elementlist(e,1)
			nodes2face(1+(e-1)*6, 2)=Elementlist(e,5)
			nodes2face(1+(e-1)*6, 3)=Elementlist(e,6)
			nodes2face(1+(e-1)*6, 4)=Elementlist(e,2)
            if (j==1) then
                facecount = facecount +1
                nodes2face(1+(e-1)*6, 5) = facecount
            else
				nodes2face(1+(e-1)*6, 5) = nodes2face(6+(e-nel_x-1)*6,5)
            endif
        enddo

        do i = 1,nel_x
            e = i + nel_x*(j-1) + nel_x*nel_y*(k-1)
			nodes2face(2+(e-1)*6,1)=Elementlist(e,1)
			nodes2face(2+(e-1)*6,2)=Elementlist(e,3)
			nodes2face(2+(e-1)*6,3)=Elementlist(e,7)
			nodes2face(2+(e-1)*6,4)=Elementlist(e,5)
            if (i == 1) then
                facecount = facecount +1
                nodes2face(2+(e-1)*6,5) = facecount
            else
                nodes2face(2+(e-1)*6,5) = nodes2face(3+(e-2)*6,5)
            endif

            nodes2face(3+(e-1)*6,1)=Elementlist(e,2)
			nodes2face(3+(e-1)*6,2)=Elementlist(e,6)
			nodes2face(3+(e-1)*6,3)=Elementlist(e,8)
			nodes2face(3+(e-1)*6,4)=Elementlist(e,4)
            facecount = facecount +1;
            nodes2face(3+(e-1)*6,5) = facecount;
        enddo

        do i = 1,nel_x
            e = i + nel_x*(j-1) + nel_x*nel_y*(k-1)
			nodes2face(4+(e-1)*6,1)=Elementlist(e,1)
			nodes2face(4+(e-1)*6,2)=Elementlist(e,2)
			nodes2face(4+(e-1)*6,3)=Elementlist(e,4)
			nodes2face(4+(e-1)*6,4)=Elementlist(e,3)
            if (k == 1) then
                facecount = facecount +1
                nodes2face(4+(e-1)*6,5) = facecount;
            else
                nodes2face(4+(e-1)*6,5) = nodes2face(5+(e-nel_x*nel_y-1)*6,5)
            endif
        enddo

        do i = 1,nel_x
            e = i + nel_x*(j-1) + nel_x*nel_y*(k-1)
			nodes2face(5+(e-1)*6,1)=Elementlist(e,5)
			nodes2face(5+(e-1)*6,2)=Elementlist(e,7)
			nodes2face(5+(e-1)*6,3)=Elementlist(e,8)
			nodes2face(5+(e-1)*6,4)=Elementlist(e,6)
            facecount = facecount +1
            nodes2face(5+(e-1)*6,5) = facecount;
        enddo

        do i= 1,nel_x
            e = i + nel_x*(j-1) + nel_x*nel_y*(k-1)
			nodes2face(6+(e-1)*6,1)=Elementlist(e,3)
			nodes2face(6+(e-1)*6,2)=Elementlist(e,4)
			nodes2face(6+(e-1)*6,3)=Elementlist(e,8)
			nodes2face(6+(e-1)*6,4)=Elementlist(e,7)
			facecount = facecount +1
            nodes2face(6+(e-1)*6,5) = facecount
        enddo
    enddo
enddo

nofaces = facecount

!do j=1,6*nel
!	print *, "nodes2face", nodes2face(j,:)
!enddo

!element2face:
!allocate(element2face(nel,6))

do e = 1,nel
    element2face(e,1) = nodes2face(3+(e-1)*6,5)
    element2face(e,2) = nodes2face(2+(e-1)*6,5)
    element2face(e,3) = nodes2face(6+(e-1)*6,5)
    element2face(e,4) = nodes2face(1+(e-1)*6,5)
    element2face(e,5) = nodes2face(5+(e-1)*6,5)
    element2face(e,6) = nodes2face(4+(e-1)*6,5)
enddo

!do j=1,nel
!	print *, "element2face", element2face(j,:)
!enddo

!face2element:
!allocate(face2element(nofaces,6))
do e = 1,nofaces
	face2element(e,:)=(/0,0,0,0,0,0/)
	!print *, "face2element", face2element(e,:)
enddo

do e = 1,nel
    do f = 1,6
		row=nodes2face(f+(e-1)*6,5)
        if (face2element(row,1)==0) then
            face2element(row,1:4) = nodes2face(f+(e-1)*6,1:4)
            face2element(row,5) = e
        else
            face2element(row,6) = e
        endif
    enddo
enddo

!do j=1,nofaces
!	print *, "face2element", face2element(j,:)
!enddo

CommonFaceNo=0
do j = 1,size(face2element,1)
	if (face2element(j,6)/= 0) then
		CommonFaceNo=CommonFaceNo+1
	endif
enddo
!print *, "No of shared faces", CommonFaceNo
end subroutine
!-----------------------------------------------
!        ***for generating our own mesh**
!-----------------------------------------------
subroutine set_BC(TopFace,RHSFace,LHSFace,BackFace,BottomFace,FrontFace,count1, &
 count2, count3, count4, count5, count6,nodes2face,element2face,face2element, nel)

integer :: TopFace(:,:), RHSFace(:,:), LHSFace(:,:), BackFace(:,:)
integer :: BottomFace(:,:), FrontFace(:,:)
integer :: count1, count2, count3, count4, count5, count6, i, nel
integer :: nodes2face(:,:), element2face(:,:),face2element(:,:)
integer, allocatable :: top_face(:)
integer, allocatable :: bottom_face(:)
integer, allocatable :: LHS_face(:)
integer, allocatable :: RHS_face(:)
integer, allocatable :: front_face(:)
integer, allocatable :: back_face(:)

allocate(top_face(count1))
allocate(bottom_face(count5))
allocate(LHS_face(count3))
allocate(RHS_face(count2))
allocate(front_face(count6))
allocate(back_face(count4))

!get face numbers
top_face(:) = find_face(nel,nodes2face,TopFace,count1)
bottom_face(:) = find_face(nel,nodes2face,BottomFace,count5)
LHS_face(:) = find_face(nel,nodes2face,LHSFace,count3)
RHS_face(:) = find_face(nel,nodes2face,RHSFace,count2)
front_face(:) = find_face(nel,nodes2face,FrontFace,count6)
back_face(:) = find_face(nel,nodes2face,BackFace,count4)

!print *, "!!", top_face(:), bottom_face(:), LHS_face(:), RHS_face(:)
!print *, "!!", front_face(:),back_face(:)

!set up dirichlet surfaces
allocate(DCInlet(count3))
allocate(DCOutlet(count2))
DCInlet(:) = LHS_face(:)
DCOutlet(:)= RHS_face(:)


!set up Neumann surfaces - normal flux in/out = neumann wall
allocate(NCWall(count1+count4+count5+count6))
do i=1,count1
	NCWall(i)=top_face(i)
enddo
do i=1,count6
	NCWall(count1+i)=front_face(i)
enddo
do i=1,count4
	NCWall(count1+count6+i)=back_face(i)
enddo
do i=1,count5
	NCWall(count1+count6+count4+i)=bottom_face(i)
enddo
!print *, "!!",NCWall(:)

allocate(NCInlet(count3))
NCInlet(:) = LHS_face(:)
!print *, "!!",NCInlet(:)

end subroutine
!------------------------------------------------
subroutine set_boundary_condition(nel, inlet,outlet,walls,n2faces, &
	DCInlet,DCOutlet, NCWall, NCInlet)
integer :: inlet(:,:),outlet(:,:),walls(:,:),n2faces(:,:)
integer :: DCInlet(:),DCOutlet(:), NCWall(:), NCInlet(:)
integer :: nel

DCInlet(:)= find_face(nel,n2faces,inlet,size(inlet,1))
DCOutlet(:)= find_face(nel, n2faces,outlet,size(outlet,1))
NCWall(:)= find_face(nel, n2faces,walls,size(walls,1))
NCInlet(:)= find_face(nel, n2faces,inlet,size(inlet,1))
!write (*,*) "NCInlet=", NCInlet(:)

end subroutine
!------------------------------------------------
!
!------------------------------------------------
function find_face(nel,n2face,surface,counting) result(DC_surface)
integer :: counting, j, i,nel , surface(:,:), n2face(:,:)
integer :: DC_surface(counting)

do i = 1,counting
    do j = 1,size(n2face,1)!6*nel
		!print *, "!!",surface(i,1), n2face(j,1)
        if ( ( surface(i,1) == n2face(j,1) .AND. surface(i,2) == n2face(j,2) ) .AND.  &
		( surface(i,3) == n2face(j,3) .AND. surface(i,4) == n2face(j,4) ) ) then
           !DC_surface(i) = n2face(j,5)
		   DC_surface(i) = j
        endif
    enddo
enddo

!do i = 1,counting
!	print *, "!!", DC_surface(i)
!enddo
end function
!------------------------------------------------
! This subroutine reads the .msh file generated
! using ICEM meshing software
!-------------------------------------------------
subroutine readmshfile(vertices,InternalFaces,INLET,OUTLET,WALLs)

character * ( 255 ) :: msh_filename
real (REAL_KIND), allocatable :: vertices(:,:)
integer, allocatable :: InternalFaces(:,:), INLET(:,:)
integer, allocatable :: OUTLET(:,:), WALLs(:,:)

!msh_filename = 'CYLINDERCourse10cells.msh'
msh_filename = 'CYLINDERCoarse_Ogrid_11by10.msh'
!get the data
	call gmsh_data_read (msh_filename, vertices,InternalFaces,&
							INLET,OUTLET,WALLs)


end subroutine
!------------------------------------------------
subroutine gmsh_data_read ( gmsh_filename,vertices, &
							InternalFaces,INLET,OUTLET,WALLs)


implicit none

character * ( * ) :: gmsh_filename
character (100) :: cur_line, crap
character*4:: hexaval1,hexaval2,hexaval3,hexaval4,hexaval5,hexaval6
integer, parameter :: line_buf_len= 1024*4
character(LEN=line_buf_len) :: InS
integer :: input, size1, status1, I, startunit
integer :: narray(32), decimaln
integer :: input_stat, ii, results, indx,Startstate, lineNum
real (REAL_KIND) :: x, y, z
integer, allocatable :: InternalFaces(:,:), INLET(:,:)
integer, allocatable :: OUTLET(:,:), WALLs(:,:)
real (REAL_KIND), allocatable :: vertices(:,:)

lineNum=0

    !call get_unit ( input )
    input=nfinput
    !write (*,*) "fale_name: ", gmsh_filename
    open ( unit = input, file = gmsh_filename, status = 'old', iostat = input_stat )
    !write (*,*) "iostat", input_stat
    if ( input_stat .ne. 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'GMSH_DATA_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open input file:' &
		 // trim ( gmsh_filename )
        stop 1
    end if

    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	lineNum=lineNum+1

	DO WHILE ( size1<6 .OR. cur_line(size1-1:size1).NE.')(' )
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
		lineNum=lineNum+1
	end Do
	!print *, "!!", cur_line
	!print *, "length of the character", size


	!*******Start reading the vertices(x,y,z)*******************************
	read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	!write(*,*) 'size1: ',size1
	I=0
	DO WHILE ( size1>6 .OR. cur_line(size1-1:size1).NE.'))' )
		I=I+1
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
		!print *, "!!", cur_line
	end Do

	rewind(input)
	allocate(vertices(I,3))
	do ii=1, lineNum
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) crap
	enddo
	do ii=1,I
		read(input, *) vertices(ii,:)
		lineNum=lineNum+1
	enddo
	!******* END reading the vertices(x,y,z)***********************************


	DO WHILE ( size1<6 .OR. cur_line(size1-1:size1).NE.')(' )
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
		lineNum=lineNum+1
	end Do
	!*******Read internal faces(nod1,node2,node3,node4,RightElement,Leftelement)********
	read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	I=0
	DO WHILE ( size1>6 .OR. cur_line(size1:size1).NE.')' )
		I=I+1
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	end Do
	rewind(input)
	allocate(InternalFaces(I,6))
	do ii=1, lineNum
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) crap
	enddo
	do ii=1,I
		read(input, *) hexaval1,hexaval2,hexaval3,hexaval4,hexaval5,hexaval6
		read(hexaval1, '(z4)') InternalFaces(ii,1)
		read(hexaval2, '(z4)') InternalFaces(ii,2)
		read(hexaval3, '(z4)') InternalFaces(ii,3)
		read(hexaval4, '(z4)') InternalFaces(ii,4)
		read(hexaval5, '(z4)') InternalFaces(ii,5)
		read(hexaval6, '(z4)') InternalFaces(ii,6)
		lineNum=lineNum+1
	enddo
	!******* End reading internal faces *************************************************
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1

    !write(*,*) 'lineNum, size1: ',lineNum,size1
	DO WHILE ( size1<6 .OR. cur_line(size1-1:size1).NE.')(' )
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
		lineNum=lineNum+1
        !write(*,*) 'lineNum, size1: ',lineNum,size1
	end Do
	!*******Read INLET faces(nod1,node2,node3,node4,RightElement,Leftelement)********
	read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	I=0
	DO WHILE ( size1>6 .OR. cur_line(size1:size1).NE.')' )
		I=I+1
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	end Do
	rewind(input)
	allocate(INLET(I,6))
	do ii=1, lineNum
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) crap
	enddo
	do ii=1,I
		read(input, *) hexaval1,hexaval2,hexaval3,hexaval4,hexaval5,hexaval6
		read(hexaval1, '(z4)') INLET(ii,1)
		read(hexaval2, '(z4)') INLET(ii,2)
		read(hexaval3, '(z4)') INLET(ii,3)
		read(hexaval4, '(z4)') INLET(ii,4)
		read(hexaval5, '(z4)') INLET(ii,5)
		read(hexaval6, '(z4)') INLET(ii,6)
		lineNum=lineNum+1
	enddo
	!******* End reading INLET faces *************************************************
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
	DO WHILE ( size1<6 .OR. cur_line(size1-1:size1).NE.')(' )
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
		lineNum=lineNum+1
	end Do
	!*******Read OUTLET faces(nod1,node2,node3,node4,RightElement,Leftelement)********
	read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	I=0
	DO WHILE ( size1>6 .OR. cur_line(size1:size1).NE.')' )
		I=I+1
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	end Do
	rewind(input)
	allocate(OUTLET(I,6))
	do ii=1, lineNum
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) crap
	enddo
	do ii=1,I
		read(input, *) hexaval1,hexaval2,hexaval3,hexaval4,hexaval5,hexaval6
		read(hexaval1, '(z4)') OUTLET(ii,1)
		read(hexaval2, '(z4)') OUTLET(ii,2)
		read(hexaval3, '(z4)') OUTLET(ii,3)
		read(hexaval4, '(z4)') OUTLET(ii,4)
		read(hexaval5, '(z4)') OUTLET(ii,5)
		read(hexaval6, '(z4)') OUTLET(ii,6)
		lineNum=lineNum+1
	enddo
	!******* End reading OUTLET faces *************************************************
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1
    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
    lineNum=lineNum+1

	DO WHILE ( size1<6 .OR. cur_line(size1-1:size1).NE.')(' )
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
		lineNum=lineNum+1
	end Do
	!*******Read WALL faces(nod1,node2,node3,node4,RightElement,Leftelement)********
	read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	I=0
	DO WHILE ( size1>6 .OR. cur_line(size1:size1).NE.')' )
		I=I+1
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	end Do
	rewind(input)
	allocate(WALLs(I,6))
	do ii=1, lineNum
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) crap

	enddo
	do ii=1,I
		read(input, *) hexaval1,hexaval2,hexaval3,hexaval4,hexaval5,hexaval6
		read(hexaval1, '(z4)') WALLs(ii,1)
		read(hexaval2, '(z4)') WALLs(ii,2)
		read(hexaval3, '(z4)') WALLs(ii,3)
		read(hexaval4, '(z4)') WALLs(ii,4)
		read(hexaval5, '(z4)') WALLs(ii,5)
		read(hexaval6, '(z4)') WALLs(ii,6)
		!print *, "WALL", WALL(ii,:)
		lineNum=lineNum+1
	enddo
	!******* End reading WALL faces *************************************************

close (input)
end subroutine

!---------------------------------------------------------------
! This subroutine reads in the mesh properties from a .k mesh
! file which is generated with ICEM meshing software.
! these properties mainly include vertices coordinates and elements
! including the vertices numbers according their ordering.
!------------------------------------------------------------------
subroutine readKfile(NODES,ELEMENT)
character * ( 255 ) :: k_filename
integer, allocatable :: ELEMENT(:,:)
real (REAL_KIND), allocatable :: NODES(:,:)

!k_filename = 'CYLINDERCourse10cells.k'
k_filename = 'CYLINDERCoarse_Ogrid_11by10.k'
!get the data
	call k_data_read(k_filename,NODES,ELEMENT)


end subroutine
!--------------------------------------------------------------
subroutine k_data_read(k_filename,NODES,ELEMENT)
character * ( * ) :: k_filename
character (100) :: cur_line, crap
integer :: input,status1,size1,lineNum
integer :: input_stat,I, ii, jj
integer, allocatable :: ELEMENT(:,:)
real (REAL_KIND), allocatable :: NODES(:,:)

lineNum=0

    !call get_unit ( input )
    input=nfinput
    open ( unit = input, file = k_filename, status = 'old', &
	 iostat = input_stat )
    ! write (*,*) "iostat", input_stat
    if ( input_stat .ne. 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'KMSH_DATA_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open input file:' &
		 // trim ( k_filename )
        stop 1
    end if

    read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	lineNum=lineNum+1
	DO WHILE ( cur_line(1:size1).NE.'*NODE' )
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
		lineNum=lineNum+1
	end Do
	!print *, "!!", cur_line
	!print *, "Line#", lineNum
	read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	lineNum=lineNum+1

	!*******Start reading the NODES(#,x,y,z)*******************************
	read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	I=0
	DO WHILE ( cur_line(size1:size1).NE.'$' )
		I=I+1
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
		!print *, "!!", cur_line
	end Do

	rewind(input)
	allocate(NODES(I,4))
	do ii=1, lineNum
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) crap
	enddo
	do ii=1,I
		read(input, *) NODES(ii,:)
		!print *, "NODES", NODES(ii,:)
		lineNum=lineNum+1
	enddo
	!******* END reading the NODEs(#,x,y,z)***********************************

	read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	lineNum=lineNum+1
	DO WHILE ( cur_line(1:size1).NE.'*ELEMENT_SOLID' )
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
		lineNum=lineNum+1
	end Do
	!print *, "!!", cur_line

	!*******Start reading the ELEMENT(node1,node2...node8)*******************************
	read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	I=0
	DO WHILE ( cur_line(1:size1).NE.'*ELEMENT_SHELL' )
		I=I+1
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
		!print *, "!!", cur_line
	end Do

	rewind(input)
	allocate(ELEMENT(I/2,8))
	do ii=1, lineNum
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) crap
	enddo
	jj=0
	do ii=1,I
	if (mod(ii,2)==0)then
		jj=jj+1
		read(input, *) ELEMENT(jj,:)
		!print *, "ELEMENT", ELEMENT(jj,:)
		lineNum=lineNum+1
	else
		read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) crap
		lineNum=lineNum+1
	endif
	enddo
	!print *, "Line#", lineNum
	!******* END reading the ELEMENTs(x,y,z)***********************************

close (input)

end subroutine
!---------------------------------------------------------------
 subroutine hex_to_decimal(hex, decimaln)
    character*(*) :: hex, decimaln
    integer :: i
    read (hex,'(z4)') i
    write (decimaln,*) i
    return
end subroutine

!---------------------------------------------------------------
! This works by replacing the specified non-number characters
! with ' '.
! It is not completely flexible - you need to know how many
! numbers (n) to read from the string with spaces.
!---------------------------------------------------------------
subroutine parse(string,narray,n)
implicit none
character*(100) :: string
integer :: narray(*)
integer :: n, loc1, loc2

do
	loc1 = index (string,'(')
	if (loc1 > 0) then
		string(loc1:loc1) = ' '
	endif
	loc2 = index (string,')')
	if (loc2 > 0) then
		string(loc2:loc2) = ' '
	endif
    if (loc1 == 0 .and. loc2 == 0) exit
enddo
write(*,*) 'new string: ',string
read(string,*)  narray(1:3)
read(string,'(z4)')  narray(4:4)
read(string,*)  narray(5:n)
end subroutine
!-------------------------------------------------
  elemental subroutine str2int(str,results,sta)
    implicit none
    ! Arguments
    character(len=*),intent(in) :: str
    integer,intent(out)         :: results
    integer,intent(out)         :: sta

    read(str,'(a)',iostat=sta)  results
  end subroutine str2int

!-------------------------------------------------
subroutine get_unit ( iunit )

!*********************************************************************
!
! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
      implicit none

      integer i
      integer iunit
      logical value

      iunit = 0

      do i = 1, 99

        if ( i .ne. 5 .and. i .ne. 6 .and. i .ne. 9 ) then

          inquire ( unit = i, opened = value, err = 10 )

          if ( .not. value ) then
            iunit = i
            return
          end if

        end if

10      continue

      end do

      return
end subroutine
!------------------------------------------
function linspace(d1,d2,n)

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
real (REAL_KIND), INTENT(IN) :: d1, d2
real (REAL_KIND) :: linspace(n)

INTEGER :: indxi


linspace(1) = d1
DO indxi= 0,n-2
   linspace(indxi+1) = d1+(DBLE(indxi)*(d2-d1))/DBLE(n-1)
END DO
linspace(n) = d2

!MATLAB
!linspace = [d1+(0:n-2)*(d2-d1)/(floor(n)-1) d2];


end function
!--------------------------------------------------
subroutine meshgrid(xgv, ygv, zgv, X, Y, Z)
  implicit none
  real,intent(in)   :: xgv(:), ygv(:), zgv(:)
  real,intent(out)  :: X(:,:,:), Y(:,:,:), Z(:,:,:)
  integer           :: sX, sY, sZ, i

  sX = size(xgv) ; sY = size(ygv) ; sZ = size(zgv)

  do i=1,sZ
    X(:,:,i) = spread( xgv, 1, sY )
    Y(:,:,i) = spread( ygv, 2, sX )
  enddo ! i
  do i=1,sX
    Z(i,:,:) = spread( zgv, 1, sY)
  enddo ! i
end subroutine
!--------------------------------------------------
subroutine Read_face(file_name, NoRow, Nocol,n2f)
character * ( * ) :: file_name
integer :: input, status1,size1
integer :: input_stat,I, ii,Nocol,NoRow
integer, allocatable :: n2f(:,:)
character (100) :: cur_line

	!call get_unit ( input )
    input=nfinput
    open ( unit = input, file = file_name, status = 'old', &
	 iostat = input_stat )

    if ( input_stat .ne. 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'FACE_DATA_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open input file "' &
		 // trim ( file_name ) // '"'
        stop 1
    end if
	allocate(n2f(NoRow,Nocol))

	!read(input, '(a)', ADVANCE='NO', iostat=status1, SIZE=size1) cur_line
	!print *, "", cur_line
    do ii=1,NoRow
		read(input, *) n2f(ii,:)
		!print *, "n2f", n2f(ii,:)
	enddo
close (input)
end subroutine
!--------------------------------------------------
subroutine Read_real(fileplace,file_name, NoRow, Nocol,n2f)

character * ( * ) :: file_name, fileplace
integer :: input, status1,size1
integer :: input_stat,I, ii,Nocol,NoRow
integer :: line_no
real (REAL_KIND), allocatable :: n2f(:,:)

IF (ALLOCATED (n2f)) DEALLOCATE (n2f)

	!call get_unit ( input )
	input=nfinput
	!write (*,*) "filename:", fileplace//file_name
    open ( unit = input, file = trim(fileplace)//trim(file_name), status = 'old', &
	 iostat = input_stat )
    !write (*,*) "iostat", input_stat
    if ( input_stat .ne. 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'REAL_READ_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open input file:' &
		 // trim ( file_name )
        stop 1
    end if
	allocate(n2f(NoRow,Nocol))

    do ii=1,NoRow
		read(input, *) n2f(ii,:)
		!print *, "n2f", n2f(ii,:)
	enddo
	!print *, "size n2f", size(n2f,1), size(n2f,2)
close (input)
end subroutine
!--------------------------------------------------
subroutine Sampling_grid(Ep,coordinates,Elemlist)

!integer :: nel, nel_x ,nel_y , nel_z, nnp
!real (REAL_KIND) :: x_min, x_max, y_min, y_max, z_min, z_max
!real (REAL_KIND) :: x_width, y_width, z_width, volume
!integer :: nsd, nen ,ndof, neq
type(element_type) :: Ep
real (REAL_KIND), allocatable :: coordinates(:,:)
integer, allocatable :: Elemlist(:,:)


Ep%nsd  = 3                             ! number of space dimensions
Ep%nel_z = 50                       ! number of elements on z axis
Ep%nel_y = 20                       ! number of elements on y axis
Ep%nel_x = 20                      ! number of elements on x axis
Ep%nel  = Ep%nel_y*Ep%nel_x*Ep%nel_z                 ! number of elements
Ep%nnp  = (Ep%nel_z+1)*(Ep%nel_y+1)*(Ep%nel_x+1)     ! number of nodal nodes
Ep%nen  = 8                             ! number of element nodes
Ep%ndof = 1                            ! degrees-of-freedom per node
Ep%neq  = Ep%nnp*Ep%ndof                ! number of equations
Ep%x_min = -200             !min x coordinates - in mm
Ep%x_max = 200          !max x coordinates - in mm
Ep%y_min = -200             !min y coordinates - in mm
Ep%y_max = 200            !max y coordinates - in mm
Ep%z_min = 0             !min z coordinates - in mm
Ep%z_max = 1000            !max z coordinates - in mm
Ep%x_width = (Ep%x_max-Ep%x_min)/Ep%nel_x
Ep%y_width = (Ep%y_max-Ep%y_min)/Ep%nel_y
Ep%z_width = (Ep%z_max-Ep%z_min)/Ep%nel_z
Ep%volume = Ep%x_width*Ep%y_width*Ep%z_width
!Jx = nel_x/(x_max-x_min)
!Jy = nel_y/(y_max-y_min)
!Jz = nel_z/(z_max-z_min)

	call  MESH_Node_3D (Ep,coordinates)
	call MESH_Element_3D(Ep,Elemlist)


end subroutine
!------------------------------------------------------
subroutine Sampling_grid_multiple(unit_x, unit_y, unit_z, Ep, Epm,&
									coordinates,Elemlist)

integer :: unit_x, unit_y, unit_z
!integer :: EP_nel_x ,EP_nel_y, EP_nel_z
!real (REAL_KIND) :: EP_x_min, EP_x_max, EP_y_min, EP_y_max, EP_z_min, EP_z_max
!integer :: nel, nel_x ,nel_y , nel_z, nnp
!real (REAL_KIND) :: x_min, x_max, y_min, y_max, z_min, z_max
!real (REAL_KIND) :: x_width, y_width, z_width, volume
!integer :: nsd, nen ,ndof, neq
type(element_type) :: Ep
type(element_type) :: Epm
real (REAL_KIND), allocatable :: coordinates(:,:)
integer, allocatable :: Elemlist(:,:)

Epm%nel_z = EP%nel_z/unit_z                       ! number of elements on z axis
Epm%nel_y = EP%nel_y/unit_y                       ! number of elements on y axis
Epm%nel_x = EP%nel_x/unit_x                       ! number of elements on x axis
!print *, "number of elements on z axis=", nel_z
!print *, "number of elements on y axis=", nel_y
!print *, "number of elements on x axis=", nel_x

Epm%nel  = Epm%nel_y*Epm%nel_x*Epm%nel_z                      ! number of elements
Epm%nnp  = (Epm%nel_z+1)*(Epm%nel_y+1)*(Epm%nel_x+1)          ! number of nodal nodes

Epm%nen  = 8                             ! number of element nodes
Epm%ndof = 1                             ! degrees-of-freedom per node
Epm%neq  = Epm%nnp*Epm%ndof                      ! number of equations

! mesh dimensions
Epm%x_min = EP%x_min                 !min x coordinates - in mm
Epm%x_max = EP%x_max          !max x coordinates - in mm
Epm%y_min = EP%y_min                 !min y coordinates - in mm
Epm%y_max = EP%y_max          !max y coordinates - in mm
Epm%z_min = EP%z_min                 !min z coordinates - in mm
Epm%z_max = EP%z_max          !max z coordinates - in mm

Epm%x_width = (Epm%x_max-Epm%x_min)/Epm%nel_x
Epm%y_width = (Epm%y_max-Epm%y_min)/Epm%nel_y
Epm%z_width = (Epm%z_max-Epm%z_min)/Epm%nel_z
Epm%volume = Epm%x_width*Epm%y_width*Epm%z_width;

	call  MESH_Node_3D (Epm, coordinates)
	call MESH_Element_3D(Epm,Elemlist)
 call mesh_cylinder_volume(EPm%volume,Elemlist,EPm%nel_x, &
			EPm%nel_y,EPm%nel_z,coordinates,EPm%vol_fraction)
!print *, "vol_fraction=", EPm%vol_fraction
end subroutine
!------------------------------------------------------
subroutine mesh_cylinder_volume(volume,Elemlist,nel_x, &
			nel_y,nel_z,coordinates,VolDist)

integer :: nel_x ,nel_y , nel_z, total_elems
integer :: Elemlist(:,:), indices(8)
integer :: i , j
real (REAL_KIND) :: x(8),x1, x0,y(8), y1, y0,z(8), z1, z0
real (REAL_KIND) :: angle_between, segmentprop, triangle_height
real (REAL_KIND) :: chord_length, segment_area,corner_triangle_area
real (REAL_KIND) :: approx_circle_area, volume, coordinates(:,:)
real (REAL_KIND), allocatable :: VolDist(:)

total_elems=nel_x*nel_y*nel_z
allocate(VolDist(total_elems))
VolDist=1.
do i=1,total_elems
   indices(:)=Elemlist(i,:)
   do j=1, 8
		x(j)=abs(coordinates(indices(j),1))
		y(j)=abs(coordinates(indices(j),2))
		z(j)=abs(coordinates(indices(j),3))
   enddo
   x1=max(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8))
   x0=min(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8))
   y1=max(y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8))
   y0=min(y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8))
   z1=max(z(1),z(2),z(3),z(4),z(5),z(6),z(7),z(8))
   z0=min(z(1),z(2),z(3),z(4),z(5),z(6),z(7),z(8))
   if ((sqrt(x1**2+y1**2)) < 200) then
       VolDist(i)=1.
   elseif ((sqrt(x1**2+y0**2)>=200).AND.(sqrt(x0**2+y1**2)>=200)) then
       angle_between=atan2(sqrt(200**2-x0**2),x0)-atan2(y0,sqrt(200**2-y0**2))
       segmentprop=abs((angle_between)/(2*Pi))
       triangle_height=200*cos(0.5*(angle_between))
       chord_length=200*2*sin(0.5*(angle_between))
       segment_area=Pi*200**2*segmentprop-0.5*triangle_height*chord_length
       corner_triangle_area=0.5*((sqrt(200**2-y0**2)-x0)*(sqrt(200**2-x0**2)-y0))
    approx_circle_area=corner_triangle_area+segment_area
    VolDist(i)=approx_circle_area*(z1-z0)/volume

   elseif ((sqrt(x1**2+y0**2)<200).AND.(sqrt(x0**2+y0**2)<200) ) then
		angle_between=atan2(sqrt(200**2-x0**2),x0)-atan2(sqrt(200**2-x1**2),x1)
       segmentprop=abs((angle_between)/(2*Pi))
       triangle_height=200*cos(0.5*(angle_between))
       chord_length=200*2*sin(0.5*(angle_between))
       segment_area=Pi*200**2*segmentprop-0.5*triangle_height*chord_length
       approx_circle_area=segment_area+(min(sqrt(200**2-x0**2), &
	   sqrt(200**2-x1**2))-y0)*(x1-x0)+abs(sqrt(200**2-x0**2)-sqrt(200**2-x1**2))*(x1-x0)*0.5
       VolDist(i)=approx_circle_area*(z1-z0)/volume
   elseif ((sqrt(x0**2+y1**2)<=200).AND.(sqrt(x0**2+y1**2)<=200) ) then
        angle_between=atan2(sqrt(200**2-y0**2),y0)-atan2(sqrt(200**2-y1**2),y1)
       segmentprop=abs((angle_between)/(2*Pi))
       triangle_height=200*cos(0.5*(angle_between))
       chord_length=200*2*sin(0.5*(angle_between))
       segment_area=Pi*200**2*segmentprop-0.5*triangle_height*chord_length
       approx_circle_area=segment_area+(min(sqrt(200**2-y0**2), &
	   sqrt(200**2-y1**2))-x0)*(y1-y0)+abs(sqrt(200**2-y0**2)-sqrt(200**2-y1**2))*(y1-y0)*0.5
       VolDist(i)=approx_circle_area*(z1-z0)/volume
   endif
enddo


end  subroutine
!-----------------------------------------------------
end module
