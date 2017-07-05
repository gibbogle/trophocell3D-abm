module sparse

use global
use Mesh_Generate
use sparsekit
use m_unista
use umf_solver


contains
!------------------------------------------
!
!------------------------------------------

subroutine READ_Bmatrix(nel)

integer :: nel
character * ( 255 ) :: fileplace, Bmatrix_filename

fileplace= "./B_matrix/"

do j = 1,nel

   	write (Bmatrix_filename,'("B_matrix_",I0,".txt")') j
	call Read_real(fileplace,Bmatrix_filename,6,6,B_MATRIX(j)%Inside_B)
enddo

end subroutine
!------------------------------------------
!
!------------------------------------------
subroutine sparseMatrices(nel,element2face,All_Surfaces,k,A, IA, JA)

integer :: nel
real(REAL_KIND) :: B(6,6)
real(REAL_KIND) :: k(:), BB(6*6*nel-CommonFaceNo), C(6)
real(REAL_KIND) :: CC(6*nel)
real(REAL_KIND) :: A(8*6*nel-CommonFaceNo)
integer :: IB(6*6*nel-CommonFaceNo), JB(6*6*nel-CommonFaceNo)!, ones(6)
integer :: IC(6*nel), JC(6*nel),All_Surfaces(:,:), element2face(:,:)
integer :: IA(8*6*nel-CommonFaceNo),JA(8*6*nel-CommonFaceNo)
!real(REAL_KIND), allocatable :: Ak(:,:)
integer :: ii,j, I(6), signum(6), DiagSig(6,6)
integer :: counter, Recount1, Ccounter, commonCount
logical :: ok1
character * ( 255 ) :: fileplace, Bmatrix_filename
!integer, parameter :: out_unit=20

!open (unit=out_unit,file="missed.txt",action="write",status="replace")


!allocate(Ak(6,6))
!Ak(1,:)= 1/(2.*6.) *(/2.,-1.,0.,0.,0.,0./)
!Ak(2,:)= 1/(2.*6.) *(/-1.,2.,0.,0.,0.,0./)
!Ak(3,:)= 1/(2.*6.) *(/0.,0.,2.,-1.,0.,0./)
!Ak(4,:)= 1/(2.*6.) *(/0.,0.,-1.,2.,0.,0./)
!Ak(5,:)= 1/(2.*6.) *(/0.,0.,0.,0.,2.,-1./)
!Ak(6,:)= 1/(2.*6.) *(/0.,0.,0.,0.,-1.,2./)

fileplace= "./B_matrix/"

Recount1=0
!ones=(/1,1,1,1,1,1/)
counter=0
Ccounter=0
commonCount=0
do j = 1,nel
   I(:) = element2face(j,:) ! diag matrix of face numbers of the element
                         ! gets face numbers in the element
	!print *, "I=", I(:)
	signum=(/1,1,1,1,1,1/)
	do ii=1,6
		if (j==All_Surfaces(I(ii),6)) then
			signum(ii)=-1 ! -1 if edge is on T-
		endif
	enddo
   !print *, "signum=", signum(:)
    !k = 1 !permeability

	! Need to times hydraulic conductivity here
	DiagSig(:,:) = diag(signum(:),size(signum(:)))

	B(:,:) = blood_viscosity/k(j)*matmul(matmul(DiagSig, B_MATRIX(j)%Inside_B),DiagSig)
	!print *, "B=", B
	C(:)=matmul(ones(6),diag(signum,size(signum(:))))
	do ii=1,6
		do jj =1,6
			!if (j==3619) then
				!print *, "!!", B(ii,jj)/= 0
				!print *, "!!", B(ii,jj)
			!endif
			if (B(ii,jj)/= 0 ) then
				if (counter > 6) then
					ok1 = .false.
					!if (j==3619) then
						!write (out_unit,*) "IB, JB=", I(ii), I(jj)
					!endif
					Recount1=findequal(IB(1:counter),JB(1:counter),I(ii),I(jj),ok1)
					!if (j==3619) then
					!	write (out_unit,*) "Recount1, ok=",Recount1,ok1
					!endif

					if ( ok1 ) then
						BB(Recount1) = B(ii,jj)+BB(Recount1)
					else
						counter=counter+1
						BB(counter) = B(ii,jj)
						IB(counter) = I(ii)
						JB(counter) = I(jj)
					endif
				else
						counter=counter+1
						BB(counter) = B(ii,jj)
						IB(counter) = I(ii)
						JB(counter) = I(jj)
				endif
			endif
		enddo
		!print *, "count=", counter
		Ccounter = Ccounter+1
		CC(Ccounter)= C(ii)
		IC(Ccounter)=I(ii)
		JC(Ccounter)=j
	enddo
enddo

!print *, "BB=", BB(1:10)
!print *, "IB=", IB(1:10)
!print *, "JB=", JB(1:10)

!do j=1,size(BB(:))
!	A(j) = BB(j)
!enddo
!do j=1,size(CC(:))
!	A(size(BB(:))+j) = CC(j)
!enddo
!do j=1,size(CC(:))
!	A(size(BB(:))+size(CC(:))+j) = CC(j)
!enddo

! Global stiffness matrix A
A(:) =[BB(:), CC(:),CC(:)]
IA(:) = [IB(:),IC(:),JC(:)+nofaces*ones(size(JC(:)))]
JA(:) = [JB(:),JC(:)+nofaces*ones(size(JC(:))),IC(:)]
!do j=1, 8*6*nel-CommonFaceNo
!    write (out_unit,*) " Matrix A=", A(j), IA(j), JA(j)
!enddo
 !close (out_unit)

end subroutine

!-----------------------------------------------
!
!-----------------------------------------------
function diag(matrix,rank) result(diagmatrix)
integer :: rank,i
integer :: matrix(rank)
real(REAL_KIND) :: diagmatrix(rank,rank)


diagmatrix(:,:) = 0
!print *, "zeros=", diagmatrix(:,:)
do i = 1, rank
   diagmatrix(i,i) = matrix(i)
end do
end function
!-----------------------------------------------
!
!-----------------------------------------------
function findequal(IB, JB,Row,Col,ok)

integer :: IB(:), JB(:), Row,Col , findequal, ii
logical :: ok
 ok = .false.
do ii=1,size(IB)
	if (IB(ii)==Row .AND. JB(ii)==Col) then
		findequal=ii
		!print *, "recount", findequal
		ok = .true.
		return
	endif
enddo

end function
!-----------------------------------------------
!
!-----------------------------------------------
function ones(Msize)

integer :: Msize, i, ones(Msize)
	do i= 1, Msize
		ones(i) = 1
	enddo
!print *, "ones=", ones(:)
end function
!-----------------------------------------------
!
!-----------------------------------------------
subroutine Compute_BodyForce(nel,F,x, FFaces)
integer :: nel
real (REAL_KIND) :: F(:), x(:)
integer :: i, tmp(nofaces+nel), FFaces(nofaces+nel-size(NCWall))

F(:) = 0
x(:) = 0
tmp(:) = 0
! DirichletINLET
do i= 1, size(DCInlet)
	F(DCInlet(i)) = outletPressure!inletPressure
enddo
! Dirichlet OUTLET
do i= 1, size(DCOutlet)
	F(DCOutlet(i)) = inletPressure!outletPressure
enddo

!Neumann Wall
do i =1,size(NCWall)
	x(NCWall(i)) = 0
enddo

do i =1,size(NCWall)
	tmp(NCWall(i)) = 1
enddo
FFaces(:) = findNONzeros(nel,tmp(:))
!print *, "FreeFaces=", FreeFaces(:)
end subroutine
!-----------------------------------------------
!
!-----------------------------------------------
function findNONzeros(nel, TMP)

integer :: nel
integer :: TMP(:), findNONzeros(nofaces+nel-size(NCWall))
integer :: i, count1

count1=0
do i = 1, size(TMP(:))
	if (TMP(i)==0) then
	count1=count1+1
		findNONzeros(count1) = i
	endif
enddo


end function

!-------------------------------------------------
!
!-------------------------------------------------
subroutine FEsolve
integer :: unit_x,unit_y,unit_z
integer :: i, j, k, e,  count_empty, count1
real (REAL_KIND) :: k_empty, cell_radius
type(element_type) :: EP
type(element_type) :: EPm
real (REAL_KIND), allocatable :: EP_coordinates(:,:)
integer, allocatable :: EP_Elemlist(:,:)
real (REAL_KIND), allocatable :: EPm_coordinates(:,:)
integer, allocatable :: EPm_Elemlist(:,:)
real (REAL_KIND), allocatable :: uniqueEPm_x(:)
real (REAL_KIND), allocatable :: uniqueEPm_z(:)
real (REAL_KIND), allocatable :: Cell_in(:,:)
real (REAL_KIND), allocatable :: Cell_online(:,:)
real (REAL_KIND), allocatable :: Cell_onedge (:,:)
real (REAL_KIND), allocatable :: Cell_oncorner(:,:)
integer :: size_Cell_in, size_Cell_online,size_Cell_onedge, size_Cell_oncorner
integer, allocatable :: cellcount_inside_elem(:)
integer, allocatable :: cellcount_online_elem(:)
integer, allocatable :: cellcount_onedge_elem(:)
integer, allocatable :: cellcount_oncorner_elem(:)
real (REAL_KIND), allocatable :: vol_ratio(:), k_conduct(:)
real (REAL_KIND), allocatable :: k_nonempty(:)
real (REAL_KIND), allocatable :: k_node_mat(:,:)
real (REAL_KIND), allocatable :: k_ref(:)
real (REAL_KIND), allocatable :: k_elem(:)
real (REAL_KIND), allocatable :: k_Cylinder_elem(:)
!real (REAL_KIND), allocatable :: CCK_list(:,:)

!Stiffness Matrix
real (REAL_KIND), allocatable :: AMatrix(:)
integer, allocatable :: IAMatrix(:)
integer, allocatable :: JAMatrix(:)
integer, allocatable :: FreeFaces(:)
!Q-P Matrix
real (REAL_KIND), allocatable :: Q_P(:)
!Body Force Matrix
real (REAL_KIND), allocatable :: BodyForce(:)
real (REAL_KIND), allocatable :: Ax(:)


!-----------------------------------------
! Hydraulic conductivity
!-----------------------------------------
	call Sampling_grid(EP, EP_coordinates,EP_Elemlist)
	!write(*,*) "Sampling fine grids=", EP%nel
unit_x=5
unit_y=5
unit_z=5
	call Sampling_grid_multiple(unit_x, unit_y, unit_z, EP, EPm,EPm_coordinates,EPm_Elemlist)
	!write(*,*) "Sampling coarse grids=", EPm%nel
	!print *, "Sampling coarse grids=", EPm_Elemlist(1,:)
	!print *, "Sampling coarse grids=", EPm_coordinates(:,1)

	!allocate(cell_list(MAX_NLIST))
	!lastID = 0
	!call random_cell_generation(kcell)

	!do i=1,kcell
	!	print *, "Cell location=", cell_list(i)%centre(:,1)
	!end do
	!print *, "number of cells=", kcell
	!print *, "Cell location=", cell_list(10)%centre(:,1)

	!allocate(unique_EPm_x(EPm%nel_x+1))
	call unique(EPm_coordinates(:,1),uniqueEPm_x)
	!call unique(EPm_coordinates(:,2),uniqueEPm_y) !uniqueEPm_x=uniqueEPm_y
	call unique(EPm_coordinates(:,3),uniqueEPm_z)
	!print *, "unique_EPm_x=", uniqueEPm_x(:)
	!write(*,*) "kcell=", nlist
	call grouping(nlist,uniqueEPm_x,uniqueEPm_z,Cell_in,Cell_online,Cell_onedge, &
		Cell_oncorner, size_Cell_in, size_Cell_online,size_Cell_onedge,size_Cell_oncorner)
		!print *, "Cell in =", Cell_in(1,:)
		!print *, "Cell online =", Cell_online(1,:)
		!print *, "size_Cell_online =", size_Cell_online
		!print *, "Cell on corner =", Cell_oncorner(1,:)
		!print *, "size Cell in =", size_Cell_in
		!print *, "size Cell onedge =", size_Cell_onedge
		!print *, "size Cell oncornet =", size_Cell_oncorner

	call cellcount(EPm,EPm_coordinates,EPm_Elemlist,Cell_in,size_Cell_in,cellcount_inside_elem)
	!write(*,*) "cellcount_inside_elem=", cellcount_inside_elem(:)
	call cellcount(EPm,EPm_coordinates,EPm_Elemlist,Cell_online,size_Cell_online,cellcount_online_elem)
	!write(*,*) "cellcount_online_elem =", cellcount_online_elem(:)
	call cellcount(EPm,EPm_coordinates,EPm_Elemlist,Cell_onedge,size_Cell_onedge,cellcount_onedge_elem)
	!write(*,*) "cellcount_onedge_elem =", cellcount_onedge_elem(:)
	call cellcount(EPm,EPm_coordinates,EPm_Elemlist,Cell_oncorner,size_Cell_oncorner,cellcount_oncorner_elem)
	!write(*,*) "cellcount_oncorner_elem =", cellcount_oncorner_elem(:)

	k_empty=10
	count_empty=0
	j=0
	IF (ALLOCATED (vol_ratio)) DEALLOCATE (vol_ratio)
	allocate(vol_ratio(EPm%nel))
	IF (ALLOCATED (k_conduct)) DEALLOCATE (k_conduct)
	allocate(k_conduct(EPm%nel))
	IF (ALLOCATED (k_nonempty)) DEALLOCATE (k_nonempty)
	allocate(k_nonempty(EPm%nel))

cell_radius = Raverage
do i=1,EPm%nel
    vol_ratio(i)= (cellcount_inside_elem(i)*(4./3.)*Pi*cell_radius**3+ & !assumed entirely within the element
        cellcount_online_elem(i)*2./3.*Pi*cell_radius**3+ & !if elements are on a line then they are shared between two elements
        cellcount_onedge_elem(i)/3.*Pi*cell_radius**3+ & !if elements are on an edge then they are shared between 4 elements
        cellcount_oncorner_elem(i)*4./(8.*3.)*Pi*cell_radius**3)&
		/(EPm%volume*EPm%vol_fraction(i)) !corner elements are shared between 8 elements, corner edge by 4

	if (vol_ratio(i)==0.0D+00) then
       k_conduct(i) = k_empty
       count_empty = count_empty+1
    else
       j=j+1;
       k_nonempty(j)=(2*cell_radius)**2/180*(1-vol_ratio(i))**3/(vol_ratio(i)**2) ! in um^2
		!print *, " k_nonempty=", k_nonempty(j)
       if (k_nonempty(j)>k_empty) then
            k_conduct(i) =k_empty
        else
            k_conduct(i) =k_nonempty(j)
        endif
     endif
enddo
	!write(*,*) " vol_ratio=", vol_ratio(:)
	!print *, " k_conduct=", k_conduct(:)

! Convert to nodal values of sampling grid
IF (ALLOCATED (k_node_mat)) DEALLOCATE (k_node_mat)
allocate(k_node_mat(EPm%nnp,2))
k_node_mat(:,:) = 0.0D+00

do e = 1,EPm%nel
    do i = 1,EPm%nen
        k_node_mat(EPm_Elemlist(e,i),2) = k_node_mat(EPm_Elemlist(e,i),2)+1
        if (k_node_mat(EPm_Elemlist(e,i),2)>1) then
            k_node_mat(EPm_Elemlist(e,i),1) = min(k_node_mat(EPm_Elemlist(e,i),1),k_conduct(e))
        !(k_node_mat(EPm_Elemlist(e,i),1)* &
		!						(k_node_mat(EPm_Elemlist(e,i),2)-1) + &
		!				k_conduct(e))/k_node_mat(EPm_Elemlist(e,i),2)
		else
            k_node_mat(EPm_Elemlist(e,i),1) = k_conduct(e)
        endif
    enddo
enddo
!print *, " k_node_mat=", k_node_mat(:,:)
IF (ALLOCATED (k_ref)) DEALLOCATE (k_ref)
allocate(k_ref(size(EP_coordinates,1)))
 call refined_k_grid(EP_coordinates,k_node_mat, EPm,EPm_coordinates,EPm_Elemlist, k_ref)
IF (ALLOCATED (k_elem)) DEALLOCATE (k_elem)
allocate(k_elem(EP%nel))
call find_k_mid(k_ref,EP,EP_Elemlist,k_elem)
!write(*,*) " k_ref=", k_ref(:)
!write(*,*) " k_elem=", k_elem(1:10)
IF (ALLOCATED (k_Cylinder_elem)) DEALLOCATE (k_Cylinder_elem)
allocate(k_Cylinder_elem(NofElements))
!allocate(CCK_list(nel,3))
 call FIND_K_Cylinder(k_elem,EP,EP_Elemlist,EP_coordinates,Centroids, &
				CylinCub_list, k_Cylinder_elem)
!write(*,*) " k_Cylinder_elem=", k_Cylinder_elem(1:20)


!------------------------------------------
! Stiff Matrix
! AMatrix = stiff matrix in the format of compressed spars row format
! with IAMatrix to save the extents of rows indices
! and JAMatrix to save the extents of columns indices
!-----------------------------------------
IF (ALLOCATED (AMatrix)) DEALLOCATE (AMatrix)
IF (ALLOCATED (IAMatrix)) DEALLOCATE (IAMatrix)
IF (ALLOCATED (JAMatrix)) DEALLOCATE (JAMatrix)
	allocate(AMatrix(8*6*NofElements-CommonFaceNo))
	allocate(IAMatrix((8*6*NofElements-CommonFaceNo)))
	allocate(JAMatrix((8*6*NofElements-CommonFaceNo)))
	call sparseMatrices(NofElements,ElToFace,AllSurfs,k_Cylinder_elem, AMatrix,IAMatrix,JAMatrix)
!------------------------------------------
! set body force and other boundary conditions (Dirichlet)
!-------------------------------------------
	allocate(BodyForce(nofaces+nel))
	allocate(Q_P(nofaces+nel))
	allocate(FreeFaces(nofaces+nel-size(NCWall)))
	allocate(Ax(nofaces+nel))
	call Compute_BodyForce(nel,BodyForce,Q_P,FreeFaces)
	call ax_st ( size(BodyForce), size(AMatrix), IAMatrix, JAMatrix, &
	AMatrix, Q_P, Ax ) !computes A*x, A is a sparse matrix
	BodyForce(:) = BodyForce(:) - Ax(:)
	!print *, "Body Force=", BodyForce(:)
	!print *, "QP=", Q_P(:)

count1=0
!print *, "", size(FreeFaces), size(IAMatrix), size(JAMatrix)
!open (unit=out_unit,file="notfound.txt",action="write",status="replace")
!write (out_unit,*) "Free Face=", FreeFaces
do i = 1,size(FreeFaces)
	do j= 1, size(IAMatrix)
		if (FreeFaces(i)==IAMatrix(j)) then
			!print *, "", FreeFaces(i), IAMatrix(j), JAMatrix(j)
			!write (out_unit,*) "i, JA=", i, JAMatrix(j)
			do k = 1, size(FreeFaces)
				if (FreeFaces(k)==JAMatrix(j)) then
				!write (out_unit,*) "k, FreeFaces=", k, FreeFaces(k)
				count1 = count1+1
				endif
			enddo
		endif
	enddo
enddo


!close(out_unit)
!print *, "", count1
!print *, "FreeFaces=", FreeFaces(1:10)

	call solver(AMatrix,IAMatrix,JAMatrix,BodyForce,Q_P,&
	FreeFaces,size(FreeFaces), count1)

	!write(*,*) "final=", Q_P(1:20)

end subroutine
!--------------------------------------------------
!---------------------------------------------------------------
! we solve to get the answers for x(Freeface) =A(Freeface,Freeface)\b(FreeFace)
!---------------------------------------------------------------
subroutine solver(A,IA,JA,BForce,QP,FFace,sizeFFace,n)

integer :: i, j, k,p, FFace(:), sizeFFace, ncc, icc(n), jao(n), row_ind(n)
integer :: count1, count2, n, ccc(n+1), iao(sizeFFace+1), MaxIter, col_ptr(sizeFFace+1)
real(REAL_KIND) :: A(:), NewA(n), acc(n), ao(n), TOLER, acsc(n)
real(REAL_KIND) :: solution(sizeFFace), residual(sizeFFace)
integer :: IA(:),JA(:), NewIA(n), NewJA(n)
integer :: rowcount, columncount
integer :: isolve_mode
real(4) :: t1, t2

real(REAL_KIND) :: BForce(:), QP(:), newBForce(sizeFFace)
real(REAL_KIND) :: newQP(sizeFFace), Ax(sizeFFace)

!integer, parameter :: out_unit=20



!open (unit=out_unit,file="rearranged A.txt",action="write",status="replace")


count1 = 0
rowcount=0
do i = 1, sizeFFace
	rowcount=rowcount+1
	columncount=0
	do j= 1, size(IA)
		if (FFace(i)==IA(j)) then
			do k = 1, sizeFFace
				if (JA(j)==FFace(k)) then
					count1 = count1+1
					NewA(count1) = A(j)
					NewIA(count1) = rowcount !IA(j)
					NewJA(count1) = k !JA(j)
					!write (out_unit,*) " New A=",NewA(count1),NewIA(count1),NewJA(count1)
				endif
			enddo
		endif
	enddo
	newBForce(i) = BForce(FFace(i))
	newQP(i) = QP(FFace(i))
enddo


   call coocsr (sizeFFace,size(NEWA(:)),NewA,NEWIA,NEWJA,ao,jao,iao )
   !write (out_unit,*) " New A=",ao
   !print *, "column indices", jao
   !print *, "row indices", iao
   !print *, "real values", ao
   call rearrange_cr (sizeFFace, iao, jao, ao )
	!print *, "column indices", jao
	!print *, "row indices", iao
	!print *, "real values", ao
	!write (out_unit,*) " real values=",ao
	!close (out_unit)

    call csrcsc (sizeFFace, 1, 1, ao, jao, iao, acsc, row_ind, col_ptr )

!----------------------------------------------------------------
! convert from 1-based to 0-based
!----------------------------------------------------------------
do j = 1, sizeFFace+1
	col_ptr(j) = col_ptr(j) - 1
enddo
do p = 1, size(NEWA(:))
	row_ind(p) = row_ind(p) - 1
enddo

isolve_mode = SOLVE_NO_ITERATION
call cpu_time(t1)
call solve_umf(isolve_mode, sizeFFace, col_ptr, row_ind, acsc, newBForce, solution)
call cpu_time(t2)
write(*,'(a,f8.1)') 'Solve time (sec): ',t2-t1

! print the residual.
call resid (sizeFFace, col_ptr, row_ind, acsc,  solution, newBForce, residual)

   newQP = solution
   !******check error of x compare to b********
   !call ax_st ( sizeFFace, size(NEWA), NEWIA, NEWJA, NEWA, newQP, Ax )
	!error(:) = Ax(:) - newBForce(:)
	!print *, "error of final result=", error
	!*******************************************

	do i = 1, sizeFFace
		QP(FFace(i))= newQP(i)
	enddo
	!print *, "final=", QP


!close (out_unit)
end subroutine solver
!-------------------------------------------------------
!*****************************************************************************
 subroutine rearrange_cr ( n, ia, ja, a )
!!! REARRANGE_CR sorts a sparse compressed row matrix.
    !    This routine guarantees that the entries in the CR matrix
    !    are properly sorted.
    !
    !    After the sorting, the entries of the matrix are rearranged in such
    !    a way that the entries of each column are listed in ascending order
    !    of their column values.
    !    Input, integer :: N, the order of the system.
    !
    !    Input, integer :: NZ_NUM, the number of nonzeros.
    !
    !    Input, integer :: IA(N+1), the compressed row indices.
    !
    !    Input/output, integer :: JA(NZ_NUM), the column indices.
    !    On output, these may have been rearranged by the sorting.
    !
    !    Input/output, real (REAL_KIND) :: A(NZ_NUM), the matrix values.  On output,
    !    the matrix values may have been moved somewhat because of the sorting.
    !
    implicit none

    integer :: n
    integer :: ia(*) !ia(n+1)
    integer :: ja(*) !ja(nz_num)
    real (REAL_KIND) :: a(*) !a(nz_num)

    integer :: i
    integer :: i4temp
    integer :: k
    integer :: l
    real (REAL_KIND) :: r8temp

    do i = 1, n

       do k = ia(i), ia(i+1) - 2
          do l = k + 1, ia(i+1) - 1

             if ( ja(l) < ja(k) ) then
                i4temp = ja(l)
                ja(l)  = ja(k)
                ja(k)  = i4temp

                r8temp = a(l)
                a(l)   = a(k)
                a(k)   = r8temp
             end if

          end do
       end do

    end do

    return
 end subroutine rearrange_cr
 !-----------------------------------------------
end module
