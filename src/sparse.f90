module sparse

use global
use Mesh_Generate
use sparsekit
use m_unista
!use umf_solver
use FEmapper

implicit none

INTERFACE
    ! The BIND(C) tells the compiler that this is an "interoperable"
    ! procedure.  The compiler adjusts the naming conventions as
    ! appropriate for the companion C processor.
	subroutine slu_solve(nprocs, n, nnz, a, asub, xa, rhs) bind(C)
    USE,INTRINSIC :: ISO_C_BINDING  ! Declares C kinds
	integer(c_int) :: nprocs, n, nnz, asub(*), xa(*)
	real(c_double) :: a(*), rhs(*)
	end subroutine slu_solve
END INTERFACE

contains

!------------------------------------------
!------------------------------------------
subroutine READ_Bmatrix(nel)

integer :: nel, j
character * ( 255 ) :: fileplace, Bmatrix_filename

fileplace= "./B_matrix/"

do j = 1,nel

   	write (Bmatrix_filename,'("B_matrix_",I0,".txt")') j
	call Read_real(fileplace,Bmatrix_filename,6,6,B_MATRIX(j)%Inside_B)
enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
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
integer :: ii, jj, j, I(6), signum(6), DiagSig(6,6)
integer :: counter, Recount1, Ccounter, commonCount
logical :: ok1
character * ( 255 ) :: fileplace, Bmatrix_filename

!open (unit=out_unit,file="missed.txt",action="write",status="replace")
!allocate(Ak(6,6))
!Ak(1,:)= 1/(2.*6.) *(/2.,-1.,0.,0.,0.,0./)
!Ak(2,:)= 1/(2.*6.) *(/-1.,2.,0.,0.,0.,0./)
!Ak(3,:)= 1/(2.*6.) *(/0.,0.,2.,-1.,0.,0./)
!Ak(4,:)= 1/(2.*6.) *(/0.,0.,-1.,2.,0.,0./)
!Ak(5,:)= 1/(2.*6.) *(/0.,0.,0.,0.,2.,-1./)
!Ak(6,:)= 1/(2.*6.) *(/0.,0.,0.,0.,-1.,2./)

!fileplace= "./B_matrix/"

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

findequal = 0
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
!-----------------------------------------------
function ones(Msize)

integer :: Msize, i, ones(Msize)
	do i= 1, Msize
		ones(i) = 1
	enddo
!print *, "ones=", ones(:)
end function

!-----------------------------------------------
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

end subroutine

!-----------------------------------------------
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
!-------------------------------------------------
subroutine FEsolve
integer :: unit_x,unit_y,unit_z
integer :: i, j, k, e,  count_empty, count1, nnz
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
real (REAL_KIND), allocatable :: Ax(:)

integer, parameter :: out_unit=20

!open (unit=out_unit,file="Amatrix.txt",action="write",status="replace")
!open (unit=out_unit,file="kconduct.txt",action="write",status="replace")
!open (unit=out_unit,file="Cellcount.txt",action="write",status="replace")
!-----------------------------------------
! Hydraulic conductivity
!-----------------------------------------
	call Sampling_grid(EP, EP_coordinates,EP_Elemlist)
	!call logger('Done Sampling Grid')
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
	!print *, "unique_EPm_z=", uniqueEPm_z(:)
	!print *, "unique_EPm_x=", uniqueEPm_x(:)
	!write(*,*) "kcell=", nlist
	!write(*,*) "kcell=", ncells
	call grouping(ncells,uniqueEPm_x,uniqueEPm_z,Cell_in,Cell_online,Cell_onedge, &
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
	!print *, "cellcount_inside_elem =", size(cellcount_inside_elem)
	call cellcount(EPm,EPm_coordinates,EPm_Elemlist,Cell_online,size_Cell_online,cellcount_online_elem)
	!write(*,*) "cellcount_online_elem =", cellcount_online_elem(:)
	call cellcount(EPm,EPm_coordinates,EPm_Elemlist,Cell_onedge,size_Cell_onedge,cellcount_onedge_elem)
	!write(*,*) "cellcount_onedge_elem =", cellcount_onedge_elem(:)
	call cellcount(EPm,EPm_coordinates,EPm_Elemlist,Cell_oncorner,size_Cell_oncorner,cellcount_oncorner_elem)
	!write(*,*) "cellcount_oncorner_elem =", cellcount_oncorner_elem(:)

!write(*,*) " Number of EPm elements=", EPm%nel
!do i=1,size(cellcount_inside_elem)
 !   write (*,*) cellcount_inside_elem(i)
!enddo


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
	if (no_cells) vol_ratio(i) = 0

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
!write(*,*) " Number of EPm elements=", EPm%nel
!do i=1,EPm%nel
 !       write (out_unit,*) k_conduct(i)
!enddo
	!write(*,*) " vol_ratio=", vol_ratio(:)
	!print *, " k_conduct=", k_conduct(:)
!call logger('Done finding k_conduct')
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
        !if (k_node_mat(EPm_Elemlist(e,i),2)>1 .AND. EPm_coordinates(EPm_Elemlist(e,i),3)<600 ) then
        !    k_node_mat(EPm_Elemlist(e,i),1) = min(k_node_mat(EPm_Elemlist(e,i),1),k_conduct(e))
        !elseif (k_node_mat(EPm_Elemlist(e,i),2)>1 .AND. EPm_coordinates(EPm_Elemlist(e,i),3)>=600 &
        !    .AND. EPm_coordinates(EPm_Elemlist(e,i),3)<700) then
        !    k_node_mat(EPm_Elemlist(e,i),1) =( (k_node_mat(EPm_Elemlist(e,i),1)* (k_node_mat(EPm_Elemlist(e,i),2)-1) + &
		!				k_conduct(e))/k_node_mat(EPm_Elemlist(e,i),2)+ &
		!				k_node_mat(EPm_Elemlist(e-16,i),1) )/2

        !!(k_node_mat(EPm_Elemlist(e,i),1)* &
		!!						(k_node_mat(EPm_Elemlist(e,i),2)-1) + &
		!!				k_conduct(e))/k_node_mat(EPm_Elemlist(e,i),2)
		!else
         !   k_node_mat(EPm_Elemlist(e,i),1) = k_conduct(e)
        !endif
    enddo
enddo
!write(*,*) " Number of EPm elements=", EPm%nnp
!do i=1,EPm%nnp
!        write (out_unit,*) k_node_mat(i,:)
!enddo
!call logger('Done finding k_node_mat')
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
!write(*,*) " Number of elements=", NofElements
!do i=1,NofElements
!        !write (out_unit,*) i
!        write (out_unit,*) " k_Cylinder=",k_Cylinder_elem(i)
!enddo
!call logger('Done FIND_K_Cylinder')

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
    !do i=1,size(AMatrix)
    !    write (out_unit,*) AMatrix(i)
    !enddo
!------------------------------------------
! set body force and other boundary conditions (Dirichlet)
!-------------------------------------------
	IF (ALLOCATED (Ax)) DEALLOCATE (Ax)
	allocate(Ax(nofaces+NofElements))
	IF (ALLOCATED (BodyForce)) DEALLOCATE (BodyForce)
	allocate(BodyForce(nofaces+NofElements))
	IF (ALLOCATED (Q_P)) DEALLOCATE (Q_P)
    allocate(Q_P(nofaces+NofElements))
    IF (ALLOCATED (FreeFaces)) DEALLOCATE (FreeFaces)
    allocate(FreeFaces(nofaces+NofElements-size(NCWall)))

	call Compute_BodyForce(NofElements,BodyForce,Q_P,FreeFaces)
	call ax_st ( size(BodyForce), size(AMatrix), IAMatrix, JAMatrix, &
	AMatrix, Q_P, Ax ) !computes A*x, A is a sparse matrix
	BodyForce(:) = BodyForce(:) - Ax(:)
	!print *, "Body Force = ", BodyForce(:)
	!print *, "QP=", Q_P(:)

count1=0
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
nnz = size(AMatrix)
write(*,*) 'size(FreeFaces), count1, nnz: ',size(FreeFaces), count1, nnz
do k = 1,nnz
	if (AMatrix(k) == 0) then
		write(*,*) 'Zero value in AMatrix at: ',k
		stop
	endif
enddo
call solver(AMatrix,IAMatrix,JAMatrix,BodyForce,Q_P,&
FreeFaces,size(FreeFaces), count1)
 !do i=1,nofaces+NofElements
!	write(nflog,'(a,3f10.1)') "final flow=", Q_P(i)
!enddo

!close (out_unit)
end subroutine

!---------------------------------------------------------------
! we solve to get the answers for x(Freeface) =A(Freeface,Freeface)\b(FreeFace)
! n = nz = number of non-zero entries
!---------------------------------------------------------------
subroutine solver(A,IA,JA,BForce,QP,FFace,sizeFFace,n)

integer :: i, j, k,p, FFace(:), sizeFFace, ncc, icc(n), jao(n), row_ind(n)
integer :: count1, count2, n, ccc(n+1), iao(sizeFFace+1), MaxIter, col_ptr(sizeFFace+1)
real(REAL_KIND) :: A(:), NewA(n), acc(n), ao(n), TOLER, acsc(n)
real(REAL_KIND) :: solution(sizeFFace), residual(sizeFFace)
real(REAL_KIND) :: out_vel1, DV_velocity, ave_outflow_veloc
integer :: IA(:),JA(:), NewIA(n), NewJA(n)
integer :: rowcount, columncount!, DV_Number(size(DCOutlet))
integer :: isolve_mode, nnz, M
real(4) :: t1, t2

real(REAL_KIND) :: BForce(:), QP(:), newBForce(sizeFFace)
real(REAL_KIND) :: newQP(sizeFFace), Ax(sizeFFace)

integer, allocatable :: Ap_save(:), Ai_save(:)
real(8), allocatable :: Ax_save(:)
logical :: use_umf_solver = .false.

allocate(Ap_save(sizeFFace+1))
allocate(Ai_save(n))
allocate(Ax_save(n))

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
					NewIA(count1) = rowcount 	!IA(j)
					NewJA(count1) = k 			!JA(j)
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
	!do i=1,size(NEWA(:))
    !    write (out_unit,*) ao(i)
    !enddo
	!close (out_unit)

    call csrcsc (sizeFFace, 1, 1, ao, jao, iao, acsc, row_ind, col_ptr )
!call logger('Done csrcsc')
!write (*,*) "nz",size(NEWA(:))
!do i=1,size(NEWA(:))
 !   write (out_unit,*) acsc(i)
!enddo

! At this point the non-zero entries are stored in compressed sparse column format (1-based arrays):
!	acsc(nnz)
!	col_ptr(M+1)
!	row_ind(nnz)
! where
!	M = order of stiffness matrix = size(NEWA)
!	nnz = count of non-zero entries
! This is also called Harwell/Boeing format (HB)

M = sizeFFace
nnz = count1
!call writeHB(M,nnz,col_ptr,row_ind,acsc,newBForce)

if (use_umf_solver) then

#if 0
	!========== start umf solver=============================
	isolve_mode = SOLVE_NO_ITERATION
	call cpu_time(t1)
	call solve_umf(isolve_mode, M, col_ptr, row_ind, acsc, newBForce, solution)
!	write(nflog,*) 'solver_harwell_umf: M, nnz: ',M, nnz, size(NewA)
!	call solve_harwell_umf(M, nnz, acsc, row_ind, col_ptr, newBForce, solution)
	call cpu_time(t2)
	write(*,'(a,f8.1)') 'Solve time (sec): ',t2-t1

	! print the residual.
	call resid1 (M, col_ptr, row_ind, acsc,  solution, newBForce, residual)
	!========== end umf solver=============================
#endif

else
	!----------------------------------------------------------------
	! convert from 1-based to 0-based
	!----------------------------------------------------------------
	do j = 1, sizeFFace+1
		col_ptr(j) = col_ptr(j) - 1
	enddo
	do p = 1, size(NEWA(:))
		row_ind(p) = row_ind(p) - 1
	enddo

	!========== start slu solver=============================
	solution= newBForce
	Ap_save = col_ptr
	Ai_save = row_ind
	Ax_save = acsc
	call cpu_time(t1)
	!n = sizeFFace
	!nnz = size(NewA)
	!write(nflog,*) 'n,nnz: ',n,nnz
	!write(nflog,*) 'col_ptr'
	!write(nflog,'(11i7)') col_ptr
	!write(nflog,*) 'row_ind'
	!write(nflog,'(11i7)') row_ind
	!write(nflog,*) 'acsc (nzvals)'
	!write(nflog,'(5e15.8)') NewA
	!write(nflog,*) 'rhs'
	!write(nflog,'(5e14.5)') solution

	write(nflog,*) 'slu_solver: M, nnz: ',M, nnz, size(NewA)
	call slu_solve(4, sizeFFace, size(NewA(:)), acsc, row_ind, col_ptr, solution);
	!!UMF call solve(isolve_mode, n, Ap, Ai, Ax, b, x)
	call cpu_time(t2)
	write(logmsg,'(a,f8.1)') 'Solve time (cpu_time/nprocs) (sec): ',(t2-t1)/4
	call logger(logmsg)
	!write(nflog,*) 'solution'
	!write(nflog,'(5e14.5)') solution

	!! print the residual.
	!write(*,*) 'After solve, residual:'
	call resid0 (sizeFFace, Ap_save, Ai_save, Ax_save, solution, newBForce, residual)
	!========== end slu solver=============================
endif

write(nflog,*) 'Solution:'
do i = 1,M
	write(nflog,'(i8,e12.3)') i,solution(i)
enddo

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
	out_vel1=0
	write(nflog2,'(i6,$)') istep
    do i = 1,size(DCInlet)
        out_vel1 = QP(DCInlet(i)) + out_vel1
        !DV_Number(i) = DCOutlet(i)  !first column stores face number
        DV_velocity = QP(DCInlet(i)) !2nd column stores velocity on the face
        write(nflog2,'(i20,$)') DCInlet(i)
        write(nflog2,'(1f20.3,$)') DV_velocity
    enddo
    write(nflog2,*)
    ave_outflow_veloc = out_vel1/i
    write(nflog2,'(a,1f15.3)') 'Average out flow: ', ave_outflow_veloc

!close (out_unit)
end subroutine solver

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
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
!-----------------------------------------------
subroutine GetVel
real (REAL_KIND) :: elem_vertices(8,3), fsum(3)
integer, allocatable :: cell_elem(:)
real (REAL_KIND), allocatable :: cell_in_reference(:,:)
real (REAL_KIND) :: coordPhy(3), coordref(3)
integer :: elem_nodes(8), i, j, k, e, II(6),signum(6), ND, Nvert
real (REAL_KIND) :: xhat, yhat, zhat, Uhat(6)
real (REAL_KIND) :: uhat_x, uhat_y, uhat_z, x_corners(8,3)
real (REAL_KIND) :: D_matrix(3,3),JJ,u_total(3)!, velocity_vec(3,1)
!real (REAL_KIND), allocatable :: u_x(:), u_y(:), u_z(:)

! transfer cell position from physical to the reference elements:
	!print *, "nlist=",nlist
	IF (ALLOCATED (cell_elem)) DEALLOCATE (cell_elem)
	allocate(cell_elem(ncells))
	call find_cell_in_CylElem(e_Z,NofElements_z,nlist,Cylinder_ELEMENT,Cylinder_vertices,layered_el,cell_elem)
	!do i=1,100!nlist
!		print *, "cell in cyl elem=", cell_elem(i)
		!write(*,'(a,i6)') "cell in cyl elem=", cell_elem(i)
	!enddo
	!call logger('Done finding cell_elem')
	ND = 3
	Nvert = 8
	IF (ALLOCATED (cell_in_reference)) DEALLOCATE (cell_in_reference)
	allocate(cell_in_reference(ncells,3))
	do i=1,ncells
        !write(nflog,'(i6,$)') cell_elem(i)
		!print *, "element=",Element(cell_elem(i),:)
		if (cell_list(i)%state == GONE_BACK .or. cell_list(i)%state == GONE_THROUGH) cycle
		elem_nodes(:) = Cylinder_ELEMENT(cell_elem(i),:)
		do j=1,8
			elem_vertices(j,1:3)=Cylinder_vertices(elem_nodes(j),1:3)
		enddo

		coordPhy(1:3) = cell_list(i)%centre(1:3,1)
		call mapper(ND,Nvert,elem_vertices,coordPhy,coordref)
		cell_in_reference(i,:) = coordref(1:3)
	enddo

	!write(nflog,*)
	!call logger('Done transfer to reference element')
!--------------------------------------------------
!it is to test if the cell position found in the reference element is correct or not
	!do i=1,10
    !    write(*,'(a,3e12.3)') "cell in original elem=",cell_list(1)%centre(1:3,1)
	!	write(*,'(a,3e12.3)') "cell in reference elem=", cell_in_reference(1,:)
	!enddo
!do j=1,ncells
!	fsum=0
!	if (cell_list(j)%state == GONE_BACK .or. cell_list(j)%state == GONE_THROUGH) cycle
 !       elem_nodes(1)=Cylinder_ELEMENT(cell_elem(j),1)
!		elem_nodes(2)=Cylinder_ELEMENT(cell_elem(j),2)
!		elem_nodes(3)=Cylinder_ELEMENT(cell_elem(j),3)
!		elem_nodes(4)=Cylinder_ELEMENT(cell_elem(j),4)
!		elem_nodes(5)=Cylinder_ELEMENT(cell_elem(j),5)
!		elem_nodes(6)=Cylinder_ELEMENT(cell_elem(j),6)
!		elem_nodes(7)=Cylinder_ELEMENT(cell_elem(j),7)
!		elem_nodes(8)=Cylinder_ELEMENT(cell_elem(j),8)
!	do i = 1,Nvert
!        fsum(1:3) = fsum(1:3) + Cylinder_vertices(elem_nodes(i),1:3)*fshape(Nvert,i,cell_in_reference(j,:))
!    enddo
!    write(*,'(a,3e12.3)') 'Error: ',fsum(1:ND) - cell_list(j)%centre(1:3,1)
! enddo
!----------------------------------------------------------------------------
!test versus Matlab

!allocate(u_x(NofElements))
!allocate(u_y(NofElements))
!allocate(u_z(NofElements))

!do e=1, NofElements
!    xhat=0
!    yhat=0
!    zhat=0
!    II(:) = ElToFace(e,:)
!    signum=(/1,1,1,1,1,1/)
!    do k=1,6
!		if (e==AllSurfs(II(k),6)) then
!			signum(k)=-1 ! -1 if edge is on T-
!		endif
!	enddo

	!write(nflog,'(a,6i6)' ) 'signum ',signum(:)
!	do i=1,6
!        Uhat(i) = Q_P(ElToFace(e,i))*signum(i)
!    enddo
!    uhat_x=(Uhat(1)-Uhat(2))*1/8+(Uhat(1)+Uhat(2))*xhat/8
!    uhat_y=(Uhat(6)-Uhat(5))*1/8+ (Uhat(5)+Uhat(6))*zhat/8
!    uhat_z=(Uhat(4)-Uhat(3))*1/8+ (Uhat(3)+Uhat(4))*yhat/8
    !write(nflog,'(a,i6)' ) 'x_corner',e
!    do j=1,8
!        x_corners(j,:)=(/Cylinder_vertices((Cylinder_ELEMENT(e,j)),1),&
!        Cylinder_vertices((Cylinder_ELEMENT(e,j)),2) , Cylinder_vertices((Cylinder_ELEMENT(e,j)),3)/)
        !write(nflog,'(3f10.2)' ) x_corners(j,:)
!    enddo
!    D_matrix=0.0D0
!    D_matrix(1,1)=(-1.)*(1.-yhat)*(1.-zhat)/(8.)*x_corners(1,1)+ &
!			(1.-yhat)*(1.-zhat)/(8.)*x_corners(2,1)+ &
!			(1.+yhat)*(1.-zhat)/(8.)*x_corners(3,1)+ &
!			(-1.)*(1.+yhat)*(1.-zhat)/(8.)*x_corners(4,1)+ &
!			(-1.)*(1.-yhat)*(1.+zhat)/(8.)*x_corners(5,1)+ &
!			(1.-yhat)*(1.+zhat)/(8.)*x_corners(6,1)+ &
!			(1.+yhat)*(1.+zhat)/(8.)*x_corners(7,1)+ &
!			(-1.)*(1.+yhat)*(1.+zhat)/(8.)*x_corners(8,1)
!	D_matrix(1,2)=(1.-xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(1,1)+ &
!            (1.+xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(2,1)+ &
!			(1.+xhat)*(1.-zhat)/(8.)*x_corners(3,1)+ &
!			(1.-xhat)*(1.-zhat)/(8.)*x_corners(4,1)+ &
!			(1.-xhat)*(-1.)*(1.+zhat)/(8.)*x_corners(5,1)+ &
!			(1.+xhat)*(-1.)*(1.+zhat)/(8.)*x_corners(6,1)+ &
!			(1.+xhat)*(1.+zhat)/(8.)*x_corners(7,1)+ &
!			(1.-xhat)*(1.+zhat)/(8.)*x_corners(8,1)
!	D_matrix(1,3)=(1-xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(1,1)+ &
!			(1.+xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(2,1)+ &
!			(1.+xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(3,1)+ &
!			(1.-xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(4,1)+ &
!			(1.-xhat)*(1.-yhat)*(1.)/(8.)*x_corners(5,1)+ &
!			(1.+xhat)*(1.-yhat)*(1.)/(8.)*x_corners(6,1)+ &
!			(1.+xhat)*(1.+yhat)*(1.)/(8.)*x_corners(7,1)+ &
!			(1.-xhat)*(1.+yhat)*(1.)/(8.)*x_corners(8,1)
!	D_matrix(2,1)=(-1.)*(1.-yhat)*(1.-zhat)/(8.)*x_corners(1,2)+&
!			(1.-yhat)*(1.-zhat)/(8.)*x_corners(2,2)+ &
!			(1.+yhat)*(1.-zhat)/(8.)*x_corners(3,2)+ &
!			(-1.)*(1.+yhat)*(1.-zhat)/(8.)*x_corners(4,2)+ &
!			(-1.)*(1.-yhat)*(1.+zhat)/(8.)*x_corners(5,2)+ &
!			(1.-yhat)*(1.+zhat)/(8.)*x_corners(6,2)+ &
!			(1.+yhat)*(1.+zhat)/(8.)*x_corners(7,2)+ &
!			(-1.)*(1.+yhat)*(1.+zhat)/(8.)*x_corners(8,2)
!	D_matrix(2,2)=(1.-xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(1,2)+ &
!			(1.+xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(2,2)+ &
!			(1.+xhat)*(1.-zhat)/(8.)*x_corners(3,2)+ &
!			(1.-xhat)*(1.-zhat)/(8.)*x_corners(4,2)+ &
!			(1.-xhat)*(-1.)*(1+zhat)/(8.)*x_corners(5,2)+ &
!			(1.+xhat)*(-1.)*(1+zhat)/(8.)*x_corners(6,2)+ &
!			(1.+xhat)*(1.+zhat)/(8.)*x_corners(7,2)+ &
!			(1.-xhat)*(1.+zhat)/(8.)*x_corners(8,2)
!	D_matrix(2,3)=(1.-xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(1,2)+ &
!			(1.+xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(2,2)+ &
!			(1.+xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(3,2)+ &
!			(1.-xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(4,2)+ &
!			(1.-xhat)*(1.-yhat)*(1.)/(8.)*x_corners(5,2)+ &
!			(1.+xhat)*(1.-yhat)*(1.)/(8.)*x_corners(6,2)+ &
!			(1.+xhat)*(1.+yhat)*(1.)/(8.)*x_corners(7,2)+ &
!			(1.-xhat)*(1.+yhat)*(1.)/(8.)*x_corners(8,2)
!	D_matrix(3,1)=(-1.)*(1.-yhat)*(1.-zhat)/(8.)*x_corners(1,3)+ &
!			(1.-yhat)*(1.-zhat)/(8.)*x_corners(2,3)+ &
!			(1.+yhat)*(1.-zhat)/(8.)*x_corners(3,3)+ &
!			(-1.)*(1.+yhat)*(1.-zhat)/(8.)*x_corners(4,3)+ &
!			(-1.)*(1.-yhat)*(1.+zhat)/(8.)*x_corners(5,3)+ &
!			(1.-yhat)*(1.+zhat)/(8.)*x_corners(6,3)+ &
!			(1.+yhat)*(1.+zhat)/(8.)*x_corners(7,3)+ &
!			(-1.)*(1.+yhat)*(1.+zhat)/(8.)*x_corners(8,3)
!	D_matrix(3,2)=(1.-xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(1,3)+ &
!			(1.+xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(2,3)+ &
!			(1.+xhat)*(1.-zhat)/(8.)*x_corners(3,3)+ &
!			(1.-xhat)*(1.-zhat)/(8.)*x_corners(4,3)+ &
!			(1.-xhat)*(-1.)*(1.+zhat)/(8.)*x_corners(5,3)+ &
!			(1.+xhat)*(-1.)*(1.+zhat)/(8.)*x_corners(6,3)+ &
!			(1.+xhat)*(1.+zhat)/(8.)*x_corners(7,3)+ &
!			(1.-xhat)*(1.+zhat)/(8.)*x_corners(8,3)
!	D_matrix(3,3)=(1.-xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(1,3)+ &
!			(1.+xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(2,3)+ &
!			(1.+xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(3,3)+ &
!			(1.-xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(4,3)+ &
!			(1.-xhat)*(1.-yhat)*(1.)/(8.)*x_corners(5,3)+ &
!			(1.+xhat)*(1.-yhat)*(1.)/(8.)*x_corners(6,3)+ &
!			(1.+xhat)*(1.+yhat)*(1.)/(8.)*x_corners(7,3)+ &
!			(1.-xhat)*(1.+yhat)*(1.)/(8.)*x_corners(8,3)

    !do i=1,3
        !write(nflog,'(3f10.3)' ) D_matrix(2,1),D_matrix(2,2),D_matrix(2,3)
    !enddo


	!velocity_vec=RESHAPE([uhat_x, uhat_y, uhat_z],(/3,1/))
	!write(nflog,'(a,f10.2)' ) 'determinant', JJ
!	write(nflog,'(a,3f10.1)' ) 'velocity ', uhat_x, uhat_y, uhat_z
	!write(nflog,'(a,3f15.2)' ) 'matmul', matmul(D_matrix(1,:), velocity_vec), matmul(D_matrix(2,:), velocity_vec) &
	!, matmul(D_matrix(3,:), velocity_vec)
!	write(nflog,'(3f10.3)' ) D_matrix(2,1),D_matrix(2,2),D_matrix(2,3)
!	write(nflog,'(a,3f15.2)' ) 'mine ', D_matrix(2,1),uhat_x,D_matrix(2,1)*uhat_x
!	write(nflog,'(a,3f15.2)' ) 'mine ', D_matrix(2,2),uhat_y,D_matrix(2,2)*uhat_y
!	write(nflog,'(a,3f15.2)' ) 'mine ', D_matrix(2,3),uhat_z,D_matrix(2,3)*uhat_z

!	u_total= matmul(D_matrix, [uhat_x, uhat_y, uhat_z])
!	JJ=FindDet(D_matrix,3)
!	u_total=u_total*1/JJ

!	u_x(e)=u_total(1)
!	u_y(e)=u_total(2)
!	u_z(e)=u_total(3)
    !write(*,*) "u_cell_z=", u_z(e)
     !write(nflog,'(a,3f10.1)' ) 'element velocity ', u_x(e),u_y(e),u_z(e)
     !write(nflog,'(a,6i6)' ) 'element to face ',ElToFace(e,:)

!enddo
!write(nflog,*)


IF (ALLOCATED (u_cell_x)) DEALLOCATE (u_cell_x)
IF (ALLOCATED (u_cell_y)) DEALLOCATE (u_cell_y)
IF (ALLOCATED (u_cell_z)) DEALLOCATE (u_cell_z)
allocate(u_cell_x(ncells))
allocate(u_cell_y(ncells))
allocate(u_cell_z(ncells))

do i=1, ncells
    if (cell_list(i)%state == GONE_BACK .or. cell_list(i)%state == GONE_THROUGH) cycle
    xhat=cell_in_reference(i,1)
    yhat=cell_in_reference(i,2)
    zhat=cell_in_reference(i,3)
    e= cell_elem(i)
    II(:) = ElToFace(e,:)
    signum=(/1,1,1,1,1,1/)
    do k=1,6
		if (e==AllSurfs(II(k),6)) then
			signum(k)=-1 ! -1 if edge is on T-
		endif
	enddo
    do j=1,6
        Uhat(j) = Q_P(ElToFace(e,j))*signum(j)
    enddo
    uhat_x=(Uhat(1)-Uhat(2))*1/8+(Uhat(1)+Uhat(2))*xhat/8
    uhat_y=(Uhat(6)-Uhat(5))*1/8+ (Uhat(5)+Uhat(6))*zhat/8
    uhat_z=(Uhat(4)-Uhat(3))*1/8+ (Uhat(3)+Uhat(4))*yhat/8
    do j=1,8
!        x_corners(j,:)=(/Cylinder_vertices((Cylinder_ELEMENT(e,j)),1),&
!        Cylinder_vertices((Cylinder_ELEMENT(e,j)),2) , Cylinder_vertices((Cylinder_ELEMENT(e,j)),3)/)
		x_corners(j,:) = Cylinder_vertices(Cylinder_ELEMENT(e,j),:)
    enddo
    D_matrix(1,1)=(-1.)*(1.-yhat)*(1.-zhat)/(8.)*x_corners(1,1)+ &
			(1.-yhat)*(1.-zhat)/(8.)*x_corners(2,1)+ &
			(1.+yhat)*(1.-zhat)/(8.)*x_corners(3,1)+ &
			(-1.)*(1.+yhat)*(1.-zhat)/(8.)*x_corners(4,1)+ &
			(-1.)*(1.-yhat)*(1.+zhat)/(8.)*x_corners(5,1)+ &
			(1.-yhat)*(1.+zhat)/(8.)*x_corners(6,1)+ &
			(1.+yhat)*(1.+zhat)/(8.)*x_corners(7,1)+ &
			(-1.)*(1.+yhat)*(1.+zhat)/(8.)*x_corners(8,1)
	D_matrix(1,2)=(1.-xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(1,1)+ &
            (1.+xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(2,1)+ &
			(1.+xhat)*(1.-zhat)/(8.)*x_corners(3,1)+ &
			(1.-xhat)*(1.-zhat)/(8.)*x_corners(4,1)+ &
			(1.-xhat)*(-1.)*(1.+zhat)/(8.)*x_corners(5,1)+ &
			(1.+xhat)*(-1.)*(1.+zhat)/(8.)*x_corners(6,1)+ &
			(1.+xhat)*(1.+zhat)/(8.)*x_corners(7,1)+ &
			(1.-xhat)*(1.+zhat)/(8.)*x_corners(8,1)
	D_matrix(1,3)=(1-xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(1,1)+ &
			(1.+xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(2,1)+ &
			(1.+xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(3,1)+ &
			(1.-xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(4,1)+ &
			(1.-xhat)*(1.-yhat)*(1.)/(8.)*x_corners(5,1)+ &
			(1.+xhat)*(1.-yhat)*(1.)/(8.)*x_corners(6,1)+ &
			(1.+xhat)*(1.+yhat)*(1.)/(8.)*x_corners(7,1)+ &
			(1.-xhat)*(1.+yhat)*(1.)/(8.)*x_corners(8,1)
	D_matrix(2,1)=(-1.)*(1.-yhat)*(1.-zhat)/(8.)*x_corners(1,2)+&
			(1.-yhat)*(1.-zhat)/(8.)*x_corners(2,2)+ &
			(1.+yhat)*(1.-zhat)/(8.)*x_corners(3,2)+ &
			(-1.)*(1.+yhat)*(1.-zhat)/(8.)*x_corners(4,2)+ &
			(-1.)*(1.-yhat)*(1.+zhat)/(8.)*x_corners(5,2)+ &
			(1.-yhat)*(1.+zhat)/(8.)*x_corners(6,2)+ &
			(1.+yhat)*(1.+zhat)/(8.)*x_corners(7,2)+ &
			(-1.)*(1.+yhat)*(1.+zhat)/(8.)*x_corners(8,2)
	D_matrix(2,2)=(1.-xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(1,2)+ &
			(1.+xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(2,2)+ &
			(1.+xhat)*(1.-zhat)/(8.)*x_corners(3,2)+ &
			(1.-xhat)*(1.-zhat)/(8.)*x_corners(4,2)+ &
			(1.-xhat)*(-1.)*(1+zhat)/(8.)*x_corners(5,2)+ &
			(1.+xhat)*(-1.)*(1+zhat)/(8.)*x_corners(6,2)+ &
			(1.+xhat)*(1.+zhat)/(8.)*x_corners(7,2)+ &
			(1.-xhat)*(1.+zhat)/(8.)*x_corners(8,2)
	D_matrix(2,3)=(1.-xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(1,2)+ &
			(1.+xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(2,2)+ &
			(1.+xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(3,2)+ &
			(1.-xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(4,2)+ &
			(1.-xhat)*(1.-yhat)*(1.)/(8.)*x_corners(5,2)+ &
			(1.+xhat)*(1.-yhat)*(1.)/(8.)*x_corners(6,2)+ &
			(1.+xhat)*(1.+yhat)*(1.)/(8.)*x_corners(7,2)+ &
			(1.-xhat)*(1.+yhat)*(1.)/(8.)*x_corners(8,2)
	D_matrix(3,1)=(-1.)*(1.-yhat)*(1.-zhat)/(8.)*x_corners(1,3)+ &
			(1.-yhat)*(1.-zhat)/(8.)*x_corners(2,3)+ &
			(1.+yhat)*(1.-zhat)/(8.)*x_corners(3,3)+ &
			(-1.)*(1.+yhat)*(1.-zhat)/(8.)*x_corners(4,3)+ &
			(-1.)*(1.-yhat)*(1.+zhat)/(8.)*x_corners(5,3)+ &
			(1.-yhat)*(1.+zhat)/(8.)*x_corners(6,3)+ &
			(1.+yhat)*(1.+zhat)/(8.)*x_corners(7,3)+ &
			(-1.)*(1.+yhat)*(1.+zhat)/(8.)*x_corners(8,3)
	D_matrix(3,2)=(1.-xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(1,3)+ &
			(1.+xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(2,3)+ &
			(1.+xhat)*(1.-zhat)/(8.)*x_corners(3,3)+ &
			(1.-xhat)*(1.-zhat)/(8.)*x_corners(4,3)+ &
			(1.-xhat)*(-1.)*(1.+zhat)/(8.)*x_corners(5,3)+ &
			(1.+xhat)*(-1.)*(1.+zhat)/(8.)*x_corners(6,3)+ &
			(1.+xhat)*(1.+zhat)/(8.)*x_corners(7,3)+ &
			(1.-xhat)*(1.+zhat)/(8.)*x_corners(8,3)
	D_matrix(3,3)=(1.-xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(1,3)+ &
			(1.+xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(2,3)+ &
			(1.+xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(3,3)+ &
			(1.-xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(4,3)+ &
			(1.-xhat)*(1.-yhat)*(1.)/(8.)*x_corners(5,3)+ &
			(1.+xhat)*(1.-yhat)*(1.)/(8.)*x_corners(6,3)+ &
			(1.+xhat)*(1.+yhat)*(1.)/(8.)*x_corners(7,3)+ &
			(1.-xhat)*(1.+yhat)*(1.)/(8.)*x_corners(8,3)


	u_total=matmul(D_matrix, [uhat_x, uhat_y, uhat_z])
	JJ=FindDet(D_matrix,3)
	u_total=(1/JJ)*u_total

	u_cell_x(i)=u_total(1)
	u_cell_y(i)=u_total(2)
	u_cell_z(i)=u_total(3)
    !write(*,*) "u_cell_z=", u_cell_z(i)
	!F_shear(i,1)=0.0002*u_cell_x(i)-0.0002
	!F_shear(i,2)=0.0002*u_cell_y(i)-0.0002
	!F_shear(i,3)=0.0002*u_cell_z(i)-0.0002
	!write(*,*) "shear force on cell=", F_shear(i,:)

	!if ( abs(u_cell_y(i))> abs(u_cell_z(i)) ) then
     !   write(nflog,'(a,3f6.1)' ) 'position ', xhat,yhat,zhat
     !   write(nflog,'(a,3f10.1)' ) 'uhat_x,uhat_y,uhat_z ', uhat_x,uhat_y,uhat_z
      !  write(nflog,'(a,i6)' ) 'element number ', e
        !write(*,'(a,3f10.1)' ) 'cell velocity ', u_cell_x(i),u_cell_y(i),u_cell_z(i)
    !else
    !    write(nflog,'(a,i6)' ) 'element number ', e
    !endif
	if (i == 18) then
		write(nflog,'(i3,i6,5f8.2,e12.3)') i, cell_elem(i),xhat,yhat,zhat,u_cell_x(i),u_cell_y(i),u_cell_z(i)
	endif
enddo

end subroutine

!-------------------------------------------------------------------------
!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!-------------------------------------------------------------------------

REAL(REAL_KIND) FUNCTION FindDet(matrix, n)
    IMPLICIT NONE
    REAL (REAL_KIND), DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL (REAL_KIND) :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO

    !Calculate determinant by finding product of diagonal elements
   FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO

END FUNCTION FindDet

!This function is written only to find the determinant of a 3*3 matrix.
!The commented function (above) does it for a general square matrices but it transfomrs
! the matrix into an upper triangle matrix
!REAL FUNCTION FindDet(matrix, n)
!IMPLICIT NONE
!    REAL (REAL_KIND), DIMENSION(n,n) :: matrix
!    INTEGER, INTENT(IN) :: n

!    if (n .ne. 3) then
!        Write(*,'(a)') 'This function only computes determinant for  a 3by3 square matrix'
!    else
!      FindDet =  matrix(1,1)*(matrix(2,2)*matrix(3,3) - matrix(3,2)*matrix(2,3)) &
!       + matrix(1,2)*(matrix(3,1)*matrix(2,3) - matrix(2,1)*matrix(3,3))  &
!       + matrix(1,3)*(matrix(2,1)*matrix(3,2) - matrix(3,1)*matrix(2,2))
 !   endif
!end function

!=======================================================================
!== resid ==============================================================
! Compute the residual, r = Ax-b, its max component and norm.
! Note that this assumes that A is zero-based.
!=======================================================================
subroutine resid0 (n, Ap, Ai, Ax, x, b, r)
integer :: n, Ap(:), Ai(:), j, i, p, k, row, col
double precision :: Ax(:), x(:), b(:), r(:), rmax, rnorm, aij

r(1:n) = -b(1:n)

k = 0
do col = 1,n
!	do p = Ap(j) + 1, Ap(j+1)
!		i = Ai(p) + 1
!		aij = Ax(p)
!		r(i) = r(i) + aij * x(j)
	do p = 1, Ap(col+1) - Ap(col)
		k = k+1
		row = Ai(k) + 1
		aij = Ax(k)
		r(row) = r(row) + aij * x(col)
	enddo
enddo
rmax = 0
rnorm = 0
do i = 1, n
	rmax = max(rmax, r(i))
	rnorm = rnorm + r(i)*r(i)
enddo
rnorm = sqrt(rnorm)
write(logmsg,'(a,2e12.3)') 'Residual (A*x-b): rmax, rnorm: ', rmax, rnorm
call logger(logmsg)
end subroutine

!== resid1 ==============================================================
! Compute the residual, r = Ax-b, its max component and norm.
! Note that this assumes that A is one-based.
!=======================================================================
subroutine resid1 (n, Ap, Ai, Ax, x, b, r)
integer :: n, Ap(:), Ai(:), j, i, p, k, row, col
double precision :: Ax(:), x(:), b(:), r(:), rmax, rnorm, aij

r(1:n) = -b(1:n)

k = 0
do col = 1,n
	do p = 1, Ap(col+1) - Ap(col)
		k = k+1
		row = Ai(k)
		aij = Ax(k)
		r(row) = r(row) + aij * x(col)
	enddo
enddo
rmax = 0
rnorm = 0
do i = 1, n
	rmax = max(rmax, r(i))
	rnorm = rnorm + r(i)*r(i)
enddo
rnorm = sqrt(rnorm)
write(logmsg,'(a,2e12.3)') 'Residual (A*x-b): rmax, rnorm: ', rmax, rnorm
call logger(logmsg)
end subroutine


!-----------------------------------------------
! GB to check velocity distribution
!-----------------------------------------------
subroutine getPointVel(p,v,e)
real (REAL_KIND) :: p(3), v(3)
integer :: e
real (REAL_KIND) :: elem_vertices(8,3), fsum(3)
integer, allocatable :: cell_elem(:)
real (REAL_KIND) :: cell_in_reference(3)
real (REAL_KIND) :: coordPhy(3), coordref(3)
integer :: elem_nodes(8), j, k, II(6), signum(6), ND, Nvert
real (REAL_KIND) :: xhat, yhat, zhat, Uhat(6)
real (REAL_KIND) :: uhat_x, uhat_y, uhat_z, x_corners(8,3)
real (REAL_KIND) :: D_matrix(3,3),JJ,u_total(3)

! transfer cell position from physical to the reference elements:
call getPointElement(p,e)
if (e == 0) then
	v = 0
	return
endif
ND = 3
Nvert = 8
elem_nodes(:) = Cylinder_ELEMENT(e,:)
elem_vertices(:,1:3) = Cylinder_vertices(elem_nodes(:),1:3)
coordPhy(:) = p
call mapper(ND,Nvert,elem_vertices,coordPhy,coordref)
cell_in_reference(:) = coordref(:)

xhat = cell_in_reference(1)
yhat = cell_in_reference(2)
zhat = cell_in_reference(3)
II(:) = ElToFace(e,:)
signum=(/1,1,1,1,1,1/)
do k = 1,6
	if (e == AllSurfs(II(k),6)) then
		signum(k) = -1 ! -1 if edge is on T-
	endif
enddo
do j = 1,6
    Uhat(j) = Q_P(ElToFace(e,j))*signum(j)
enddo
uhat_x = (Uhat(1)-Uhat(2))*1/8+(Uhat(1)+Uhat(2))*xhat/8
uhat_y = (Uhat(6)-Uhat(5))*1/8+ (Uhat(5)+Uhat(6))*zhat/8
uhat_z = (Uhat(4)-Uhat(3))*1/8+ (Uhat(3)+Uhat(4))*yhat/8
do j = 1,8
	x_corners(j,:) = Cylinder_vertices(Cylinder_ELEMENT(e,j),:)
enddo
D_matrix(1,1) = (-1.)*(1.-yhat)*(1.-zhat)/(8.)*x_corners(1,1)+ &
		(1.-yhat)*(1.-zhat)/(8.)*x_corners(2,1)+ &
		(1.+yhat)*(1.-zhat)/(8.)*x_corners(3,1)+ &
		(-1.)*(1.+yhat)*(1.-zhat)/(8.)*x_corners(4,1)+ &
		(-1.)*(1.-yhat)*(1.+zhat)/(8.)*x_corners(5,1)+ &
		(1.-yhat)*(1.+zhat)/(8.)*x_corners(6,1)+ &
		(1.+yhat)*(1.+zhat)/(8.)*x_corners(7,1)+ &
		(-1.)*(1.+yhat)*(1.+zhat)/(8.)*x_corners(8,1)
D_matrix(1,2) = (1.-xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(1,1)+ &
        (1.+xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(2,1)+ &
		(1.+xhat)*(1.-zhat)/(8.)*x_corners(3,1)+ &
		(1.-xhat)*(1.-zhat)/(8.)*x_corners(4,1)+ &
		(1.-xhat)*(-1.)*(1.+zhat)/(8.)*x_corners(5,1)+ &
		(1.+xhat)*(-1.)*(1.+zhat)/(8.)*x_corners(6,1)+ &
		(1.+xhat)*(1.+zhat)/(8.)*x_corners(7,1)+ &
		(1.-xhat)*(1.+zhat)/(8.)*x_corners(8,1)
D_matrix(1,3) = (1-xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(1,1)+ &
		(1.+xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(2,1)+ &
		(1.+xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(3,1)+ &
		(1.-xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(4,1)+ &
		(1.-xhat)*(1.-yhat)*(1.)/(8.)*x_corners(5,1)+ &
		(1.+xhat)*(1.-yhat)*(1.)/(8.)*x_corners(6,1)+ &
		(1.+xhat)*(1.+yhat)*(1.)/(8.)*x_corners(7,1)+ &
		(1.-xhat)*(1.+yhat)*(1.)/(8.)*x_corners(8,1)
D_matrix(2,1) = (-1.)*(1.-yhat)*(1.-zhat)/(8.)*x_corners(1,2)+&
		(1.-yhat)*(1.-zhat)/(8.)*x_corners(2,2)+ &
		(1.+yhat)*(1.-zhat)/(8.)*x_corners(3,2)+ &
		(-1.)*(1.+yhat)*(1.-zhat)/(8.)*x_corners(4,2)+ &
		(-1.)*(1.-yhat)*(1.+zhat)/(8.)*x_corners(5,2)+ &
		(1.-yhat)*(1.+zhat)/(8.)*x_corners(6,2)+ &
		(1.+yhat)*(1.+zhat)/(8.)*x_corners(7,2)+ &
		(-1.)*(1.+yhat)*(1.+zhat)/(8.)*x_corners(8,2)
D_matrix(2,2) = (1.-xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(1,2)+ &
		(1.+xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(2,2)+ &
		(1.+xhat)*(1.-zhat)/(8.)*x_corners(3,2)+ &
		(1.-xhat)*(1.-zhat)/(8.)*x_corners(4,2)+ &
		(1.-xhat)*(-1.)*(1+zhat)/(8.)*x_corners(5,2)+ &
		(1.+xhat)*(-1.)*(1+zhat)/(8.)*x_corners(6,2)+ &
		(1.+xhat)*(1.+zhat)/(8.)*x_corners(7,2)+ &
		(1.-xhat)*(1.+zhat)/(8.)*x_corners(8,2)
D_matrix(2,3) = (1.-xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(1,2)+ &
		(1.+xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(2,2)+ &
		(1.+xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(3,2)+ &
		(1.-xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(4,2)+ &
		(1.-xhat)*(1.-yhat)*(1.)/(8.)*x_corners(5,2)+ &
		(1.+xhat)*(1.-yhat)*(1.)/(8.)*x_corners(6,2)+ &
		(1.+xhat)*(1.+yhat)*(1.)/(8.)*x_corners(7,2)+ &
		(1.-xhat)*(1.+yhat)*(1.)/(8.)*x_corners(8,2)
D_matrix(3,1) = (-1.)*(1.-yhat)*(1.-zhat)/(8.)*x_corners(1,3)+ &
		(1.-yhat)*(1.-zhat)/(8.)*x_corners(2,3)+ &
		(1.+yhat)*(1.-zhat)/(8.)*x_corners(3,3)+ &
		(-1.)*(1.+yhat)*(1.-zhat)/(8.)*x_corners(4,3)+ &
		(-1.)*(1.-yhat)*(1.+zhat)/(8.)*x_corners(5,3)+ &
		(1.-yhat)*(1.+zhat)/(8.)*x_corners(6,3)+ &
		(1.+yhat)*(1.+zhat)/(8.)*x_corners(7,3)+ &
		(-1.)*(1.+yhat)*(1.+zhat)/(8.)*x_corners(8,3)
D_matrix(3,2) = (1.-xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(1,3)+ &
		(1.+xhat)*(-1.)*(1.-zhat)/(8.)*x_corners(2,3)+ &
		(1.+xhat)*(1.-zhat)/(8.)*x_corners(3,3)+ &
		(1.-xhat)*(1.-zhat)/(8.)*x_corners(4,3)+ &
		(1.-xhat)*(-1.)*(1.+zhat)/(8.)*x_corners(5,3)+ &
		(1.+xhat)*(-1.)*(1.+zhat)/(8.)*x_corners(6,3)+ &
		(1.+xhat)*(1.+zhat)/(8.)*x_corners(7,3)+ &
		(1.-xhat)*(1.+zhat)/(8.)*x_corners(8,3)
D_matrix(3,3) = (1.-xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(1,3)+ &
		(1.+xhat)*(1.-yhat)*(-1.)/(8.)*x_corners(2,3)+ &
		(1.+xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(3,3)+ &
		(1.-xhat)*(1.+yhat)*(-1.)/(8.)*x_corners(4,3)+ &
		(1.-xhat)*(1.-yhat)*(1.)/(8.)*x_corners(5,3)+ &
		(1.+xhat)*(1.-yhat)*(1.)/(8.)*x_corners(6,3)+ &
		(1.+xhat)*(1.+yhat)*(1.)/(8.)*x_corners(7,3)+ &
		(1.-xhat)*(1.+yhat)*(1.)/(8.)*x_corners(8,3)

u_total = matmul(D_matrix, [uhat_x, uhat_y, uhat_z])
JJ = FindDet(D_matrix,3)
u_total = (1/JJ)*u_total
v = u_total

end subroutine

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
subroutine test_getPointVel
real(REAL_KIND) :: p(3), v(3), p0(3), dxyz(8,3), delta
integer :: i, k, e, node(8), nd, ix, iy, iz
logical :: use_simple = .true.

delta = 0.01
i = 0
do ix = -1,1,2
do iy = -1,1,2
do iz = -1,1,2
	i = i+1
	dxyz(i,:) = delta*[ix, iy, iz]
enddo
enddo
enddo
if (use_simple) then
	do
		write(*,*) 'Enter p (x,y,z):'
		read(*,*) p
		if (p(3) == 0) stop
		call getPointVel(p,v,e)
		if (e == 0) then
			write(*,*) 'Outside mesh, e = 0'
			cycle
		endif
		write(*,'(a,3e12.3)') 'Vel: ',v

		node = cylinder_element(e,:)
		write(*,'(a,i6,a,8i6)') 'p is in element: ',e,' with nodes: ',node
		do k = 1,8
			nd = node(k)
			p0 = cylinder_vertices(nd,:)
			write(*,'(a,3i8,3f10.2)') 'element, k, node, p0: ',e,k,nd,p0
		enddo

	enddo
else
	p = [-180, 0, 900]
	call getPointElement(p,e)
	node = cylinder_element(e,:)
	write(*,'(a,3f10.2)') 'point p: ',p
	write(*,'(a,i6,a,8i6)') 'p is in element: ',e,' with nodes: ',node
	do k = 1,8
		nd = node(k)
		p0 = cylinder_vertices(nd,:)
		write(*,'(a,3i8,3f10.2)') 'element, k, node, p0: ',e,k,nd,p0
		do i = 1,8
			p = p0 + dxyz(i,:)
			call getPointVel(p,v,e)
			write(*,'(a,i8,3f10.2,3x,3f10.1)') 'e,p,v: ',e,p,v
		enddo
	enddo
endif
stop

end subroutine
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
subroutine makeVelDist(x0,y0,z0)
real(REAL_KIND) :: x0, y0, z0

real(REAL_KIND) :: x, y, z, r, p(3), delta, vv(3)
real(REAL_KIND), allocatable :: v(:,:,:)
integer :: ix, iy, iz, i, e

write(*,*) 'makeVelDist: writing to velocity.dat'
delta = 0.5
allocate(v(-200:200,-200:200,3))
do ix = -200,200
	write(*,*) 'ix: ',ix
	do iy = -200,200
		x = ix
		y = iy
		r = sqrt(x*x + y*y)
		if (r > 200 - delta) then
			v(ix,iy,:) = 0
		else
			p = [x, y, z0]
			call getPointVel(p,v(ix,iy,:),e)
			if (e == 0) then
				v(ix,iy,:) = 0
			endif
		endif
	enddo
enddo
open(nfvel,file='velocity.dat',status='replace')
write(nfvel,'(a,f8.1)') 'z = ',z0
do i = 1,3
	if (i == 1) then
		write(nfvel,*) 'v_x'
	elseif (i == 2) then
		write(nfvel,*) 'v_y'
	elseif (i == 3) then
		write(nfvel,*) 'v_z'
	endif
	do ix = -200,200
		write(nfvel,'(401f9.1)') v(ix,:,i)
	enddo
enddo
deallocate(v)

write(nfvel,*)
x = 0
y = 0
write(nfvel,'(a,2f8.1,a)') 'velocity at: ',x0,y0,' over range of z'
do iz = 0,1000
	z = iz
	p = [x0, y0, z]
	call getPointVel(p,vv,e)
	if (e == 0) then
		vv = 0
	endif
	write(nfvel,'(f8.1,2x,3f9.1)') z,vv
enddo
close(nfvel)
stop

end subroutine

#if 0
!---------------------------------------------------------------------------------
! Stiffness matrix is MxM, with ne non-zero entries.
! A, row_ind, col_ptr is CSC (compressed spare column) format of stiffness matrix
! MR1 = Maximum number of equations = 21000
!---------------------------------------------------------------------------------
subroutine solve_harwell_umf(M, ne, A, row_ind, col_ptr, rhs, solution)
real(REAL_KIND) :: A(:), rhs(:), solution(:)
integer :: M, ne, row_ind(:), col_ptr(:)

!PARAMETER (LINDEX=2500000,LVALUE=7500000,NNZC=700000)
!REAL(4) SEQ(LVALUE),UWORK(4*MR1),X(MR1)
real(4), allocatable :: SEQ(:), UWORK(:), X(:), R1(:)
real(4) :: CNTL(10),RINFO(20)	!,RTM(MR1)

INTEGER INFO(40),ICNTL(20),KEEP(20)
!integer ICOL(MR1+1), IROW(NNZC),IWORK(LINDEX)	!IDST(30),JDST(30),
integer, allocatable :: IWORK(:)

integer :: row, col, k, j, MR1, LINDEX, LVALUE, I

MR1 = 2*M
LINDEX = 10*ne + 2*M + 1	! 3*ne+2*n+1
LVALUE = 200*ne		! 2*ne
!NSZF = n
write(*,*) 'solve_harwell_umf: M, ne: ',M,ne
write(*,*) 'LINDEX, LVALUE: ',LINDEX, LVALUE

allocate(SEQ(LVALUE))
allocate(UWORK(4*MR1))
allocate(X(MR1))
allocate(R1(MR1))
allocate(IWORK(LINDEX))

CALL UMS2IN(ICNTL,CNTL,KEEP)

ICNTL(4)=0
ICNTL(6)=1
ICNTL(7)=32

! Need to convert from CSC to sparse triplet format
k = 0
col = 1
do
	do j = 1,col_ptr(col+1) - col_ptr(col)
		k = k+1
		SEQ(k) = A(k)
		row = row_ind(k)
		IWORK(k) = row
		IWORK(k+ne) = col
	enddo
	col = col+1
	if (col > M) exit
enddo

!do k = 1,ne
!	write(nflog,'(2i8,e12.3)') iwork(k),iwork(k+ne),seq(k)
!enddo

do i = 1,M
	R1(i) = rhs(i)		! convert real(8) to real(4)
enddo

write(*,*) 'call UMS2FA'
CALL UMS2FA(M, ne, 0, .FALSE., LVALUE, LINDEX, SEQ, IWORK, KEEP, CNTL, ICNTL, INFO, RINFO)

WRITE(*,*) 'LINDEX - DIM: ",LINDEX," REQ: ',INFO(18)
WRITE(*,*) 'LVALUE - DIM: ",LVALUE," REQ: ',INFO(20)

IF(INFO(1) .NE. 0) THEN
	WRITE(*,*) 'ERROR IN FACTORIZATION: ',INFO(1)
	STOP
ENDIF

write(*,*) 'call UMS2SO'
CALL UMS2SO(M, 0, .FALSE., LVALUE, LINDEX, SEQ, IWORK, KEEP, R1, X, UWORK, CNTL, ICNTL, INFO, RINFO)

IF(INFO(1) .NE. 0) THEN
	WRITE(*,*) 'ERROR IN SOLVE: ',INFO(1)
	STOP
ENDIF

write(*,*) 'solution:'
DO I = 1,M
!	R1(I) = X(I)
	solution(I) = X(I)
	write(nflog,'(i6,e12.3)') I,solution(I)
ENDDO

deallocate(SEQ)
deallocate(UWORK)
deallocate(X)
deallocate(R1)
deallocate(IWORK)

end subroutine

#endif

!--------------------------------------------------------------------------------------
! Save the CSC format stiffness matrix and the rhs.
!--------------------------------------------------------------------------------------
subroutine writeHB(M,nnz,col_ptr,row_ind,A,rhs)
integer :: M, nnz, col_ptr(:), row_ind(:)
real(REAL_KIND) :: A(:), rhs(:)
integer, parameter :: nfhb=30, nfrhs=31

open(nfhb,file='HB.dat',status='replace')
write(nfhb,*) 'M, nnz: ',M,nnz
write(nfhb,'(11i7)') col_ptr
write(nfhb,'(11i7)') row_ind
write(nfhb,'(5e16.8)') A
close(nfhb)

open(nfrhs,file='RHS.dat',status='replace')
write(nfrhs,'(5e16.8)') rhs
close(nfrhs)

stop
end subroutine


end module
