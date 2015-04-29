! IL module specifications
! Each module requires a parameter named IL2_NP, IL4_NP, etc
! conveying the number of real parameters needed to characterize the
! state of a cognate T cell with respect to this cytokine.
! Data exchanges with the module then occur via a section of
! the type(cog_type) structure which is a subarray of IL2_NP
! REALs (for example).
! Initialization of the parameters for a cytokine is by reading
! data from a file.  The name of the file holding all the parameter
! data is passed as the argument to IL2_setparameters() (for example).

!------------------------------------------------------------------------------
! IL2 module
! Implements Fallon-Lauffenburger (2000) model of IL2/CD25 receptor dynamics.
! Matlab code: Lauff_a.m and Lauff_deriv_a.m
! Note the error in Fallon-Lauffenburger in the expression for rate of change
! of Li - the scale factor M_pM to convert molar (M) to picomolar (pM) is missing.
! Added to the F-L formulation is a model for the rate of synthesis of IL2,
! as a function of S_TCR and S_Rec.  Note that IL2 production rate (as
! computed by IL2_rate()) is separate from the balance of IL2 binding and
! release from the receptor and recycling from the endosomes.
! To prevent the uptake of IL2 (binding to the receptor) from exceeding the
! amount of IL2 available at the site, the apparent IL2 conc must be limited.
! C_IL2 must be set so that Kf*C_IL2*state(IL2_Rec_Surf)*dt does not exceed
! the number of available molecules, which is IL2conc*(V/L_um3)*(NAvo/M_pM)
! One way to achieve this is to compute the concentration of IL-2 as another
! ODE variable.
!
! xcase = .true. allows the concentration in the site to vary over the timestep,
! thus avoiding problems when the concentration is low.
! The vector xstate() is state() with an extra element to hold the IL-2 conc.
!------------------------------------------------------------------------------
module IL2

use rkf45

implicit none

private

integer, public, parameter :: IL2_NP = 6
integer, public, parameter :: IL2_Rec_Surf     = 1
integer, public, parameter :: IL2_Complex_Surf = 2
integer, public, parameter :: IL2_Rec_Int      = 3
integer, public, parameter :: IL2_Complex_Int  = 4
integer, public, parameter :: IL2_Ligand_Int   = 5
integer, public, parameter :: IL2_Store        = 6
real, public, parameter :: IL2_CONSTIT_RATE = 0.0001     ! constitutive rate of IL2 production mols/min/um^3 (0.0001)
real, public, parameter :: IL2_DECAY_RATE = 0.05            !0.05  ! TESTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real, public, parameter :: NAvo = 6.022e23          ! Avogadro's number
real, public, parameter :: M_pM = 1.0e12            ! convert M -> pM
real, public, parameter :: L_um3 = 1.0e15           ! convert L -> um^3

real, parameter :: Kr = 0.0138;             ! dissociation rate (min^-1)
real, parameter :: Kf = Kr/11.1;            ! association rate constant (pM^-1 min^-1)
real, parameter :: Kre = 8*Kr;              ! dissociation rate constant, endosome (min^-1)
real, parameter :: Kfe = Kre/1000;          ! association rate constant, endosome (pM^-1 min^-1)
real, parameter :: Kt = 0.007;              ! constitutive internalization rate constant (min^-1)
!real, parameter :: Ksyn = 0.0011;          ! induced receptor synthesis rate (min^-1) (NOT USED)
real, parameter :: Ke = 0.04;               ! internalization rate constant (min^-1)
real, parameter :: Kx = 0.15;               ! recycling rate constant (min^-1)
real, parameter :: Kh = 0.035;              ! degradation rate constant (min^-1)
real, parameter :: Kstim = 0.2;            ! CD25 stimulation rate constant (min^-1)
real, parameter :: Tdecay = 6;              ! integrated CD25 stimulation decay time constant (hour)
real, parameter :: Kdecay = 1/(60*Tdecay);  ! integrated CD25 stimulation decay rate (min^-1)
real, parameter :: Vs = 10.5;               ! constitutive receptor synthesis rate (min^-1) (10.5)
real, parameter :: Vendo = 1.0e-14;         ! total endosomal volume (L/cell)

integer, parameter :: CD8=2
integer, parameter :: IL2_HILL_N = 4        ! should be 6?
real, parameter :: IL2_HILL_C = 1000
real, parameter :: IL2_MAX_RATE = 200     ! max IL2 synthesis rate /min
integer, parameter :: CD25_HILL_N = 2
real, parameter :: CD25_HILL_C = 10000
real, parameter :: CD25_MAX_RATE = 10     ! max IL2 synthesis rate /min
real, parameter :: IL2_PROD_CD8_CD4 = 0.2

integer :: celltype, ID
real :: S_TCR, C_IL2, delt, vol
logical :: producing

public :: IL2_setparameters, IL2_update, IL2_init_state, IL2_divide, IL2_deriv
contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine IL2_setparameters(filename)
character(LEN=*) :: filename
end subroutine

!---------------------------------------------------------------------
! ctype = type of T cell
! t = start time (min)
! TCRstim = integrated TCR stimulation
! On entry state(:) holds the array of (REAL) state variables for a cell,
! and statep(:) holds the arrays of time derivatives of state variables
! On return the arrays hold the updated values.
! first = true on first call for this cell (setup), false otherwise
! IL2prod = true if IL-2/CD25 production, false otherwise
! IL2conc = local concentration of the cytokine (input)
! dt = time step (min)
! mrate = nett rate of IL2 secretion (secreted - consumed) (molecules/min)
!
! The results from this code agree very closely with those from Matlab,
! (Lauff_a.m) when xcase=.false.
!---------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE:  Need to differentiate between CD4 and CD8 T cells, since CD8 cells DO NOT MAKE IL-2 !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  IL-2 is made by Th1 CD4 cells, not by Th2 CD4 or CD8                                !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine IL2_update(cogID,ctype,t,TCRstim,state,statep,first,IL2prod,IL2conc,Vc,dt,mrate,dbgflag)
integer :: cogID, ctype
real, intent(INOUT) :: state(:), statep(:), IL2conc
real, intent(IN) :: t, TCRstim, dt, Vc
logical :: first, IL2prod, dbgflag
logical :: xcase = .true.
real :: mrate
real :: tstart, tend, relerr, abserr, an
real :: save_state(IL2_NP), save_statep(IL2_NP)
real :: xstate(IL2_NP+1),xstatep(IL2_NP+1)
integer :: flag, k	! nfout = 11

!if (cogID >= 53) then
!    write(*,'(i4,7f8.2)') cogID,TCRstim,state
!endif
celltype = ctype
S_TCR = TCRstim
C_IL2 = IL2conc
ID = cogID
delt = dt
vol = Vc
producing = IL2prod
tstart = t
tend = tstart + dt
if (dbgflag) then
    save_state = state
    save_statep = statep
endif

!if (producing) then
!    dIL2_dt = IL2_rate(state(IL2_Store))
!else
!    dIL2_dt = 0
!endif
!a = Kf*state(IL2_Rec_Surf)
!b = Kr*state(IL2_Complex_Surf) + Kx*state(IL2_Ligand_Int)*Vendo*NAvo/M_pM + dIL2_dt
!C_IL2max = 0.99*(IL2conc*Vc*NAvo/(dt*M_pM*L_um3) + b)/a
!C_IL2 = min(C_IL2max,IL2conc)
!mrate = -a*C_IL2 + b

!write(*,*) mrate,IL2conc,C_IL2max,C_IL2

!mrate = (-Kf*C_IL2*state(IL2_Rec_Surf) + Kr*state(IL2_Complex_Surf) &
!    + Kx*state(IL2_Ligand_Int)*Vendo*NAvo/M_pM + dIL2_dt)   ! Molecules.min^-1 [x M_pM.min/L/Navo to get pM]
!if (mrate < 0) then
!    C_IL2max = IL2conc*(Vc/L_um3)*(NAvo/M_pM)/(-mrate*dt)
!    C_IL2max = (-mrate/Vc)*dt*M_pM*L_um3/Navo
!    C_IL2 = min(IL2conc,C_IL2max)
!    write(*,*) mrate,IL2conc,C_IL2max,C_IL2
!else
!    C_IL2 = IL2conc
!endif

if (xcase) then
    xstate(1:IL2_NP) = state(1:IL2_NP)
    xstate(IL2_NP+1) = IL2conc              ! add ambient IL-2 as an extra state variable
    xstatep(1:IL2_NP) = statep(1:IL2_NP)
    xstatep(IL2_NP+1) = 0
endif
if (first) then
    if (xcase) then
        call xIL2_deriv(tstart,xstate,xstatep)
    else
        call IL2_deriv(tstart,state,statep)
    endif
    flag = 1
else
    flag = 2
endif

abserr = sqrt ( epsilon ( abserr ) )
relerr = sqrt ( epsilon ( relerr ) )

if (dbgflag) write(*,*) 'call r4_rkf45'
if (xcase) then
    call r4_rkf45 ( xIL2_deriv, IL2_NP+1, xstate, xstatep, tstart, tend, relerr, abserr, flag )
else
    call r4_rkf45 ( IL2_deriv, IL2_NP, state, statep, tstart, tend, relerr, abserr, flag )
endif
if (flag /= 2) then
    write(*,*) 'r4_rkf45: bad flag: ',cogID,flag,first
    write(*,*) 'IL2: ',IL2conc
    write(*,'(6f8.2)') save_state,save_statep
    if (xcase) then
        write(*,'(7f8.2)') xstate,xstatep
    else
        write(*,'(6f8.2)') state,statep
    endif
    stop
endif
if (dbgflag) write(*,*) 'did r4_rkf45: flag: ',flag

!mrate = 0   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!if (mrate < 0) then
!    write(*,'(i3,5f8.4)') ctype,mrate,-Kf*C_IL2*state(IL2_Rec_Surf),Kr*state(IL2_Complex_Surf),Kx*state(IL2_Ligand_Int)*Vendo*NAvo/M_pM,dIL2_dt
!endif
!dtsub = dt/nsub
!
!Rs = state(IL2_Rec_Surf)       ! surface free receptors
!Cs = state(IL2_Complex_Surf)    ! surface complex
!Ri = state(IL2_Rec_Int)        ! endosomal free receptor
!Ci = state(IL2_Complex_Int)     ! endosomal complex
!Li = state(IL2_Ligand_Int)      ! free endosomal IL2
!S_Rec = state(IL2_Store)       ! integrated CD25 signalling
!
!do k = 1,nsub
!
!    if (producing) then
!        dIL2_dt = IL2_rate(ctype,S_Rec,S_TCR)
!        dCD25_dt = CD25_rate(ctype,S_Rec,S_TCR)
!    else
!        dIL2_dt = 0
!        dCD25_dt = 0
!    endif
!    dRs_dt = -Kf*C_IL2*Rs + Kr*Cs - Kt*Rs + Vs + dCD25_dt    ! Rs
!    dCs_dt = Kf*C_IL2*Rs - (Kr + Ke)*Cs                      ! Cs
!    dRi_dt = -Kfe*Li*Ri + Kre*Ci + Kt*Rs - Kh*Ri             ! Ri
!    dCi_dt = Kfe*Li*Ri - (Kre + Kh)*Ci + Ke*Cs               ! Ci
!    dLi_dt = (-Kfe*Li*Ri + Kre*Ci)*M_pM/(Vendo*NAvo) - Kx*Li ! Li
!    dS_dt  = Kstim*Cs - Kdecay*S_Rec                        ! S_Rec
!
!    Rs     = Rs     + dtsub*dRs_dt
!    Cs     = Cs     + dtsub*dCs_dt
!    Ri     = Ri     + dtsub*dRi_dt
!    Ci     = Ci     + dtsub*dCi_dt
!    Li     = Li     + dtsub*dLi_dt
!    S_Rec = S_Rec + dtsub*dS_dt
!
!enddo
!
!state(IL2_Rec_Surf)    = Rs
!state(IL2_Complex_Surf) = Cs
!state(IL2_Rec_Int)     = Ri
!state(IL2_Complex_Int)  = Ci
!state(IL2_Ligand_Int)   = Li
!state(IL2_Store)        = S_Rec
!
!mrate = (-Kf*C_IL2*Rs + Kr*Cs + Kx*Li*Vendo*NAvo/M_pM + dIL2_dt)   ! Molecules.min^-1 (x M_pM.min/L/Navo to get pM)
!mrate = dLi_dt

if (dbgflag) then
    write(*,'(a,L2,3f8.2)') 'IL2_update:            ',producing, C_IL2, S_TCR
    write(*,'(a,6f8.2)') 'Rs, Cs, Ri, Ci, Li, S: ',state(1:6)
    write(*,'(a,6f8.2)') 'rates:                 ',statep(1:6)
    write(*,'(a,5f8.2)') 'IL2 rate:              ',mrate
endif
if (xcase) then
    state(1:IL2_NP) = xstate(1:IL2_NP)
    IL2conc = xstate(IL2_NP+1)
    statep(1:IL2_NP) = xstatep(1:IL2_NP)
    mrate = xstatep(IL2_NP+1)    ! actually rate of increase of IL-2 conc.
endif
do k = 1,6
    if (state(k) < 0) then
        write(*,'(a,i3,6f10.4)') 'cogID, t, S_TCR, IL2conc, mrate: ',cogID,t,S_TCR,IL2conc,mrate, &
            IL2_rate(), CD25_rate(state(IL2_NP))
        write(*,*) 'Bad state: ',state

        an = S_TCR**IL2_HILL_N
        write(*,'(i3,4e12.3)') IL2_HILL_N,S_TCR,an,an/(IL2_HILL_C**IL2_HILL_N + an), &
            IL2_MAX_RATE*an/(IL2_HILL_C**IL2_HILL_N + an)
        stop
    endif
enddo
end subroutine

!--------------------------------------------------------------------------
! Initialize the state vector
!--------------------------------------------------------------------------
subroutine IL2_init_state(state,statep)
real, intent(INOUT) :: state(:), statep(:)

statep = 0
state = 0
state(IL2_Rec_Surf) = 1500
state(IL2_Complex_Surf) = 1
state(IL2_Rec_Int) = 300
state(IL2_Complex_Int) = 1
state(IL2_Ligand_Int) = 1
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine IL2_deriv(t,y,yp)
real :: t, y(*), yp(*)
real :: Li, Rs, Ri, Cs, Ci, S_Rec, dIL2_dt, dCD25_dt

Rs = y(IL2_Rec_Surf)       ! surface free receptors
Cs = y(IL2_Complex_Surf)    ! surface complex
Ri = y(IL2_Rec_Int)        ! endosomal free receptor
Ci = y(IL2_Complex_Int)     ! endosomal complex
Li = y(IL2_Ligand_Int)      ! free endosomal IL2
S_Rec = y(IL2_Store)       ! integrated CD25 signalling

if (producing) then
    dIL2_dt = IL2_rate()
    dCD25_dt = CD25_rate(S_Rec)
else
    dIL2_dt = 0
    dCD25_dt = 0
endif

yp(IL2_Rec_Surf)    = -Kf*C_IL2*Rs + Kr*Cs - Kt*Rs + dCD25_dt + Vs    ! Rs
yp(IL2_Complex_Surf) = Kf*C_IL2*Rs - (Kr + Ke)*Cs                      ! Cs
yp(IL2_Rec_Int)     = -Kfe*Li*Ri + Kre*Ci + Kt*Rs - Kh*Ri             ! Ri
yp(IL2_Complex_Int)  = Kfe*Li*Ri - (Kre + Kh)*Ci + Ke*Cs               ! Ci
yp(IL2_Ligand_Int)   = (-Kfe*Li*Ri + Kre*Ci)*M_pM/(Vendo*NAvo) - Kx*Li ! Li
yp(IL2_Store)        = Kstim*Cs - Kdecay*S_Rec                        ! S_Rec

end subroutine

!--------------------------------------------------------------------------
! Computes mrate
!--------------------------------------------------------------------------
subroutine xIL2_deriv(t,y,yp)
real :: t, y(*), yp(*)
real :: Li, Rs, Ri, Cs, Ci, S_Rec, dIL2_dt, dCD25_dt, mrate

Rs = y(IL2_Rec_Surf)       ! surface free receptors
Cs = y(IL2_Complex_Surf)    ! surface complex
Ri = y(IL2_Rec_Int)        ! endosomal free receptor
Ci = y(IL2_Complex_Int)     ! endosomal complex
Li = y(IL2_Ligand_Int)      ! free endosomal IL2
S_Rec = y(IL2_Store)       ! integrated CD25 signalling
C_IL2 = y(IL2_NP+1)         ! IL2 conc

if (producing) then
    dIL2_dt = IL2_rate()
    dCD25_dt = CD25_rate(S_Rec)
else
    dIL2_dt = 0
    dCD25_dt = 0
endif

yp(IL2_Rec_Surf)    = -Kf*C_IL2*Rs + Kr*Cs - Kt*Rs + dCD25_dt + Vs    ! Rs
yp(IL2_Complex_Surf) = Kf*C_IL2*Rs - (Kr + Ke)*Cs                      ! Cs
yp(IL2_Rec_Int)     = -Kfe*Li*Ri + Kre*Ci + Kt*Rs - Kh*Ri             ! Ri
yp(IL2_Complex_Int)  = Kfe*Li*Ri - (Kre + Kh)*Ci + Ke*Cs               ! Ci
yp(IL2_Ligand_Int)   = (-Kfe*Li*Ri + Kre*Ci)*M_pM/(Vendo*NAvo) - Kx*Li ! Li
yp(IL2_Store)        = Kstim*Cs - Kdecay*S_Rec                        ! S_Rec
mrate = -Kf*C_IL2*Rs + Kr*Cs + Kx*Li*Vendo*NAvo/M_pM + dIL2_dt
    ! Molecules.min^-1 [x M_pM.min/L/Navo to get pM]
yp(IL2_NP+1) = mrate*M_pM*L_um3/(vol*Navo)      !!!!!!!!!!!!  Need to check this !!!!!!!!!!

!if (ID == 2) then
!    write(*,'(a,5e12.3)') 'ID: ', Cs, dIL2_dt, C_IL2, yp(IL2_store), y(IL2_store)
!endif
end subroutine

!--------------------------------------------------------------------------
! Rate of synthesis of CD25 receptor molecules.
! S_Rec = integrated CD25 stimulation
! S_TCR = integrated TCR stimulation
!--------------------------------------------------------------------------
real function CD25_rate(S_Rec)
real :: S_Rec
!real :: n,an,afactor
real :: Smax, k, rmax, c
integer :: n

Smax = 10000
!k = 0.01
k = CD25_MAX_RATE/Smax
if (S_TCR < Smax) then
    rmax = k*S_TCR
else
    rmax = k*Smax
endif
c = CD25_HILL_C
n = CD25_HILL_N
CD25_rate = rmax*S_Rec**n/(c**n + S_Rec**n)
if (celltype == CD8) then
	CD25_rate = CD25_rate*IL2_PROD_CD8_CD4
endif
!write(*,*) 'S_Rec, CD25_rate: ',S_Rec,CD25_rate

!n = IL2_HILL_N
!an = S_TCR**n
!afactor = an/(IL2_HILL_C**n + an)
!CD25_rate = Vs + CD25_MAX_RATE*afactor
!if (ctype == CD8_CELL) then
!	CD25_rate = CD25_rate*IL2_PROD_CD8_CD4
!endif
end function

!--------------------------------------------------------------------------
! Rate of synthesis of IL2 molecules.
! S_Rec = integrated CD25 stimulation
! S_TCR = integrated TCR stimulation
! Only CD4 cells synthesize IL2
!--------------------------------------------------------------------------
real function IL2_rate()
real :: n,an,afactor

if (celltype == CD8) then
	IL2_rate = 0
	return
endif
n = IL2_HILL_N
an = S_TCR**n
afactor = an/(IL2_HILL_C**n + an)
IL2_rate = IL2_MAX_RATE*afactor

end function

!---------------------------------------------------------------------
! Handles the change in the state vector that occurs on cell division.
! On entry state1() is the state vector of the dividing cell.
! On return state1() holds the new state for the dividing cell,
! and state2() holds the state of the other progeny cell
!---------------------------------------------------------------------
subroutine IL2_divide(state1,state2,statep1,statep2)
real :: state1(:), state2(:), statep1(:), statep2(:)
state2 = state1/2
state1 = state2
statep2 = statep1/2
statep1 = statep2
end subroutine

end module IL2

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! IL7 module
module IL7
implicit none
integer, parameter :: IL7_NP = 5
integer, parameter :: IL7_Rfree = 1
integer, parameter :: IL7_Rcomp = 2
integer, parameter :: IL7_Rfree_Int = 3
integer, parameter :: IL7_Rcomp_Int = 4
integer, parameter :: IL7_Store = 5

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine IL7_setparameters(filename)
character(LEN=*) :: filename
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine IL7_update(val)
real :: val(:)
real :: Y
Y = val(IL7_Store)
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine IL7_divide(state1,state2,statep1,statep2)
real :: state1(:), state2(:), statep1(:), statep2(:)
state2 = state1/2
state1 = state2
statep2 = statep1/2
statep1 = statep2
end subroutine

end module IL7

module IL_dummy
integer, parameter :: IL4_NP = 0
integer, parameter :: IL9_NP = 0
integer, parameter :: IL15_NP = 0
integer, parameter :: IL21_NP = 0
contains

subroutine IL4_setparameters(filename)
character(LEN=*) :: filename
end subroutine
subroutine IL4_update(val)
real :: val(:)
end subroutine

subroutine IL4_divide(state1,state2,statep1,statep2)
real :: state1(:), state2(:), statep1(:), statep2(:)
state2 = state1/2
state1 = state2
statep2 = statep1/2
statep1 = statep2
end subroutine

subroutine IL9_setparameters(filename)
character(LEN=*) :: filename
end subroutine

subroutine IL9_update(val)
real :: val(:)
end subroutine

subroutine IL9_divide(state1,state2,statep1,statep2)
real :: state1(:), state2(:), statep1(:), statep2(:)
state2 = state1/2
state1 = state2
statep2 = statep1/2
statep1 = statep2
end subroutine

subroutine IL15_setparameters(filename)
character(LEN=*) :: filename
end subroutine

subroutine IL15_update(val)
real :: val(:)
end subroutine

subroutine IL15_divide(state1,state2,statep1,statep2)
real :: state1(:), state2(:), statep1(:), statep2(:)
state2 = state1/2
state1 = state2
statep2 = statep1/2
statep1 = statep2
end subroutine

subroutine IL21_setparameters(filename)
character(LEN=*) :: filename
end subroutine
subroutine IL21_update(val)
real :: val(:)
end subroutine

subroutine IL21_divide(state1,state2,statep1,statep2)
real :: state1(:), state2(:), statep1(:), statep2(:)
state2 = state1/2
state1 = state2
statep2 = statep1/2
statep1 = statep2
end subroutine

end module IL_dummy

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! CD69 module
module CD69
implicit none

private
real :: K1_S1PR1 != 0.01
real :: K2_S1PR1 != 0.05
real :: K1_CD69 != 0.04
real :: K2_CD69 != 0.01

public :: CD69_setparameters, S1PR1_update

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine CD69_setparameters(ks1,ks2,kc1,kc2)
real :: ks1, ks2, kc1, kc2

K1_S1PR1 = ks1
K2_S1PR1 = ks2
K1_CD69 = kc1
K2_CD69 = kc2
end subroutine

!---------------------------------------------------------------------
! TCR stimulation (stimrate) drives CD69 expression (by K1_CD69, limited
! by (1-CD69) factor), which also decays (by K2_CD69).
! The S1PR1 (S1P receptor level) tends to grow (at K1_S1PR1, limited by
! (1-S1PR1) factor), and is held down by CD69 (rate constant K2_S1PR1).
! When a T cell enters the LN from the blood the S1PR1 level is low
! (because of high S1P in the blood) and CD69 is low.
! In the case of a non-cognate cell, CD69 stays low, but S1PR1 rises over
! a period of hours, making the cell susceptible to chemotaxis towards
! the exits.
! In the case of cognate cells, TCR stimulation drives CD69 up, and
! the high CD69 level keeps S1PR1 expression low, making the cell
! insensitive to the chemotactic influence of S1P, and thus reducing
! the probability of cell exit.
! NEEDED:
! Need a model for the decrease in CD69 with time, or (more likely) 
! number of divisions.  How does continuing TCR signalling combine
! with that?
!---------------------------------------------------------------------
subroutine S1PR1_update(CD69,S1PR1,stimrate,dt)
real :: CD69, S1PR1, stimrate, dt

!zp(1) = KC1*(1-CD69)*TCR - KC2*CD69;
!zp(2) = KS1*(1-S1PR1) - KS2*CD69*S1PR1;

CD69 = CD69 + (K1_CD69*(1-CD69)*stimrate - K2_CD69*CD69)*dt
CD69 = max(CD69,0.0)
CD69 = min(CD69,1.0)
S1PR1 = S1PR1 + (K1_S1PR1*(1-S1PR1) - K2_S1PR1*CD69*S1PR1)*dt
S1PR1 = max(S1PR1,0.0)
S1PR1 = min(S1PR1,1.0)
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine CD69_update1(CD69,stimrate,dt)
real :: CD69, stimrate, dt

!CD69 = CD69 + (K_CD69*stimrate - Kd_CD69*CD69)*dt
!CD69 = max(CD69,0.0)
end subroutine


end module
