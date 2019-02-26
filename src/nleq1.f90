module nleq1_mod

implicit none

interface
	subroutine func(N,X,F,IFAIL)
	integer :: N, IFAIL
	real(8) :: X(N), F(N)
	end subroutine
end interface

interface
	subroutine jacfunc(N,LDJAC,X,DFDX,IFAIL)
	integer :: N, LDJAC, IFAIL
	real(8) :: X(N), DFDX(LDJAC,N)
	end subroutine
end interface

contains

!     ------------------------------------------------------------
!
!   Title
!
!     Numerical solution of nonlinear (NL) equations (EQ)
!     especially designed for numerically sensitive problems.
!     (E)asy-to-use driver routine for NLEQ1.
!
!   Written by        U. Nowak, L. Weimann
!   Purpose           Solution of systems of highly nonlinear equations
!   Method            Damped affine invariant Newton method
!              (see references below)
!   Category          F2a. - Systems of nonlinear equations
!   Keywords          Nonlinear equations, Newton methods
!   Version           2.3
!   Revision          November 1991
!   Latest Change     November 1991
!   Library           CodeLib
!   Code              Fortran 77, Double Precision
!   Environment       Standard Fortran 77 environment on PC's,
!              workstations and hosts.
!   Copyright     (c) Konrad-Zuse-Zentrum fuer
!              Informationstechnik Berlin (ZIB)
!              Takustrasse 7, D-14195 Berlin-Dahlem
!              phone : + 49/30/84185-0
!              fax   : + 49/30/84185-125
!   Contact           Lutz Weimann
!              ZIB, Division Scientific Computing,
!                   Department Scientific Software
!              phone : + 49/30/84185-185
!              fax   : + 49/30/84185-107
!              e-mail: weimann@zib.de
!
!     References:
!
!     /1/ P. Deuflhard:
!         Newton Methods for Nonlinear Problems. -
!         Affine Invariance and Adaptive Algorithms.
!         Series Computational Mathematics 35, Springer (2004)
!
!     /2/ U. Nowak, L. Weimann:
!         A Family of Newton Codes for Systems of Highly Nonlinear
!         Equations - Algorithm, Implementation, Application.
!         ZIB, Technical Report TR 90-10 (December 1990)
!
!  ---------------------------------------------------------------
!
!* Licence
!    You may use or modify this code for your own non commercial
!    purposes for an unlimited time.
!    In any case you should not deliver this code without a special
!    permission of ZIB.
!    In case you intend to use the code commercially, we oblige you
!    to sign an according licence agreement with ZIB.
!
!* Warranty
!    This code has been tested up to a certain level. Defects and
!    weaknesses, which may be included in the code, do not establish
!    any warranties by ZIB. ZIB does not take over any liabilities
!    which may follow from aquisition or application of this code.
!
!* Software status
!    This code is under care of ZIB and belongs to ZIB software class 1.
!
!     ------------------------------------------------------------
!
!     Summary:
!     ========
!     Damped Newton-algorithm for systems of highly nonlinear
!     equations - damping strategy due to Ref. (1).
!
!     (The iteration is done by subroutine N1INT actually. NLEQ1E
!      calls the standard interface driver NLEQ1, which itself does
!      some house keeping and builds up workspace.)
!
!     Jacobian approximation by numerical differences.
!
!     The numerical solution of the arising linear equations is
!     done by means of the subroutines *GEFA and *GESL ( Gauss-
!     algorithm with column-pivoting and row-interchange;
!     replace '*' by 'S' or 'D' for single or double precision
!     version respectively).
!
!     ------------------------------------------------------------
!
!     External subroutine to be supplied by the user
!     ==============================================
!
!     FCN(N,X,F,IFAIL) Ext    Function subroutine - the problem function
!                             (The routine must be named exactly FCN !)
!       N              Int    Number of vector components (input)
!                             Must not be altered!
!       X(N)           Dble   Vector of unknowns (input)
!                             Must not be altered!
!       F(N)           Dble   Vector of function values (output)
!       IFAIL          Int    FCN evaluation-failure indicator. (output)
!                             On input:  Has always value 0 (zero).
!                             On output: Indicates failure of FCN eval-
!                                uation, if having a nonzero value.
!                             If <0: NLEQ1E will be terminated with
!                                    IFAIL returned via IERR.
!                             If =1: A new trial Newton iterate will be
!                                    computed, with the damping factor
!                                    reduced to it's half.
!                             If =2: A new trial Newton iterate will
!                                    computed, with the damping factor
!                                    reduced by a reduction factor,
!                                    which must be output through F(1)
!                                    by FCN, and it's value must be
!                                    >0 and < 1.
!
!     Parameters list description
!     ===========================
!
!     Name   Type    In/Out Meaning
!
!     N      Int     In     Number of unknowns ( N .LE. 50 )
!     X(N)   Dble    In     Initial estimate of the solution
!             Out    Solution values ( or final values,
!                           respectively )
!     RTOL   Dble    In     Required relative precision of
!                           solution components -
!                           RTOL.GE.EPMACH*TEN*N
!             Out    Finally achieved accuracy
!     IERR   Int     Out    Return value parameter
!                           < 0 Termination forced by user function FCN
!                               due to IFAIL<0 on output, IERR is set
!                               to IFAIL
!                           = 0 successfull completion of the iteration,
!                               solution has been computed
!                           > 0 see list of error/warning messages below
!
!    Error and warning messages:
!    ===========================
!
!      1    Termination, since Jacobian matrix became singular
!      2    Termination after 100 iterations
!      3    Termination, since damping factor became to small
!      4    Warning: Superlinear or quadratic convergence slowed down
!           near the solution.
!           Iteration has been stopped therefore with an approximation
!           of the solution not such accurate as requested by RTOL,
!           because possibly the RTOL requirement may be too stringent
!           (i.e. the nonlinear problem is ill-conditioned)
!      5    Warning: Iteration stopped with termination criterion
!           (using RTOL as requested precision) satisfied, but no
!           superlinear or quadratic convergence has been indicated yet.
!           Therefore, possibly the error estimate for the solution may
!           not match good enough the really achieved accuracy.
!     20    Bad input value to parameter N; 1 .LE. N .LE. 50 required
!     21    Nonpositive value for RTOL supplied
!     82    Termination, because user routine FCN returned with IFAIL>0
!
!     Note   : in case of failure:
!        -    use better initial guess
!        -    or refine model
!        -    or use non-standard options and/or analytical Jacobian
!             via the standard interface NLEQ1
!
!     Machine dependent constants used:
!     =================================
!
!     DOUBLE PRECISION EPMACH  in  N1PCHK, N1INT
!     DOUBLE PRECISION GREAT   in  N1PCHK
!     DOUBLE PRECISION SMALL   in  N1PCHK, N1INT, N1SCAL
!
!     Subroutines called: NLEQ1
!
!     ------------------------------------------------------------
!     End Prologue
      SUBROUTINE NLEQ1E(FCN,DUMMY,N,X,RTOL,IERR)
!      use FEmapper
!     Begin Prologue NLEQ1E
      INTEGER N
      DOUBLE PRECISION X(N), RTOL
      INTEGER IERR

      INTEGER NMAX, LIOPT, LIWK, LRWK, LUPRT
      PARAMETER ( NMAX=50 )	! Gib was 50
      PARAMETER ( LIOPT=50, LIWK=NMAX+50, LRWK=(NMAX+13)*NMAX+60 )
      PARAMETER ( LUPRT=6 )
!      EXTERNAL NLEQ1	! FCN
      procedure(func) :: FCN
      DOUBLE PRECISION XSCAL(NMAX)
      INTEGER IOPT(LIOPT), IWK(LIWK)
      DOUBLE PRECISION RWK(LRWK)
      INTEGER I, NIW, NRW
      CHARACTER CHGDAT*20, PRODCT*8
      procedure(jacfunc) :: DUMMY
! added this
      logical :: suppress_error_messages
!
!     Version: 2.3               Latest change:
!     -----------------------------------------
!
      DATA      CHGDAT      /'November 15, 1991   '/
      DATA      PRODCT      /'NLEQ1E  '/

! added this
      suppress_error_messages = .false.

!     Begin
      NIW = N+50
      NRW = (N+13)*N+60
!     Checking dimensional parameter N
      IF ( N.LT.1 .OR. N.GT.NMAX ) THEN
        WRITE(LUPRT,1001) NMAX, N
1001    FORMAT(/,' Error: Bad input to parameter N supplied', &
              /,8X,'choose 1 .LE. N .LE. ',I3, &
                    ' , your input is: N = ',I5)
        IERR = 20
        RETURN
      ENDIF
      DO 10 I=1,LIOPT
        IOPT(I) = 0
10    CONTINUE
      DO 20 I=1,NIW
        IWK(I) = 0
20    CONTINUE
      DO 30 I=1,NRW
        RWK(I) = 0.0D0
30    CONTINUE
      DO 40 I=1,N
        XSCAL(I) = 0.0D0
40    CONTINUE
!     Print errors, warnings, monitor and time monitor
!     to standard output
      IOPT(11) = 1	!3
      if (suppress_error_messages) then
		  IOPT(11) = 0
	endif
      IOPT(12) = LUPRT
      IOPT(13) = 0	!3
      IOPT(14) = LUPRT
      IOPT(19) = 1
      IOPT(20) = LUPRT
      IOPT(31) = 3	! degree of nonlinearity 1..4
!     Maximum number of Newton iterations
      IWK(31) = 100

      CALL NLEQ1(N,FCN,DUMMY,X,XSCAL,RTOL,IOPT,IERR,LIWK,IWK,LRWK,RWK)
      IF (IERR.EQ.82 .AND. IWK(23).LT.0) IERR=IWK(23)
!     End of subroutine NLEQ1E
      RETURN
      END subroutine


      SUBROUTINE NLEQ1(N,FCN,JAC,X,XSCAL,RTOL,IOPT,IERR,LIWK,IWK,LRWK,RWK)
!     Begin Prologue NLEQ1
      INTEGER N
!      EXTERNAL FCN,JAC
      procedure(func) :: FCN
      procedure(jacfunc) :: JAC
      DOUBLE PRECISION X(N),XSCAL(N)
      DOUBLE PRECISION RTOL
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
!     ------------------------------------------------------------
!
!   Title
!
!     Numerical solution of nonlinear (NL) equations (EQ)
!     especially designed for numerically sensitive problems.
!
!   Written by        U. Nowak, L. Weimann
!   Purpose           Solution of systems of highly nonlinear equations
!   Method            Damped affine invariant Newton method
!              (see references below)
!   Category          F2a. - Systems of nonlinear equations
!   Keywords          Nonlinear equations, Newton methods
!   Version           2.3
!   Revision          September 1991
!   Latest Change     July 2000
!   Library           CodeLib
!   Code              Fortran 77, Double Precision
!   Environment       Standard Fortran 77 environment on PC's,
!              workstations and hosts.
!   Copyright     (c) Konrad-Zuse-Zentrum fuer
!              Informationstechnik Berlin (ZIB)
!              Takustrasse 7, D-14195 Berlin-Dahlem
!              phone : + 49/30/84185-0
!              fax   : + 49/30/84185-125
!   Contact           Lutz Weimann
!              ZIB, Division Scientific Computing,
!                   Department Scientific Software
!              phone : + 49/30/84185-185
!              fax   : + 49/30/84185-107
!              e-mail: weimann@zib.de
!
!     References:
!
!     /1/ P. Deuflhard:
!         Newton Methods for Nonlinear Problems. -
!         Affine Invariance and Adaptive Algorithms.
!         Series Computational Mathematics 35, Springer (2004)
!
!     /2/ U. Nowak, L. Weimann:
!         A Family of Newton Codes for Systems of Highly Nonlinear
!         Equations - Algorithm, Implementation, Application.
!         ZIB, Technical Report TR 90-10 (December 1990)
!
!  ---------------------------------------------------------------
!
!* Licence
!    You may use or modify this code for your own non commercial
!    purposes for an unlimited time.
!    In any case you should not deliver this code without a special
!    permission of ZIB.
!    In case you intend to use the code commercially, we oblige you
!    to sign an according licence agreement with ZIB.
!
!* Warranty
!    This code has been tested up to a certain level. Defects and
!    weaknesses, which may be included in the code, do not establish
!    any warranties by ZIB. ZIB does not take over any liabilities
!    which may follow from acquisition or application of this code.
!
!* Software status
!    This code is under care of ZIB and belongs to ZIB software class 1.
!
!     ------------------------------------------------------------
!
!     Summary:
!     ========
!     Damped Newton-algorithm for systems of highly nonlinear
!     equations - damping strategy due to Ref. (1).
!
!     (The iteration is done by subroutine N1INT currently. NLEQ1
!      itself does some house keeping and builds up workspace.)
!
!     Jacobian approximation by numerical differences or user
!     supplied subroutine JAC.
!
!     The numerical solution of the arising linear equations is
!     done by means of the subroutines *GEFA and *GESL ( Gauss-
!     algorithm with column-pivoting and row-interchange ) in the
!     dense matrix case, or by the subroutines *GBFA and *GBSL in
!     the band matrix case from LINPACK (replace '*' by 'S' or 'D'
!     for single or double precision version respectively).
!     For special purposes these routines may be substituted.
!
!     This is a driver routine for the core solver N1INT.
!
!     ------------------------------------------------------------
!
!     Parameters list description (* marks inout parameters)
!     ======================================================
!
!     External subroutines (to be supplied by the user)
!     =================================================
!
!     (Caution: Arguments declared as (input) must not
!        be altered by the user subroutines ! )
!
!     FCN(N,X,F,IFAIL) Ext    Function subroutine
!       N              Int    Number of vector components (input)
!       X(N)           Dble   Vector of unknowns (input)
!       F(N)           Dble   Vector of function values (output)
!       IFAIL          Int    FCN evaluation-failure indicator. (output)
!                             On input:  Has always value 0 (zero).
!                             On output: Indicates failure of FCN eval-
!                                uation, if having a value <= 2.
!                             If <0: NLEQ1 will be terminated with
!                                    error code = 82, and IFAIL stored
!                                    to IWK(23).
!                             If =1: A new trial Newton iterate will
!                                    computed, with the damping factor
!                                    reduced to it's half.
!                             If =2: A new trial Newton iterate will
!                                    computed, with the damping factor
!                                    reduced by a reduct. factor, which
!                                    must be output through F(1) by FCN,
!                                    and it's value must be >0 and < 1.
!                             Note, that if IFAIL = 1 or 2, additional
!                             conditions concerning the damping factor,
!                             e.g. the minimum damping factor or the
!                             bounded damping strategy may also influ-
!                             ence the value of the reduced damping
!                             factor.
!
!     JAC(N,LDJAC,X,DFDX,IFAIL)
!                Ext    Jacobian matrix subroutine
!       N                 Int    Number of vector components (input)
!       LDJA!             Int    Leading dimension of Jacobian array
!                                (input)
!       X(N)              Dble   Vector of unknowns (input)
!       DFDX(LDJAC,N)     Dble   DFDX(i,k): partial derivative of
!                                I-th component of FCN with respect
!                                to X(k) (output)
!       IFAIL             Int    JAC evaluation-failure indicator.
!                                (output)
!                                Has always value 0 (zero) on input.
!                                Indicates failure of JAC evaluation
!                                and causes termination of NLEQ1,
!                                if set to a negative value on output
!
!
!     Input parameters of NLEQ1
!     =========================
!
!     N              Int    Number of unknowns
!    * X(N)           Dble   Initial estimate of the solution
!    * XSCAL(N)       Dble   User scaling (lower threshold) of the
!                           iteration vector X(N)
!    * RTOL           Dble   Required relative precision of
!                           solution components -
!                           RTOL.GE.EPMACH*TEN*N
!    * IOPT(50)       Int    Array of run-time options. Set to zero
!                           to get default values (details see below)
!
!     Output parameters of NLEQ1
!     ==========================
!
!    * X(N)           Dble   Solution values ( or final values,
!                           respectively )
!    * XSCAL(N)       Dble   After return with IERR.GE.0, it contains
!                           the latest internal scaling vector used
!                           After return with IERR.EQ.-1 in onestep-
!                           mode it contains a possibly adapted
!                           (as described below) user scaling vector:
!                           If (XSCAL(I).LT. SMALL) XSCAL(I) = SMALL ,
!                           If (XSCAL(I).GT. GREAT) XSCAL(I) = GREAT .
!                           For SMALL and GREAT, see section machine
!                           constants below  and regard note 1.
!    * RTOL           Dble   Finally achieved (relative) accuracy.
!                           The estimated absolute error of component i
!                           of x_out is approximately given by
!                             abs_err(i) = RTOL * XSCAL_out(i) ,
!                           where (approximately)
!                             XSCAL_out(i) =
!                                max(abs(X_out(i)),XSCAL_in(i)).
!     IERR           Int    Return value parameter
!                           =-1 sucessfull completion of one iteration
!                               step, subsequent iterations are needed
!                               to get a solution. (stepwise mode only)
!                           = 0 successfull completion of the iteration,
!                               solution has been computed
!                           > 0 see list of error messages below
!
!     Note 1.
!        The machine dependent values SMALL and EPMACH are
!        gained from calls of the ZIBCONST.
!
!     Workspace parameters of NLEQ1
!     =============================
!
!     LIWK           Int    Declared dimension of integer workspace.
!                           Required minimum (for standard linear system
!                           solver) : N+50
!  *  IWK(LIWK)      Int    Integer Workspace
!     LRWK           Int    Declared dimension of real workspace.
!                           Required minimum (for standard linear system
!                           solver and Jacobian computed by numerical
!                           approximation - if the Jacobian is computed
!                           by a user subroutine JAC, decrease the
!                           expressions noted below by N):
!                           for full case Jacobian: (N+NBROY+13)*N+61
!                           for a band-matrix Jacobian:
!                              (2*ML+MU+NBROY+14)*N+61 with
!                           ML = lower bandwidth , MU = upper bandwidth
!                           NBROY = Maximum number of Broyden steps
!                           (Default: if Broyden steps are enabled, e.g.
!                                                IOPT(32)=1            -
!                                       NBROY=N (full Jacobian),
!                                            =ML+MU+1 (band Jacobian),
!                                       but at least NBROY=10
!                                     else (if IOPT(32)=0) -
!                                       NBROY=0 ;
!                            see equally named IOPT and IWK-fields below)
!    * RWK(LRWK)      Dble   Real Workspace
!
!     Note 2a.  A test on sufficient workspace is made. If this
!        test fails, IERR is set to 10 and an error-message
!        is issued from which the minimum of required
!        workspace size can be obtained.
!
!     Note 2b.  The first 50 elements of IWK and RWK are partially
!        used as input for internal algorithm parameters (for
!        details, see below). In order to set the default values
!        of these parameters, the fields must be set to zero.
!        Therefore, it's recommanded always to initialize the
!        first 50 elements of both workspaces to zero.
!
!    Options IOPT:
!    =============
!
!     Pos. Name   Default  Meaning
!
!       1  QSUCC  0        =0 (.FALSE.) initial call:
!                             NLEQ1 is not yet initialized, i.e. this is
!                             the first call for this nonlinear system.
!                             At successfull return with MODE=1,
!                             QSUCC is set to 1.
!                   =1 (.TRUE.) successive call:
!                             NLEQ1 is initialized already and is now
!                             called to perform one or more following
!                             Newton-iteration steps.
!                             ATTENTION:
!                                Don't destroy the contents of
!                                IOPT(i) for 1 <= i <= 50 ,
!                                IWK(j)  for 1 <= j < NIWKFR and
!                                RWK(k)  for 1 <= k < NRWKFR.
!                                (Nevertheless, some of the options, e.g.
!                                 FCMIN, SIGMA, MPR..., can be modified
!                                 before successive calls.)
!       2  MODE   0        =0 Standard mode initial call:
!                             Return when the required accuracy for the
!                             iteration vector is reached. User defined
!                             parameters are evaluated and checked.
!                             Standard mode successive call:
!                             If NLEQ1 was called previously with MODE=1,
!                             it performs all remaining iteration steps.
!                   =1 Stepwise mode:
!                             Return after one Newton iteration step.
!       3  JACGEN 0        Method of Jacobian generation
!                   =0 Standard method is JACGEN=2
!                   =1 User supplied subroutine JAC will be
!                             called to generate Jacobian matrix
!                   =2 Jacobian approximation by numerical
!                             differentation (see subroutines N1JAC
!                             and N1JACB)
!                   =3 Jacobian approximation by numerical
!                             differentation with feedback control
!                             (see subroutines N1JCF and N1JCFB)
!       4  MSTOR  0        =0 The Jacobian A is a dense matrix
!                   =1 A is a band matrix
!       5                  Reserved
!       6  ML     0        Lower bandwidth of A (excluding the
!                   diagonal);
!                   IOPT(6) ignored, if IOPT(4).NE.1
!       7  MU     0        Upper bandwidth of A (excluding the
!                   diagonal);
!                   IOPT(7) ignored, if IOPT(4).NE.1
!       8                  Reserved
!       9  ISCAL  0        Determines how to scale the iterate-vector:
!                   =0 The user supplied scaling vector XSCAL is
!                             used as a (componentwise) lower threshold
!                             of the current scaling vector
!                   =1 The vector XSCAL is always used as the
!                             current scaling vector
!      10                  Reserved
!      11  MPRERR 0        Print error messages
!                   =0 No output
!                   =1 Error messages
!                   =2 Warnings additionally
!                   =3 Informal messages additionally
!      12  LUERR  6        Logical unit number for error messages
!      13  MPRMON 0        Print iteration Monitor
!                   =0 No output
!                   =1 Standard output
!                   =2 Summary iteration monitor additionally
!                   =3 Detailed iteration monitor additionally
!                   =4,5,6 Outputs with increasing level addi-
!                             tional increasing information for code
!                             testing purposes. Level 6 produces
!                             in general extremely large output!
!      14  LUMON  6        Logical unit number for iteration monitor
!      15  MPRSOL 0        Print solutions
!                   =0 No output
!                   =1 Initial values and solution values
!                   =2 Intermediate iterates additionally
!      16  LUSOL  6        Logical unit number for solutions
!      17..18              Reserved
!      19  MPRTIM 0        Output level for the time monitor
!                   = 0 : no time measurement and no output
!                   = 1 : time measurement will be done and
!                                summary output will be written -
!                                regard note 5a.
!      20  LUTIM  6        Logical output unit for time monitor
!      21..30              Reserved
!      31  NONLIN 3        Problem type specification
!                   =1 Linear problem
!                             Warning: If specified, no check will be
!                             done, if the problem is really linear, and
!                             NLEQ1 terminates unconditionally after one
!                             Newton-iteration step.
!                   =2 Mildly nonlinear problem
!                   =3 Highly nonlinear problem
!                   =4 Extremely nonlinear problem
!      32  QRANK1 0        =0 (.FALSE.) Rank-1 updates by Broyden-
!                             approximation are inhibited.
!                   =1 (.TRUE.) Rank-1 updates by Broyden-
!                             approximation are allowed.
!      33  QORDI  0        =0 (.FALSE.) Standard program mode
!                   =1 (.TRUE.)  Special program mode:
!                             Ordinary Newton iteration is done, e.g.:
!                             No damping strategy and no monotonicity
!                             test is applied
!      34  QSIMPL 0        =0 (.FALSE.) Standard program mode
!                   =1 (.TRUE.)  Special program mode:
!                             Simplified Newton iteration is done, e.g.:
!                             The Jacobian computed at the starting
!                             point is fixed for all subsequent
!                             iteration steps, and
!                             no damping strategy and no monotonicity
!                             test is applied.
!      35  QNSCAL 0        Inhibit automatic row scaling:
!                   =0 (.FALSE.) Automatic row scaling of
!                             the linear system is activ:
!                             Rows i=1,...,N will be divided by
!                             max j=1,...,N (abs(a(i,j)))
!                   =1 (.TRUE.) No row scaling of the linear
!                             system. Recommended only for well row-
!                             scaled nonlinear systems.
!      36..37              Reserved
!      38  IBDAMP          Bounded damping strategy switch:
!                   =0 The default switch takes place, dependent
!                             on the setting of NONLIN (=IOPT(31)):
!                             NONLIN = 0,1,2,3 -> IBDAMP = off ,
!                             NONLIN = 4 -> IBDAMP = on
!                   =1 means always IBDAMP = on
!                   =2 means always IBDAMP = off
!      39  IORMON          Convergence order monitor
!                   =0 Standard option is IORMON=2
!                   =1 Convergence order is not checked,
!                             the iteration will be always proceeded
!                             until the solution has the required
!                             precision RTOL (or some error condition
!                             occured)
!                   =2 Use additional 'weak stop' criterion:
!                             Convergence order is monitored
!                             and termination due to slowdown of the
!                             convergence may occur.
!                   =3 Use additional 'hard stop' criterion:
!                             Convergence order is monitored
!                             and termination due to superlinear
!                             convergence slowdown may occur.
!                   In case of termination due to convergence
!                   slowdown, the warning code IERR=4 will be
!                   set.
!                   In cases, where the Newton iteration con-
!                   verges but superlinear convergence order has
!                   never been detected, the warning code IERR=5
!                   is returned.
!      40..45              Reserved
!      46..50              User options (see note 5b)
!
!     Note 3:
!         If NLEQ1 terminates with IERR=2 (maximum iterations)
!         or  IERR=3 (small damping factor), you may try to continue
!         the iteration by increasing NITMAX or decreasing FCMIN
!         (see RWK) and setting QSUCC to 1.
!
!     Note 4 : Storage of user supplied banded Jacobian
!        In the band matrix case, the following lines may build
!        up the analytic Jacobian A;
!        Here AFL denotes the quadratic matrix A in dense form,
!        and ABD the rectangular matrix A in banded form :
!
!            ML = IOPT(6)
!            MU = IOPT(7)
!            MH = MU+1
!            DO 20 J = 1,N
!              I1 = MAX0(1,J-MU)
!              I2 = MIN0(N,J+ML)
!              DO 10 I = I1,I2
!                K = I-J+MH
!                ABD(K,J) = AFL(I,J)
!           10        CONTINUE
!           20      CONTINUE
!
!         The total number of rows needed in  ABD  is  ML+MU+1 .
!         The  MU by MU  upper left triangle and the
!         ML by ML  lower right triangle are not referenced.
!
!     Note 5a:
!        The integrated time monitor calls the machine dependent
!        subroutine SECOND to get the current time stamp in form
!        of a real number (Single precision). As delivered, this
!        subroutine always return 0.0 as time stamp value. Refer
!        to the compiler- or library manual of the FORTRAN compiler
!        which you currently use to find out how to get the current
!        time stamp on your machine.
!
!     Note 5b:
!         The user options may be interpreted by the user replacable
!         routines N1SOUT, N1FACT, N1SOLV - the distributed version
!         of N1SOUT currently uses IOPT(46) as follows:
!         0 = standard plotdata output (may be postprocessed by a user-
!             written graphical program)
!         1 = plotdata output is suitable as input to the graphical
!             package GRAZIL (based on GKS), which has been developed
!             at ZIB.
!
!
!    Optional INTEGER input/output in IWK:
!    =======================================
!
!     Pos. Name          Meaning
!
!      1   NITER  IN/OUT Number of Newton-iterations
!      2                 reserved
!      3   NCORR  IN/OUT Number of corrector steps
!      4   NFCN   IN/OUT Number of FCN-evaluations
!      5   NJA!    IN/OUT Number of Jacobian generations or
!                 JAC-calls
!      6                 reserved
!      7                 reserved
!      8   NFCNJ  IN/OUT Number of FCN-evaluations for Jacobian
!                 approximation
!      9   NREJR1 IN/OUT Number of rejected Newton iteration steps
!                 done with a rank-1 approximated Jacobian
!     10..11             Reserved
!     12   IDCODE IN/OUT Output: The 8 decimal digits program identi-
!                 fication number ppppvvvv, consisting of the
!                 program code pppp and the version code vvvv.
!                 Input: If containing a negative number,
!                 it will only be overwritten by the identi-
!                 fication number, immediately followed by
!                 a return to the calling program.
!     13..15             Reserved
!     16   NIWKFR OUT    First element of IWK which is free to be used
!                 as workspace between Newton iteration steps
!                 for standard linear solvers: 51
!     17   NRWKFR OUT    First element of RWK which is free to be used
!                 as workspace between Newton iteration steps.
!                 For standard linear solvers and numerically
!                 approximated Jacobian computed by one of the
!                 expressions:
!                 (N+7+NBROY)*N+61        for a full Jacobian
!                 (2*ML+MU+8+NBROY)*N+61  for a banded Jacobian
!                 If the Jacobian is computed by a user routine
!                 JAC, subtract N in both expressions.
!     18   LIWKA  OUT    Length of IWK currently required
!     19   LRWKA  OUT    Length of RWK currently required
!     20..22             Reserved
!     23   IFAIL  OUT    Set in case of failure of N1FACT (IERR=80),
!                 N1SOLV (IERR=81), FCN (IERR=82) or JAC(IERR=83)
!                 to the nonzero IFAIL value returned by the
!                 routine indicating the failure .
!     24   ICONV  OUT    Current status of of the convergence monitor
!                 (only if convergence order monitor is on -
!                  see IORMON(=IOPT(39)))
!                 =0: No convergence indicated yet
!                 =1: Damping factor is 1.0d0
!                 =2: Superlinear convergence in progress
!                 =3: Quadratic convergence in progress
!     25..30             Reserved
!     31   NITMAX IN     Maximum number of permitted iteration
!                 steps (default: 50)
!     32                 Reserved
!     33   NEW    IN/OUT Count of consecutive rank-1 updates
!     34..35             Reserved
!     36   NBROY  IN     Maximum number of possible consecutive
!                 iterative Broyden steps. The total real
!                 workspace needed (RWK) depends on this value
!                 (see LRWK above).
!                 Default is N (see parameter N),
!                 if MSTOR=0 (=IOPT(4)),
!                 and ML+MU+1 (=IOPT(6)+IOPT(7)+1), if MSTOR=1
!                 (but minimum is always 10) -
!                 provided that Broyden is allowed.
!                 If Broyden is inhibited, NBROY is always set to
!                 zero.
!     37..50             Reserved
!
!    Optional REAL input/output in RWK:
!    ====================================
!
!     Pos. Name          Meaning
!
!      1..16             Reserved
!     17   CONV   OUT    The achieved relative accuracy after the
!                 current step
!     18   SUMX   OUT    Natural level (not Normx of printouts)
!                 of the current iterate, i.e. Sum(DX(i)**2),
!                 where DX = scaled Newton correction.
!     19   DLEVF  OUT    Standard level (not Normf of printouts)
!                 of the current iterate, i.e. Norm2(F(X)),
!                 where F =  nonlinear problem function.
!     20   FCBND  IN     Bounded damping strategy restriction factor
!                 (Default is 10)
!     21   FCSTRT IN     Damping factor for first Newton iteration -
!                 overrides option NONLIN, if set (see note 6)
!     22   FCMIN  IN     Minimal allowed damping factor (see note 6)
!     23   SIGMA  IN     Broyden-approximation decision parameter
!                 Required choice: SIGMA.GE.1. Increasing this
!                 parameter make it less probable that the algo-
!                 rith performs rank-1 updates.
!                 Rank-1 updates are inhibited, if
!                 SIGMA.GT.1/FCMIN is set. (see note 6)
!     24   SIGMA2 IN     Decision parameter about increasing damping
!                 factor to corrector if predictor is small.
!                 Required choice: SIGMA2.GE.1. Increasing this
!                 parameter make it less probable that the algo-
!                 rith performs rank-1 updates.
!     25                 Reserved
!     26   AJDEL  IN     Jacobian approximation without feedback:
!                 Relative pertubation for components
!                 (Default: sqrt(epmach*10), epmach: relative
!                  machine precision)
!     27   AJMIN  IN     Jacobian approximation without feedback:
!                 Threshold value (Default: 0.0d0)
!                   The absolute pertubation for component k is
!                   computed by
!                   DELX := AJDEL*max(abs(Xk),AJMIN)
!     28  ETADIF  IN     Jacobian approximation with feedback:
!                 Target value for relative pertubation ETA of X
!                 (Default: 1.0d-6)
!     29  ETAINI  IN     Jacobian approximation with feedback:
!                 Initial value for denominator differences
!                 (Default: 1.0d-6)
!     30..50             Reserved
!
!     Note 6:
!       The default values of the internal parameters may be obtained
!       from the monitor output with at least IOPT field MPRMON set to 2
!       and by initializing the corresponding RWK-fields to zero.
!
!    Error and warning messages:
!    ===========================
!
!      1    Termination, since jacobian matrix became singular
!      2    Termination after NITMAX iterations ( as indicated by
!           input parameter NITMAX=IWK(31) )
!      3    Termination, since damping factor became to small
!      4    Warning: Superlinear or quadratic convergence slowed down
!           near the solution.
!           Iteration has been stopped therefore with an approximation
!           of the solution not such accurate as requested by RTOL,
!           because possibly the RTOL requirement may be too stringent
!           (i.e. the nonlinear problem is ill-conditioned)
!      5    Warning: Iteration stopped with termination criterion
!           (using RTOL as requested precision) satisfied, but no
!           superlinear or quadratic convergence has been indicated yet.
!           Therefore, possibly the error estimate for the solution may
!           not match good enough the really achieved accuracy.
!     10    Integer or real workspace too small
!     20    Bad input to dimensional parameter N
!     21    Nonpositive value for RTOL supplied
!     22    Negative scaling value via vector XSCAL supplied
!     30    One or more fields specified in IOPT are invalid
!           (for more information, see error-printout)
!     80    Error signalled by linear solver routine N1FACT,
!           for more detailed information see IFAIL-value
!           stored to IWK(23)
!     81    Error signalled by linear solver routine N1SOLV,
!           for more detailed information see IFAIL-value
!           stored to IWK(23)
!           (not used by standard routine N1SOLV)
!     82    Error signalled by user routine FCN (Nonzero value
!           returned via IFAIL-flag; stored to IWK(23) )
!     83    Error signalled by user routine JAC (Nonzero value
!           returned via IFAIL-flag; stored to IWK(23) )
!
!     Note 7 : in case of failure:
!        -    use non-standard options
!        -    or turn to Newton-algorithm with rank strategy
!        -    use another initial guess
!        -    or reformulate model
!        -    or apply continuation techniques
!
!     Machine dependent constants used:
!     =================================
!
!     DOUBLE PRECISION EPMACH  in  N1PCHK, N1INT
!     DOUBLE PRECISION GREAT   in  N1PCHK
!     DOUBLE PRECISION SMALL   in  N1PCHK, N1INT, N1SCAL
!
!     Subroutines called: N1PCHK, N1INT
!
!     ------------------------------------------------------------
!     End Prologue
!
!     Summary of changes:
!     ===================
!
!     2.2.1  91, May 30     Time monitor included
!     2.2.2  91, May 30     Bounded damping strategy implemented
!     2.2.3  91, June 19    AJDEL, AJMIN as RWK-options for JACGEN.EQ.2,
!                           ETADIF, ETAINI as RWK-opts. for JACGEN.EQ.3
!                           FCN-count changed for anal. Jacobian
!     2.2.4  91, August  9  Convergence order monitor included
!     2.2.5  91, August 13  Standard Broyden updates replaced by
!                           iterative Broyden
!     2.2.6  91, Sept.  16  Damping factor reduction by FCN-fail imple-
!                           mented
!     2.3    91, Dec.   20  New Release for CodeLib
!            92, March  11  Level of accepted simplified correction
!                           stored to RWK(IRWKI+4)
!            00, July   12  RTOL output-value bug fixed
!
!     ------------------------------------------------------------
!
!     PARAMETER (IRWKI=xx, LRWKI=yy)
!     IRWKI: Start position of internally used RWK part
!     LRWKI: Length of internally used RWK part
!     (current values see parameter statement below)
!
!     INTEGER L4,L5,L51,L6,L61,L62,L63,L7,L71,L8,L9,L10,L11,L12,L121,
!             L13,L14,L20
!     Starting positions in RWK of formal array parameters of internal
!     routine N1INT (dynamically determined in driver routine NLEQ1,
!     dependent on N and options setting)
!
!     Further RWK positions (only internally used)
!
!     Position  Name     Meaning
!
!     IRWKI     FCKEEP   Damping factor of previous successfull iter.
!     IRWKI+1   FCA      Previous damping factor
!     IRWKI+2   FCPRI    A priori estimate of damping factor
!     IRWKI+3   DMYCOR   Number My of latest corrector damping factor
!                 (kept for use in rank-1 decision criterium)
!     IRWKI+4   SUMXS    natural level of accepted simplified correction
!     IRWKI+(5..LRWKI-1) Free
!
!     Internal arrays stored in RWK (see routine N1INT for descriptions)
!
!     Position  Array         Type   Remarks
!
!     L4        A(M1,N)       Perm   M1=N (full mode) or
!                                    M1=2*IOPT(6)+IOPT(7)+1 (band mode)
!     L41       DXSAVE(N,NBROY)
!                             Perm   NBROY=IWK(36) (Default: N or 0)
!     L5        DX(N)         Perm
!     L51       DXQ(N)        Perm
!     L6        XA(N)         Perm
!     L61       F(N)          Perm
!     L62       FW(N)         Perm
!     L63       XWA(N)        Perm
!     L7        FA(N)         Perm
!     L71       ETA(N)        Perm   Only used for JACGEN=IOPT(3)=3
!     L9        XW(N)         Temp
!     L11       DXQA(N)       Temp
!     L12       T1(N)         Temp
!     L121      T2(N)         Temp
!     L13       T3(N)         Temp
!     L14                     Temp   Start position of array workspace
!                                    needed for linear solver
!
!
!      EXTERNAL N1INT
      INTRINSIC DBLE
      INTEGER IRWKI, LRWKI
      PARAMETER (IRWKI=51, LRWKI=10)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER NITMAX,LUERR,LUMON,LUSOL,MPRERR,MPRMON,MPRSOL, &
      MSTOR,M1,M2,NRWKFR,NRFRIN,NRW,NIWKFR,NIFRIN,NIW,NONLIN,JACGEN
      INTEGER L4,L41,L5,L51,L6,L61,L62,L63,L7,L71,L8,L9,L11,L12,L121, &
      L13,L14,L20
      DOUBLE PRECISION FC,FCMIN,PERCI,PERCR
      LOGICAL QINIMO,QRANK1,QFCSTR,QSUCC,QBDAMP,QSIMPL
      CHARACTER CHGDAT*20, PRODCT*8
!     Which version ?
      LOGICAL QVCHK
      INTEGER IVER, MPRTIM, LUTIM, ML, MU, NBROY
      PARAMETER( IVER=21112301 )
!
!     Version: 2.3               Latest change:
!     -----------------------------------------
!
      DATA      CHGDAT      /'July 12, 2000       '/
      DATA      PRODCT      /'NLEQ1   '/
!     Begin
      IERR = 0
      QVCHK = IWK(12).LT.0
      IWK(12) = IVER
      IF (QVCHK) RETURN
!        Print error messages?
      MPRERR = IOPT(11)
      LUERR = IOPT(12)
      IF (LUERR .EQ. 0) THEN
        LUERR = 6
        IOPT(12)=LUERR
      ENDIF
!        Print iteration monitor?
      MPRMON = IOPT(13)
      LUMON = IOPT(14)
      IF (LUMON .LE. 0 .OR. LUMON .GT. 99) THEN
        LUMON = 6
        IOPT(14)=LUMON
      ENDIF
!        Print intermediate solutions?
      MPRSOL = IOPT(15)
      LUSOL = IOPT(16)
      IF (LUSOL .EQ. 0) THEN
        LUSOL = 6
        IOPT(16)=LUSOL
      ENDIF
!        Print time summary statistics?
      MPRTIM = IOPT(19)
      LUTIM = IOPT(20)
      IF (LUTIM .EQ. 0) THEN
        LUTIM = 6
        IOPT(20)=LUTIM
      ENDIF
      QSUCC = IOPT(1).EQ.1
      QINIMO = MPRMON.GE.1.AND..NOT.QSUCC
!     Print NLEQ1 heading lines
      IF(QINIMO)THEN
10000   FORMAT('   N L E Q 1  *****  V e r s i o n  ', &
               '2 . 3 ***',//,1X,'Newton-Method ', &
               'for the solution of nonlinear systems',//)
        WRITE(LUMON,10000)
      ENDIF
!     Check input parameters and options
      CALL N1PCHK(N,X,XSCAL,RTOL,IOPT,IERR,LIWK,IWK,LRWK,RWK)
!     Exit, if any parameter error was detected till here
      IF (IERR.NE.0) RETURN
!
      MSTOR=IOPT(4)
      IF (MSTOR.EQ.0) THEN
        M1=N
        M2=N
      ELSE IF (MSTOR.EQ.1) THEN
        ML=IOPT(6)
        MU=IOPT(7)
        M1=2*ML+MU+1
        M2=ML+MU+1
      ENDIF
      JACGEN=IOPT(3)
      IF (JACGEN.EQ.0) JACGEN=2
      IOPT(3)=JACGEN
      QRANK1=IOPT(32).EQ.1
      QSIMPL=IOPT(34).EQ.1
      IF (QRANK1) THEN
        NBROY=IWK(36)
        IF (NBROY.EQ.0) NBROY=MAX(M2,10)
        IWK(36)=NBROY
      ELSE
        NBROY=0
      ENDIF
!     WorkSpace: RWK
      L4=IRWKI+LRWKI
      L41=L4+M1*N
      L5=L41+NBROY*N
      L51=L5+N
      L6=L51+N
      L61=L6+N
      L62=L61+N
      L63=L62+N
      L7=L63+N
      L71=L7+N
      IF (JACGEN.NE.3) THEN
        L8=L71
      ELSE
        L8=L71+N
      ENDIF
      NRWKFR = L8
      L9=L8
      L11=L9+N
      L12=L11+N
      L121=L12+N
      L13=L121+N
      L14=L13+N
      NRW=L14-1
!     End WorkSpace at NRW
!     WorkSpace: IWK
      L20=51
      NIWKFR = L20
      IF (QRANK1.OR.QSIMPL) NIWKFR = NIWKFR+N
      NIW=L20-1
!     End WorkSpace at NIW
      IWK(16) = NIW+1
      IWK(17) = NRW+1
      NIFRIN = NIW+1
      NRFRIN = NRW+1

      IF(NRW.GT.LRWK.OR.NIW.GT.LIWK)THEN
        IERR=10
      ELSE
        IF(QINIMO)THEN
          PERCR = DBLE(NRW)/DBLE(LRWK)*100.0D0
          PERCI = DBLE(NIW)/DBLE(LIWK)*100.0D0
!         Print statistics concerning workspace usage
10050     FORMAT(' Real    Workspace declared as ',I9, &
          ' is used up to ',I9,' (',F5.1,' percent)',//, &
          ' Integer Workspace declared as ',I9, &
          ' is used up to ',I9,' (',F5.1,' percent)',//)
          WRITE(LUMON,10050)LRWK,NRW,PERCR,LIWK,NIW,PERCI
        ENDIF
        IF(QINIMO)THEN
10051     FORMAT(/,' N =',I4,//,' Prescribed relative ', &
          'precision',D10.2,/)
          WRITE(LUMON,10051)N,RTOL
10052     FORMAT(' The Jacobian is supplied by ',A)
          IF (JACGEN.EQ.1) THEN
            WRITE(LUMON,10052) 'a user subroutine'
          ELSE IF (JACGEN.EQ.2) THEN
             WRITE(LUMON,10052) &
              'numerical differentiation (without feedback strategy)'
          ELSE IF (JACGEN.EQ.3) THEN
             WRITE(LUMON,10052) &
              'numerical differentiation (feedback strategy included)'
          ENDIF
10055     FORMAT(' The Jacobian will be stored in ',A,' mode')
          IF (MSTOR.EQ.0) THEN
            WRITE(LUMON,10055) 'full'
          ELSE IF (MSTOR.EQ.1) THEN
            WRITE(LUMON,10055) 'banded'
10056       FORMAT(' Lower bandwidth : ',I3,'   Upper bandwidth : ',I3)
            WRITE(LUMON,10056) ML,MU
          ENDIF
10057     FORMAT(' Automatic row scaling of the Jacobian is ',A,/)
          IF (IOPT(35).EQ.1) THEN
            WRITE(LUMON,10057) 'inhibited'
          ELSE
            WRITE(LUMON,10057) 'allowed'
          ENDIF
        ENDIF
        NONLIN=IOPT(31)
        IF (IOPT(38).EQ.0) QBDAMP = NONLIN.EQ.4
        IF (IOPT(38).EQ.1) QBDAMP = .TRUE.
        IF (IOPT(38).EQ.2) QBDAMP = .FALSE.
        IF (QBDAMP) THEN
          IF (RWK(20).LT.ONE) RWK(20) = TEN
        ENDIF
10064   FORMAT(' Rank-1 updates are ',A)
        IF (QINIMO) THEN
          IF (QRANK1) THEN
            WRITE(LUMON,10064) 'allowed'
          ELSE
            WRITE(LUMON,10064) 'inhibited'
          ENDIF
10065     FORMAT(' Problem is specified as being ',A)
          IF (NONLIN.EQ.1) THEN
            WRITE(LUMON,10065) 'linear'
          ELSE IF (NONLIN.EQ.2) THEN
            WRITE(LUMON,10065) 'mildly nonlinear'
          ELSE IF (NONLIN.EQ.3) THEN
            WRITE(LUMON,10065) 'highly nonlinear'
          ELSE IF (NONLIN.EQ.4) THEN
            WRITE(LUMON,10065) 'extremely nonlinear'
          ENDIF
10066     FORMAT(' Bounded damping strategy is ',A,:,/, &
                 ' Bounding factor is ',D10.3)
          IF (QBDAMP) THEN
            WRITE(LUMON,10066) 'active', RWK(20)
          ELSE
            WRITE(LUMON,10066) 'off'
          ENDIF
10067     FORMAT(' Special mode: ',A,' Newton iteration will be done')
          IF (IOPT(33).EQ.1) WRITE(LUMON,10067) 'Ordinary'
          IF (IOPT(34).EQ.1) WRITE(LUMON,10067) 'Simplified'
        ENDIF
!       Maximum permitted number of iteration steps
        NITMAX=IWK(31)
        IF (NITMAX.LE.0) NITMAX=50
        IWK(31)=NITMAX
10068   FORMAT(' Maximum permitted number of iteration steps : ', &
               I6)
        IF (QINIMO) WRITE(LUMON,10068) NITMAX
!       Initial damping factor for highly nonlinear problems
        QFCSTR=RWK(21).GT.ZERO
        IF (.NOT.QFCSTR) THEN
          RWK(21)=1.0D-2
          IF (NONLIN.EQ.4) RWK(21)=1.0D-4
        ENDIF
!       Minimal permitted damping factor
        IF (RWK(22).LE.ZERO) THEN
          RWK(22)=1.0D-4
          IF (NONLIN.EQ.4) RWK(22)=1.0D-8
        ENDIF
        FCMIN=RWK(22)
!       Rank1 decision parameter SIGMA
        IF (RWK(23).LT.ONE) RWK(23)=3.0D0
        IF (.NOT.QRANK1) RWK(23)=10.0D0/FCMIN
!       Decision parameter about increasing too small predictor
!       to greater corrector value
        IF (RWK(24).LT.ONE) RWK(24)=10.0D0/FCMIN
!       Starting value of damping factor (FCMIN.LE.FC.LE.1.0)
        IF(NONLIN.LE.2.AND..NOT.QFCSTR)THEN
!         for linear or mildly nonlinear problems
          FC = ONE
        ELSE
!         for highly or extremely nonlinear problems
          FC = RWK(21)
        ENDIF
!       Simplied Newton iteration implies ordinary Newton it. mode
        IF (IOPT(34).EQ.1) IOPT(33)=1
!       If ordinary Newton iteration, factor is always one
        IF (IOPT(33).EQ.1) FC=1.0D0
        RWK(21)=FC
        IF (MPRMON.GE.2.AND..NOT.QSUCC) THEN
10069     FORMAT(//,' Internal parameters:',//, &
            ' Starting value for damping factor FCSTART = ',D9.2,/, &
            ' Minimum allowed damping factor FCMIN = ',D9.2,/, &
            ' Rank-1 updates decision parameter SIGMA = ',D9.2)
          WRITE(LUMON,10069) RWK(21),FCMIN,RWK(23)
        ENDIF
!       Store lengths of currently required workspaces
        IWK(18) = NIFRIN-1
        IWK(19) = NRFRIN-1

!       Initialize and start time measurements monitor

        IF ( IOPT(1).EQ.0 .AND. MPRTIM.NE.0 ) THEN
          CALL MONINI (' NLEQ1',LUTIM)
          CALL MONDEF (0,'NLEQ1')
          CALL MONDEF (1,'FCN')
          CALL MONDEF (2,'Jacobi')
          CALL MONDEF (3,'Lin-Fact')
          CALL MONDEF (4,'Lin-Sol')
          CALL MONDEF (5,'Output')
          CALL MONSTR (IERR)
        ENDIF

        IERR=-1
!       If IERR is unmodified on exit, successive steps are required
!       to complete the Newton iteration
        IF (NBROY.EQ.0) NBROY=1
        CALL N1INT(N,FCN,JAC,X,XSCAL,RTOL,NITMAX,NONLIN,IOPT,IERR, &
        LRWK,RWK,NRFRIN,LIWK,IWK,NIFRIN, &
        M1,M2,NBROY, &
        RWK(L4),RWK(L41),RWK(L5),RWK(L51),RWK(L6),RWK(L63),RWK(L61), &
        RWK(L7), &
        RWK(L71),RWK(L9),RWK(L62),RWK(L11),RWK(L12),RWK(L121),RWK(L13), &
        RWK(21),RWK(22),RWK(23),RWK(24),RWK(IRWKI+1),RWK(IRWKI), &
        RWK(IRWKI+2),RWK(IRWKI+3),RWK(17),RWK(18),RWK(IRWKI+4),RWK(19), &
        MSTOR,MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,IWK(1),IWK(3), &
        IWK(4),IWK(5),IWK(8),IWK(9),IWK(33),IWK(24),QBDAMP)

        IF (MPRTIM.NE.0.AND.IERR.NE.-1.AND.IERR.NE.10) THEN
          CALL MONHLT
!          CALL MONPRT		! Gib suppressed this
        ENDIF

!       Free workspaces, so far not used between steps
        IWK(16) = NIWKFR
        IWK(17) = NRWKFR
      ENDIF
!     Print statistics
      IF (MPRMON.GE.1.AND.IERR.NE.-1.AND.IERR.NE.10) THEN
10080   FORMAT(/, '   ******  Statistics * ', A8, ' *******', /, &
                  '   ***  Newton iterations : ', I7,'  ***', /, &
                  '   ***  Corrector steps   : ', I7,'  ***', /, &
                  '   ***  Rejected rk-1 st. : ', I7,'  ***', /, &
                  '   ***  Jacobian eval.    : ', I7,'  ***', /, &
                  '   ***  Function eval.    : ', I7,'  ***', /, &
                  '   ***  ...  for Jacobian : ', I7,'  ***', /, &
                  '   *************************************', /)
        WRITE (LUMON,10080) PRODCT,IWK(1),IWK(3),IWK(9),IWK(5), &
        IWK(4),IWK(8)
      ENDIF
!     Print workspace requirements, if insufficient
      IF (IERR.EQ.10) THEN
10090   FORMAT(//,' ',20('*'),'Workspace Error',20('*'))
        IF (MPRERR.GE.1) WRITE(LUERR,10090)
        IF(NRW.GT.LRWK)THEN
10091     FORMAT(/,' Real Workspace dimensioned as',1X,I9, &
          1X,'must be enlarged at least up to ',I9,//)
          IF (MPRERR.GE.1) WRITE(LUERR,10091)LRWK,NRFRIN-1
        ENDIF
        IF(NIW.GT.LIWK)THEN
10092     FORMAT(/,' Integer Workspace dimensioned as ', &
          I9,' must be enlarged at least up ', &
          'to ',I9,//)
          IF (MPRERR.GE.1) WRITE(LUERR,10092)LIWK,NIFRIN-1
        ENDIF
      ENDIF
!     End of subroutine NLEQ1
      RETURN
      END subroutine

      SUBROUTINE N1PCHK(N,X,XSCAL,RTOL,IOPT,IERR,LIWK,IWK,LRWK,RWK)
!     Begin Prologue N1PCHK
      INTEGER N
      DOUBLE PRECISION X(N),XSCAL(N)
      DOUBLE PRECISION RTOL
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
!     ------------------------------------------------------------
!
!     Summary :
!
!     N 1 P C H K : Checking of input parameters and options
!            for NLEQ1.
!
!     Parameters:
!     ===========
!
!     See parameter description in driver routine.
!
!     Subroutines called: ZIBCONST
!
!     Machine dependent constants used:
!     =================================
!
!     EPMACH = relative machine precision
!     GREAT = squareroot of maxreal divided by 10
!     SMALL = squareroot of "smallest positive machine number
!             divided by relative machine precision"
      DOUBLE PRECISION EPMACH,GREAT,SMALL
!
!     ------------------------------------------------------------
!     End Prologue
!
      INTRINSIC DBLE
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
!
      integer, PARAMETER :: NUMOPT=50
      INTEGER MSTOR,IOPTL(NUMOPT),IOPTU(NUMOPT)
      DOUBLE PRECISION TOLMIN,TOLMAX,DEFSCL
      integer :: MPRERR, LUERR, NONLIN, I
!
      DATA IOPTL /0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1, &
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                  0,0,0,0,0,0,0,0,0,0, &
                  -9999999,-9999999,-9999999,-9999999,-9999999/
      DATA IOPTU /1,1,3,1,0,9999999,9999999,0,1,0,3,99,6,99,3,99,0,0,1, &
                  99,0,0,0,0,0,0,0,0,0,0,4,1,1,1,1, &
                  0,0,2,3,0,0,0,0,0,0, &
                  9999999,9999999,9999999,9999999,9999999/
!
      CALL ZIBCONST(EPMACH,SMALL)
      GREAT  = 1.0D0/SMALL
      IERR = 0
!        Print error messages?
      MPRERR = IOPT(11)
      LUERR = IOPT(12)
      IF (LUERR .LE. 0 .OR. LUERR .GT. 99) THEN
        LUERR = 6
        IOPT(12)=LUERR
      ENDIF

!     Checking dimensional parameter N
      IF ( N.LE.0 ) THEN
        IF (MPRERR.GE.1)  WRITE(LUERR,10011) N
10011   FORMAT(/,' Error: Bad input to dimensional parameter N supplied' &
               ,/,8X,'choose N positive, your input is: N = ',I5)
        IERR = 20
      ENDIF

!     Problem type specification by user
      NONLIN=IOPT(31)
      IF (NONLIN.EQ.0) NONLIN=3
      IOPT(31)=NONLIN

!     Checking and conditional adaption of the user-prescribed RTOL
      IF (RTOL.LE.ZERO) THEN
        IF (MPRERR.GE.1)  &
            WRITE(LUERR,'(/,A)') ' Error: Nonpositive RTOL supplied'
        IERR = 21
      ELSE
        TOLMIN = EPMACH*TEN*DBLE(N)
        IF(RTOL.LT.TOLMIN) THEN
          RTOL = TOLMIN
          IF (MPRERR.GE.2)  &
            WRITE(LUERR,10012) 'increased ','smallest',RTOL
        ENDIF
        TOLMAX = 1.0D-1
        IF(RTOL.GT.TOLMAX) THEN
          RTOL = TOLMAX
          IF (MPRERR.GE.2)  &
            WRITE(LUERR,10012) 'decreased ','largest',RTOL
        ENDIF
10012   FORMAT(/,' Warning: User prescribed RTOL ',A,'to ', &
               'reasonable ',A,' value RTOL = ',D11.2)
      ENDIF
!
!     Test user prescribed accuracy and scaling on proper values
      IF (N.LE.0) RETURN
      IF (NONLIN.GE.3) THEN
        DEFSCL = RTOL
      ELSE
        DEFSCL = ONE
      ENDIF
      DO 10 I=1,N
        IF (XSCAL(I).LT.ZERO) THEN
          IF (MPRERR.GE.1) THEN
            WRITE(LUERR,10013) I
10013       FORMAT(/,' Error: Negative value in XSCAL(',I5,') supplied')
          ENDIF
          IERR = 22
        ENDIF
        IF (XSCAL(I).EQ.ZERO) XSCAL(I) = DEFSCL
        IF ( XSCAL(I).GT.ZERO .AND. XSCAL(I).LT.SMALL ) THEN
          IF (MPRERR.GE.2) THEN
            WRITE(LUERR,10014) I,XSCAL(I),SMALL
10014       FORMAT(/,' Warning: XSCAL(',I5,') = ',D9.2,' too small, ', &
                   'increased to',D9.2)
          ENDIF
          XSCAL(I) = SMALL
        ENDIF
        IF (XSCAL(I).GT.GREAT) THEN
          IF (MPRERR.GE.2) THEN
            WRITE(LUERR,10015) I,XSCAL(I),GREAT
10015       FORMAT(/,' Warning: XSCAL(',I5,') = ',D9.2,' too big, ', &
                   'decreased to',D9.2)
          ENDIF
          XSCAL(I) = GREAT
        ENDIF
10    CONTINUE
!     Special dependence on Jacobian's storage mode for ML and MU
      MSTOR = IOPT(4)
      IF (MSTOR.EQ.0) THEN
        IOPTU(6)=0
        IOPTU(7)=0
      ELSE IF (MSTOR.EQ.1) THEN
        IOPTU(6)=N-1
        IOPTU(7)=N-1
      ENDIF
!     Checks options
      DO 20 I=1,30
        IF (IOPT(I).LT.IOPTL(I) .OR. IOPT(I).GT.IOPTU(I)) THEN
          IERR=30
          IF (MPRERR.GE.1) THEN
            WRITE(LUERR,20001) I,IOPT(I),IOPTL(I),IOPTU(I)
20001       FORMAT(' Invalid option specified: IOPT(',I2,')=',I12,';', &
                   /,3X,'range of permitted values is ',I8,' to ',I8)
          ENDIF
        ENDIF
20    CONTINUE
!     End of subroutine N1PCHK
      RETURN
      END subroutine

      SUBROUTINE N1INT(N,FCN,JAC,X,XSCAL,RTOL,NITMAX,NONLIN,IOPT,IERR, &
      LRWK,RWK,NRWKFR,LIWK,IWK,NIWKFR,M1,M2,NBROY, &
      A,DXSAVE,DX,DXQ,XA,XWA,F,FA,ETA,XW,FW,DXQA,T1,T2,T3,FC,FCMIN, &
      SIGMA,SIGMA2,FCA,FCKEEP,FCPRI,DMYCOR,CONV,SUMX,SUMXS,DLEVF,MSTOR, &
      MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,NITER,NCORR,NFCN,NJAC, &
      NFCNJ,NREJR1,NEW,ICONV,QBDAMP)
!     Begin Prologue N1INT
      INTEGER N
!      EXTERNAL FCN,JAC
      procedure(func) :: FCN
      procedure(jacfunc) :: JAC
      DOUBLE PRECISION X(N),XSCAL(N)
      DOUBLE PRECISION RTOL
      INTEGER NITMAX,NONLIN
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER NRWKFR,LIWK
      INTEGER IWK(LIWK)
      INTEGER NIWKFR,M1,M2,NBROY
      DOUBLE PRECISION A(M1,N),DXSAVE(N,NBROY)
      DOUBLE PRECISION DX(N),DXQ(N),XA(N),XWA(N),F(N),FA(N),ETA(N)
      DOUBLE PRECISION XW(N),FW(N),DXQA(N),T1(N),T2(N),T3(N)
      DOUBLE PRECISION FC,FCMIN,SIGMA,SIGMA2,FCA,FCKEEP,CONV,SUMX,SUMXS, &
                       DLEVF,FCPRI,DMYCOR
      INTEGER MSTOR,MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,NITER, &
      NCORR,NFCN,NJAC,NFCNJ,NREJR1,NEW,ICONV
      LOGICAL QBDAMP
!     ------------------------------------------------------------
!
!     Summary :
!
!     N 1 I N T : Core routine for NLEQ1 .
!     Damped Newton-algorithm for systems of highly nonlinear
!     equations especially designed for numerically sensitive
!     problems.
!
!     Parameters:
!     ===========
!
!       N,FCN,JAC,X,XSCAL,RTOL
!                  See parameter description in driver routine
!
!       NITMAX      Int    Maximum number of allowed iterations
!       NONLIN      Int    Problem type specification
!                   (see IOPT-field NONLIN)
!       IOPT        Int    See parameter description in driver routine
!       IERR        Int    See parameter description in driver routine
!       LRWK        Int    Length of real workspace
!       RWK(LRWK)   Dble   Real workspace array
!       NRWKFR      Int    First free position of RWK on exit
!       LIWK        Int    Length of integer workspace
!       IWK(LIWK)   Int    Integer workspace array
!       NIWKFR      Int    First free position of IWK on exit
!       M1          Int    Leading dimension of Jacobian array A
!                   for full case Jacobian: N
!                   for banded Jacobian: 2*ML+MU+1;
!                   ML, MU see IOPT-description in driver routine
!       M2          Int    for full case Jacobian: N
!                   for banded Jacobian: ML+MU+1
!       NBROY       Int    Maximum number of possible consecutive
!                   iterative Broyden steps. (See IWK(36))
!       A(M1,N)     Dble   Holds the Jacobian matrix (decomposed form
!                   after call of linear decomposition
!                   routine)
!       DXSAVE(X,NBROY)
!            Dble   Used to save the quasi Newton corrections of
!                   all previously done consecutive Broyden
!                   steps.
!       DX(N)       Dble   Current Newton correction
!       DXQ(N)      Dble   Simplified Newton correction J(k-1)*X(k)
!       XA(N)       Dble   Previous Newton iterate
!       XWA(N)      Dble   Scaling factors used for latest decomposed
!                   Jacobian for column scaling - may differ
!                   from XW, if Broyden updates are performed
!       F(N)        Dble   Function (FCN) value of current iterate
!       FA(N)       Dble   Function (FCN) value of previous iterate
!       ETA(N)      Dble   Jacobian approximation: updated scaled
!                   denominators
!       XW(N)       Dble   Scaling factors for iteration vector
!       FW(N)       Dble   Scaling factors for rows of the system
!       DXQA(N)     Dble   Previous Newton correction
!       T1(N)       Dble   Workspace for linear solvers and internal
!                   subroutines
!       T2(N)       Dble   Workspace array for internal subroutines
!       T3(N)       Dble   Workspace array for internal subroutines
!       F!          Dble   Current Newton iteration damping factor.
!       FCMIN       Dble   Minimum permitted damping factor. If
!                   FC becomes smaller than this value, one
!                   of the following may occur:
!                   a.    Recomputation of the Jacobian
!                                matrix by means of difference
!                                approximation (instead of Rank1
!                                update), if Rank1 - update
!                                previously was used
!                   b.    Fail exit otherwise
!       SIGMA       Dble   Decision parameter for rank1-updates
!       SIGMA2      Dble   Decision parameter for damping factor
!                   increasing to corrector value
!       FCA         Dble   Previous Newton iteration damping factor.
!       FCKEEP      Dble   Keeps the damping factor as it is at start
!                   of iteration step.
!       CONV        Dble   Scaled maximum norm of the Newton-
!                   correction. Passed to RWK-field on output.
!       SUMX        Dble   Square of the natural level (see equal-
!                   named RWK-output field)
!       SUMXS       Dble   Square of the "simplified" natural level
!                   (see equal-named RWK-internal field)
!       DLEVF       Dble   The standard level (see equal-
!                   named RWK-output field)
!       MSTOR       Int    see description of IOPT-field MSTOR
!       MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,
!       NITER,NCORR,NFCN,NJAC,NFCNJ,NREJR1,NEW :
!                   See description of equal named IWK-fields
!                   in the driver subroutine
!       QBDAMP      Logic  Flag, that indicates, whether bounded damping
!                   strategy is active:
!                   .true.  = bounded damping strategy is active
!                   .false. = normal damping strategy is active
!
!     Internal double variables
!     =========================
!
!       AJDEL    See RWK(26) (num. diff. without feedback)
!       AJMIN    See RWK(27) (num. diff. without feedback)
!       CONVA    Holds the previous value of CONV .
!       DMUE     Temporary value used during computation of damping
!         factors predictor.
!       EPDIFF   sqrt(10*epmach) (num. diff. with feedback)
!       ETADIF   See description of RWK(28) (num. diff. with feedback)
!       ETAINI   Initial value for all ETA-components (num. diff. fb.)
!       ETAMAX   Maximum allowed pertubation (num. diff. with feedback)
!       ETAMIN   Minimum allowed pertubation (num. diff. with feedback)
!       FCDNM    Used to compute the denominator of the damping
!         factor FC during computation of it's predictor,
!         corrector and aposteriori estimate (in the case of
!         performing a Rank1 update) .
!       FCK2     Aposteriori estimate of FC.
!       FCH      Temporary storage of new corrector FC.
!       FCMIN2   FCMIN**2 . Used for FC-predictor computation.
!       FCNUMP   Gets the numerator of the predictor formula for FC.
!       FCNMP2   Temporary used for predictor numerator computation.
!       FCNUMK   Gets the numerator of the corrector computation
!         of FC .
!       SUMXA    Natural level of the previous iterate.
!       TH       Temporary variable used during corrector- and
!         aposteriori computations of FC.
!
!     Internal integer variables
!     ==========================
!
!     IFAIL      Gets the return value from subroutines called from
!         N1INT (N1FACT, N1SOLV, FCN, JAC)
!     ISCAL      Holds the scaling option from the IOPT-field ISCAL
!     MODE       Matrix storage mode (see IOPT-field MODE)
!     NRED       Count of successive corrector steps
!     NILUSE     Gets the amount of IWK used by the linear solver
!     NRLUSE     Gets the amount of RWK used by the linear solver
!     NIWLA      Index of first element of IWK provided to the
!         linear solver
!     NRWLA      Index of first element of RWK provided to the
!         linear solver
!     LIWL       Holds the maximum amount of integer workspace
!         available to the linear solver
!     LRWL       Holds the maximum amount of real workspace
!         available to the linear solver
!
!
!     Internal logical variables
!     ==========================
!
!     QGENJ      Jacobian updating technique flag:
!         =.TRUE.  : Call of analytical subroutine JAC or
!                           numerical differentiation
!         =.FALSE. : rank1- (Broyden-) update
!     QINIS!     Iterate initial-scaling flag:
!         =.TRUE.  : at first call of N1SCAL
!         =.FALSE. : at successive calls of N1SCAL
!     QSUC!      See description of IOPT-field QSUCC.
!     QJCRFR     Jacobian refresh flag:
!         set to .TRUE. if damping factor gets too small
!         and Jacobian was computed by rank1-update.
!         Indicates, that the Jacobian needs to be recomputed
!         by subroutine JAC or numerical differentation.
!     QLINIT     Initialization state of linear solver workspace:
!         =.FALSE. : Not yet initialized
!         =.TRUE.  : Initialized - N1FACT has been called at
!                           least one time.
!     QSCALE     Holds the value of .NOT.QNSCAL. See description
!         of IOPT-field QNSCAL.
!
!     Subroutines called:
!     ===================
!
!       N1FACT, N1SOLV, N1JAC,  N1JACB, N1JCF,  N1JCFB, N1LVLS,
!       N1SCRF, N1SCRB, N1SOUT, N1PRV1, N1PRV2, N1SCAL,
!       MONON,  MONOFF
!
!     Functions called:
!     =================
!
!       ZIBCONST, WNORM
!
!     Machine constants used
!     ======================
!
      DOUBLE PRECISION EPMACH,SMALL
!
!     ------------------------------------------------------------
!     End Prologue
!      EXTERNAL N1FACT, N1SOLV, N1JAC,  N1JACB, N1JCF,  N1JCFB, N1LVLS,
!     $         N1SCRF, N1SCRB, N1SOUT, N1PRV1, N1PRV2, N1SCAL,
!     $         MONON,  MONOFF, WNORM
      INTRINSIC DSQRT,DMIN1,MAX0,MIN0
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=10.0D0)
      INTEGER IFAIL,ILOOP,ISCAL,K,MODE,NRED,NILUSE,NRLUSE,NIWLA, &
      NRWLA,LIWL,LRWL,L1,L2,JACGEN
      DOUBLE PRECISION AJDEL,AJMIN,ALFA1,ALFA2,ALFA,BETA,CONVA, &
      DLEVXA,DMYPRI,DXANRM,DXNRM,EPDIFF,ETAMIN,ETAMAX, &
      ETAINI,ETADIF,FCDNM,FCMIN2,FCNUMP,FCCOR,FCNMP2,FCH,FCBND, &
      FCBH,FCK2,FCNUMK,FCREDU,DLEVFN,SUMXA,SUM1,SUM2,S1,TH,RSMALL, &
      APREC
!  WNORM
      LOGICAL QGENJ,QINISC,QSUCC,QJCRFR,QLINIT,QSCALE,QNEXT,QREP, &
              QRANK1,QMIXIO,QORDI,QSIMPL,QLU
!WEI
      INTRINSIC DLOG
      DOUBLE PRECISION CLIN0,CLIN1,CALPHA,CALPHK,ALPHAE,ALPHAK,ALPHAA, &
                       SUMXA0,SUMXA1,SUMXA2,SUMXTE,FCMON,DLOG
      INTEGER IORMON, MPRTIM, ML, MU, L3, IDUMMY, I, MODEFI
      LOGICAL QMSTOP
      SAVE CLIN0,CLIN1,CALPHA,ALPHAE,ALPHAK,ALPHAA,SUMXA0,SUMXA1,SUMXA2, &
           QMSTOP

      CALL ZIBCONST(EPMACH,SMALL)
!     Begin
!       ----------------------------------------------------------
!       1 Initialization
!       ----------------------------------------------------------
!       1.1 Control-flags and -integers
        QSUCC = IOPT(1).EQ.1
        QSCALE = .NOT. IOPT(35).EQ.1
        QORDI  = IOPT(33).EQ.1
        QSIMPL = IOPT(34).EQ.1
        QRANK1 = IOPT(32).EQ.1
        IORMON = IOPT(39)
        IF (IORMON.EQ.0) IORMON=2
        ISCAL = IOPT(9)
        MODE = IOPT(2)
        JACGEN = IOPT(3)
        QMIXIO = LUMON.EQ.LUSOL .AND. MPRMON.NE.0 .AND. MPRSOL.NE.0
        QLU    = .NOT. QSIMPL
        MPRTIM = IOPT(19)
        
        FCMON = 999	! Gib added this to prevent the runtime error of FCMON undefined
!       ----------------------------------------------------------
!       1.2 Derivated dimensional parameters
        IF (MSTOR.EQ.0) THEN
          ML=0
        ELSE IF (MSTOR.EQ.1) THEN
          ML=M1-M2
          MU=M2-1-ML
        ENDIF
!       ----------------------------------------------------------
!       1.3 Derivated internal parameters
        FCMIN2 = FCMIN*FCMIN
        RSMALL = DSQRT(TEN*RTOL)
!       ----------------------------------------------------------
!       1.4 Adaption of input parameters, if necessary
        IF(FC.LT.FCMIN) FC = FCMIN
        IF(FC.GT.ONE) FC = ONE
!       ----------------------------------------------------------
!       1.5 Initial preparations
        QJCRFR = .FALSE.
        QLINIT = .FALSE.
        IFAIL = 0
        FCBND = ZERO
        IF (QBDAMP) FCBND = RWK(20)
!       ----------------------------------------------------------
!       1.5.1 Numerical differentiation related initializations
        IF (JACGEN.EQ.2) THEN
          AJDEL = RWK(26)
          IF (AJDEL.LE.SMALL) AJDEL = DSQRT(EPMACH*TEN)
          AJMIN = RWK(27)
        ELSE IF (JACGEN.EQ.3) THEN
          ETADIF = RWK(28)
          IF (ETADIF .LE. SMALL) ETADIF = 1.0D-6
          ETAINI = RWK(29)
          IF (ETAINI .LE. SMALL) ETAINI = 1.0D-6
          EPDIFF = DSQRT(EPMACH*TEN)
          ETAMAX = DSQRT(EPDIFF)
          ETAMIN = EPDIFF*ETAMAX
        ENDIF
!       ----------------------------------------------------------
!       1.5.2 Miscellaneous preparations of first iteration step
        IF (.NOT.QSUCC) THEN
          NITER = 0
          NCORR = 0
          NREJR1 = 0
          NFCN = 0
          NJAC = 0
          NFCNJ = 0
          QGENJ = .TRUE.
          QINISC = .TRUE.
          FCKEEP = FC
          FCA = FC
          FCPRI = FC
          FCK2 = FC
          CONV = ZERO
          IF (JACGEN.EQ.3) THEN
            DO 1520 L1=1,N
              ETA(L1)=ETAINI
1520        CONTINUE
          ENDIF
          DO 1521 L1=1,N
            XA(L1)=X(L1)
1521      CONTINUE
!WEI
          ICONV = 0
          ALPHAE = ZERO
          SUMXA1 = ZERO
          SUMXA0 = ZERO
          CLIN0  = ZERO
          QMSTOP = .FALSE.
!         ------------------------------------------------------
!         1.6 Print monitor header
          IF(MPRMON.GE.2 .AND. .NOT.QMIXIO)THEN
16003       FORMAT(///,2X,66('*'))
            WRITE(LUMON,16003)
16004       FORMAT(/,8X,'It',7X,'Normf ',10X,'Normx ',8X, &
                   'Damp.Fct.',3X,'New')
            WRITE(LUMON,16004)
          ENDIF
!         --------------------------------------------------------
!         1.7 Startup step
!         --------------------------------------------------------
!         1.7.1 Computation of the residual vector
          IF (MPRTIM.NE.0) CALL MONON(1)
          CALL FCN(N,X,F,IFAIL)
          IF (MPRTIM.NE.0) CALL MONOFF(1)
          NFCN = NFCN+1
!     Exit, if ...
          IF (IFAIL.NE.0) THEN
            IERR = 82
            GOTO 4299
          ENDIF
        ELSE
          QINISC = .FALSE.
        ENDIF
!
!       Main iteration loop
!       ===================
!
!       Repeat
2       CONTINUE
!         --------------------------------------------------------
!         2 Startup of iteration step
          IF (.NOT.QJCRFR) THEN
!           ------------------------------------------------------
!           2.1 Scaling of variables X(N)
            CALL N1SCAL(N,X,XA,XSCAL,XW,ISCAL,QINISC,IOPT,LRWK,RWK)
            QINISC = .FALSE.
            IF(NITER.NE.0)THEN
!             ----------------------------------------------------
!             2.2 Aposteriori estimate of damping factor
              DO 2200 L1=1,N
                DXQA(L1)=DXQ(L1)
2200          CONTINUE
              IF (.NOT.QORDI) THEN
                FCNUMP = ZERO
                DO 2201 L1=1,N
                  FCNUMP=FCNUMP+(DX(L1)/XW(L1))**2
2201            CONTINUE
                TH = FC-ONE
                FCDNM = ZERO
                DO 2202 L1=1,N
                  FCDNM=FCDNM+((DXQA(L1)+TH*DX(L1))/XW(L1))**2
2202            CONTINUE
!        --------------------------------------------------
!        2.2.2 Decision criterion for Jacobian updating
!              technique:
!              QGENJ.EQ..TRUE. numerical differentation,
!              QGENJ.EQ..FALSE. rank1 updating
                QGENJ = .TRUE.
                IF (FC.EQ.FCPRI) THEN
                  QGENJ = FC.LT.ONE.OR.FCA.LT.ONE.OR.DMYCOR.LE.FC*SIGMA &
                          .OR. .NOT.QRANK1 .OR. NEW+2.GT.NBROY
                  FCA = FC
                ELSE
                  DMYCOR = FCA*FCA*HALF*DSQRT(FCNUMP/FCDNM)
                  IF (NONLIN.LE.3) THEN
                    FCCOR = DMIN1(ONE,DMYCOR)
                  ELSE
                    FCCOR = DMIN1(ONE,HALF*DMYCOR)
                  ENDIF
                  FCA = DMAX1(DMIN1(FC,FCCOR),FCMIN)
!$Test-begin
                  IF (MPRMON.GE.5) THEN
                    WRITE(LUMON,22201) FCCOR, FC, DMYCOR, FCNUMP, FCDNM
22201               FORMAT (/, ' +++ aposteriori estimate +++', /, &
                          ' FCCOR  = ', D18.10, '  F!     = ', D18.10, /, &
                          ' DMYCOR = ', D18.10, '  FCNUMP = ', D18.10, /, &
                          ' FCDNM  = ', D18.10, /, &
                          ' ++++++++++++++++++++++++++++', /)
                  ENDIF
!$Test-end
                ENDIF
                FCK2 = FCA
!        ------------------------------------------------------
!        2.2.1 Computation of the numerator of damping
!              factor predictor
                FCNMP2 = ZERO
                DO 221 L1=1,N
                  FCNMP2=FCNMP2+(DXQA(L1)/XW(L1))**2
221             CONTINUE
                FCNUMP = FCNUMP*FCNMP2
              ENDIF
            ENDIF
          ENDIF
          QJCRFR =.FALSE.
!         --------------------------------------------------------
!         2.3 Jacobian matrix (stored to array A(M1,N))
!         --------------------------------------------------------
!         2.3.1 Jacobian generation by routine JAC or
!        difference approximation (If QGENJ.EQ..TRUE.)
!        - or -
!        Rank-1 update of Jacobian (If QGENJ.EQ..FALSE.)
          IF (QGENJ .AND. (.NOT.QSIMPL .OR. NITER.EQ.0)) THEN
            NEW = 0
            IF (JACGEN.EQ.1) THEN
               IF (MPRTIM.NE.0) CALL MONON(2)
               CALL JAC(N,M1,X,A,IFAIL)
               IF (MPRTIM.NE.0) CALL MONOFF(2)
            ELSE
              IF (MSTOR.EQ.0) THEN
                IF (MPRTIM.NE.0) CALL MONON(2)
                IF (JACGEN.EQ.3)  &
                  CALL N1JCF(FCN,N,N,X,F,A,XW,ETA,ETAMIN,ETAMAX, &
                             ETADIF,CONV,NFCNJ,T1,IFAIL)
                IF (JACGEN.EQ.2)  &
                  CALL N1JAC(FCN, N, N, X, F, A, XW, AJDEL, AJMIN, &
                             NFCNJ, T1, IFAIL)
                IF (MPRTIM.NE.0) CALL MONOFF(2)
              ELSE IF (MSTOR.EQ.1) THEN
                IF (MPRTIM.NE.0) CALL MONON(2)
                IF (JACGEN.EQ.3)  &
                  CALL N1JCFB(FCN,N,M1,ML,X,F,A,XW,ETA,ETAMIN, &
                              ETAMAX,ETADIF,CONV,NFCNJ,T1,T2,T3,IFAIL)
                IF (JACGEN.EQ.2)  &
                  CALL N1JACB (FCN, N, M1, ML, X, F, A, XW,  AJDEL, &
                               AJMIN, NFCNJ, T1, T2, T3, IFAIL)
                IF (MPRTIM.NE.0) CALL MONOFF(2)
              ENDIF
            ENDIF
            NJAC = NJAC + 1
!     Exit, If ...
            IF (JACGEN.EQ.1 .AND. IFAIL.LT.0) THEN
              IERR = 83
              GOTO 4299
            ENDIF
            IF (JACGEN.NE.1 .AND. IFAIL.NE.0) THEN
              IERR = 82
              GOTO 4299
            ENDIF
          ELSE IF (.NOT.QSIMPL) THEN
            NEW = NEW+1
          ENDIF
          IF ( NEW.EQ.0 .AND. (QLU.OR.NITER.EQ.0) ) THEN
!             ------------------------------------------------------
!             2.3.2.1 Save scaling values
              DO 2321 L1=1,N
                XWA(L1) = XW(L1)
2321          CONTINUE
!             ------------------------------------------------------
!             2.3.2.2 Prepare Jac. for use by band-solver DGBFA/DGBSL
              IF (MSTOR.EQ.1) THEN
                DO 2322 L1=1,N
                  DO 2323 L2=M2,1,-1
                    A(L2+ML,L1)=A(L2,L1)
2323              CONTINUE
2322            CONTINUE
              ENDIF
!             ------------------------------------------------------
!             2.4 Prepare solution of the linear system
!             ------------------------------------------------------
!             2.4.1 internal column scaling of matrix A
              IF (MSTOR.EQ.0) THEN
                DO 2410 K=1,N
                  S1 =-XW(K)
                  DO 2412 L1=1,N
                    A(L1,K)=A(L1,K)*S1
2412              CONTINUE
2410            CONTINUE
              ELSE IF (MSTOR.EQ.1) THEN
                DO 2413 K=1,N
                  L2=MAX0(1+M2-K,ML+1)
                  L3=MIN0(N+M2-K,M1)
                  S1 =-XW(K)
                  DO 2415 L1=L2,L3
                    A(L1,K)=A(L1,K)*S1
2415              CONTINUE
2413            CONTINUE
              ENDIF
!             ----------------------------------------------------
!             2.4.2 Row scaling of matrix A
              IF (QSCALE) THEN
                IF (MSTOR.EQ.0) THEN
                  CALL N1SCRF(N,N,A,FW)
                ELSE IF (MSTOR.EQ.1) THEN
                  CALL N1SCRB(N,M1,ML,MU,A,FW)
                ENDIF
              ELSE
                DO 242 K=1,N
                  FW(K)=ONE
242             CONTINUE
              ENDIF
          ENDIF
!         --------------------------------------------------------
!         2.4.3 Save and scale values of F(N)
          DO 243 L1=1,N
            FA(L1)=F(L1)
            T1(L1)=F(L1)*FW(L1)
243       CONTINUE
!         --------------------------------------------------------
!         3 Central part of iteration step
!         --------------------------------------------------------
!         3.1 Solution of the linear system
!         --------------------------------------------------------
!         3.1.1 Decomposition of (N,N)-matrix A
          IF (.NOT.QLINIT) THEN
            NIWLA = IWK(18)+1
            NRWLA = IWK(19)+1
            LIWL = LIWK-NIWLA+1
            LRWL = LRWK-NRWLA+1
          ENDIF
          IF ( NEW.EQ.0 .AND. (QLU.OR.NITER.EQ.0)) THEN
            IF (MPRTIM.NE.0) CALL MONON(3)
            CALL N1FACT(N,M1,ML,MU,A,IOPT,IFAIL,LIWL,IWK(NIWLA), &
                        NILUSE,LRWL,RWK(NRWLA),NRLUSE)
            IF (MPRTIM.NE.0) CALL MONOFF(3)
            IF (.NOT.QLINIT) THEN
              NIWKFR = NIWKFR+NILUSE
              NRWKFR = NRWKFR+NRLUSE
!             Store lengths of currently required workspaces
              IWK(18) = NIWKFR-1
              IWK(19) = NRWKFR-1
            ENDIF
!       Exit Repeat If ...
            IF(IFAIL.NE.0) THEN
              IF (IFAIL.EQ.1) THEN
                IERR = 1
              ELSE
                IERR = 80
              ENDIF
              GOTO 4299
            ENDIF
          ENDIF
          QLINIT = .TRUE.
!         --------------------------------------------------------
!         3.1.2 Solution of linear (N,N)-system
          IF(NEW.EQ.0) THEN
            IF (MPRTIM.NE.0) CALL MONON(4)
            CALL N1SOLV(N,M1,ML,MU,A,T1,IOPT,IFAIL,LIWL,IWK(NIWLA), &
                       IDUMMY,LRWL,RWK(NRWLA),IDUMMY)
            IF (MPRTIM.NE.0) CALL MONOFF(4)
!     Exit Repeat If ...
            IF(IFAIL.NE.0)  THEN
              IERR = 81
              GOTO 4299
            ENDIF
          ELSE
            ALFA1=ZERO
            ALFA2=ZERO
            DO 3121 I=1,N
              ALFA1=ALFA1+DX(I)*DXQ(I)/XW(I)**2
              ALFA2=ALFA2+DX(I)**2/XW(I)**2
3121        CONTINUE
            ALFA=ALFA1/ALFA2
            BETA=ONE-ALFA
            DO 3122 I=1,N
              T1(I)=(DXQ(I)+(FCA-ONE)*ALFA*DX(I))/BETA
3122        CONTINUE
            IF(NEW.EQ.1) THEN
              DO 3123 I=1,N
                DXSAVE(I,1)=DX(I)
3123          CONTINUE
            ENDIF
            DO 3124 I=1,N
              DXSAVE(I,NEW+1)=T1(I)
              DX(I)=T1(I)
              T1(I)=T1(I)/XW(I)
3124        CONTINUE
          ENDIF
!         --------------------------------------------------------
!         3.2 Evaluation of scaled natural level function SUMX
!             scaled maximum error norm CONV
!             evaluation of (scaled) standard level function
!             DLEVF ( DLEVF only, if MPRMON.GE.2 )
!             and computation of ordinary Newton corrections
!             DX(N)
          IF (.NOT. QSIMPL) THEN
            CALL N1LVLS(N,T1,XW,F,DX,CONV,SUMX,DLEVF,MPRMON,NEW.EQ.0)
          ELSE
            CALL N1LVLS(N,T1,XWA,F,DX,CONV,SUMX,DLEVF,MPRMON,NEW.EQ.0)
          ENDIF
          DO 32 L1=1,N
            XA(L1)=X(L1)
32        CONTINUE
          SUMXA = SUMX
          DLEVXA = DSQRT(SUMXA/DBLE(FLOAT(N)))
          CONVA = CONV
          DXANRM = WNORM(N,DX,XW)
!         --------------------------------------------------------
!         3.3 A - priori estimate of damping factor FC
          IF(NITER.NE.0.AND.NONLIN.NE.1.AND.NEW.EQ.0.AND. &
             .NOT. QORDI)THEN
!           ------------------------------------------------------
!           3.3.1 Computation of the denominator of a-priori
!          estimate
            FCDNM = ZERO
            DO 331 L1=1,N
              FCDNM=FCDNM+((DX(L1)-DXQA(L1))/XW(L1))**2
331         CONTINUE
            FCDNM = FCDNM*SUMX
!           ------------------------------------------------------
!           3.3.2 New damping factor
            IF(FCDNM.GT.FCNUMP*FCMIN2 .OR. &
              (NONLIN.EQ.4 .AND. FCA**2*FCNUMP .LT. 4.0D0*FCDNM)) THEN
              DMYPRI = FCA*DSQRT(FCNUMP/FCDNM)
              FCPRI = DMIN1(DMYPRI,ONE)
              IF (NONLIN.EQ.4) FCPRI = DMIN1(HALF*DMYPRI,ONE)
            ELSE
              FCPRI = ONE
!$Test-begin
              DMYPRI = -1.0D0
!$Test-end
            ENDIF
!$Test-begin
            IF (MPRMON.GE.5) THEN
              WRITE(LUMON,33201) FCPRI, FC, FCA, DMYPRI, FCNUMP, FCDNM
33201         FORMAT (/, ' +++ apriori estimate +++', /, &
                      ' FCPRI  = ', D18.10, '  F!     = ', D18.10, /, &
                      ' FCA    = ', D18.10, '  DMYPRI = ', D18.10, /, &
                      ' FCNUMP = ', D18.10, '  FCDNM  = ', D18.10, /, &
                         ' ++++++++++++++++++++++++', /)
            ENDIF
!$Test-end
            FC = DMAX1(FCPRI,FCMIN)
            IF (QBDAMP) THEN
              FCBH = FCA*FCBND
              IF (FC.GT.FCBH) THEN
                FC = FCBH
                IF (MPRMON.GE.4) &
                  WRITE(LUMON,*) ' *** incr. rest. act. (a prio) ***'
              ENDIF
              FCBH = FCA/FCBND
              IF (FC.LT.FCBH) THEN
                FC = FCBH
                IF (MPRMON.GE.4) &
                  WRITE(LUMON,*) ' *** decr. rest. act. (a prio) ***'
              ENDIF
            ENDIF
          ENDIF
!WEI
          IF (IORMON.GE.2) THEN
            SUMXA2=SUMXA1
            SUMXA1=SUMXA0
            SUMXA0=DLEVXA
            IF (SUMXA0.EQ.ZERO) SUMXA0=SMALL
!           Check convergence rates (linear and superlinear)
!           ICONV : Convergence indicator
!            =0: No convergence indicated yet
!            =1: Damping factor is 1.0d0
!            =2: Superlinear convergence detected (alpha >=1.2)
!            =3: Quadratic convergence detected (alpha > 1.8)
            FCMON = DMIN1(FC,FCMON)
            IF (FCMON.LT.ONE) THEN
              ICONV = 0
              ALPHAE = ZERO
            ENDIF
            IF (FCMON.EQ.ONE .AND. ICONV.EQ.0) ICONV=1
            IF (NITER.GE.1) THEN
              CLIN1 = CLIN0
              CLIN0 = SUMXA0/SUMXA1
            ENDIF
            IF (ICONV.GE.1.AND.NITER.GE.2) THEN
              ALPHAK = ALPHAE
              ALPHAE = ZERO
              IF (CLIN1.LE.0.95D0) ALPHAE = DLOG(CLIN0)/DLOG(CLIN1)
              IF (ALPHAK.NE.ZERO) ALPHAK =0.5D0*(ALPHAE+ALPHAK)
              ALPHAA = DMIN1(ALPHAK,ALPHAE)
              CALPHK = CALPHA
              CALPHA = ZERO
              IF (ALPHAE.NE.ZERO) CALPHA = SUMXA1/SUMXA2**ALPHAE
              SUMXTE = DSQRT(CALPHA*CALPHK)*SUMXA1**ALPHAK-SUMXA0
              IF (ALPHAA.GE.1.2D0 .AND. ICONV.EQ.1) ICONV = 2
              IF (ALPHAA.GT.1.8D0) ICONV = 3
              IF (MPRMON.GE.4)  WRITE(LUMON,32001) ICONV, ALPHAE,  &
                                  CALPHA, CLIN0, ALPHAK, SUMXTE
32001         FORMAT(' ** ICONV: ',I1,'  ALPHA: ',D9.2, &
                     '  CONST-ALPHA: ',D9.2,'  CONST-LIN: ',D9.2,' **', &
                     /,' **',11X,'ALPHA-POST: ',D9.2,' CHECK: ',D9.2, &
                     25X,'**')
              IF ( ICONV.GE.2 .AND. ALPHAA.LT.0.9D0 ) THEN
                 IF (IORMON.EQ.3) THEN
                   IERR = 4
                   GOTO 4299
                 ELSE
                   QMSTOP = .TRUE.
                 ENDIF
              ENDIF
            ENDIF
          ENDIF
          FCMON = FC

!         --------------------------------------------------------
!         3.4 Save natural level for later computations of
!             corrector and print iterate
          FCNUMK = SUMX
          IF (MPRMON.GE.2) THEN
            IF (MPRTIM.NE.0) CALL MONON(5)
            CALL N1PRV1(DLEVF,DLEVXA,FCKEEP,NITER,NEW,MPRMON,LUMON, &
                        QMIXIO)
            IF (MPRTIM.NE.0) CALL MONOFF(5)
          ENDIF
          NRED = 0
          QNEXT = .FALSE.
          QREP  = .FALSE.
!         QREP = ITER .GT. ITMAX   or  QREP = ITER .GT. 0
!
!         Damping-factor reduction loop
!         ================================
!         DO (Until)
34        CONTINUE
!           ------------------------------------------------------
!           3.5 Preliminary new iterate
            DO 35 L1=1,N
              X(L1)=XA(L1)+DX(L1)*FC
35          CONTINUE
!           -----------------------------------------------------
!           3.5.2 Exit, if problem is specified as being linear
!     Exit Repeat If ...
            IF( NONLIN.EQ.1 )THEN
              IERR = 0
              GOTO 4299
            ENDIF
!           ------------------------------------------------------
!           3.6.1 Computation of the residual vector
            IF (MPRTIM.NE.0) CALL MONON(1)
            CALL FCN(N,X,F,IFAIL)
            IF (MPRTIM.NE.0) CALL MONOFF(1)
            NFCN = NFCN+1
!     Exit, If ...
            IF(IFAIL.LT.0)THEN
              IERR = 82
              GOTO 4299
            ENDIF
            IF(IFAIL.EQ.1 .OR. IFAIL.EQ.2) THEN
              IF (QORDI) THEN
                IERR = 82
                GOTO 4299
              ENDIF
              IF (IFAIL.EQ.1) THEN
                FCREDU = HALF
              ELSE
                FCREDU = F(1)
!     Exit, If ...
                IF (FCREDU.LE.0 .OR. FCREDU.GE.1) THEN
                  IERR = 82
                  GOTO 4299
                ENDIF
              ENDIF
              IF (MPRMON.GE.2) THEN
36101           FORMAT(8X,I2,' FCN could not be evaluated  ', &
                       8X,F7.5,4X,I2)
                WRITE(LUMON,36101)NITER,FC,NEW
              ENDIF
              FCH = FC
              FC = FCREDU*FC
              IF (FCH.GT.FCMIN) FC = DMAX1(FC,FCMIN)
              IF (QBDAMP) THEN
                FCBH = FCH/FCBND
                IF (FC.LT.FCBH) THEN
                  FC = FCBH
                  IF (MPRMON.GE.4) WRITE(LUMON,*) &
                     ' *** decr. rest. act. (FCN redu.) ***'
                ENDIF
              ENDIF
              IF (FC.LT.FCMIN) THEN
                IERR = 3
                GOTO 4299
              ENDIF
!     Break DO (Until) ...
              GOTO 3109
            ENDIF
            IF (QORDI) THEN
!             -------------------------------------------------------
!             3.6.2 Convergence test for ordinary Newton iteration
!     Exit Repeat If ...
              IF( DXANRM.LE.RTOL )THEN
                IERR = 0
                GOTO 4299
              ENDIF
            ELSE
              DO 361 L1=1,N
               T1(L1)=F(L1)*FW(L1)
361           CONTINUE
!             --------------------------------------------------
!             3.6.3 Solution of linear (N,N)-system
              IF (MPRTIM.NE.0) CALL MONON(4)
              CALL N1SOLV(N,M1,ML,MU,A,T1,IOPT,IFAIL,LIWL,IWK(NIWLA),  &
                         IDUMMY,LRWL,RWK(NRWLA),IDUMMY)
              IF (MPRTIM.NE.0) CALL MONOFF(4)
!     Exit Repeat If ...
              IF(IFAIL.NE.0)  THEN
                IERR = 81
                GOTO 4299
              ENDIF
              IF(NEW.GT.0) THEN
                DO 3630 I=1,N
                  DXQ(I) = T1(I)*XWA(I)
3630            CONTINUE
                DO 363 ILOOP=1,NEW
                  SUM1=ZERO
                  SUM2=ZERO
                  DO 3631 I=1,N
                    SUM1=SUM1+(DXQ(I)*DXSAVE(I,ILOOP))/ XW(I)**2
                    SUM2=SUM2+(DXSAVE(I,ILOOP)/XW(I))**2
3631              CONTINUE
                  BETA=SUM1/SUM2
                  DO 3632 I=1,N
                    DXQ(I)=DXQ(I)+BETA*DXSAVE(I,ILOOP+1)
                    T1(I) = DXQ(I)/XW(I)
3632              CONTINUE
363             CONTINUE
              ENDIF
!             ----------------------------------------------------
!             3.6.4 Evaluation of scaled natural level function
!            SUMX
!            scaled maximum error norm CONV and evaluation
!            of (scaled) standard level function DLEVFN
              IF (.NOT. QSIMPL) THEN
                CALL N1LVLS(N,T1,XW,F,DXQ,CONV,SUMX,DLEVFN,MPRMON, &
                            NEW.EQ.0)
              ELSE
                CALL N1LVLS(N,T1,XWA,F,DXQ,CONV,SUMX,DLEVFN,MPRMON, &
                            NEW.EQ.0)
              ENDIF
              DXNRM = WNORM(N,DXQ,XW)
!             -----------------------------------------------------
!             3.6.5 Convergence test
!     Exit Repeat If ...
              IF ( DXNRM.LE.RTOL .AND. DXANRM.LE.RSMALL .AND.  &
                  FC.EQ.ONE ) THEN
                IERR = 0
                GOTO 4299
              ENDIF
!
              FCA = FC
!             ----------------------------------------------------
!             3.6.6 Evaluation of reduced damping factor
              TH = FCA-ONE
              FCDNM = ZERO
              DO 39 L1=1,N
                FCDNM=FCDNM+((DXQ(L1)+TH*DX(L1))/XW(L1))**2
39            CONTINUE
              IF (FCDNM.NE.ZERO) THEN
                DMYCOR = FCA*FCA*HALF*DSQRT(FCNUMK/FCDNM)
              ELSE
                DMYCOR = 1.0D+35
              ENDIF
              IF (NONLIN.LE.3) THEN
                FCCOR = DMIN1(ONE,DMYCOR)
              ELSE
                FCCOR = DMIN1(ONE,HALF*DMYCOR)
              ENDIF
!$Test-begin
              IF (MPRMON.GE.5) THEN
                WRITE(LUMON,39001) FCCOR, FC, DMYCOR, FCNUMK,FCDNM, FCA
39001           FORMAT (/, ' +++ corrector computation +++', /, &
                      ' FCCOR  = ', D18.10, '  F!     = ', D18.10, /, &
                      ' DMYCOR = ', D18.10, '  FCNUMK = ', D18.10, /, &
                      ' FCDNM  = ', D18.10, '  FCA    = ', D18.10, /, &
                         ' +++++++++++++++++++++++++++++', /)
              ENDIF
!$Test-end
            ENDIF
!           ------------------------------------------------------
!           3.7 Natural monotonicity test
            IF(SUMX.GT.SUMXA .AND. .NOT.QORDI)THEN
!             ----------------------------------------------------
!             3.8 Output of iterate
              IF(MPRMON.GE.3) THEN
                IF (MPRTIM.NE.0) CALL MONON(5)
                CALL N1PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC, &
                            NITER,MPRMON,LUMON,QMIXIO,'*')
                IF (MPRTIM.NE.0) CALL MONOFF(5)
              ENDIF
              IF (QMSTOP) THEN
                IERR = 4
                GOTO 4299
              ENDIF
              FCH = DMIN1(FCCOR,HALF*FC)
              IF (FC.GT.FCMIN) THEN
                FC=DMAX1(FCH,FCMIN)
              ELSE
                FC=FCH
              ENDIF
              IF (QBDAMP) THEN
                FCBH = FCA/FCBND
                IF (FC.LT.FCBH) THEN
                  FC = FCBH
                  IF (MPRMON.GE.4) &
                    WRITE(LUMON,*) ' *** decr. rest. act. (a post) ***'
                ENDIF
              ENDIF
!WEI
              FCMON = FC
!$Test-begin
                IF (MPRMON.GE.5) THEN
                  WRITE(LUMON,39002) FC
39002             FORMAT (/, ' +++ corrector setting 1 +++', /, &
                          ' F!     = ', D18.10, /, &
                             ' +++++++++++++++++++++++++++', /)
                ENDIF
!$Test-end
              QREP = .TRUE.
              NCORR = NCORR+1
              NRED = NRED+1
!             ----------------------------------------------------
!             3.10 If damping factor is too small:
!           Refresh Jacobian,if current Jacobian was computed
!           by a Rank1-update, else fail exit
              QJCRFR  = FC.LT.FCMIN.OR.NEW.GT.0.AND.NRED.GT.1
!     Exit Repeat If ...
              IF(QJCRFR.AND.NEW.EQ.0)THEN
                IERR = 3
                GOTO 4299
              ENDIF
            ELSE
              IF (.NOT.QORDI.AND..NOT.QREP .AND. FCCOR.GT.SIGMA2*FC) &
              THEN
                IF(MPRMON.GE.3) THEN
                  IF (MPRTIM.NE.0) CALL MONON(5)
                  CALL N1PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC, &
                              NITER,MPRMON,LUMON,QMIXIO,'+')
                  IF (MPRTIM.NE.0) CALL MONOFF(5)
                ENDIF
                FC = FCCOR
!$Test-begin
                IF (MPRMON.GE.5) THEN
                  WRITE(LUMON,39003) FC
39003             FORMAT (/, ' +++ corrector setting 2 +++', /, &
                          ' F!     = ', D18.10, /, &
                          ' +++++++++++++++++++++++++++', /)
                ENDIF
!$Test-end
                QREP = .TRUE.
              ELSE
                QNEXT = .TRUE.
              ENDIF
            ENDIF
3109      CONTINUE
          IF(.NOT.(QNEXT.OR.QJCRFR)) GOTO  34
!         UNTIL ( expression - negated above)
!         End of damping-factor reduction loop
!         =======================================
          IF(QJCRFR)THEN
!           ------------------------------------------------------
!           3.11 Restore former values for repeting iteration
!         step
            NREJR1 = NREJR1+1
            DO 3111 L1=1,N
              X(L1)=XA(L1)
3111        CONTINUE
            DO 3112 L1=1,N
              F(L1)=FA(L1)
3112        CONTINUE
            IF(MPRMON.GE.2)THEN
31130           FORMAT(8X,I2,' Not accepted damping factor ', &
                       8X,F7.5,4X,I2)
                WRITE(LUMON,31130)NITER,FC,NEW
            ENDIF
            FC = FCKEEP
            FCA = FCK2
            IF(NITER.EQ.0)THEN
              FC = FCMIN
            ENDIF
            QGENJ = .TRUE.
          ELSE
!           ------------------------------------------------------
!           4 Preparations to start the following iteration step
!           ------------------------------------------------------
!           4.1 Print values
            IF(MPRMON.GE.3 .AND. .NOT.QORDI) THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL N1PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC,NITER+1, &
                          MPRMON,LUMON,QMIXIO,'*')
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ENDIF
!           Print the natural level of the current iterate and return
!           it in one-step mode
            SUMXS = SUMX
            SUMX = SUMXA
            IF(MPRSOL.GE.2.AND.NITER.NE.0) THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL N1SOUT(N,XA,2,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ELSE IF(MPRSOL.GE.1.AND.NITER.EQ.0)THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL N1SOUT(N,XA,1,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ENDIF
            NITER = NITER+1
            DLEVF = DLEVFN
!     Exit Repeat If ...
            IF(NITER.GE.NITMAX)THEN
              IERR = 2
              GOTO 4299
            ENDIF
            FCKEEP = FC
!           ------------------------------------------------------
!           4.2 Return, if in one-step mode
!
! Exit Subroutine If ...
            IF (MODE.EQ.1) THEN
              IWK(18)=NIWLA-1
              IWK(19)=NRWLA-1
              IOPT(1)=1
              RETURN
            ENDIF
          ENDIF
        GOTO 2
!       End Repeat
4299    CONTINUE
!       End of main iteration loop
!       ==========================
!       ----------------------------------------------------------
!       9 Exits
!       ----------------------------------------------------------
!       9.1 Solution exit
        APREC = -1.0D0

        IF(IERR.EQ.0 .OR. IERR.EQ.4)THEN
          IF (NONLIN.NE.1) THEN
            IF (.NOT.QORDI) THEN
              IF ( IERR.EQ.0 ) THEN
                APREC = DSQRT(SUMX/DBLE(FLOAT(N)))
                DO 91 L1=1,N
                  X(L1)=X(L1)+DXQ(L1)
91              CONTINUE
              ELSE
                APREC = DSQRT(SUMXA/DBLE(FLOAT(N)))
                IF (ALPHAA.GT.ZERO .AND. IORMON.EQ.3) THEN
                  DO 92 L1=1,N
                    X(L1)=X(L1)+DX(L1)
92                CONTINUE
                ENDIF
              ENDIF
!             Print final monitor output
              IF(MPRMON.GE.2) THEN
                IF (IERR.EQ.0) THEN
                  IF (MPRTIM.NE.0) CALL MONON(5)
                  CALL N1PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC, &
                              NITER+1,MPRMON,LUMON,QMIXIO,'*')
                  IF (MPRTIM.NE.0) CALL MONOFF(5)
                ELSE IF (IORMON.EQ.3) THEN
                  IF (MPRTIM.NE.0) CALL MONON(5)
                  CALL N1PRV1(DLEVFN,DSQRT(SUMXA/DBLE(FLOAT(N))),FC, &
                              NITER,NEW,MPRMON,LUMON,QMIXIO)
                  IF (MPRTIM.NE.0) CALL MONOFF(5)
                ENDIF
              ENDIF
              IF (  IORMON.GE.2 ) THEN
                IF ( ICONV.LE.1 .AND. ALPHAE .NE. ZERO  &
                                .AND. ALPHAK .NE. ZERO ) IERR = 5
              ENDIF
            ELSE
!           IF (QORDI) THEN
              APREC = DSQRT(SUMXA/DBLE(FLOAT(N)))
            ENDIF
            IF(MPRMON.GE.1) THEN
91001         FORMAT(///' Solution of nonlinear system ', &
              'of equations obtained within ',I3, &
              ' iteration steps',//,' Achieved relative accuracy',D10.3)
              IF (QORDI .OR. IERR.EQ.4) THEN
                WRITE(LUMON,91001) NITER,APREC
              ELSE
                WRITE(LUMON,91001) NITER+1,APREC
              ENDIF
            ENDIF
          ELSE
            IF(MPRMON.GE.1) THEN
91002         FORMAT(///' Solution of linear system ', &
              'of equations obtained by NLEQ1',//,' No estimate ', &
              'available for the achieved relative accuracy')
                WRITE(LUMON,91002)
            ENDIF
          ENDIF
        ENDIF
!       ----------------------------------------------------------
!       9.2 Fail exit messages
!       ----------------------------------------------------------
!       9.2.1 Termination, since jacobian matrix became singular
        IF(IERR.EQ.1.AND.MPRERR.GE.1)THEN
92101     FORMAT(/,' Iteration terminates due to ', &
          'singular jacobian matrix',/)
          WRITE(LUERR,92101)
        ENDIF
!       ----------------------------------------------------------
!       9.2.2 Termination after more than NITMAX iterations
        IF(IERR.EQ.2.AND.MPRERR.GE.1)THEN
92201     FORMAT(/,' Iteration terminates after NITMAX ', &
          '=',I3,'  Iteration steps')
          WRITE(LUERR,92201)NITMAX
        ENDIF
!       ----------------------------------------------------------
!       9.2.3 Damping factor FC became too small
        IF(IERR.EQ.3.AND.MPRERR.GE.1)THEN
92301     FORMAT(/,' Damping factor has become too ', &
          'small: lambda =',D10.3,2X,/)
          WRITE(LUERR,92301)FC
        ENDIF
!WEI
!       ----------------------------------------------------------
!       9.2.4.1 Superlinear convergence slowed down
        IF(IERR.EQ.4.AND.MPRERR.GE.1)THEN
92401     FORMAT(/,' Warning: Monotonicity test failed after ',A, &
                 ' convergence was already checked;',/, &
          ' RTOL requirement may be too stringent',/)
92402     FORMAT(/,' Warning: ',A,' convergence slowed down;',/, &
          ' RTOL requirement may be too stringent',/)
          IF (QMSTOP) THEN
            IF (ICONV.EQ.2) WRITE(LUERR,92401) 'superlinear'
            IF (ICONV.EQ.3) WRITE(LUERR,92401) 'quadratic'
          ELSE
            IF (ICONV.EQ.2) WRITE(LUERR,92402) 'superlinear'
            IF (ICONV.EQ.3) WRITE(LUERR,92402) 'quadratic'
          ENDIF
        ENDIF
!       ----------------------------------------------------------
!       9.2.4.2 Convergence criterion satisfied before superlinear
!        convergence has been established
! This suppresses spurious warnings in FEmapper
        IF(IERR.EQ.5.AND.MPRERR.GT.1)THEN
!        IF(IERR.EQ.5.AND.MPRERR.GE.1)THEN
          write(*,*) 'MPRERR: ',MPRERR
92410     FORMAT(/,' Warning: No quadratic or superlinear convergence ', &
                 'established yet',/, &
                 10X,'your solution may perhaps may be less accurate ', &
                 /,10X,'as indicated by the standard error estimate')
          WRITE(LUERR,92410)
        ENDIF
!       ----------------------------------------------------------
!       9.2.5 Error exit due to linear solver routine N1FACT
        IF(IERR.EQ.80.AND.MPRERR.GE.1)THEN
92501     FORMAT(/,' Error ',I5,' signalled by linear solver N1FACT')
          WRITE(LUERR,92501) IFAIL
        ENDIF
!       ----------------------------------------------------------
!       9.2.6 Error exit due to linear solver routine N1SOLV
        IF(IERR.EQ.81.AND.MPRERR.GE.1)THEN
92601     FORMAT(/,' Error ',I5,' signalled by linear solver N1SOLV')
          WRITE(LUERR,92601) IFAIL
        ENDIF
!       ----------------------------------------------------------
!       9.2.7 Error exit due to fail of user function FCN
        IF(IERR.EQ.82.AND.MPRERR.GE.1)THEN
92701     FORMAT(/,' Error ',I5,' signalled by user function FCN')
          WRITE(LUERR,92701) IFAIL
        ENDIF
!       ----------------------------------------------------------
!       9.2.8 Error exit due to fail of user function JAC
        IF(IERR.EQ.83.AND.MPRERR.GE.1)THEN
92801     FORMAT(/,' Error ',I5,' signalled by user function JAC')
          WRITE(LUERR,92801) IFAIL
        ENDIF
        IF(IERR.GE.80.AND.IERR.LE.83) IWK(23) = IFAIL
        IF ((IERR.EQ.82.OR.IERR.EQ.83).AND.NITER.LE.1.AND.MPRERR.GE.1) &
        THEN
          WRITE (LUERR,92810)
92810     FORMAT(' Try to find a better initial guess for the solution')
        ENDIF
!       ----------------------------------------------------------
!       9.3 Common exit
        IF (MPRERR.GE.3.AND.IERR.NE.0.AND.IERR.NE.4.AND.NONLIN.NE.1) &
          THEN
93100     FORMAT(/,'    Achieved relative accuracy',D10.3,2X)
          WRITE(LUERR,93100)CONVA
          APREC = CONVA
        ENDIF
        RTOL = APREC
        SUMX = SUMXA
        IF(MPRSOL.GE.2.AND.NITER.NE.0) THEN
           MODE=2
           IF (QORDI) MODE=3
           IF (MPRTIM.NE.0) CALL MONON(5)
           CALL N1SOUT(N,XA,MODE,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
           IF (MPRTIM.NE.0) CALL MONOFF(5)
        ELSE IF(MPRSOL.GE.1.AND.NITER.EQ.0)THEN
           IF (MPRTIM.NE.0) CALL MONON(5)
           CALL N1SOUT(N,XA,1,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
           IF (MPRTIM.NE.0) CALL MONOFF(5)
        ENDIF
        IF (.NOT.QORDI) THEN
          IF (IERR.NE.4) NITER = NITER+1
          DLEVF = DLEVFN
          IF(MPRSOL.GE.1)THEN
!           Print Solution or final iteration vector
            IF(IERR.EQ.0)THEN
               MODEFI = 3
            ELSE
               MODEFI = 4
            ENDIF
            IF (MPRTIM.NE.0) CALL MONON(5)
            CALL N1SOUT(N,X,MODEFI,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
            IF (MPRTIM.NE.0) CALL MONOFF(5)
          ENDIF
        ENDIF
!       Return the latest internal scaling to XSCAL
        DO 93 I=1,N
          XSCAL(I)=XW(I)
93      CONTINUE
!       End of exits
!       End of subroutine N1INT
      RETURN
      END subroutine

      SUBROUTINE N1SCAL(N,X,XA,XSCAL,XW,ISCAL,QINISC,IOPT,LRWK,RWK)
!     Begin Prologue SCALE
      INTEGER N
      DOUBLE PRECISION X(N),XSCAL(N),XA(N),XW(N)
      INTEGER ISCAL
      LOGICAL QINISC
      INTEGER IOPT(50),LRWK
      DOUBLE PRECISION RWK(LRWK)
!     ------------------------------------------------------------
!
!     Summary :
!
!     S C A L E : To be used in connection with NLEQ1 .
!       Computation of the internal scaling vector XW used for the
!       Jacobian matrix, the iterate vector and it's related
!       vectors - especially for the solution of the linear system
!       and the computations of norms to avoid numerical overflow.
!
!     Input parameters
!     ================
!
!     N         Int     Number of unknowns
!     X(N)      Dble    Current iterate
!     XA(N)     Dble    Previous iterate
!     XSCAL(N)  Dble    User scaling passed from parameter XSCAL
!                of interface routine NLEQ1
!     ISCAL     Int     Option ISCAL passed from IOPT-field
!                (for details see description of IOPT-fields)
!     QINIS!    Logical = .TRUE.  : Initial scaling
!                = .FALSE. : Subsequent scaling
!     IOPT(50)  Int     Options array passed from NLEQ1 parameter list
!     LRWK      Int     Length of real workspace
!     RWK(LRWK) Dble    Real workspace (see description above)
!
!     Output parameters
!     =================
!
!     XW(N)     Dble   Scaling vector computed by this routine
!               All components must be positive. The follow-
!               ing relationship between the original vector
!               X and the scaled vector XSCAL holds:
!               XSCAL(I) = X(I)/XW(I) for I=1,...N
!
!     Subroutines called: ZIBCONST
!
!     Machine constants used
!     ======================
!
      DOUBLE PRECISION EPMACH, SMALL
!
!     ------------------------------------------------------------
!     End Prologue
      INTRINSIC DABS,DMAX1
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      INTEGER MPRMON,LUMON, L1
      CALL ZIBCONST(EPMACH, SMALL)
!     Begin
      DO 1 L1=1,N
        IF (ISCAL.EQ.1) THEN
          XW(L1) = XSCAL(L1)
        ELSE
          XW(L1)=DMAX1(XSCAL(L1),(DABS(X(L1))+DABS(XA(L1)))*HALF,SMALL)
        ENDIF
1     CONTINUE
!$Test-Begin
      MPRMON = IOPT(13)
      IF (MPRMON.GE.6) THEN
        LUMON = IOPT(14)
        WRITE(LUMON,*) ' '
        WRITE(LUMON,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        WRITE(LUMON,*) '      X-components   Scaling-components    '
        WRITE(LUMON,10) (X(L1), XW(L1), L1=1,N)
10      FORMAT('  ',D18.10,'  ',D18.10)
        WRITE(LUMON,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        WRITE(LUMON,*) ' '
      ENDIF
!$Test-End
!     End of subroutine N1SCAL
      RETURN
      END subroutine

      SUBROUTINE N1SCRF(M,N,A,FW)
!     Begin Prologue SCROWF
      INTEGER M,N
      DOUBLE PRECISION A(M,N),FW(M)
!     ------------------------------------------------------------
!
!     Summary :
!
!     S C R O W F : Row Scaling of a (M,N)-matrix in full storage
!            mode
!
!     Input parameters (* marks inout parameters)
!     ===========================================
!
!       M           Int    Number of rows of the matrix
!       N           Int    Number of columns of the matrix
!     * A(M,N)      Dble   Matrix to be scaled
!
!     Output parameters
!     =================
!
!       FW(M)       Dble   Row scaling factors - FW(i) contains
!                   the factor by which the i-th row of A
!                   has been multiplied
!
!     ------------------------------------------------------------
!     End Prologue
      INTRINSIC DABS
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER J,K
      DOUBLE PRECISION S1,S2
!     Begin
      DO 1 K=1,M
        S1=ZERO
        DO 2 J=1,N
          S2=DABS(A(K,J))
          IF (S2.GT.S1) S1=S2
2       CONTINUE
        IF (S1.GT.ZERO) THEN
          S1=ONE/S1
          FW(K)=S1
          DO 3 J=1,N
            A(K,J)=A(K,J)*S1
3         CONTINUE
        ELSE
          FW(K)=ONE
        ENDIF
1     CONTINUE
!     End of subroutine N1SCRF
      RETURN
      END subroutine

      SUBROUTINE N1SCRB(N,LDA,ML,MU,A,FW)
!     Begin Prologue SCROWB
      INTEGER N,LDA,ML,MU
      DOUBLE PRECISION A(LDA,N),FW(N)
!     ------------------------------------------------------------
!
!     Summary :
!
!     S C R O W B : Row Scaling of a (N,N)-matrix in band storage
!            mode
!
!     Input parameters (* marks inout parameters)
!     ===========================================
!
!       N           Int    Number of rows and columns of the matrix
!       LDA         Int    Leading dimension of the matrix array
!       ML          Int    Lower bandwidth of the matrix (see IOPT)
!       MU          Int    Upper bandwidth of the matrix (see IOPT)
!     * A(LDA,N)    Dble   Matrix to be scaled (see Note in routine
!                   DGBFA for details on storage technique)
!
!     Output parameters
!     =================
!
!       FW(N)       Dble   Row scaling factors - FW(i) contains
!                   the factor by which the i-th row of
!                   the matrix has been multiplied
!
!     ------------------------------------------------------------
!     End Prologue
      INTRINSIC DABS,MAX0,MIN0
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER K,K1,L1,L2,L3,M2
      DOUBLE PRECISION S1,S2
!     Begin
      M2=ML+MU+1
      DO 1 K=1,N
        S1=ZERO
        L2=MAX0(1,K-ML)
        L3=MIN0(N,K+MU)
        K1=M2+K
        DO 2 L1=L2,L3
          S2=DABS(A(K1-L1,L1))
          IF (S2.GT.S1) S1=S2
2       CONTINUE
        IF (S1.GT.ZERO) THEN
          S1=ONE/S1
          FW(K)=S1
          DO 3 L1=L2,L3
            A(K1-L1,L1)=A(K1-L1,L1)*S1
3         CONTINUE
        ELSE
          FW(K)=ONE
        ENDIF
1     CONTINUE
!     End of subroutine N1SCRB
      RETURN
      END subroutine

      SUBROUTINE N1FACT(N,LDA,ML,MU,A,IOPT,IFAIL,LIWK,IWK,LAIWK,LRWK, &
      RWK,LARWK)
!     Begin Prologue FACT
      INTEGER N,LDA,ML,MU
      DOUBLE PRECISION A(LDA,N)
      INTEGER IOPT(50)
      INTEGER IFAIL
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LAIWK,LRWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER LARWK
      integer :: MSTOR, MPRERR, LUERR
!     ------------------------------------------------------------
!
!     Summary :
!
!     F A C T : Call linear algebra subprogram for factorization of
!        a (N,N)-matrix
!
!     Input parameters (* marks inout parameters)
!     ===========================================
!
!     N             Int    Order of the linear system
!     LDA           Int    Leading dimension of the matrix array A
!     ML            Int    Lower bandwidth of the matrix (only for
!                   banded systems)
!     MU            Int    Upper bandwidth of the matrix (only for
!                   banded systems)
!    * A(LDA,N)      Dble   Matrix storage. See main routine NLEQ1,
!                   note 4 for how to store a banded Jacobian
!     IOPT(50)      Int    Option vector passed from NLEQ1
!
!     Output parameters
!     =================
!
!     IFAIL         Int    Error indicator returned by this routine:
!                   = 0 matrix decomposition successfull
!                   = 1 decomposition failed -
!                              matrix is numerically singular
!                   =10 supplied (integer) workspace too small
!
!     Workspace parameters
!     ====================
!
!     LIWK          Int    Length of integer workspace passed to this
!                   routine (In)
!     IWK(LIWK)     Int    Integer Workspace supplied for this routine
!     LAIWK         Int    Length of integer Workspace used by this
!                   routine (out)
!     LRWK          Int    Length of real workspace passed to this
!                   routine (In)
!     RWK(LRWK)     Dble   Real Workspace supplied for this routine
!     LARWK         Int    Length of real Workspace used by this
!                   routine (out)
!
!     Subroutines called:  DGEFA, DGBFA
!
!     ------------------------------------------------------------
!     End Prologue
!      EXTERNAL DGEFA, DGBFA
!     Begin
      LAIWK = N
      LARWK = 0
!     ( Real workspace starts at RWK(NRWKFR) )
!     IF (LIWK.GE.LAIWK.AND.LRWK.GE.LARWK) THEN
      IF (LIWK.GE.LAIWK) THEN
        MSTOR = IOPT(4)
        IF (MSTOR.EQ.0) THEN
          CALL DGEFA(A,LDA,N,IWK,IFAIL)
        ELSE IF (MSTOR.EQ.1) THEN
          CALL DGBFA(A,LDA,N,ML,MU,IWK,IFAIL)
        ENDIF
        IF(IFAIL.NE.0)THEN
          IFAIL = 1
        ENDIF
      ELSE
        IFAIL = 10
        MPRERR=IOPT(11)
        LUERR=IOPT(12)
10001   FORMAT(/,' Insuffient workspace for linear solver,', &
               ' at least needed more needed : ',/, &
               ' ',A,' workspace : ',I4)
        IF (LIWK.LT.LAIWK.AND.MPRERR.GT.0)  &
          WRITE(LUERR,10001) 'Integer',LAIWK-LIWK
      ENDIF
      RETURN
      END subroutine

      SUBROUTINE N1SOLV(N,LDA,ML,MU,A,B,IOPT,IFAIL,LIWK,IWK,LAIWK, &
      LRWK,RWK,LARWK)
!     Begin Prologue SOLVE
      INTEGER N,LDA,ML,MU
      DOUBLE PRECISION A(LDA,N)
      DOUBLE PRECISION B(N)
      INTEGER IOPT(50)
      INTEGER IFAIL
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK,LAIWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER LARWK
      integer :: MSTOR
!     ------------------------------------------------------------
!
!     Summary :
!
!     S O L V E : Call linear algebra subprogram for solution of
!          the linear system A*Z = B
!
!     Parameters
!     ==========
!
!     N,LDA,ML,MU,A,IOPT,IFAIL,LIWK,IWK,LAIWK,LRWK,RWK,LARWK :
!                 See description for subroutine N1FACT.
!     B          Dble    In:  Right hand side of the linear system
!                 Out: Solution of the linear system
!
!     Subroutines called: DGESL, DGBSL
!
!     ------------------------------------------------------------
!     End Prologue
!      EXTERNAL DGESL, DGBSL
!     Begin
      MSTOR = IOPT(4)
      IF (MSTOR.EQ.0) THEN
        CALL DGESL(A,LDA,N,IWK,B,0)
      ELSE IF (MSTOR.EQ.1) THEN
        CALL DGBSL(A,LDA,N,ML,MU,IWK,B,0)
      ENDIF
      IFAIL = 0
      RETURN
      END subroutine

      SUBROUTINE N1LVLS(N,DX1,XW,F,DXQ,CONV,SUMX,DLEVF,MPRMON,QDSCAL)
!     Begin Prologue LEVELS
      INTEGER N,MPRMON
      DOUBLE PRECISION DX1(N),XW(N),F(N),DXQ(N)
      DOUBLE PRECISION CONV,SUMX,DLEVF
      LOGICAL QDSCAL
!
!     ------------------------------------------------------------
!
!     Summary :
!
!     L E V E L S : To be used in connection with NLEQ1 .
!     provides descaled solution, error norm and level functions
!
!     Input parameters (* marks inout parameters)
!     ===========================================
!
!       N              Int    Number of parameters to be estimated
!       DX1(N)         Dble   array containing the scaled Newton
!                             correction
!       XW(N)          Dble   Array containing the scaling values
!       F(N)           Dble   Array containing the residuum
!
!     Output parameters
!     =================
!
!       DXQ(N)         Dble   Array containing the descaled Newton
!                             correction
!       CONV           Dble   Scaled maximum norm of the Newton
!                             correction
!       SUMX           Dble   Scaled natural level function value
!       DLEVF          Dble   Standard level function value (only
!                             if needed for print)
!       MPRMON         Int    Print information parameter (see
!                             driver routine NLEQ1 )
!       QDSCAL         Logic  .TRUE., if descaling of DX1 required,
!                             else .FALSE.
!
!     ------------------------------------------------------------
!     End Prologue
      INTRINSIC DABS
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER L1
      DOUBLE PRECISION S1
!     Begin
      IF (QDSCAL) THEN
!       ------------------------------------------------------------
!       1.2 Descaling of solution DX1 ( stored to DXQ )
        DO 12 L1=1,N
          DXQ(L1)=DX1(L1)*XW(L1)
12      CONTINUE
      ENDIF
!     ------------------------------------------------------------
!     2 Evaluation of scaled natural level function SUMX and
!       scaled maximum error norm CONV
      CONV = ZERO
      DO 20 L1=1,N
        S1 = DABS(DX1(L1))
        IF(S1.GT.CONV) CONV=S1
20    CONTINUE
      SUMX = ZERO
      DO 21 L1=1,N
        SUMX = SUMX+DX1(L1)**2
21    CONTINUE
!     ------------------------------------------------------------
!     3 Evaluation of (scaled) standard level function DLEVF (only
!       if needed for print)
      IF(MPRMON.GE.2)THEN
        DLEVF = ZERO
        DO 3 L1=1,N
          DLEVF = DLEVF+F(L1)**2
3       CONTINUE
        DLEVF = DSQRT(DLEVF/DBLE(FLOAT(N)))
      ENDIF
!     End of subroutine N1LVLS
      RETURN
      END subroutine

      SUBROUTINE N1JAC (FCN, N, LDA, X, FX, A, YSCAL, AJDEL, AJMIN, &
                        NFCN, FU, IFAIL)
!* Begin Prologue N1JAC
!      EXTERNAL FCN
      procedure(func) :: FCN
      INTEGER N, LDA
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), AJDEL, AJMIN
      INTEGER NFCN
      DOUBLE PRECISION FU(N)
      INTEGER IFAIL
!
!  ---------------------------------------------------------------------
!
!* Title
!
!  Evaluation of a dense Jacobian matrix using finite difference
!  approximation adapted for use in nonlinear systems solver NLEQ1
!
!* Environment       Fortran 77
!             Double Precision
!             Sun 3/60, Sun OS
!* Latest Revision   May 1990
!
!
!* Parameter list description
!  --------------------------
!
!* External subroutines (to be supplied by the user)
!  -------------------------------------------------
!
!  FCN        Ext     FCN (N, X, FX, IFAIL)
!              Subroutine in order to provide right-hand
!              side of first-order differential equations
!    N        Int     Number of rows and columns of the Jacobian
!    X(N)     Dble    The current scaled iterates
!    FX(N)    Dble    Array containing FCN(X)
!    IFAIL    Int     Return code
!              Whenever a negative value is returned by FCN
!              routine N1JAC is terminated immediately.
!
!
!* Input parameters (* marks inout parameters)
!  ----------------
!
!  N          Int     Number of rows and columns of the Jacobian
!  LDA        Int     Leading Dimension of array A
!  X(N)       Dble    Array containing the current scaled
!              iterate
!  FX(N)      Dble    Array containing FCN(X)
!  YSCAL(N)   Dble    Array containing the scaling factors
!  AJDEL      Dble    Perturbation of component k: abs(Y(k))*AJDEL
!  AJMIN      Dble    Minimum perturbation is AJMIN*AJDEL
!  NFCN       Int  *  FCN - evaluation count
!
!* Output parameters (* marks inout parameters)
!  -----------------
!
!  A(LDA,N)   Dble    Array to contain the approximated
!              Jacobian matrix ( dF(i)/dx(j)in A(i,j))
!  NFCN       Int  *  FCN - evaluation count adjusted
!  IFAIL      Int     Return code non-zero if Jacobian could not
!              be computed.
!
!* Workspace parameters
!  --------------------
!
!  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
!              the numerator differences
!
!* Called
!  ------
!
      INTRINSIC DABS, DMAX1, DSIGN
!  ---------------------------------------------------------------------
!
!* End Prologue
!
!* Local variables
!  ---------------
!
      INTEGER I, K
      DOUBLE PRECISION U, W
!
!* Begin
!
      IFAIL = 0
      DO 1 K = 1,N
         W = X(K)
         U = DSIGN(DMAX1(DABS(X(K)),AJMIN,YSCAL(K))*AJDEL, X(K))
         X(K) = W + U

         CALL FCN (N, X, FU, IFAIL)
         NFCN = NFCN + 1
         IF (IFAIL .NE. 0) GOTO 99

         X(K) = W
         DO 11 I = 1,N
            A(I,K) = (FU(I) - FX(I)) / U
 11      CONTINUE
 1    CONTINUE

99    CONTINUE
      RETURN

!
!* End of N1JAC
!
      END subroutine

      SUBROUTINE N1JACB (FCN, N, LDA, ML, X, FX, A, YSCAL, AJDEL, AJMIN, &
           NFCN, FU, U, W, IFAIL)
!* Begin Prologue N1JACB
!      EXTERNAL FCN
      procedure(func) :: FCN
      INTEGER N, LDA, ML
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), AJDEL, AJMIN
      INTEGER NFCN
      DOUBLE PRECISION FU(N), U(N), W(N)
      INTEGER IFAIL
!
!  ---------------------------------------------------------------------
!
!* Title
!
!  Evaluation of a banded Jacobian matrix using finite difference
!  approximation adapted for use in nonlinear systems solver NLEQ1
!
!* Environment       Fortran 77
!             Double Precision
!             Sun 3/60, Sun OS
!* Latest Revision   May 1990
!
!
!* Parameter list description
!  --------------------------
!
!* External subroutines (to be supplied by the user)
!  -------------------------------------------------
!
!  FCN        Ext     FCN (N, X, FX, IFAIL)
!              Subroutine in order to provide right-hand
!              side of first-order differential equations
!    N        Int     Number of rows and columns of the Jacobian
!    X(N)     Dble    The current scaled iterates
!    FX(N)    Dble    Array containing FCN(X)
!    IFAIL    Int     Return code
!              Whenever a negative value is returned by FCN
!              routine N1JACB is terminated immediately.
!
!
!* Input parameters (* marks inout parameters)
!  ----------------
!
!  N          Int     Number of rows and columns of the Jacobian
!  LDA        Int     Leading Dimension of array A
!  ML         Int     Lower bandwidth of the Jacobian matrix
!  X(N)       Dble    Array containing the current scaled
!              iterate
!  FX(N)      Dble    Array containing FCN(X)
!  YSCAL(N)   Dble    Array containing the scaling factors
!  AJDEL      Dble    Perturbation of component k: abs(Y(k))*AJDEL
!  AJMIN      Dble    Minimum perturbation is AJMIN*AJDEL
!  NFCN       Int  *  FCN - evaluation count
!
!* Output parameters (* marks inout parameters)
!  -----------------
!
!  A(LDA,N)   Dble    Array to contain the approximated
!              Jacobian matrix ( dF(i)/dx(j)in A(i,j))
!  NFCN       Int  *  FCN - evaluation count adjusted
!  IFAIL      Int     Return code non-zero if Jacobian could not
!              be computed.
!
!* Workspace parameters
!  --------------------
!
!  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
!              the numerator differences
!  U(N)       Dble    Array to contain dx(i)
!  W(N)       Dble    Array to save original components of X
!
!* Called
!  ------
!
      INTRINSIC DABS, DMAX1, MAX0, MIN0
!
!* Constants
!  ---------
!
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0D0)
!
!  ---------------------------------------------------------------------
!
!* End Prologue
!
!* Local variables
!  ---------------
!
      INTEGER I, I1, I2, JJ, K, MH, MU, LDAB
!
!* Begin
!
      IFAIL = 0
      MU = LDA - 2*ML - 1
      LDAB = ML + MU + 1
      DO 1 I = 1,LDAB
         DO 11 K = 1,N
            A(I,K) = ZERO
 11      CONTINUE
 1    CONTINUE

      DO 2 JJ = 1,LDAB
         DO 21 K = JJ,N,LDAB
            W(K) = X(K)
            U(K) = DSIGN(DMAX1(DABS(X(K)),AJMIN,YSCAL(K))*AJDEL, X(K))
            X(K) = W(K) + U(K)
 21      CONTINUE

         CALL FCN (N, X, FU, IFAIL)
         NFCN = NFCN + 1
         IF (IFAIL .NE. 0) GOTO 99

         DO 22 K = JJ,N,LDAB
            X(K) = W(K)
            I1 = MAX0 (1, K - MU)
            I2 = MIN0 (N, K + ML)
            MH = MU + 1 - K
            DO 221 I = I1,I2
               A(MH+I,K) = (FU(I) - FX(I)) / U(K)
 221        CONTINUE
 22      CONTINUE

 2    CONTINUE

99    CONTINUE
      RETURN

!
!* End of N1JACB
!
      END subroutine

      SUBROUTINE N1JCF (FCN, N, LDA, X, FX, A, YSCAL, ETA, ETAMIN, &
           ETAMAX, ETADIF, CONV, NFCN, FU, IFAIL)
!* Begin Prologue N1JCF
!      EXTERNAL FCN
      procedure(func) :: FCN
      INTEGER N, LDA
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), ETA(N), &
           ETAMIN, ETAMAX, ETADIF, CONV
      INTEGER NFCN
      DOUBLE PRECISION FU(N)
      INTEGER IFAIL
!
!  ---------------------------------------------------------------------
!
!* Title
!
!  Approximation of dense Jacobian matrix for nonlinear systems solver
!  NLEQ1 with feed-back control of discretization and rounding errors
!
!* Environment       Fortran 77
!             Double Precision
!             Sun 3/60, Sun OS
!* Latest Revision   May 1990
!
!
!* Parameter list description
!  --------------------------
!
!* External subroutines (to be supplied by the user)
!  -------------------------------------------------
!
!  FCN        Ext     FCN (N, X, FX, IFAIL)
!              Subroutine in order to provide right-hand
!              side of first-order differential equations
!    N        Int     Number of rows and columns of the Jacobian
!    X(N)     Dble    The current scaled iterates
!    FX(N)    Dble    Array containing FCN(X)
!    IFAIL    Int     Return code
!              Whenever a negative value is returned by FCN
!              routine N1JCF is terminated immediately.
!
!
!* Input parameters (* marks inout parameters)
!  ----------------
!
!  N          Int     Number of rows and columns of the Jacobian
!  LDA        Int     Leading dimension of A (LDA .GE. N)
!  X(N)       Dble    Array containing the current scaled
!              iterate
!  FX(N)      Dble    Array containing FCN(X)
!  YSCAL(N)   Dble    Array containing the scaling factors
!  ETA(N)     Dble *  Array containing the scaled denominator
!              differences
!  ETAMIN     Dble    Minimum allowed scaled denominator
!  ETAMAX     Dble    Maximum allowed scaled denominator
!  ETADIF     Dble    DSQRT (1.1*EPMACH)
!              EPMACH = machine precision
!  CONV       Dble    Maximum norm of last (unrelaxed) Newton correction
!  NFCN       Int  *  FCN - evaluation count
!
!* Output parameters (* marks inout parameters)
!  -----------------
!
!  A(LDA,N)   Dble    Array to contain the approximated
!              Jacobian matrix ( dF(i)/dx(j)in A(i,j))
!  ETA(N)     Dble *  Scaled denominator differences adjusted
!  NFCN       Int  *  FCN - evaluation count adjusted
!  IFAIL      Int     Return code non-zero if Jacobian could not
!              be computed.
!
!* Workspace parameters
!  --------------------
!
!  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
!              the numerator differences
!
!* Called
!  ------
!
      INTRINSIC DABS, DMAX1, DMIN1, DSIGN, DSQRT
!
!* Constants
!  ---------

      DOUBLE PRECISION SMALL2, ZERO
      PARAMETER (SMALL2 = 0.1D0, &
                 ZERO   = 0.0D0)
!
!  ---------------------------------------------------------------------
!
!* End Prologue
!
!* Local variables
!  ---------------
!
      INTEGER I, K, IS
      DOUBLE PRECISION FHI, HG, U, SUMD, W
      LOGICAL QFINE
!
!* Begin
!
      DO 1 K = 1,N
         IS = 0
!        DO (Until)
 11         CONTINUE
            W = X(K)
            U = DSIGN (ETA(K)*YSCAL(K), X(K))
            X(K) = W + U
            CALL FCN (N, X, FU, IFAIL)
            NFCN = NFCN + 1
!           Exit, If ...
            IF (IFAIL .NE. 0) GOTO 99
            X(K) = W
            SUMD = ZERO
            DO 111 I = 1,N
               HG = DMAX1 (DABS (FX(I)), DABS (FU(I)))
               FHI = FU(I) - FX(I)
               IF (HG .NE. ZERO) SUMD = SUMD + (FHI/HG)**2
               A(I,K) = FHI / U
 111        CONTINUE
            SUMD = DSQRT (SUMD / DBLE(N))
            QFINE = .TRUE.
            IF (SUMD .NE. ZERO .AND. IS .EQ. 0)THEN
               ETA(K) = DMIN1 (ETAMAX, &
                    DMAX1 (ETAMIN, DSQRT (ETADIF / SUMD)*ETA(K)))
               IS = 1
               QFINE = CONV .LT. SMALL2 .OR. SUMD .GE. ETAMIN
            ENDIF
            IF (.NOT.(QFINE)) GOTO  11
!        UNTIL ( expression - negated above)
 1    CONTINUE
!
!     Exit from DO-loop
 99   CONTINUE
!
      RETURN
!
!* End of subroutine N1JCF
!
      END subroutine

      SUBROUTINE N1JCFB (FCN, N, LDA, ML, X, FX, A, YSCAL, ETA, &
           ETAMIN, ETAMAX, ETADIF, CONV, NFCN, FU, U, W, IFAIL)
!* Begin Prologue N1JCFB
!      EXTERNAL FCN
      procedure(func) :: FCN
      INTEGER N, LDA, ML
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), ETA(N), &
           ETAMIN, ETAMAX, ETADIF, CONV
      INTEGER NFCN
      DOUBLE PRECISION FU(N), U(N), W(N)
      INTEGER IFAIL
!
!     ---------------------------------------------------------------------
!
!* Title
!
!  Approximation of banded Jacobian matrix for nonlinear systems solver
!  NLEQ1 with feed-back control of discretization and rounding errors
!
!* Environment       Fortran 77
!             Double Precision
!             Sun 3/60, Sun OS
!* Latest Revision   May 1990
!
!
!* Parameter list description
!  --------------------------
!
!* External subroutines (to be supplied by the user)
!  -------------------------------------------------
!
!  FCN        Ext     FCN (N, X, FX, IFAIL)
!              Subroutine in order to provide right-hand
!              side of first-order differential equations
!    N        Int     Number of rows and columns of the Jacobian
!    X(N)     Dble    The current scaled iterates
!    FX(N)    Dble    Array containing FCN(X)
!    IFAIL    Int     Return code
!              Whenever a negative value is returned by FCN
!              routine N1JCFB is terminated immediately.
!
!
!* Input parameters (* marks inout parameters)
!  ----------------
!
!  N          Int     Number of rows and columns of the Jacobian
!  LDA        Int     Leading dimension of A (LDA .GE. ML+MU+1)
!  ML         Int     Lower bandwidth of Jacobian matrix
!  X(N)       Dble    Array containing the current scaled
!              iterate
!  FX(N)      Dble    Array containing FCN(X)
!  YSCAL(N)   Dble    Array containing the scaling factors
!  ETA(N)     Dble *  Array containing the scaled denominator
!              differences
!  ETAMIN     Dble    Minimum allowed scaled denominator
!  ETAMAX     Dble    Maximum allowed scaled denominator
!  ETADIF     Dble    DSQRT (1.1*EPMACH)
!              EPMACH = machine precision
!  CONV       Dble    Maximum norm of last (unrelaxed) Newton correction
!  NFCN       Int  *  FCN - evaluation count
!
!* Output parameters (* marks inout parameters)
!  -----------------
!
!  A(LDA,N)   Dble    Array to contain the approximated
!              Jacobian matrix ( dF(i)/dx(j)in A(i,j))
!  ETA(N)     Dble *  Scaled denominator differences adjusted
!  NFCN       Int  *  FCN - evaluation count adjusted
!  IFAIL      Int     Return code non-zero if Jacobian could not
!              be computed.
!
!* Workspace parameters
!  --------------------
!
!  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
!              the numerator differences
!  U(N)       Dble    Array to contain dx(i)
!  W(N)       Dble    Array to save original components of X
!
!* Called
!  ------
!
      INTRINSIC DABS, DMAX1, DMIN1, DSIGN, DSQRT, MAX0, MIN0
!
!* Constants
!  ---------
!
      DOUBLE PRECISION SMALL2, ZERO
      PARAMETER (SMALL2 = 0.1D0, &
                 ZERO   = 0.0D0)
!
!  ---------------------------------------------------------------------
!
!* End Prologue
!
!* Local variables
!  ---------------
!
      INTEGER I, IS, I1, I2, JJ, K, MH, MU, LDAB
      DOUBLE PRECISION FHI, HG, SUMD
      LOGICAL QFINE
!
!* Begin
!
      MU = LDA - 2*ML - 1
      LDAB = ML + MU + 1
      DO 1 I = 1,LDAB
         DO 11 K = 1,N
            A(I,K) = ZERO
 11      CONTINUE
 1    CONTINUE
!
      DO 2 JJ = 1,LDAB
         IS = 0
!        DO (Until)
 21         CONTINUE
            DO 211  K = JJ,N,LDAB
               W(K) = X(K)
               U(K) = DSIGN (ETA(K)*YSCAL(K), X(K))
               X(K) = W(K) + U(K)
 211        CONTINUE
!
            CALL FCN (N, X, FU, IFAIL)
            NFCN = NFCN + 1
!           Exit, If ...
            IF (IFAIL .NE. 0) GOTO 99
!
            DO 212 K = JJ,N,LDAB
               X(K) = W(K)
               SUMD = ZERO
               I1 = MAX0 (1, K - MU)
               I2 = MIN0 (N, K + ML)
               MH = MU + 1 - K
               DO 213 I = I1,I2
                  HG = DMAX1 (DABS(FX(I)), DABS(FU(I)))
                  FHI = FU(I) - FX(I)
                  IF (HG .NE. ZERO) SUMD = SUMD + (FHI/HG)**2
                  A(MH+I, K) = FHI / U(K)
 213           CONTINUE
               SUMD = DSQRT (SUMD / DBLE(N))
               QFINE = .TRUE.
               IF (SUMD .NE. ZERO .AND. IS .EQ. 0) THEN
                  ETA(K) = DMIN1 (ETAMAX, &
                       DMAX1 (ETAMIN, DSQRT (ETADIF / SUMD)*ETA(K)))
                  IS = 1
                  QFINE = CONV .LT. SMALL2 .OR. SUMD .GE. ETAMIN
               ENDIF
 212        CONTINUE
            IF (.NOT.(QFINE)) GOTO 21
!        UNTIL ( expression - negated above)
 2    CONTINUE

!     Exit from DO-loop
 99   CONTINUE

      RETURN
!
!* End of subroutine N1JCFB
!
      END subroutine

      SUBROUTINE N1PRV1(DLEVF,DLEVX,FC,NITER,NEW,MPRMON,LUMON,QMIXIO)
!     Begin Prologue N1PRV1
      DOUBLE PRECISION DLEVF,DLEVX,FC
      INTEGER NITER,NEW,MPRMON,LUMON
      LOGICAL QMIXIO
!     ------------------------------------------------------------
!
!     Summary :
!
!     N 1 P R V 1 : Printing of intermediate values (Type 1 routine)
!
!     Parameters
!     ==========
!
!     DLEVF, DLEVX   See descr. of internal double variables of N1INT
!     FC,NITER,NEW,MPRMON,LUMON
!           See parameter descr. of subroutine N1INT
!     QMIXIO Logical  = .TRUE.  , if LUMON.EQ.LUSOL
!              = .FALSE. , if LUMON.NE.LUSOL
!
!     ------------------------------------------------------------
!     End Prologue
!     Print Standard - and natural level
      IF(QMIXIO)THEN
1       FORMAT(2X,66('*'))
        WRITE(LUMON,1)
2       FORMAT(8X,'It',7X,'Normf ',10X,'Normx ',20X,'New')
        IF (MPRMON.GE.3) WRITE(LUMON,2)
3       FORMAT(8X,'It',7X,'Normf ',10X,'Normx ',8X,'Damp.Fct.',3X,'New')
        IF (MPRMON.EQ.2) WRITE(LUMON,3)
      ENDIF
4     FORMAT(6X,I4,5X,D10.3,2X,4X,D10.3,17X,I2)
      IF (MPRMON.GE.3.OR.NITER.EQ.0)  &
        WRITE(LUMON,4) NITER,DLEVF,DLEVX,NEW
5     FORMAT(6X,I4,5X,D10.3,6X,D10.3,6X,F7.5,4X,I2)
      IF (MPRMON.EQ.2.AND.NITER.NE.0)  &
        WRITE(LUMON,5) NITER,DLEVF,DLEVX,FC,NEW
      IF(QMIXIO)THEN
6       FORMAT(2X,66('*'))
        WRITE(LUMON,6)
      ENDIF
!     End of subroutine N1PRV1
      RETURN
      END subroutine

      SUBROUTINE N1PRV2(DLEVF,DLEVX,FC,NITER,MPRMON,LUMON,QMIXIO, &
                        CMARK)
!     Begin Prologue N1PRV2
      DOUBLE PRECISION DLEVF,DLEVX,FC
      INTEGER NITER,MPRMON,LUMON
      LOGICAL QMIXIO
      CHARACTER*1 CMARK
!     ------------------------------------------------------------
!
!     Summary :
!
!     N 1 P R V 2 : Printing of intermediate values (Type 2 routine)
!
!     Parameters
!     ==========
!
!     DLEVF, DLEVX   See descr. of internal double variables of N2INT
!     FC,NITER,MPRMON,LUMON
!           See parameter descr. of subroutine N2INT
!     QMIXIO Logical  = .TRUE.  , if LUMON.EQ.LUSOL
!              = .FALSE. , if LUMON.NE.LUSOL
!     CMARK Char*1    Marker character to be printed before DLEVX
!
!     ------------------------------------------------------------
!     End Prologue
!     Print Standard - and natural level, and damping
!     factor
      IF(QMIXIO)THEN
1       FORMAT(2X,66('*'))
        WRITE(LUMON,1)
2       FORMAT(8X,'It',7X,'Normf ',10X,'Normx ',8X,'Damp.Fct.')
        WRITE(LUMON,2)
      ENDIF
3     FORMAT(6X,I4,5X,D10.3,4X,A1,1X,D10.3,2X,4X,F7.5)
      WRITE(LUMON,3)NITER,DLEVF,CMARK,DLEVX,FC
      IF(QMIXIO)THEN
4       FORMAT(2X,66('*'))
        WRITE(LUMON,4)
      ENDIF
!     End of subroutine N1PRV2
      RETURN
      END subroutine

      SUBROUTINE N1SOUT(N,X,MODE,IOPT,RWK,NRW,IWK,NIW,MPRINT,LUOUT)
!     Begin Prologue SOLOUT
      INTEGER N
      DOUBLE PRECISION X(N)
      INTEGER NRW
      INTEGER MODE
      INTEGER IOPT(50)
      DOUBLE PRECISION RWK(NRW)
      INTEGER NIW
      INTEGER IWK(NIW)
      INTEGER MPRINT,LUOUT
!     ------------------------------------------------------------
!
!     Summary :
!
!     S O L O U T : Printing of iterate (user customizable routine)
!
!     Input parameters
!     ================
!
!     N         Int Number of equations/unknowns
!     X(N)   Double iterate vector
!     MODE          =1 This routine is called before the first
!               Newton iteration step
!            =2 This routine is called with an intermedi-
!               ate iterate X(N)
!            =3 This is the last call with the solution
!               vector X(N)
!            =4 This is the last call with the final, but
!               not solution vector X(N)
!     IOPT(50)  Int The option array as passed to the driver
!            routine(elements 46 to 50 may be used
!            for user options)
!     MPRINT    Int Solution print level
!            (see description of IOPT-field MPRINT)
!     LUOUT     Int the solution print unit
!            (see description of see IOPT-field LUSOL)
!
!
!     Workspace parameters
!     ====================
!
!     NRW, RWK, NIW, IWK    see description in driver routine
!
!     Use of IOPT by this routine
!     ===========================
!
!     Field 46:       =0 Standard output
!              =1 GRAZIL suitable output
!
!     ------------------------------------------------------------
!     End Prologue
      LOGICAL QGRAZ,QNORM
      integer :: L1, I
!     Begin
      QNORM = IOPT(46).EQ.0
      QGRAZ = IOPT(46).EQ.1
      IF(QNORM) THEN
1        FORMAT('  ',A,' data:',/)
         IF (MODE.EQ.1) THEN
101        FORMAT('  Start data:',/,'  N =',I5,//, &
                  '  Format: iteration-number, (x(i),i=1,...N), ', &
                  'Normf , Normx ',/)
           WRITE(LUOUT,101) N
           WRITE(LUOUT,1) 'Initial'
         ELSE IF (MODE.EQ.3) THEN
           WRITE(LUOUT,1) 'Solution'
         ELSE IF (MODE.EQ.4) THEN
           WRITE(LUOUT,1) 'Final'
         ENDIF
2        FORMAT(' ',I5)
!        WRITE          NITER
         WRITE(LUOUT,2) IWK(1)
3        FORMAT((12X,3(D18.10,1X)))
         WRITE(LUOUT,3)(X(L1),L1=1,N)
!        WRITE          DLEVF,   DLEVX
         WRITE(LUOUT,3) RWK(19),DSQRT(RWK(18)/DBLE(FLOAT(N)))
         IF(MODE.EQ.1.AND.MPRINT.GE.2) THEN
           WRITE(LUOUT,1) 'Intermediate'
         ELSE IF(MODE.GE.3) THEN
           WRITE(LUOUT,1) 'End'
         ENDIF
      ENDIF
      IF(QGRAZ) THEN
        IF(MODE.EQ.1) THEN
10        FORMAT('&name com',I3.3,:,255(7(', com',I3.3,:),/))
          WRITE(LUOUT,10)(I,I=1,N+2)
15        FORMAT('&def  com',I3.3,:,255(7(', com',I3.3,:),/))
          WRITE(LUOUT,15)(I,I=1,N+2)
16        FORMAT(6X,': X=1, Y=',I3)
          WRITE(LUOUT,16) N+2
        ENDIF
20      FORMAT('&data ',I5)
!        WRITE          NITER
        WRITE(LUOUT,20) IWK(1)
21      FORMAT((6X,4(D18.10)))
        WRITE(LUOUT,21)(X(L1),L1=1,N)
!        WRITE          DLEVF,   DLEVX
        WRITE(LUOUT,21) RWK(19),DSQRT(RWK(18)/DBLE(FLOAT(N)))
        IF(MODE.GE.3) THEN
30        FORMAT('&wktype 3111',/,'&atext x ''iter''')
          WRITE(LUOUT,30)
35        FORMAT('&vars = com',I3.3,/,'&atext y ''x',I3,'''', &
                 /,'&run')
          WRITE(LUOUT,35) (I,I,I=1,N)
36        FORMAT('&vars = com',I3.3,/,'&atext y ''',A,'''', &
                 /,'&run')
          WRITE(LUOUT,36) N+1,'Normf ',N+2,'Normx '
!39       FORMAT('&stop')
!         WRITE(LUOUT,39)
        ENDIF
      ENDIF
!     End of subroutine N1SOUT
      RETURN
      END subroutine
!     End package

!--------------------------------------------------------------------
! utils
! All NLEQ1 files packaged together
!--------------------------------------------------------------------

      SUBROUTINE ZIBCONST(EPMACH,SMALL)
      DOUBLE PRECISION EPMACH,SMALL
!
!*********************************************************************
!  Set Approximations to machine constants.
!  This routine is machine dependent.
!*********************************************************************
!
!  Output parameters
!    EPMACH    DOUBLE     Relative machine precision
!    SMALL     DOUBLE     SQRT of smallest positive machine number
!
!  The proposed values should work on Intel compatible cpus,
!  PowerPCs and Sun sparcs
!
!*********************************************************************

      EPMACH = 1.0D-17
      SMALL  = 1.0D-150
      RETURN
      END subroutine

      SUBROUTINE MON
!
!* Begin Prologue MON

!  ---------------------------------------------------------------------
!
!* Title
!
!  Time monitor to discover cpu time consumption of interesting parts
!  of a code.
!
!* Written by        U. Nowak, U. Poehle, L. Weimann
!* Purpose           Cpu time measuring
!* Category          ...
!* File              mon.f
!* Version           2.1
!* Latest Change     96/03/04 (1.3)
!* Library           CodeLib
!* Code              Fortran 77
!             Single precision
!
!  ---------------------------------------------------------------------
!
!* Summary
!  -------
!
!  MON provides different entries
!    - to initialize the time monitor,
!    - to start and to stop the measuring,
!    - to start and to stop a number of stop-watches which may be nested,
!    - to print a summary table,
!    - to define text labels for the summary table, and
!    - to store the vector of average cpu times
!
!  Errors detected during the measuring will disable the time monitor
!  but will not affect the calling program.
!
!  ---------------------------------------------------------------------
!
!* Parameter list description
!  --------------------------
!
!* Input parameters
!  ----------------
!
!  TEXTH      Char    Text to identify the measuring (up to 31
!              characters will be printed in the headline of the
!              summary table)
!              Used by MONINI
!
!  IOUNIT     Int     Output unit for error messages and summary table
!              Used by MONINI
!
!  INDX       Int     Number of a range of source code lines (a distinct
!              part of the algorithm normally) for which the
!              cpu time will be measured by a kind of stop-watch
!              Used by MONDEF, MONON, and MONOFF
!
!  NAMEH      Char    Text to label the summary lines (up to 17
!              characters will be printed)
!              Used by MONDEF
!
!* Output parameters
!  -----------------
!
!  IRET       Int     Return code of ZIBSEC
!              Used by MONSTR
!
!  AVER       Dble    Vector of average cpu time for all measurings
!              Used by MONGET
!
!
!* End Prologue
!  ------------
!
!
!* Constants
!  ---------
!
!  MAXTAB     Int     Maximum number of stop-watches
!  MNEST      Int     Maximum number of nested measurings
!
!
!* Local Variables
!  ---------------
!
!  ASE!       Real    Array of averages per call
!  INDACT     Int     Array for indices of nested measurings
!  NCALL      Int     Array of call counts
!  PC1        Real    Array of per cent values with respect to all
!              measurings
!  PC2        Real    Array of per cent values with respect to the sum
!              of all stop-watches
!  QDISAB     Log     Time monitor disabled
!  QON        Log     Array reflecting if an index is active
!  QSTART     Log     Time monitor started
!  SEC        Real    Array of cpu time measurings

      INTEGER MAXTAB, MNEST

      PARAMETER (MAXTAB = 25, MNEST=20)

      INTEGER INDACT(MNEST), NCALL(0:MAXTAB)
      REAL ASEC(0:MAXTAB), PC1(0:MAXTAB), PC2(0:MAXTAB), SEC(0:MAXTAB)
      DOUBLE PRECISION AVER(0:MAXTAB)
      LOGICAL QDISAB, QON(0:MAXTAB), QSTART
      CHARACTER NAME(0:MAXTAB)*17, TEXT*30
      CHARACTER*(*) NAMEH, TEXTH

      integer :: IONCNT, MAXIND, MONI, INFO, IOUNIT, INDX, IRET
      integer :: I, IFAIL
      real :: CPTIM, SUM, SUM0

      SAVE ASEC, CPTIM, INDACT, IONCNT, MAXIND, MONI, NAME, NCALL, PC1, &
           PC2, QDISAB, QON, QSTART, SEC, TEXT

      DATA CPTIM /0.0E0/, MONI /6/, INFO /1/, IONCNT /-1/, &
           QDISAB /.FALSE./, QSTART /.FALSE./

      RETURN


      ENTRY MONINI (TEXTH, IOUNIT)
!     Initialize monitor.
!     Has to be called first. May be called again after MONHLT.

      TEXT = TEXTH
      MONI = IOUNIT

      IF (IONCNT .GT. 0 .AND. .NOT. QDISAB) GOTO 1070

      MAXIND = 0
      IONCNT = 0
      QDISAB = .FALSE.

      DO 1000 I = 0,MAXTAB
         SEC(I) = 0.
         ASEC(I) = 0.
         NCALL(I) = 0
         QON(I) = .FALSE.
         WRITE (NAME(I), '(A, I2)') 'Part ', I
 1000 CONTINUE
      NAME(0) = 'General'

      DO 1010 I = 1,MNEST
         INDACT(I) = 0
 1010 CONTINUE

      RETURN

      ENTRY MONDEF(INDX, NAMEH)
!     Define one monitor entry.
!     May be called at any time before MONPRT.

      IF (QDISAB) RETURN
      IF (INDX.LT.0 .OR. INDX.GT.MAXTAB) GOTO 1080

      NAME(INDX) = NAMEH

      RETURN


      ENTRY MONSTR (IRET)
!     Start monitor measurings.
!     Has to be called once after initialization.

      IF (QDISAB) RETURN
      IF (IONCNT .LT. 0) GOTO 1090
      IF (IONCNT .GT. 0) GOTO 1100
      IF (QON(0)) GOTO 1110

      IFAIL = 0
      CALL ZIBSEC (CPTIM, IFAIL)

      IF (IFAIL .EQ. 0) THEN
!        Switch on general stop-watch
         SEC(0) = -CPTIM
         QON(0) = .TRUE.
         IONCNT = 1
         QSTART = .TRUE.
      ENDIF
      IRET = IFAIL

      RETURN

      ENTRY MONON (INDX)
!     Start one measuring.
!     A running stop-watch will be deactivated until the new stop-watch
!     stops.

      IF (.NOT. QSTART) RETURN
      IF (QDISAB) RETURN
      IF (IONCNT .LT. 1) GOTO 1120
      IF (INDX .GT. MAXTAB .OR. INDX .LE. 0) GOTO 1130
      IF (QON(INDX)) GOTO 1140

      MAXIND = MAX(MAXIND, INDX)
      CALL ZIBSEC (CPTIM, IFAIL)

!     Hold actual stop-watch
      SEC(INDACT(IONCNT)) = SEC(INDACT(IONCNT)) + CPTIM

!     Switch on new stop-watch
      IF (INFO .GT. 1) WRITE (MONI,*) ' Enter ', NAME(INDX), SEC(INDX)

      IONCNT = IONCNT + 1
      IF (IONCNT .GT. MNEST) GOTO 1150

      INDACT(IONCNT) = INDX
      SEC(INDX) = SEC(INDX) - CPTIM
      QON(INDX) = .TRUE.

      RETURN


      ENTRY MONOFF (INDX)
!     Stop one measuring.
!     May be called for the stop-watch started most recently.
!     The stop-watch deactivated most recently will be activated.

      IF (.NOT. QSTART) RETURN
      IF (QDISAB) RETURN
      IF (INDX .GT. MAXTAB .OR. INDX .LE. 0) GOTO 1160
      IF (INDACT(IONCNT) .NE. INDX) GOTO 1170

      CALL ZIBSEC (CPTIM, IFAIL)

!     Switch off actual stop-watch
      QON(INDX) = .FALSE.
      SEC(INDX) = SEC(INDX) + CPTIM
      NCALL(INDX) = NCALL(INDX) + 1
      IONCNT = IONCNT - 1
      IF (INFO .GT. 1) WRITE (MONI,*) ' Exit ', NAME(INDX), SEC(INDX)

!     Continue previous stop-watch
      SEC(INDACT(IONCNT)) = SEC(INDACT(IONCNT)) - CPTIM

      RETURN


      ENTRY MONHLT
!     Terminate monitor.
!     Stops all active stop-watches.

      IF (.NOT. QSTART) RETURN
      IF (QDISAB) RETURN

      CALL ZIBSEC (CPTIM, IFAIL)

      DO 1020 I=IONCNT,1,-1
         QON(INDACT(I)) = .FALSE.
         SEC(INDACT(I)) = SEC(INDACT(I)) + CPTIM
         NCALL(INDACT(I)) = NCALL(INDACT(I)) + 1
 1020 CONTINUE

      IONCNT = 0
!     This means that the time monitor has to be started by MONSTR
!     before another measuring.

      RETURN



      ENTRY MONPRT

!     Print statistics.
!     May be called after MONHLT only.

      IF (IONCNT .GT. 0) GOTO 1180
      IF (.NOT. QSTART) GOTO 1190

      SUM = 1.E-10
      DO 1030 I = 1,MAXIND
         SUM = SUM + SEC(I)
         IF (NCALL(I) .GT. 0) ASEC(I) = SEC(I)/FLOAT(NCALL(I))
 1030 CONTINUE
      SUM0 = SUM + SEC(0)
      IF (NCALL(0) .GT. 0) ASEC(0) = SEC(0)/FLOAT(NCALL(0))

      DO 1040 I = 1,MAXIND
         PC1(I) = 100.*SEC(I)/SUM0
         PC2(I) = 100.*SEC(I)/SUM
 1040 CONTINUE
      PC1(0) = 100.*SEC(0)/SUM0

      WRITE (MONI,9000)
      WRITE (MONI,9010)
      WRITE (MONI,9020) ' '
 9000 FORMAT (///)
 9010 FORMAT (1X, 77('#'))
 9020 FORMAT (' #   ', A, T75,  '   #')

      IF (QDISAB) THEN
         WRITE (MONI,9020) ' '
         WRITE (MONI,9020) &
              'Warning  The following results may be misleading'
         WRITE (MONI,9020) &
              'because an error occured and disabled the time monitor'
      ENDIF

      WRITE (MONI,9020) ' '
      WRITE (MONI,9030) 'Results from time monitor program for:', TEXT
 9030 FORMAT (' #   ', A40, 1X, A31, '#')

      WRITE (MONI,9020) ' '
      WRITE (MONI,9040) 'Total time:', SUM0, 'Sum of parts:', SUM
 9040 FORMAT (' #   ', A11, F11.3, 5X, A13, F11.3, 21X, '#')

      WRITE (MONI,9020) ' '
      WRITE (MONI,9050)
 9050 FORMAT (' #   ', 2X, 'Name', 14X, 'Calls', 7X, 'Time', 4X, &
           'Av-time', 4X, '% Total', 6X, '% Sum   #')

      WRITE (MONI,9060) NAME(0), NCALL(0), SEC(0), ASEC(0), PC1(0)
 9060 FORMAT (' #   ', A17, I8, F11.3, F11.4, F11.2, 14X, '#')

      DO 1050 I = 1,MAXIND
         WRITE (MONI,9070) NAME(I), NCALL(I), SEC(I), ASEC(I), PC1(I), &
              PC2(I)
 1050 CONTINUE

 9070 FORMAT (' #   ', A17, I8, F11.3, F11.4, F11.2, F11.2, 3X, '#')


      WRITE (MONI,9020) ' '
      WRITE (MONI,9010)
      WRITE (MONI,9000)

      RETURN


      ENTRY MONGET(AVER)
!     Store average cpu times per call.
!     May be called at any time after MONPRT.

      IF (.NOT. QSTART) RETURN
      IF (QDISAB) RETURN

      DO 1060 I=0,MAXIND
         AVER(I) = DBLE(ASEC(I))
 1060 CONTINUE

      RETURN
!
!  Error exits
!
 1070 CONTINUE
      WRITE (MONI,9080) 'MONINI', 'Time monitor is running already.'
      GOTO 1200

 1080 CONTINUE
      WRITE (MONI,9090) 'MONDEF', 'Index out of range', INDX
      GOTO 1200

 1090 CONTINUE
      WRITE (MONI,9080) 'MONSTR', &
           'Time monitor has to be initialized by MONINI first.'
      GOTO 1200

 1100 CONTINUE
      WRITE (MONI,9080) 'MONSTR', 'Time monitor is running already.'
      GOTO 1200

 1110 CONTINUE
      WRITE (MONI,9080) 'MONSTR', &
           'Time monitor has been started already.'
      GOTO 1200

 1120 CONTINUE
      WRITE (MONI,9080) 'MONON', 'Time monitor is not yet started.'
      GOTO 1200

 1130 CONTINUE
      WRITE (MONI,9090) 'MONON', 'Index out of range', INDX
      GOTO 1200

 1140 CONTINUE
      WRITE (MONI,9090) 'MONON', &
           'Measuring is running already for this INDX', INDX
      GOTO 1200

 1150 CONTINUE
      WRITE (MONI,9100) 'MONON', 'Nesting is too deep.', &
           'The following indices are active', (INDACT(I),I=0,IONCNT)
      GOTO 1200

 1160 CONTINUE
      WRITE (MONI,9090) 'MONOFF', 'Index out of range', INDX
      GOTO 1200

 1170 CONTINUE
      WRITE (MONI,9110) 'MONOFF', 'Measuring ', INDX, &
           'cannot be stopped.', &
           'The following indices are active', (INDACT(I),I=0,IONCNT)
      GOTO 1200

 1180 CONTINUE
      WRITE (MONI,9080) 'MONPRT', 'Time monitor is still running.'
      GOTO 1200

 1190 CONTINUE
      WRITE (MONI,9080) 'MONPRT', 'Time monitor was not started.'
      GOTO 1200

 1200 QDISAB = .TRUE.
      RETURN

 9080 FORMAT (/, ' ++ Error in subroutine ', A, ' ++',/, 4X, A)

 9090 FORMAT (/, ' ++ Error in subroutine ', A, ' ++',/, 4X, A, &
           ' (', I6, ').')

 9100 FORMAT (/, ' ++ Error in subroutine ', A, ' ++', 4X, A,/, 4X, A, &
           (10I4))

 9110 FORMAT (/, ' ++ Error in subroutine ', A, ' ++', 4X, A, I3, 1X, A, &
           /, 4X, A, (10I4))
!
!  End subroutine monitor
!
      END subroutine

      SUBROUTINE ZIBSEC (CPTIM, IFAIL)

!*********************************************************************
!  Set CPTIM to cpu time in seconds.
!  This routine is machine dependent.
!*********************************************************************
!
!  Output parameters
!    CPTIM    REAL        Cpu time in seconds
!    IFAIL    INTEGER     Errorcode
!
!  CSUN code works on sun and MacOSX (g77).
!  CXLF code works on MacOSX (xlf)
!  CVIS code works with Visual Fortran (Windows)
!
!*********************************************************************

      REAL CPTIM
      INTEGER IFAIL
	real tim

!SUN      REAL RTIME(2)
!XLF      REAL(4) etime_
!XLF      TYPE TB_TYPE
!XLF      SEQUENCE
!XLF      REAL(4) USRTIME
!XLF      REAL(4) SYSTIME
!XLF      END TYPE
!XLF      TYPE (TB_TYPE) ETIME_STRUCT

      CPTIM = 0
      IFAIL = 0

!SUN      CPTIM = ETIME(RTIME)
!SUN      CPTIM = RTIME(1)

!XLF      CPTIM = etime_(ETIME_STRUCT)

      CALL CPU_TIME(TIM)	! CVIS
	CPTIM = TIM

      RETURN
      END subroutine

      DOUBLE PRECISION FUNCTION WNORM(N,Z,XW)
      INTEGER N
      DOUBLE PRECISION Z(N), XW(N)
!     ------------------------------------------------------------
!
!     Summary :
!
!     E N O R M : Return the norm to be used in exit (termination)
!          criteria
!
!     Input parameters
!     ================
!
!     N         Int Number of equations/unknowns
!     Z(N)     Dble  The vector, of which the norm is to be computed
!     XW(N)    Dble  The scaling values of Z(N)
!
!     Output
!     ======
!
!     WNORM(N,Z,XW)  Dble  The mean square root norm of Z(N) subject
!                   to the scaling values in XW(N):
!                   = Sqrt( Sum(1,...N)((Z(I)/XW(I))**2) / N )
!
!     ------------------------------------------------------------
!     End Prologue
      INTEGER I
      DOUBLE PRECISION S
!     Begin
      S = 0.0D0
      DO 10 I=1,N
        S = S + ( Z(I)/XW(I) ) ** 2
10    CONTINUE
      WNORM = DSQRT( S / DBLE(FLOAT(N)) )
!     End of function WNORM
      RETURN
      END function

! Linear algebra

      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(N),INFO
      DOUBLE PRECISION A(LDA,N)
!
!     dgefa factors a double precision matrix by gaussian elimination.
!
!     dgefa is usually called by dgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
!
!     on entry
!
!        a       double precision(lda, n)
!         the matrix to be factored.
!
!        lda     integer
!         the leading dimension of the array  a .
!
!        n       integer
!         the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!         which were used to obtain it.
!         the factorization can be written  a = l*u  where
!         l  is a product of permutation and unit lower
!         triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!         an integer vector of pivot indices.
!
!        info    integer
!         = 0  normal value.
!         = k  if  u(k,k) .eq. 0.0 .  this is not an error
!              condition for this subroutine, but it does
!              indicate that dgesl or dgedi will divide by zero
!              if called.  use  rcond  in dgeco for a reliable
!              indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,dscal,idamax
!
!     internal variables
!
      DOUBLE PRECISION T
!      INTEGER IDAMAX,J,K,KP1,L,NM1
      INTEGER J,K,KP1,L,NM1

!
!     gaussian elimination with partial pivoting
!
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
!
!        find l = pivot index
!
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
!
!        zero pivot implies this column already triangularized
!
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
!
!           interchange if necessary
!
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
!
!           compute multipliers
!
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
!
!           row elimination with column indexing
!
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END subroutine

      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
      INTEGER LDA,N,IPVT(N),JOB
      DOUBLE PRECISION A(LDA,N),B(N)
!
!     dgesl solves the double precision system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by dgeco or dgefa.
!
!     on entry
!
!        a       double precision(lda, n)
!         the output from dgeco or dgefa.
!
!        lda     integer
!         the leading dimension of the array  a .
!
!        n       integer
!         the order of the matrix  a .
!
!        ipvt    integer(n)
!         the pivot vector from dgeco or dgefa.
!
!        b       double precision(n)
!         the right hand side vector.
!
!        job     integer
!         = 0         to solve  a*x = b ,
!         = nonzero   to solve  trans(a)*x = b  where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if dgeco has set rcond .gt. 0.0
!        or dgefa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call dgesl(a,lda,n,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,ddot
!
!     internal variables
!
!      DOUBLE PRECISION DDOT,T
      DOUBLE PRECISION T
      INTEGER K,KB,L,NM1

      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
!
!        job = 0 , solve  a * x = b
!        first solve  l*y = b
!
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
!
!        now solve  u*x = y
!
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE

!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b

         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
!
!        now solve trans(l)*x = y
!
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END subroutine

      SUBROUTINE DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
      INTEGER LDA,N,ML,MU,IPVT(N),INFO
      DOUBLE PRECISION ABD(LDA,N)
!
!     dgbfa factors a double precision band matrix by elimination.
!
!     dgbfa is usually called by dgbco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!
!     on entry
!
!        abd     double precision(lda, n)
!         contains the matrix in band storage.  the columns
!         of the matrix are stored in the columns of  abd  and
!         the diagonals of the matrix are stored in rows
!         ml+1 through 2*ml+mu+1 of  abd .
!         see the comments below for details.
!
!        lda     integer
!         the leading dimension of the array  abd .
!         lda must be .ge. 2*ml + mu + 1 .
!
!        n       integer
!         the order of the original matrix.
!
!        ml      integer
!         number of diagonals below the main diagonal.
!         0 .le. ml .lt. n .
!
!        mu      integer
!         number of diagonals above the main diagonal.
!         0 .le. mu .lt. n .
!         more efficient if  ml .le. mu .
!     on return
!
!        abd     an upper triangular matrix in band storage and
!         the multipliers which were used to obtain it.
!         the factorization can be written  a = l*u  where
!         l  is a product of permutation and unit lower
!         triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!         an integer vector of pivot indices.
!
!        info    integer
!         = 0  normal value.
!         = k  if  u(k,k) .eq. 0.0 .  this is not an error
!              condition for this subroutine, but it does
!              indicate that dgbsl will divide by zero if
!              called.  use  rcond  in dgbco for a reliable
!              indication of singularity.
!
!     band storage
!
!           if  a  is a band matrix, the following program segment
!           will set up the input.
!
!            ml = (band width below the diagonal)
!            mu = (band width above the diagonal)
!            m = ml + mu + 1
!            do 20 j = 1, n
!               i1 = max0(1, j-mu)
!               i2 = min0(n, j+ml)
!               do 10 i = i1, i2
!                  k = i - j + m
!                  abd(k,j) = a(i,j)
!         10    continue
!         20 continue
!
!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
!           in addition, the first  ml  rows in  abd  are used for
!           elements generated during the triangularization.
!           the total number of rows needed in  abd  is  2*ml+mu+1 .
!           the  ml+mu by ml+mu  upper left triangle and the
!           ml by ml  lower right triangle are not referenced.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,dscal,idamax
!     fortran max0,min0
!
!     internal variables
!
      DOUBLE PRECISION T
!      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
      INTEGER I,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1

      M = ML + MU + 1
      INFO = 0
!
!     zero initial fill-in columns
!
      J0 = MU + 2
      J1 = MIN0(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0D0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
!
!     gaussian elimination with partial pivoting
!
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1
!
!        zero next fill-in column
!
         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0D0
  40       CONTINUE
   50    CONTINUE
!
!        find l = pivot index
!
         LM = MIN0(ML,N-K)
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M
!
!        zero pivot implies this column already triangularized
!
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100
!
!           interchange if necessary
!
            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE
!
!           compute multipliers
!
            T = -1.0D0/ABD(M,K)
            CALL DSCAL(LM,T,ABD(M+1,K),1)
!
!           row elimination with column indexing
!
            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N
      RETURN
      END subroutine

      SUBROUTINE DGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)
      INTEGER LDA,N,ML,MU,IPVT(N),JOB
      DOUBLE PRECISION ABD(LDA,N),B(N)
!
!     dgbsl solves the double precision band system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by dgbco or dgbfa.
!
!     on entry
!
!        abd     double precision(lda, n)
!         the output from dgbco or dgbfa.
!
!        lda     integer
!         the leading dimension of the array  abd .
!
!        n       integer
!         the order of the original matrix.
!
!        ml      integer
!         number of diagonals below the main diagonal.
!
!        mu      integer
!         number of diagonals above the main diagonal.
!
!        ipvt    integer(n)
!         the pivot vector from dgbco or dgbfa.
!
!        b       double precision(n)
!         the right hand side vector.
!
!        job     integer
!         = 0         to solve  a*x = b ,
!         = nonzero   to solve  trans(a)*x = b , where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if dgbco has set rcond .gt. 0.0
!        or dgbfa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,ddot
!     fortran min0
!
!     internal variables
!
!      DOUBLE PRECISION DDOT,T
      DOUBLE PRECISION T
      INTEGER K,KB,L,LA,LB,LM,M,NM1

      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
!
!        job = 0 , solve  a * x = b
!        first solve l*y = b
!
         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN0(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
!
!        now solve  u*x = y
!
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
!
!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b
!
         DO 60 K = 1, N
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
!
!        now solve trans(l)*x = y
!
         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN0(ML,N-K)
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END subroutine

      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
      DOUBLE PRECISION DX(*),DY(*),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N

      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        code for both increments equal to 1
!
!        clean-up loop
!
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END subroutine

      SUBROUTINE  DSCAL(N,DA,DX,INCX)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!
      DOUBLE PRECISION DA,DX(*)
      INTEGER I,INCX,M,MP1,N,NINCX
!
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
!
!        code for increment not equal to 1
!
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
!
!        code for increment equal to 1
!

!        clean-up loop

   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END subroutine

      INTEGER FUNCTION IDAMAX(N,DX,INCX)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!
      DOUBLE PRECISION DX(*),DMAX
      INTEGER I,INCX,IX,N

      IDAMAX = 0
      IF( N .LT. 1 ) RETURN
      IDAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
!
!        code for increment not equal to 1
!
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5
         IDAMAX = I
         DMAX = DABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
!
!        code for increment equal to 1
!
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
         IF(DABS(DX(I)).LE.DMAX) GO TO 30
         IDAMAX = I
         DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END function

      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
      DOUBLE PRECISION DX(*),DY(*),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N

      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
!
!        code for both increments equal to 1

!
!        clean-up loop
!
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) + &
         DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END function

end module
