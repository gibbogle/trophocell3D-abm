! Modified for parallel use (OMP)

module rkc_90

implicit none

integer, parameter :: nflogg=12
double precision :: rkc_sprad

public :: rkc_comm, rkc_sprad

type rkc_comm
    private
    integer :: nfe,nsteps,naccpt,nrejct,nfesig,maxm
end type


contains

      subroutine rkc(comm,neqn,f,y,t,tend,rtol,atol,info,work,idid,icase)
!--------------------------------------------------------------------------
!
!  ABSTRACT:  RKC integrates initial value problems for systems of first 
!  order ordinary differential equations.  It is based on a family of
!  explicit Runge-Kutta-Chebyshev formulas of order two.  The stability
!  of members of the family increases quadratically in the number of 
!  stages m. An estimate of the spectral radius is used at each step to 
!  select the smallest m resulting in a stable integration. RKC is 
!  appropriate for the solution to modest accuracy of mildly stiff problems
!  with eigenvalues of Jacobians that are close to the negative real axis.  
!  For such problems it has the advantages of explicit one-step methods and
!  very low storage. If it should turn out that RKC is using m far beyond 
!  100, the problem is not mildly stiff and alternative methods should be 
!  considered.  Answers can be obtained cheaply anywhere in the interval 
!  of integration by means of a continuous extension evaluated in the 
!  subroutine RKCINT.
!
!  The initial value problems arising from semi-discretization of 
!  diffusion-dominated parabolic partial differential equations and of 
!  reaction-diffusion equations, especially in two and three spatial 
!  variables, exemplify the problems for which RKC was designed.  Two
!  example programs, ExA and ExB, are provided that show how to use RKC.
!
!---------------------------------------------------------------------------
!  USAGE:  RKC integrates a system of NEQN first order ordinary differential
!  equations specified by a subroutine F from T to TEND.  The initial values
!  at T are input in Y(*).  On all returns from RKC, Y(*) is an approximate 
!  solution at T.  In the computation of Y(*), the local error has been 
!  controlled at each step to satisfy a relative error tolerance RTOL and 
!  absolute error tolerances ATOL(*).  The array INFO(*) specifies the way 
!  the problem is to be solved.  WORK(*) is a work array. IDID reports 
!  success or the reason the computation has been terminated.
!
!  FIRST CALL TO RKC  
!
!  You must provide storage in your calling program for the arrays in the 
!  call list -- Y(NEQN), INFO(4), WORK(8+5*NEQN).  If INFO(2) = 0, you can
!  reduce the storage for the work array to WORK(8+4*NEQN).  ATOL may be 
!  a scalar or an array.  If it is an array, you must provide storage for 
!  ATOL(NEQN).  You must declare F in an external statement, supply the 
!  subroutine F and the function SPCRAD, and initialize the following 
!  quantities:
!
!    NEQN:  The number of differential equations. Integer.
!          
!    T:     The initial point of the integration. Double precision.
!           Must be a variable.
!
!    TEND:  The end of the interval of integration.  Double precision.
!           TEND may be less than T.
!          
!    Y(*):  The initial value of the solution.  Double precision array
!           of length NEQN.
!
!    F:     The name of a subroutine for evaluating the differential
!           equation.  It must have the form
!
!             subroutine f(neqn,t,y,dy)
!             integer          neqn
!             double precision t,y(neqn),dy(neqn)
!             dy(1)    = ...     
!             ...
!             dy(neqn) = ...
!             return
!             end
!
!  RTOL,
!  ATOL(*):  At each step of the integration the local error is controlled
!            so that its RMS norm is no larger than tolerances RTOL, ATOL(*).
!            RTOL is a double precision scalar. ATOL(*) is either a double 
!            precision scalar or a double precision array of length NEQN.
!            RKC is designed for the solution of problems to modest accuracy.
!            Because it is based on a method of order 2, it is relatively
!            expensive to achieve high accuracy.
!
!            RTOL is a relative error tolerance.  You must ask for some
!            relative accuracy, but you cannot ask for too much for the
!            precision available.  Accordingly, it is required that
!            0.1 >= RTOL >= 10*uround. (See below for the machine and
!            precision dependent quantity uround.) 
!
!            ATOL is an absolute error tolerance that can be either a 
!            scalar or an array.  When it is an array, the tolerances are
!            applied to corresponding components of the solution and when
!            it is a scalar, it is applied to all components.  A scalar
!            tolerance is reasonable only when all solution components are
!            scaled to be of comparable size.  A scalar tolerance saves a
!            useful amount of storage and is convenient.  Use INFO(*) to
!            tell RKC whether ATOL is a scalar or an array.
!            
!            The absolute error tolerances ATOL(*) must satisfy ATOL(i) >= 0
!            for i = 1,...,NEQN.  ATOL(j)= 0 specifies a pure relative error 
!            test on component j of the solution, so it is an error if this
!            component vanishes in the course of the integration.
!
!            If all is going well, reducing the tolerances by a factor of 
!            0.1 will reduce the error in the computed solution by a factor
!            of roughly 0.2.
!
!  INFO(*)   Integer array of length 4 that specifies how the problem
!            is to be solved.
!
!  INFO(1):  RKC integrates the initial value problem from T to TEND. 
!            This is done by computing approximate solutions at points 
!            chosen automatically throughout [T, TEND].  Ordinarily RKC
!            returns at each step with an approximate solution. These
!            approximations show how y behaves throughout the interval.
!            The subroutine RKCINT can be used to obtain answers anywhere
!            in the span of a step very inexpensively. This makes it 
!            possible to obtain answers at specific points in [T, TEND]
!            and to obtain many answers very cheaply when attempting to
!            locating where some function of the solution has a zero
!            (event location).  Sometimes you will be interested only in
!            a solution at TEND, so you can suppress the returns at each
!            step along the way if you wish.
!
!  INFO(1)  = 0 Return after each step on the way to TEND with a
!               solution Y(*) at the output value of T.
!
!           = 1 Compute a solution Y(*) at TEND only.
!
!  INFO(2):  RKC needs an estimate of the spectral radius of the Jacobian.
!            You must provide a function that must be called SPCRAD and 
!            have the form
!
!              double precision function spcrad(neqn,t,y)
!              integer          neqn
!              double precision t,y(neqn)
!
!              spcrad = < expression depending on info(2) >
!               
!              return
!              end
!          
!            You can provide a dummy function and let RKC compute the 
!            estimate. Sometimes it is convenient for you to compute in 
!            SPCRAD a reasonably close upper bound on the spectral radius,
!            using, e.g., Gershgorin's theorem.  This may be faster and/or 
!            more reliable than having RKC compute one.
! 
!  INFO(2)  = 0  RKC is to compute the estimate internally.  
!                Assign any value to SPCRAD.
!
!           = 1  SPCRAD returns an upper bound on the spectral
!                radius of the Jacobian of f at (t,y). 
!
!  INFO(3):  If you know that the Jacobian is constant, you should say so.
!
!  INFO(3)  = 0  The Jacobian may not be constant.
!
!           = 1  The Jacobian is constant.
!
!  INFO(4):  You must tell RKC whether ATOL is a scalar or an array.
!
!  INFO(4)  = 0  ATOL is a double precision scalar.
!
!           = 1  ATOL is a double precision array of length NEQN.
!             
!  WORK(*):  Work array.  Double precision array of length at least 
!            8 + 5*NEQN if INFO(2) = 0 and otherwise, 8 + 4*NEQN.
!
!  IDID:     Set IDID = 0 to initialize the integration.
!
!  icase:    Added by Gib.  Conveys the chemo that is being solved for.
!
!  RETURNS FROM RKC
!
!  T:        The integration has advanced to T.
!
!  Y(*):     The solution at T.
!
!  IDID:     The value of IDID reports what happened.
!
!                          SUCCESS
!
!    IDID     = 1 T = TEND, so the integration is complete.
!
!             = 2 Took a step to the output value of T.  To continue on
!                 towards TEND, just call RKC again.   WARNING:  Do not 
!                 alter any argument between calls.
! 
!                 The last step, HLAST, is returned as WORK(1). RKCINT 
!                 can be used to approximate the solution anywhere in 
!                 [T-HLAST, T] very inexpensively using data in WORK(*).
!
!                 The work can be monitored by inspecting data in RKCDID.
!
!                          FAILURE
!
!             = 3 Improper error control: For some j, ATOL(j) = 0
!                 and Y(j) = 0.
!
!             = 4 Unable to achieve the desired accuracy with the
!                 precision available.  A severe lack of smoothness in
!                 the solution y(t) or the function f(t,y) is likely.
!
!             = 5 Invalid input parameters:  NEQN <= 0, RTOL > 0.1, 
!                 RTOL < 10*UROUND, or ATOL(i) < 0 for some i.
!
!             = 6 The method used by RKC to estimate the spectral
!                 radius of the Jacobian failed to converge.
!
!  RKCDID is a labelled common block that communicates statistics
!         about the integration process:
!         common /rkcdid/   nfe,nsteps,naccpt,nrejct,nfesig,maxm
!
!         The integer counters are:
!
!        NFE      number of evaluations of F used
!                   to integrate the initial value problem
!        NSTEPS   number of integration steps
!        NACCPT   number of accepted steps
!        NREJCT   number of rejected steps
!        NFESIG   number of evaluations of F used 
!                   to estimate the spectral radius
!        MAXM     maximum number of stages used
!
!        This data can be used to monitor the work and terminate a run
!        that proves to be unacceptably expensive.  Also, if MAXM should
!        be far beyond 100, the problem is too expensive for RKC and
!        alternative methods should be considered.
!  
!--------------------------------------------------------------------------
!
!   CAUTION: MACHINE/PRECISION ISSUES
!
!     UROUND (the machine precision) is the smallest number such that 
!     1 + UROUND > 1, where 1 is a floating point number in the working 
!     precision. UROUND is set in a parameter statement in RKC. Its 
!     value depends on both the precision and the machine used, so it
!     must be set appropriately.  UROUND is the only constant in RKC
!     that depends on the precision.
!
!     This version of RKC is written in double precision. It can be changed
!     to single precision by replacing DOUBLE PRECISION in the declarations
!     by REAL and changing the type of the floating point constants set in
!     PARAMETER statements from double precision to real.
!            
!--------------------------------------------------------------------------
!
!  Authors: B.P. Sommeijer and J.G. Verwer
!           Centre for Mathematics and Computer Science (CWI)
!           Kruislaan 413
!           1098 SJ  Amsterdam
!           The Netherlands
!           e-mail: bsom@cwi.nl
!         
!           L.F. Shampine
!           Mathematics Department
!           Southern Methodist University
!           Dallas, Texas 75275-0156
!           USA
!           e-mail: lshampin@mail.smu.edu
!
!  Details of the methods used and the performance of RKC can be
!  found in 
!
!         B.P. Sommeijer, L.F. Shampine and J.G. Verwer
!         RKC: an Explicit Solver for Parabolic PDEs.
!              Technical Report MAS-R9715, CWI, Amsterdam, 1997
!
!  This source code for RKC and some examples, as well as the
!  reference solution to the second example can also be obtained
!  by anonymous ftp from the address ftp://ftp.cwi.nl/pub/bsom/rkc
!------------------------------------------------------------------
      type(rkc_comm) :: comm
      integer          neqn,info(*),idid, icase
      double precision y(neqn),t,tend,rtol,atol(*),work(*)

!*********************************************************************
!  uround is set here for IEEE double precision arithmetic.      
      double precision uround
      parameter       (uround=2.22d-16)
!*********************************************************************

      double precision zero,rmax,rmin
      parameter       (zero=0d0,rmax=0.1d0,rmin=10d0*uround)
      integer          i,ptr1,ptr2,ptr3,ptr4
      logical          array,valid
!      save
      integer          nfe,nsteps,naccpt,nrejct,nfesig,maxm
      common /rkcdid/  nfe,nsteps,naccpt,nrejct,nfesig,maxm
      external         f

      if(idid .eq. 0) then
!----------------------
!  Test the input data.
!----------------------  
         array = info(4) .eq. 1
         valid = neqn .gt. 0
         if((rtol .gt. rmax) .or. (rtol .lt. rmin)) valid = .false.
         if(atol(1) .lt. zero) valid = .false.
         if(array) then
           do i = 2, neqn
             if(atol(i) .lt. zero) valid = .false.
           enddo
         endif      
         if(.not. valid) then
           idid = 5
           return
         endif
!-----------------------------------
!  Initialize counters and pointers.
!-----------------------------------
!         nfe = 0
!         nsteps = 0
!         naccpt = 0
!         nrejct = 0
!         nfesig = 0
!         maxm = 0
         comm%nfe = 0
         comm%nsteps = 0
         comm%naccpt = 0
         comm%nrejct = 0
         comm%nfesig = 0
         comm%maxm = 0
!-----------------------------------------------------------
!  work(*) contains information needed for interpolation,
!  continuation after a return, and working storage. Items
!  relevant here are:
!
!  The last step taken, hlast, is work(1).
!  The current t is work(2).
!  The number of equations, neqn, is work(3).
!  The unit roundoff, uround, is work(4).
!  The square root of uround, sqrtu, is work(5).
!  The maximum step size, hmax, is work(6).
!  The base address for the solution is ptr1 = nint(work(7)).
!  The solution at t starts at ptr1.
!  The derivative of the solution at t starts at ptr2.
!  The solution at t-hlast starts at ptr3.
!  The derivative of the solution at t-hlast starts at ptr4.
!  The estimated dominant eigenvector starts at ptr4 + neqn.
!------------------------------------------------------------
        work(2) = t
        work(3) = neqn
        work(4) = uround
        work(5) = sqrt(uround)
        ptr1 = 8
        work(7) = ptr1
        ptr2 = ptr1 + neqn
        ptr3 = ptr2 + neqn
        ptr4 = ptr3 + neqn
      elseif(idid .ne. 2) then
        write(*,*) ' RKC was called with an illegal value of IDID.'
        stop
      endif
           
      call rkclow(comm,neqn,t,tend,y,f,info,rtol,atol,work, &
                  work(ptr1),work(ptr2),work(ptr3),work(ptr4),idid,icase)
      end subroutine
      
!----------------------------------------------------------------------
!  RKC is an interface to RKCLOW where the actual solution takes place.
!----------------------------------------------------------------------     
      subroutine rkclow(comm,neqn,t,tend,y,f,info,rtol,atol,work, &
                       yn,fn,vtemp1,vtemp2,idid,icase)
      type(rkc_comm) :: comm
      integer          neqn,info(*),idid, icase
      double precision t,tend,y(*),rtol,atol(*),work(*),yn(*),fn(*),vtemp1(*),vtemp2(*)
      external         f
     
      double precision one,onep1,onep54,p1,p4,p8,ten,zero,one3rd,two3rd
      parameter       (one=1d0,onep1=1.1d0,onep54=1.54d0,p1=0.1d0,p4=0.4d0,p8=0.8d0,ten=10d0, &
                       zero=0d0,one3rd=1d0/3d0,two3rd=2d0/3d0)
      integer          i,m,mmax,nstsig
      double precision absh,est,err,errold,fac,h,hmax,hmin,hold,sprad,tdir,temp1,temp2,	 & !spcrad,
                       uround,wt,ylast,yplast,at
      logical          array,last,newspc,jacatt
!      save
      integer          nfe,nsteps,naccpt,nrejct,nfesig,maxm
      common /rkcdid/  nfe,nsteps,naccpt,nrejct,nfesig,maxm
      
      interface 
		double precision function spcrad(neqn,t,y)
		integer :: neqn
		double precision :: t, y(neqn)
        end function
      end interface

!---------------------------------
!    Initialize on the first call.
!---------------------------------
      if(idid .eq. 0) then
        array = info(4) .eq. 1
        uround = work(4)
        mmax = nint(sqrt(rtol/(10d0*uround)))
        mmax = max(mmax,2)
        newspc = .true.
        jacatt = .false.
        nstsig = 0
!!$omp parallel do default(shared)
        do i = 1, neqn
          yn(i) = y(i)
        enddo
!!$omp end parallel do
        call f(neqn,t,yn,fn,icase)
!        nfe = nfe + 1
        comm%nfe = comm%nfe + 1
        tdir = sign(one,tend - t)
        hmax = abs(tend - t)
        work(6) = hmax
        hmin = ten*uround*max(abs(t),hmax)
        rkc_sprad = 0
      endif
!------------------------------------      
!  Start of loop for taking one step.  
!------------------------------------
      errold = 0
      hold = 0
20    continue      
!----------------------------------------------
!  Estimate the spectral radius of the Jacobian 
!  when newspc = .true..  A convergence failure 
!  in rkcrho is reported by idid = 6.
!----------------------------------------------
      if(newspc) then
        if(info(2) .eq. 1) then
          sprad = spcrad(neqn,t,yn)
        else
          call rkcrho(comm,neqn,t,f,yn,fn,vtemp1,vtemp2,work,sprad,idid,icase)
          
          rkc_sprad = max(sprad,rkc_sprad)
    
          if(idid .eq. 6) return
        endif
        jacatt = .true.
      endif
!-------------------------------
!  Compute an initial step size.
!-------------------------------        
      if(nsteps .eq. 0) then
        absh = hmax
        if(sprad*absh .gt. one) absh = one/sprad
        absh = max(absh,hmin)
!!$omp parallel do default(shared)
        do i = 1,neqn
          vtemp1(i) = yn(i) + absh*fn(i)
        enddo 
!!$omp end parallel do
        call f(neqn,t+absh,vtemp1,vtemp2,icase)
!        nfe = nfe + 1
        comm%nfe = comm%nfe + 1
        est = zero
        at = atol(1)
        do i = 1,neqn
          if(array) at = atol(i)
          wt = at + rtol*abs(yn(i))
          if(wt .eq. zero) then
            idid = 3
            return
          endif
          est = est + ((vtemp2(i) - fn(i))/wt)**2
        enddo         
        est = absh*sqrt(est/neqn)
        if(p1*absh .lt. hmax*sqrt(est)) then
          absh = max(p1*absh/sqrt(est), hmin)
        else 
          absh = hmax
        endif
      endif        
!------------------------------------------------------------
!  Adjust the step size and determine the number of stages m.
!------------------------------------------------------------      
      last = .false.
      if(onep1*absh .ge. abs(tend - t)) then
        absh = abs(tend - t)
        last = .true.
      endif
      m = 1 + int(sqrt(onep54*absh*sprad + one))
!----------------------------------------------------------
!  Limit m to mmax to control the growth of roundoff error.
!----------------------------------------------------------
      if(m .gt. mmax) then
        m = mmax
        absh = (m**2 - 1)/(onep54*sprad)
        last = .false.
      endif
!      maxm = max(m,maxm)
      comm%maxm = max(m,comm%maxm)
!--------------------------------------------
!  A tentative solution at t+h is returned in
!  y and its slope is evaluated in vtemp1(*).
!--------------------------------------------
      h = tdir*absh
      hmin = ten*uround*max(abs(t),abs(t + h))
      call step(neqn,f,t,yn,fn,h,m,y,vtemp1,vtemp2,icase)
      call f(neqn,t+h,y,vtemp1,icase)
!      nfe = nfe + m
      comm%nfe = comm%nfe + m
!      nsteps = nsteps + 1
      comm%nsteps = comm%nsteps + 1
!-------------------------------------------------------------
!  Estimate the local error and compute its weighted RMS norm.
!-------------------------------------------------------------
      err = zero
      at = atol(1)
      do i = 1, neqn
        if(array) at = atol(i)
        wt = at + rtol*max(abs(y(i)),abs(yn(i)))
        if(wt .eq. zero) then
          idid = 3
          return
        endif
        est = p8*(yn(i) - y(i)) + p4*h*(fn(i) + vtemp1(i))
        err = err + (est/wt)**2
      enddo
      err = sqrt(err/neqn)

      if(err .gt. one) then
!-------------------
!  Step is rejected.
!-------------------
!        nrejct = nrejct + 1
        comm%nrejct = comm%nrejct + 1
        absh = p8*absh/(err**one3rd)
        if(absh .lt. hmin) then
          idid = 4
          return
        else
          newspc = .not. jacatt
          goto 20
        endif
      endif
!-------------------
!  Step is accepted.
!-------------------
!      naccpt = naccpt + 1
      comm%naccpt = comm%naccpt + 1
      t = t + h
      jacatt = info(3) .eq. 1
      nstsig = mod(nstsig+1,25)
      newspc = .false.
      if((info(2) .eq. 1) .or. (nstsig .eq. 0)) newspc = .not. jacatt
!------------------------------------------------------
!  Update the data for interpolation stored in work(*).
!------------------------------------------------------      
      work(1) = h
      work(2) = t
!!$omp parallel do default(shared) private(ylast, yplast)
      do i = 1, neqn
         ylast = yn(i)
         yplast = fn(i)
         yn(i) = y(i)
         fn(i) = vtemp1(i)
         vtemp1(i) = ylast
         vtemp2(i) = yplast
      enddo
!!$omp end parallel do
      fac = ten
      if(naccpt .eq. 1) then
        temp2 = err**one3rd
        if(p8 .lt. fac*temp2) fac = p8/temp2
      else
        temp1 = p8*absh*errold**one3rd
        temp2 = abs(hold)*err**two3rd
        if(temp1 .lt. fac*temp2) fac = temp1/temp2
      endif
      absh = max(p1,fac)*absh
      absh = max(hmin,min(hmax,absh))
      errold = err
      hold = h
      h = tdir*absh
      if(last) then
        idid = 1
        return
      elseif(info(1) .eq. 0) then
        idid = 2
        return
      else
        goto 20
      endif
      end subroutine

      subroutine step(neqn,f,t,yn,fn,h,m,y,yjm1,yjm2,icase)
!--------------------------------------------------
!  Take a step of size H from T to T+H to get Y(*).
!--------------------------------------------------
      integer          neqn,m,icase
      double precision t,yn(neqn),fn(neqn),h,y(neqn),yjm1(neqn),yjm2(neqn)
      external         f
     
      double precision one,two,four,c13,zero
      parameter       (one=1d0,two=2d0,four=4d0,c13=13d0,zero=0d0)
      integer          i,j
      double precision ajm1,arg,bj,bjm1,bjm2,dzj,dzjm1,dzjm2,d2zj,d2zjm1,d2zjm2,mu,mus,nu, &
                       temp1,temp2,thj,thjm1,thjm2,w0,w1,zj,zjm1,zjm2

      w0 = one + two/(c13*m**2)
      temp1 = w0**2 - one
      temp2 = sqrt(temp1)
      arg = m*log(w0 + temp2)
      w1 = sinh(arg)*temp1 / (cosh(arg)*m*temp2 - w0*sinh(arg))
      bjm1 = one/(two*w0)**2
      bjm2 = bjm1
!---------------------------
!  Evaluate the first stage.
!---------------------------
      mus = w1*bjm1
!!$omp parallel do default(shared)
      do i = 1, neqn
        yjm1(i) = yn(i) + h*mus*fn(i)
        yjm2(i) = yn(i)
      enddo
!!$omp end parallel do
      thjm2  = zero
      thjm1  = mus
      zjm1   = w0
      zjm2   = one
      dzjm1  = one
      dzjm2  = zero
      d2zjm1 = zero
      d2zjm2 = zero
!------------------------------
!  Evaluate stages j = 2,...,m.
!------------------------------
      do j = 2, m
        zj   =   two*w0*zjm1 - zjm2
        dzj  =   two*w0*dzjm1 - dzjm2 + two*zjm1
        d2zj =   two*w0*d2zjm1 - d2zjm2 + four*dzjm1
        bj   =   d2zj/dzj**2
        ajm1 =   one - zjm1*bjm1
        mu   =   two*w0*bj/bjm1
        nu   = - bj/bjm2
        mus  =   mu*w1/w0
!---------------------------------------------
!  Use the y array for temporary storage here.         
!---------------------------------------------
        call f(neqn,t + h*thjm1,yjm1,y,icase)
!!$omp parallel do default(shared)
        do i = 1, neqn
          y(i) = mu*yjm1(i) + nu*yjm2(i) + (one - mu - nu)*yn(i) + h*mus*(y(i) - ajm1*fn(i))
        enddo
!!$omp end parallel do
        thj = mu*thjm1 + nu*thjm2 + mus*(one - ajm1)
!------------------------------------
!  Shift the data for the next stage.
!------------------------------------
        if(j .lt. m) then
!!$omp parallel do default(shared)
          do i = 1, neqn
            yjm2(i) = yjm1(i)
            yjm1(i) = y(i)
          enddo
!!$omp end parallel do
          thjm2  = thjm1
          thjm1  = thj
          bjm2   = bjm1
          bjm1   = bj
          zjm2   = zjm1
          zjm1   = zj
          dzjm2  = dzjm1
          dzjm1  = dzj
          d2zjm2 = d2zjm1
          d2zjm1 = d2zj
        endif
      enddo
      end subroutine

      subroutine rkcint(work,arg,yarg)
!-------------------------------------------------------------------------
!  RKCINT is used to compute approximate solutions at specific t and to 
!  compute cheaply the large number of approximations that may be needed
!  for plotting or locating when events occur.
!
!  After a step to T, RKC provides HLAST, the step just taken, in WORK(1).
!  In other entries of WORK(*) it provides the data needed to interpolate 
!  anywhere in [T-HLAST, T]. YARG(*), the approximate solution at t = ARG 
!  computed by interpolation in RKCINT has the same order of accuracy as 
!  the Y(*) computed directly by RKC.
!  
!  INPUT:
!
!    WORK(*)   Double precision array returned by RKC.
!
!    ARG       The point at which a solution is desired. Double precision.
!
!  OUTPUT:
!
!    YARG(*)   The approximate solution at t = ARG.  Double precision
!              array of length neqn.
!--------------------------------------------------------------------------
      double precision work(*),arg,yarg(*)

      double precision one,two,three
      parameter       (one=1d0,two=2d0,three=3d0)
      integer          i,neqn,ptr1,ptr2,ptr3,ptr4
      double precision a1,a2,b1,b2,s,hlast,t,tlast
     
!---------------------------------------------------------------------
!  The data needed for interpolation are stored in work(*) as follows:
!
!  The last step taken, hlast, is work(1).
!  The current t is work(2).
!  The number of equations, neqn, is work(3).
!  The base address for the solution is ptr1 = nint(work(7))
!  The solution at t starts at ptr1.
!  The derivative of the solution at t starts at ptr2.
!  The solution at t-hlast starts at ptr3.
!  The derivative of the solution at t-hlast starts at ptr4.
!---------------------------------------------------------------------
      hlast = work(1)
      t = work(2)
      tlast = t - hlast
      neqn = nint(work(3))
      ptr1 = nint(work(7))
      ptr2 = ptr1 + neqn
      ptr3 = ptr2 + neqn
      ptr4 = ptr3 + neqn

      s  = (arg - tlast)/hlast
      a1 = (one + two*s)*(s - one)**2
      a2 = (three - two*s)*s**2
      b1 = hlast*s*(s - one)**2
      b2 = hlast*(s - one)*s**2

!!$omp parallel do default(shared)
      do i = 1, neqn
        yarg(i) = a1*work(ptr3+i-1) + a2*work(ptr1+i-1) + b1*work(ptr4+i-1) + b2*work(ptr2+i-1)
      enddo
!!$omp end parallel do
      end subroutine

      subroutine rkcrho(comm,neqn,t,f,yn,fn,v,fv,work,sprad,idid,icase)
!---------------------------------------------------------------
!  RKCRHO attempts to compute a close upper bound, SPRAD, on
!  the spectral radius of the Jacobian matrix using a nonlinear
!  power method.  A convergence failure is reported by IDID = 6.
! Gib added sprad_factor
!---------------------------------------------------------------
      type(rkc_comm) :: comm
      integer          neqn,idid,icase
      double precision t,yn(neqn),fn(neqn),v(neqn),fv(neqn),work(*),sprad
      external         f
     
      integer          itmax
      parameter       (itmax=50)
      double precision zero,one,onep2,p01
      parameter       (zero=0d0,one=1d0,onep2=1.2d0,p01=0.01d0)
      integer          i,iter,index,ptr5
      double precision uround,sqrtu,ynrm,sigma,sigmal,dynrm,dfnrm,vnrm,small
      integer          nfe,nsteps,naccpt,nrejct,nfesig,maxm
	  real(8), parameter :: sprad_factor = 1
      common /rkcdid/  nfe,nsteps,naccpt,nrejct,nfesig,maxm

      uround = work(4)
      sqrtu = work(5)
!------------------------------------------------------------
!  hmax = work(6).  sprad smaller than small = 1/hmax are not
!  interesting because they do not constrain the step size.  
!------------------------------------------------------------
      small = one/work(6)
!---------------------------------------------------------
!  The initial slope is used as guess when nsteps = 0 and
!  thereafter the last computed eigenvector.  Some care
!  is needed to deal with special cases. Approximations to
!  the eigenvector are normalized so that their Euclidean
!  norm has the constant value dynrm.
!---------------------------------------------------------
      ptr5 = nint(work(7)) + 4*neqn
      if(nsteps .eq. 0) then
        do i = 1,neqn
          v(i) = fn(i)
        enddo
      else
        do i = 1,neqn          
          v(i) = work(ptr5+i-1)
        enddo
      endif
      ynrm = zero
      vnrm = zero
      do i = 1,neqn
        ynrm = ynrm + yn(i)**2     
        vnrm = vnrm + v(i)**2
      enddo
      ynrm = sqrt(ynrm)
      vnrm = sqrt(vnrm)
      if(ynrm .ne. zero .and. vnrm .ne. zero) then
        dynrm = ynrm*sqrtu
        do i = 1,neqn
          v(i) = yn(i) + v(i)*(dynrm/vnrm)
        enddo          
      elseif(ynrm .ne. zero) then
        dynrm = ynrm*sqrtu
        do i = 1, neqn
          v(i) = yn(i) + yn(i)*sqrtu
        enddo
      elseif(vnrm .ne. zero) then
        dynrm = uround
        do i = 1,neqn
          v(i) = v(i)*(dynrm/vnrm)
        enddo
      else          
        dynrm = uround
        do i = 1,neqn
          v(i) = dynrm
        enddo          
      endif
!--------------------------------------------
!  Now iterate with a nonlinear power method.
!--------------------------------------------
      sigma = zero
      do iter = 1, itmax
        call f(neqn,t,v,fv,icase)
!        nfesig = nfesig + 1
        comm%nfesig = comm%nfesig + 1
        dfnrm = zero
        do i = 1, neqn
          dfnrm = dfnrm + (fv(i) - fn(i))**2
        enddo   
        dfnrm = sqrt(dfnrm)
        sigmal = sigma
        sigma = dfnrm/dynrm
!----------------------------------------------------------
!  sprad is a little bigger than the estimate sigma of the 
!  spectral radius, so is more likely to be an upper bound.
!----------------------------------------------------------
        sprad = onep2*sigma
        sprad = sprad * sprad_factor	! Gib
        if(iter >= 2 .and. abs(sigma - sigmal) <= max(sigma,small)*p01) then
          do i = 1,neqn          
            work(ptr5+i-1) = v(i) - yn(i)
          enddo
          return
        endif
!--------------------------------------
!  The next v(*) is the change in f
!  scaled so that norm(v - yn) = dynrm.
!--------------------------------------
        if(dfnrm .ne. zero) then
          do i = 1,neqn
            v(i) = yn(i) + (fv(i) - fn(i))*(dynrm/dfnrm)
          enddo           
        else
!-------------------------------------------------------
!  The new v(*) degenerated to yn(*)--"randomly" perturb 
!  current approximation to the eigenvector by changing 
!  the sign of one component.
!-------------------------------------------------------
          index = 1 + mod(iter,neqn)
          v(index) = yn(index) - (v(index) - yn(index))
        endif
      enddo
!-------------------------------------------
!  Set flag to report a convergence failure.
!-------------------------------------------
      idid = 6
      end subroutine

end module