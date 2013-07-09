// Copyright (C) 2010 - 2013
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2013-03-27 td
// last changed:
//

/*
C
C  ---------------------------------------------------------------------
C
C* Title
C
C    Integrator for stiff systems of autonomous ordinary differential
C    equations.
C
C* Written by        P. Deuflhard, U. Nowak, U. Poehle
C* Purpose           Solution of systems of initial value problems
C* Method            Semi-implicit mid-point rule with
C                    h**2-extrapolation
C* Category          i1a2a. - System of stiff first order differential
C                             equations
C* Keywords          extrapolation, ODE, mid-point rule, stiff
C* Version           1.1 , July 1989
C* Latest Change     February 1991
C* Library           CodeLib
C* Code              Fortran 77
C                    Double Precision
C* Environment       Standard version for FORTRAN77 environments on
C                    PCs, workstations, and hosts
C* Copyright     (c) Konrad-Zuse-Zentrum fuer Informationstechnik
C                    Berlin (ZIB)
C                    Takustrasse 7, D-14195 Berlin-Dahlem
C                    phone : + 49/30/84185-0
C                    fax   : + 49/30/84185-125
C* Contact           Uwe Poehle
C                    ZIB, Scientific Software Group
C                    phone : + 49/30/84185-241
C                    fax   : + 49/30/84185-107
C                    e-mail: poehle@zib.de
C
C  ---------------------------------------------------------------------
C
C* Licence
C  -------
C
C  You may use or modify this code for your own non-commercial
C  purposes for an unlimited time.
C  In any case you should not deliver this code without a special
C  permission of ZIB.
C  In case you intend to use the code commercially, we oblige you
C  to sign an according licence agreement with ZIB.
C
C
C* Warranty
C  --------
C
C  This code has been tested up to a certain level. Defects and
C  weaknesses, which may be included in the code, do not establish
C  any warranties by ZIB. ZIB does not take over any liabilities
C  which may follow from aquisition or application of this code.
C
C
C* Software status
C  ---------------
C
C  This code is under care of ZIB and belongs to ZIB software
C  class I.
C
C
C  ---------------------------------------------------------------------
C
C* Short doc:
C
C  Jacobian approximation by numerical differences
C  or user supplied subroutine M1JAC (see constant QJACUS)
C
C  The numerical solution of the arising linear equations is done by
C  means of the subroutines DGEFA and DGESL.  For special purposes these
C  routines may be substituted.
C
C
C* References:
C
C /1/ P. Deuflhard:
C     A Semi-Implicit Midpoint Rule for Stiff Systems of Ordinary
C     Differential Equations
C     Num. Math. 41, 373 - 398 (1983) .
C
C /2/ P. Deuflhard:
C     Order and Stepsize Control in Extrapolation Methods
C     Numer. Math. 41, 399-422 (1983)
C
C /3/ P. Deuflhard:
C     Uniqueness Theorems for Stiff ODE Initial Value Problems
C     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin,
C     Preprint SC-87-3 (1987)
C
C
C* External subroutine: (to be supplied by the user)
C
C    FCN           EXT  Subroutine FCN(N,T,Y,DY,IFAIL)
C                       Right-hand side of first-order
C                       differential equations
C                       N      Number of first-order ODE's
C                       T      Actual position
C                       Y(N)   Values at T
C                       DY(N)  Derivatives at T
C                       IFAIL  Error return code
C
C    SOUT          EXT  Subroutine SOUT(N,TOLD,T,Y)
C                       Callback for intermediate solution points
C                       Only called if KFLAG = 4
C                       N      Number of first-order ODEs
C                       TOLD   Position before
C                       T      Actual intermediate position
C                       Y(N)   Values at T
C
C
C* Parameters: (* marks transient parameters)
C
C    N         I   IN   Number of ODE'S
C  * T         D   IN   Starting point of integration
C                       (T .LE. TEND)
C                  OUT  Achieved final point of integration
C  * Y         D   IN   Array of initial values Y(1),...,Y(N)
C                  OUT  Array of final values
C    TEND      D   IN   Prescribed final point of integration
C    TOL       D   IN   Prescribed relative precision (.GT.0)
C    HMAX      D   IN   Maximum permitted stepsize
C  * H         D   IN   Initial stepsize guess
C                  OUT  Stepsize proposal for next integration step
C                       (H .EQ. 0. ,if METAN1 fails to proceed)
C  * KFLAG     I   IN   Print parameter
C                        0   No output
C                        1   Integration monitor
C                        2   Intermediate solution points  T, Y(I),I=1,N
C                        3   Integration monitor and intermediate points
C                        4   Intermediate points via callback routine SOUT
C                        5   Monitor and use of callback routine SOUT
C                  OUT  Error flag
C                       .GE. 0  Successful integration
C                               (KFLAG not altered internally)
C                       -1   TEND .LT. T
C                       -2   More than NSTMAX basic integration steps
C                            per interval have been performed
C                       -3   More than JRMAX stepsize reductions
C                            occurred per basic integration step
C                       -4   Stepsize proposal for next basic
C                            integration too small
C
*/

typedef void (*FcnMetan)(int* n,
                           double* t, double* y, double* dy,
                           int* ifail);

typedef void (*SoutMetan)(int* n,
                            double* tOld, double* t, double* y, double* dy);

// original METAN API

extern void metan1_(
            int*       n,
            FcnMetan   fcn,
            double*   t0,
            double*   y0,
            double*   tEnd,
            double*   tol,
            double*   hMax,
            double*   h,
            int*      kFlag,
            SoutMetan  sout);

