// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-02-08 td
// last changed:
//

/*
c
c-----------------------------------------------------------------------
c
c     Extrapolation integrator  for the solution of  linearly-implicit
c     differential-algebraic systems of the form
c
c          B (t,y) * y' (t) = f (t,y)
c
c     with B a (n,n)-matrix of rank less or equal n.
c
c-----------------------------------------------------------------------
c
c     Copyright (C) 2000 - 2004,
c     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
c     ALL RIGHTS RESERVED
c
c     Written by:
c
c     R. Ehrig, U. Nowak,
c     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
c     Takustrasse 7
c     D-14195 Berlin-Dahlem
c
c     phone : +49-30-84185-0
c     fax   : +49-30-84185-125
c     e-mail: ehrig@zib.de, nowak@zib.de
c     URL   : http://www.zib.de
c
c-----------------------------------------------------------------------
c
c            *************************************************
c            **                                             **
c            **  This is version 4.3A of January, 31, 2004  **
c            **                                             **
c            *************************************************
c
c-----------------------------------------------------------------------
c
c     Overview of current versions:
c
c     4.3A   Non-sparse dense  or  banded  Jacobians. Direct  solvers:
c            LAPACK routines. Based on BLAS and LAPACK routines.
c
c     4.3B   Sparse Jacobians. Direct  solvers: MA28 routines from the
c            Harwell subroutine library, iterative  solvers: GMRES and
c            BICGSTAB  with a  variable ILU  preconditioner. Based  on
c            BLAS and LAPACK routines.
c
c     Versions with other solver routines are available on request.
c
c-----------------------------------------------------------------------
c
c     NOTICE: "The LIMEX  program may be used SOLELY  for educational,
c     research, and benchmarking  purposes by no-profit organizations.
c     Commercial and other organizations  may make use of LIMEX SOLELY
c     for benchmarking  purposes only. LIMEX may  be modified by or on
c     behalf of the user for  such use but  at no time  shall LIMEX or
c     any such  modified version  of LIMEX become the  property of the
c     user. LIMEX  is provided  without warranty  of any  kind, either
c     expressed  or implied.  Neither the  authors nor their employers
c     shall be liable  for any direct or consequential  loss or damage
c     whatsoever  arising out  of the  use or  misuse of  LIMEX by the
c     user. LIMEX  must not  be sold. You  may make  copies  of LIMEX,
c     but  this NOTICE  and the  Copyright notice  must appear  in all
c     copies. Any  other  use  of LIMEX  requires written  permission.
c     Your use of LIMEX is an implicit agreement to these conditions."
c
c     LIMEX also is  available for commercial use. Contact the authors
c     who will provide details of price and condition of use.
c
c-----------------------------------------------------------------------
c
c     Description
c
c     LIMEX  solves linearly-implicit  differential-algebraic  systems
c     (DAEs) of the form
c
c        B (t,y) * y' (t) = f (t,y).
c
c     If n is the size of the system, B is a (n,n)-matrix of rank less
c     or equal n. General  conditions for  the applicability  of LIMEX
c     are a regular  matrix pencil B + h A with A the Jacobian of  the
c     residual of the DAE and an index of the DAE less or equal 1. The
c     discretization  of LIMEX  is based  on the  elementary  linearly
c     implicit Euler discretization
c
c        ( B(t,y(k)) - h J ) ( y(k+1) - y(k) ) = h f(t(k+1),y(k))
c
c     where J is the (approximate) Jacobian of the residual
c
c         d              |
c        -- ( f - B y' ) |      .
c        dy              |t=t(0)
c
c     Combined  with  extrapolation  this one-step  method permits  an
c     adaptive control of stepsize and order. Within the extrapolation
c     process  one  computes  for  a  basic stepsize H  approximations
c     T(j,1) for y(t(0)+H) using  the discretization  above with step-
c     sizes  h(j) = H / j, j = 1, ..., j(max). Then the  extrapolation
c     tableau recursively defines higher order approximations T(j,k):
c
c                            T(j,k-1) - T(j-1,k-1)
c        T(j,k) = T(j,k-1) + ---------------------, k = 1, ..., j.
c                            j / ( j - k + 1 ) - 1
c
c     As error estimates the subdiagonal differences T(j,j) - T(j,j-1)
c     are taken.
c
c     The efficiency of LIMEX mainly depends on the performance of the
c     evaluation of the Jacobian  and in particular on the solution of
c     the  linear systems. Therefore  different versions  of LIMEX are
c     available with  different methods to create the  Jacobian and to
c     solve the linear systems.
c
c     Jacobian generation: numerical or analytical (user supplied)
c
c     Linear systems solvers: direct for dense matrices, LIMEX version
c     4.3A, direct or iterative methods (GMRES or BICGSTAB) for sparse
c     matrices, LIMEX version 4.3B. The  code  for dense  matrices de-
c     termines the optimal matrix representation, full or banded, from
c     the  ratio between  the system  size and the number of lower and
c     upper diagonals in the Jacobian.
c
c     Throughout  the code a local error control is implemented, which
c     requires that the  local error of a component y(i) is  less than
c     rTol * abs ( y(i) ) + aTol(i). This approach enables  to specify
c     more or less sensitive components of the solution vector.
c
c     The code has an additional option to compute  consistent initial
c     values (CIVs). For this  task an  extrapolated Newton  iteration
c     is used to  determine correct  start values at t = t(0) for  the
c     algebraic variables (which are known only  approximately in many
c     applications) and the derivatives  of the differential variables
c     (which are unknown in most applications), which  satisfy exactly
c     the  DAE. The algebraic  variables are these  components of  the
c     solution  vector, for  which  the  corresponding columns  in the
c     matrix B(t,y) contain zeros only. All other  components will  be
c     accounted  as  differential  variables,  whose  values  are  not
c     changed during the CIV computation. However, in many cases where
c     only the  right derivatives at  the start of the integration are
c     unknown, LIMEX solves the DAE even without CIV determination.
c
c     Furthermore a dense output option is availabe. With this  option
c     one may compute solutions of the DAE on a dense set of points in
c     the integration interval with nearly the same accuracy as at the
c     automatically  selected integration  points. The method bases on
c     Hermitian interpolation and considers all intermediate solutions
c     computed in  the extrapolation process. The dense  output option
c     may be combined with the  one-step mode, i.e. the return  to the
c     calling  program after  every  integration  step, and  also with
c     continuation  calls. These options  enable to obtain efficiently
c     solutions  on an  arbitrary user defined  set of values  for the
c     independent variable.
c
c     If information  about the  structure of the  Jacobian  matrix is
c     needed, the  matrix in a specific  integration step is available
c     in a PostScript plot.
c
c     Since in many  applications some components of the solution have
c     to be nonnegative, due to their physical interpretation, one may
c     set for these  components that only positive values ( => 0 ) are
c     allowed. If in the extrapolation negative values  are occurring,
c     the stepsize is reduced by a heuristic factor.
c
c-----------------------------------------------------------------------
c
c     References:
c
c     1) P. Deuflhard, U. Nowak:
c        Extrapolation integrators for quasilinear implicit ODEs.
c        In: P. Deuflhard, B. Engquist (eds.): Large  scale scientific
c        computing, Birkhaeuser, Prog. Sci. Comp. 7, pp. 37-50 (1987)
c
c     2) P. Deuflhard, E. Hairer, J. Zugck:
c        One step and extrapolation methods for differential-algebraic
c        systems.
c        Num. Math. 51, pp. 501-516 (1987)
c
c     3) P. Deuflhard:
c        Recent progress in extrapolation methods for ODEs.
c        SIAM Rev. 27, pp. 505-535 (1985)
c
c     4) R. Ehrig, U. Nowak, L. Oeverdieck, P. Deuflhard:
c        Advanced extrapolation  methods for large  scale differential
c        algebraic problems.
c        In: High  Performance Scientific  and Engineering  Computing,
c        H.-J. Bungartz,  F. Durst, and  Chr. Zenger  (eds.),  Lecture
c        Notes  in  Computational  Science  and  Engineering, Springer
c        Vol 8, pp. 233-244 (1999)
c
c-----------------------------------------------------------------------
c
c     Installation notes:
c
c     Requires the BLAS (Basic Linear Algebra Subprograms) library and
c     the LAPACK  library. Ideally, one  should  use  vendor-optimized
c     BLAS and LAPACK routines. If these are not available, one should
c     use the routines collected in source form in the file
c
c        'LIMEX4_3_Auxiliaries.f'
c
c     This  file is  part  of  the  LIMEX4_3A distribution. The  whole
c     Fortran BLAS and LAPACK libraries are  available from the netlib
c     repository  at http://cm.bell-labs.com/netlib/master/readme.html
c     or one of its mirrors.
c
c     This code is written almost  purely in ANSI Fortran 77 standard,
c     only the following non-standard features are used:
c
c        do ...  end do  constructions (without statement labels)
c
c        implicit none   to use only explicitly declared variables
c
c        include '....'  for inserting code
c
c        symbolic names may be longer than 6 characters
c
c        lower case and underscore characters used
c
c     All dimension statements in LIMEX  refer to the size definitions
c     in the include file 'LIMEX4_3_Size_Definitions.h'. This file has
c     to consist of four parameter statements of the following form:
c
c-----------------------------------------------------------------------
c
c     Define the  maximum number of  equations. This statement sets an
c     upper limit for the number of equations defining the DAE.
c
c     parameter ( Max_Nr_of_Equations    =
c
c-----------------------------------------------------------------------
c
c     Define the maximum number  of non-zero entries in the  Jacobian.
c     This statement defines an upper limit for the number of non-zero
c     entries  in  the  Jacobian of  the  residual  of  the  DAE. This
c     parameter is not used  in the LIMEX version 4.3A, it is included
c     here  only for consistency of the include files of the different
c     LIMEX versions.
c
c     parameter ( Max_Non_Zeros_Jacobian =
c
c-----------------------------------------------------------------------
c
c     Define the  maximum  number of non-zero entries in the left-hand
c     side matrix B.
c
c     parameter ( Max_Non_Zeros_B       =
c
c-----------------------------------------------------------------------
c
c     Define the  maximum number of  lower diagonals in  the Jacobian,
c     i.e. the lower  half-bandwith of  this matrix. The value of this
c     parameter  must  be  between 0 (no lower  diagonals) and n (full
c     lower triangular submatrix possible). The actual number of upper
c     diagonals is defined  in LIMEX by the control parameter Iopt(8).
c     This parameter is not used in the LIMEX version 4.3B.
c
c     parameter ( Max_Lower_Diagonals    =
c
c-----------------------------------------------------------------------
c
c     Define the  maximum number of  upper diagonals in  the Jacobian,
c     i.e. the upper  half-bandwith of  this matrix. The value of this
c     parameter  must  be  between 0 (no upper  diagonals) and n (full
c     upper triangular submatrix possible). The actual number of upper
c     diagonals is defined in LIMEX by the control  parameter Iopt(9).
c     This parameter is not used in the LIMEX version 4.3B.
c
c     parameter ( Max_Upper_Diagonals    =
c
c-----------------------------------------------------------------------
c
c     Define the maximum number of vectors which are  available within
c     GMRES or  BICGSTAB. For GMRES  this statement  implies  an upper
c     limit for the dimension of the Krylov  subspaces. This dimension
c     is  defined by  Iopt(20), see  below. GMRES  then requires  that
c     Max_It_Vectors => Iopt(20) + 4 (this  relation  is  checked). If
c     the Krylov subspace is exhausted without  attaining the required
c     accuracy, a  restart is done. Reasonable values for Iopt(20) are
c     in many cases 10 or 20, then Max_It_Vectors should  be set to 15
c     resp. 25. The greater  the dimension  of the Krylov subspace is,
c     the  faster GMRES  converges  (in theory), but  needs a  rapidly
c     increasing number of floating point  operations and is also less
c     stable. If  BICGSTAB is  used, Max_It_Vectors must be 6 at least
c     (this  condition is  checked). If  no iterative  method will  be
c     used, Max_It_Vectors may  be set equal 1. This parameter  is not
c     used in  the LIMEX version 4.3A, it  is included  here only  for
c     consistency  of   the  include  files  of  the  different  LIMEX
c     versions.
c
c     parameter ( Max_It_Vectors         =
c
c-----------------------------------------------------------------------
c
c     Arguments
c
c-----------------------------------------------------------------------
c
c     All restrictions  are  checked, but  no further  control of  the
c     consistency of the input data will be performed.
c
c     n          Integer variable, must be set by caller on entry, not
c                modified. Size of the  differential algebraic system.
c                Restriction: Max_Nr_of_Equations => n => 1.
c
c     Fcn        Name of an  external subroutine computing  the values
c                of  the  right-hand  side  function  f(t,y)  and  the
c                entries in the matrix B (see below).
c
c     Jacobian   Name of an  external subroutine computing  analytical
c                derivatives of the residual f(t,y) - B(t,y) * y' (see
c                below).
c
c     t_Begin    Real variable, must be set by caller on  entry. Inde-
c                pendent variable, starting point  of the integration.
c                On exit, the  value  of the  independent value, up to
c                LIMEX has solved the DAE successfully.
c
c     t_End      Real variable, must  be set  by caller on  entry, not
c                modified. Independent  variable, final  point  of the
c                integration.
c
c     y          Real  array  of  size  n, must  be set  by caller  on
c                entry. The  initial  values  of the  solution at  the
c                starting  point  t_Begin. On  exit,  y  contains  the
c                solution  at  the  final  point t_End, resp. at  this
c                value of the independent variable, up  to which LIMEX
c                has solved the DAE successfully.
c
c     ys         Real array of size n. Derivatives  of the solution at
c                t_Begin. Only the values  of ys for  the differential
c                variables  are needed. If the  values for ys are  not
c                known, set ys to  zero. If only estimates are  known,
c                these  should be  used. However, the  option  Iopt(6)
c                for the  computation of consistent initial values can
c                be used to  start always with  correct values for ys.
c                On exit, ys contains  the derivatives at t_End, resp.
c                at this  value  of  the  independent  variable, up to
c                LIMEX has solved the DAE successfully.
c
c     rTol       Real  array  of size 1 at least for Iopt(11) = 0 or n
c                for Iopt(11) = 1. Values  must be  set by  caller  on
c                entry,  not  modified.  Relative  error   tolerances,
c                interpreted as a scalar or a vector  depending on the
c                value of Iopt(11). The code keeps the local  error of
c                y(i)  below  rTol(i)*abs(y(i))+aTol(i).
c
c                Restriction: rTol(*) > 0.
c
c     aTol       Real array  of size 1 at least  for Iopt(11) = 0 or n
c                for  Iopt(11) = 1. Values  must be set by  caller  on
c                entry,  not  modified.  Absolute  error   tolerances,
c                interpreted as a scalar or a vector  depending on the
c                value of Iopt(11).
c
c                Restriction: aTol(*) => 0.
c
c                Remark: For any  component y(i), which could  be zero
c                during the integration  one has to set aTol(i) > 0 to
c                prevent  divisions  by  zero.  In  many  well  scaled
c                problems a reasonable setting is aTol(i) = rTol(i).
c
c     h          Real  variable, must  be  set  by  caller  on  entry.
c                Initial  stepsize  guess, if h  is  set  to zero, the
c                program puts  h = h_Min. This  variable is  currently
c                set  in a  parameter statement to 1.d0-4, but  may be
c                changed. On exit h is the  estimated optimal stepsize
c                for the next  integration step. If an  error occurred
c                h is set to zero.
c
c     Iopt       Integer  array  of  size 30, values  from  Iopt(1) to
c                Iopt(18) must be set by caller on  entry. Integration
c                control parameters.
c
c                Iopt(1)  Integration monitoring, not modified
c                         = 0 : no output
c                         = 1 : standard output
c                         = 2 : additional integration monitor
c
c                         Restriction: 0 <= Iopt(1) <= 2.
c
c                Iopt(2)  The  unit  number  for  monitor  output, not
c                         modified. If Iopt(2) = 0, the  default value
c                         for the unit number is 6.
c
c                         Restriction: 0 <= Iopt(2) if Iopt(1) > 0.
c
c                Iopt(3)  Solution output, not modified
c                         = 0 : no solution output
c                         = 1 : initial and final solution values
c                         = 2 : additional  solution  values at inter-
c                               mediate points
c
c                         Restriction: 0 <= Iopt(3) <= 2.
c
c                Iopt(4)  The  unit  number for  solution  output, not
c                         modified. If Iopt(4) = 0, the  default value
c                         for the unit number is 6.
c
c                         Restriction: 0 <= Iopt(4) if Iopt(3) > 0.
c
c                Iopt(5)  Singular or nonsingular  matrix B, not modi-
c                         fied
c                         = 0 : matrix B may be singular
c                         = 1 : matrix B is nonsingular
c                         For theoretical  reasons Iopt(5) = 0 reduces
c                         the  maximum  order  to 5, even  if  in most
c                         cases  LIMEX  also  runs  successfully  with
c                         Iopt(5) = 1  and  singular  matrices B(t,y).
c
c                         Restriction: 0 <= Iopt(5) <= 1.
c
c                Iopt(6)  Determination  of consistent  initial values
c                         (CIVs), not modified
c                         = 0 : no determination of CIVs
c                         = 1 : determination of CIVs
c                         If Iopt(6) = 1, the code computes before the
c                         integration  consistent initial  values. For
c                         this  reason, new  values for  the algebraic
c                         variables  and  for the  derivatives  of the
c                         differential  variables are  computed, which
c                         satisfy the DAE up to the required accuracy.
c
c                         Restriction: 0 <= Iopt(6) <= 1.
c
c                Iopt(7)  Numerical  or analytical  generation  of the
c                         Jacobian matrix, not modified
c                         = 0 : Numerical  difference approximation of
c                               Jacobian of the residual
c                         = 1 : Analytic Jacobian supplied by the user
c                               subroutine Jacobian
c
c                         Restriction: 0 <= Iopt(7) <= 1.
c
c                Iopt(8)  Lower bandwidth of the  Jacobian matrix, not
c                         modified. If  the Jacobian is a band matrix,
c                         the computing time needed for the evaluation
c                         of this matrix can be substantially reduced.
c                         If Iopt(8) => n, the whole lower  triangular
c                         submatrix will be computed. If Iopt(8) (and/
c                         or Iopt(9)) is less  than n, the Jacobian is
c                         a band  matrix. Then LIMEX uses  specialized
c                         subroutines  for the factorization  and  the
c                         solution of linear systems, which are faster
c                         than  the corresponding  algorithms for full
c                         matrices. Also the storage needed is reduced
c                         if  the parameters  Max_Upper_Diagonals  and
c                         Max_Lower_Diagonals are  set properly in the
c                         include file LIMEX4_3_Size_Definitions.h.
c
c                         Restriction:
c                         0 <= Iopt(8) <= Max_Lower_Diagonals.
c
c                Iopt(9)  Upper bandwidth of the  Jacobian matrix, not
c                         modified. If  the Jacobian is a band matrix,
c                         the computing time needed for the evaluation
c                         of this matrix can be substantially reduced.
c                         If Iopt(9) => n, the whole upper  triangular
c                         submatrix will be computed. If Iopt(9) (and/
c                         or Iopt(8)) is less  than n, the Jacobian is
c                         a band matrix. Read then the hints above for
c                         Iopt(8).
c
c                         Restriction:
c                         0 <= Iopt(9) <= Max_Upper_Diagonals.
c
c                Iopt(10) Determines whether the code tries to use the
c                         Jacobians in more than one integration step,
c                         not modified. This can lead to a significant
c                         increase of the overall performance.
c                         = 0 : no reuse of the Jacobian
c                         = 1 : reuse of the  Jacobian in the follwing
c                               integration steps
c                         If Iopt(10) = 1, the  code  computes the  so
c                         called contractivity factor, an estimate for
c                         the current linearity of  the DAE. Depending
c                         from an  upper  bound for  this  factor, the
c                         current Jacobian  is used also  in the  next
c                         integration step. The upper bound is defined
c                         by ThMin. This  variable  is set  within the
c                         parameter statement  below and should always
c                         selected between 0 (no  reuse the  Jacobian)
c                         and 1 (try almost in every step to reuse the
c                         Jacobian). A recommended value for the first
c                         computations  is 2.0d-2. If  in the  current
c                         step a  new Jacobian  is computed, the  step
c                         number is marked in the monitor  output with
c                         an asterisk '*'. If  the stepsize is reduced
c                         due to  an unsatisfying  convergence  within
c                         the extrapolation, the true current Jacobian
c                         is recomputed if necessary.
c
c                         Restriction: 0 <= Iopt(10) <= 1.
c
c                Iopt(11) Switch  for error  tolerances, not modified.
c                         = 0 : both rTol and aTol are scalars
c                         = 1 : both rTol and aTol are vectors
c
c                         Restriction: 0 <= Iopt(11) <= 1.
c
c                Iopt(12) Switch for the one step mode, not  modified.
c                         If Iopt(12) is equal 1, LIMEX returns during
c                         the integration to the calling program after
c                         each integration  step. If  Iopt(12) = 2 and
c                         the  dense  output  option is  active , i.e.
c                         Iopt(13) > 0, the code  returns only  on the
c                         specified dense  output points (see  below).
c                         If LIMEX returns during the integration, the
c                         following  parameters of  LIMEX contain  the
c                         current  values and  may be used  within the
c                         calling program:
c
c                         t_Begin: the current value of t (independent
c                                  variable)
c                         y      : the current solution
c
c                         To  assure  a correct  continuation  of  the
c                         integration, no parameter of LIMEX should be
c                         changed  before the  next call  of LIMEX. If
c                         Iopt(12) = 0, LIMEX does  not return  to the
c                         calling program  until the whole integration
c                         is finished.
c
c                         Restriction: 0 <= Iopt(12) <= 2.
c
c                Iopt(13) Dense output option, not modified.
c                         = 0 : no dense output
c                         = 1 : dense  output  on  equidistant  points
c                               within the whole integration interval.
c                               The number of such points is specified
c                               by Iopt(14).
c                         = 2 : dense  output  on  equidistant  points
c                               within  every  integration  step.  The
c                               number of such points  is specified by
c                               Iopt(14).
c                         = 3 : output  of the solution  at the end of
c                               every integration  interval and  dense
c                               output on some additional points, such
c                               the distance  of two  output points is
c                               always  less or  equal t_Max. tMax  is
c                               specified by Ropt(2).
c
c                         If  one sets  Iopt(13) > 0 and Iopt(12) = 2,
c                         LIMEX returns to  the calling program at all
c                         dense output points.
c
c                         The solution at such points is  computed via
c                         Hermite interpolation. The size of the error
c                         of the interpolated solution may be slightly
c                         larger than  specified by  the rTol and aTol
c                         values.
c
c                         Restriction: 0 <= Iopt(13) <= 3.
c
c                Iopt(14) The number of  equidistant output points, if
c                         the dense output  option  Iopt(13) = 1 or 2,
c                         not  modified. If  Iopt(13) = 1, Iopt(14) is
c                         the number of output points within the whole
c                         integration  interval. If Iopt(14) = 0 or 1,
c                         no dense output will be provided.
c
c                         Example: if  t_Begin = 0 and  t_End = 1 with
c                         Iopt(13) = 1 and  Iopt(14) = 6  the solution
c                         at  the  points 0.0, 0.2, 0.4, 0.6, 0.8, 1.0
c                         will be provided.
c
c                         If Iopt(13) = 2, Iopt(14) is  the number  of
c                         output points within every integration step.
c                         If  Iopt(14) = 0, no  dense  output will  be
c                         provided.  If  Iopt(14) = 1 or 2, the  dense
c                         output consists of the solutions at the ends
c                         of the integration intervals.
c
c                         Example: if the current integration interval
c                         is [0,1], with Iopt(13) = 2 and Iopt(14) = 6
c                         the solution at 0.0, 0.2, 0.4, 0.6, 0.8, 1.0
c                         will be computed in this interval.
c
c                         Restriction: 0 <= Iopt(14)  if  Iopt(13) = 1
c                                      or 2.
c
c                Iopt(15) The unit number  for dense solution  output,
c                         not  modified.  Iopt(15) = 0  supresses  the
c                         dense output, which is then available in the
c                         calling program by setting Iopt(12) = 2.
c
c                         Restriction: 0 <= Iopt(15) if Iopt(13) > 0.
c
c                Iopt(16) Type  of call, may be modified. Iopt(16) = 0
c                         specifies  the  current  call  as an initial
c                         call, i.e. the  first  call  of LIMEX  for a
c                         given  problem. Thus  LIMEX performs  then a
c                         full data check and initialization. Further-
c                         more the solution  output  and, if required,
c                         the dense  output will  be printed  also  at
c                         t_Begin. Iopt(16) = 1 specifies a successive
c                         call, i.e. the integration of the last LIMEX
c                         call will  be continued. Then  the  solution
c                         and/or the dense  output will not be printed
c                         at t_Begin. To assure a correct continuation
c                         of the  further integration, no parameter of
c                         LIMEX except t_End should be changed between
c                         the LIMEX calls. If the code returns with no
c                         error, Iopt(16) on return will be always set
c                         to 1, thus enabling the continuation  of the
c                         integration.
c
c                         Restriction: 0 <= Iopt(16) <= 1.
c
c                Iopt(17) This option determines the behavior of LIMEX
c                         at  t_End, not  modified. With  Iopt(17) = 0
c                         the code integrates  exactly up to t_End. If
c                         Iopt(17) = 1, then  LIMEX may  have  already
c                         computed internally solutions for t > t_End.
c                         This option enables to  get most efficiently
c                         solutions  on an arbitrary user defined  set
c                         of values of the independent  variable t. If
c                         f(t,y)  or B(t,y) is  defined only  until an
c                         upper  limit, Ropt(3) may  be used to assure
c                         that t is always less or equal this limit.
c
c                         Restriction: 0 <= Iopt(17) <= 1.
c
c                Iopt(18) By  means of  this option one may generate a
c                         PostScript plot of the  Jacobian  matrix. If
c                         Iopt(18) = i > 0, in the  integration step i
c                         the  file 'Jacobian.ps'  will  be  generated
c                         representing the  non-zero structure  of the
c                         current  Jacobian. If in  this  step no  new
c                         Jacobian is computed, since Iopt(10) = 1, no
c                         PostScript plot will be  generated. The sub-
c                         routine JacPlot internally  uses the logical
c                         unit-number 99. If Iopt(18) = 0 no plot will
c                         be generated, if Iopt(18) = -1 the structure
c                         of the initial Jacobian (step 0) is plotted.
c
c                         Restriction: -1 <= Iopt(18).
c
c                Iopt(19) - Iopt(23)
c                         Parameters  used  only in  the LIMEX version
c                         4.3B (sparse Jacobians).
c
c                Iopt(24) On return the number of function evaluations
c                         during the integration without those for the
c                         computation of the Jacobian.
c
c                Iopt(25) On return the number of function evaluations
c                         needed for the  numerical computation of the
c                         Jacobian.
c
c                Iopt(26) On return the number of LU decompositions by
c                         the LAPACK routines dgetrf or dgbtrf.
c
c                Iopt(27) On return  the number of  back-substitutions
c                         by the LAPACK routines dgetrs or dgbtrs.
c
c                Iopt(28) On return the number of integration steps.
c
c                Iopt(29) On  return  the  number  of Jacobian  matrix
c                         evaluations.
c
c                Iopt(30) Parameter  used  only in  the  LIMEX version
c                         4.3B (sparse Jacobians).
c
c     Ropt       Real array of size 5, Ropt(1) to Ropt(3) must  be set
c                by caller on entry. Integration control parameter.
c
c                Ropt(1)  The maximum allowed  stepsize, not modified.
c                         If Ropt(1) = 0, the maximum allowed stepsize
c                         is set to the default value t_End - t_Begin.
c
c                         Restriction: 0 <= Ropt(1).
c
c                Ropt(2)  The  maximal  distance of  two dense  output
c                         points, not  modified. Used only if Iopt(13)
c                         is equal 3.
c
c                         Restriction: 0 < Ropt(2) if Iopt(13) = 3.
c
c                Ropt(3)  The upper limit for the independent variable
c                         t, not modified. Used only  if Iopt(17) = 1.
c                         If Ropt(3)  is less  than t_End the value of
c                         Ropt(3) will be ignored  and no upper  limit
c                         will be considered.
c
c                Ropt(4), Ropt(5)
c                         Parameters used  only in  the  LIMEX version
c                         4.3B (sparse Jacobians).
c
c     IPos       Integer array of size n, values must be set by caller
c                on entry, not modified. If IPos(i) = 1, the values of
c                the corresponding  solution component will be checked
c                to be  nonnegative. If  during  the  extrapolation  a
c                negative value occurs, the stepsize is reduced by the
c                heuristic factor RedPos, which  is set in a parameter
c                statement to 0.5, but  may  be changed. If IPos(i) is
c                not equal 1, no check will be performed.
c
c     IFail      Integer array of size 3, need not to be set by caller
c                on entry. Error indicator field.
c
c                IFail(1) Internal LIMEX error code. Zero if  no error
c                         occurred. If  Iopt(1) < 0, for each error an
c                         explanatory message is output on the current
c                         monitor  unit. The error  codes -1 until -31
c                         indicate  wrong or  inconsistent input data,
c                         the other error codes indicate errors during
c                         the integration.
c
c                Error-code IFail(1)            Description
c
c                         -1            Max_Nr_of_Equations  <  n   or
c                                       n < 1
c
c                         -2            rTol(i) <= 0 for some index i
c
c                         -3            aTol(i) < 0 for some index i
c
c                         -4            Iopt(1) < 0 or Iopt(1) > 2
c
c                         -5            Iopt(2) < 0 and Iopt(1) > 0
c
c                         -6            Iopt(3) < 0 or Iopt(3) > 2
c
c                         -7            Iopt(4) < 0 and Iopt(3) > 0
c
c                         -8            Iopt(5) < 0 or Iopt(5) > 1
c
c                         -9            Iopt(6) < 0 or Iopt(6) > 1
c
c                        -10            Iopt(7) < 0 or Iopt(7) > 1
c
c                        -11            Max_Lower_Diagonals <  Iopt(8)
c                                       or Iopt(8) < 0
c
c                        -12            Max_Upper_Diagonals <  Iopt(9)
c                                       or Iopt(9) < 0
c
c                        -13            Iopt(10) < 0 or Iopt(10) > 1
c
c                        -14            Iopt(11) < 0 or Iopt(11) > 1
c
c                        -15            Iopt(12) < 0 or Iopt(12) > 1
c
c                        -16            Iopt(13) < 0 or Iopt(13) > 3
c
c                        -17            Iopt(14) < 0 and  Iopt(13) = 1
c                                       or 2
c
c                        -18            Iopt(15) < 0 and Iopt(13) > 0
c
c                        -19            Iopt(16) < 0  or  Iopt(16) > 1
c
c                        -20            Iopt(17) < 0  or  Iopt(17) > 1
c
c                        -21            Iopt(18) < - 1
c
c                        -27            Ropt(1) < 0
c
c                        -28            Ropt(2) <= 0 and Iopt(13) = 3
c
c                        -32            An  initial value of y is less
c                                       than zero, but  by the setting
c                                       of IPos it should be => 0.
c
c                        -33            The routine Fcn  returns  with
c                                       an   error ( FcnInfo < 0 )  in
c                                       the first call  of Fcn in  the
c                                       current integration interval.
c
c                        -34            More than Max_Non_Zeros_B non-
c                                       zero entries in  matrix B(t,y)
c                                       defined.
c
c                        -35            The routine Fcn  returns  with
c                                       an  error (FcnInfo < 0) in the
c                                       numerical  evaluation  of  the
c                                       Jacobian.
c
c                        -36            The  routine Jacobian  returns
c                                       with an error (JacInfo < 0).
c
c                        -39            Internal error  in the  LAPACK
c                                       routines dgetrf or dgbtrf, see
c                                       IFail(2) and the LAPACK User's
c                                       Guide.
c
c                        -43            Problem not solvable with this
c                                       LIMEX version,  most  probable
c                                       reason: index of the DAE > 1.
c
c                        -44            Problem not solvable with this
c                                       LIMEX version,  most  probable
c                                       reason:  initial  values   not
c                                       consistent  or  index  of  the
c                                       DAE > 1.
c
c                        -45            A value of the solution vector
c                                       y within  the CIV  computation
c                                       is negative, but  should be by
c                                       => 0 by the settings of IPos.
c
c                        -46            More  stepsize reductions than
c                                       allowed  due  to failing  con-
c                                       vergence in the extrapolation.
c                                       This  is  controlled  by   the
c                                       parameter Max_Step_Red_Ex.
c
c                        -47            Singular matrix pencil B - hA,
c                                       problem   not  solvable   with
c                                       LIMEX.
c
c                        -48            More  integration  steps  than
c                                       allowed. This is controlled by
c                                       the parameter Max_Int_Steps.
c
c                        -49            More internal Newton steps for
c                                       the computation  of consistent
c                                       initial  values than  allowed.
c                                       This  is  controlled  by   the
c                                       parameter Max_Step_CIV.
c
c                        -50            Stepsize  too  small, in  most
c                                       cases  due  to  many  stepsize
c                                       reductions.
c
c                IFail(2) Internal error code of the subroutine dgetrf
c                         or  dgbtrf (IFail(1) = - 39). See the LAPACK
c                         User's Guide for  an explanation. Zero if no
c                         error occurred.
c
c                IFail(3) This  error code  signals whether  the error
c                         occurred  during the  computation of CIVs or
c                         not. IFail(3) = 0: error during integration,
c                         IFail(3) = 1: error  during CIV computation.
c                         For  both IFail(3) = 0 or 1 the error  codes
c                         IFail(1)  and  IFail(2) specify the  type of
c                         the error.
c
c     Remark: In  addition to  these errors  LIMEX may  output on  the
c             current monitor unit the following warning messages:
c
c             1. Warning: component xxx does  not have  an  asymptotic
c                         expansion in the first integration step
c
c                This messsage  occurrs if for  some components of the
c                solution  no asymptotic expansion in  powers of h can
c                be found. The reason  for this can be either an index
c                of the DAE which is greater than 1 or (in most cases)
c                the lack of consistent initial  values. Thus, if this
c                message occurs  very frequently, the  computation  of
c                consistent initial  values by setting Iopt(6) = 1 may
c                prevent the  warning. However, in  many  cases  where
c                this  warning is  printed, the  integration  proceeds
c                without problems.
c
c             2. Warning: switch  to direct  solution of  higher order
c                         approximations
c
c                This warning  occurrs, if during  the integration the
c                higher order approximations can not be computed using
c                a Richardson iteration, if the direct solution of the
c                linear  systems is  selected. If an  iterative method
c                is  selected, this message  indicates that the higher
c                order  solutions can  not be obtained with  the first
c                order matrix ILU factorization as  preconditioner. In
c                both  cases, an appropriate  direct approach  will be
c                tried, but this can  be more time  consuming. If  the
c                problem is well posed, then this warning  indicates a
c                strongly varying matrix B(t,y). If the linear systems
c                are solved with direct methods, the maximum number of
c                iterations in the  Richardson  process  is defined by
c                the variable Max_Iter_Rich. The value of the limit is
c                defined in a parameter  statement and may  be changed
c                if necessary.
c
c-----------------------------------------------------------------------
c
c     Subroutines and functions called
c
c-----------------------------------------------------------------------
c
c     Subroutines to be supplied by the user as externals:
c
c     Fcn ( n, nz, t, y, f, b, ir, ic, FcnInfo )
c     ------------------------------------------
c
c     This subroutine  computes the right-hand side of the DAE and the
c     non-zero entries in the left-hand side matrix B(t,y).
c
c     Parameters:
c
c     n          Integer variable, size  of the DAE system, should not
c                be changed.
c
c     nz         Integer  variable, must  be  set  in the  subroutine.
c                Number of non-zero entries in the matrix B(t,y).
c
c     t          Real  variable, should not be  changed. Current value
c                of the independent variable.
c
c     y          Real array of size n, should not  be changed. Current
c                values of  the solution  vector, for  which f(t,y) is
c                required.
c
c     f          Real array of  size n, must be set in the subroutine.
c                Values of f(t,y) at the current values of t and y.
c
c     b          Real array of size nz, must be set in the subroutine.
c                Values of the non-zero entries in  the matrix B(t,y).
c
c     ir         Integer  array of  size nz, must  be set  in the sub-
c                routine. Row  indices of  the non-zero entries in the
c                matrix B(t,y).
c
c     ic         Integer  array of  size nz, must  be set  in the sub-
c                routine. Column  indices  of the  non-zero entries in
c                the matrix B(t,y).
c
c     FcnInfo    Integer variable. FcnInfo indicates  the type  of the
c                call  of Fcn. If FcnInfo = 0, the  call is  the first
c                call of  Fcn  in  the  current  integration  step. If
c                FcnInfo = 1 this  call is  within  the  extrapolation
c                process. FcnInfo = 2 indicates  a call of  Fcn within
c                the  numerical  evaluation of  the Jacobian. In  some
c                application  it may be  reasonable for efficiency, to
c                compute f and/or B during the numerical evaluation of
c                the Jacobian only with  limited accuracy. Furthermore
c                FcnInfo on return is  used as error  indicator. If an
c                error in the subroutine  Fcn occurs, set FcnInfo < 0.
c                Then LIMEX reduces the stepsize by the factor RedFcn.
c                This  variable   is  currently  set  in  a  parameter
c                statement  to 0.5, but  may  be  changed. If  such an
c                error  occurs during  the numerical evaluation of the
c                Jacobian, the program exits with an error message. If
c                no error occurs, set FcnInfo = 0.
c
c     Remarks:   The order  of the  entries in  matrix B is arbitrary.
c                If the indices  ir and ic are constant, they  must be
c                defined only once in the first call of Fcn.
c
c
c     Jacobian ( n, t, y, ys, Jac, LDJac, ml, mu, Full_or_Band,
c     ---------------------------------------------------------
c                JacInfo )
c                ---------
c
c     This  subroutine computes  the Jacobian  matrix, i.e. analytical
c     derivatives of the residual f(t,y) - B(t,y) * y'.
c
c     Parameters:
c
c     n          Integer variable, size  of the DAE system, should not
c                be changed.
c
c     t          Real  variable, should not be  changed. Current value
c                of the independent variable.
c
c     y          Real array of size n, should not  be changed. Current
c                values of  the solution  vector, for  which f(t,y) is
c                required.
c
c     ys         Real array of size n, should not  be changed. Current
c                values of the derivatives y'.
c
c     Jac        Real matrix of  dimension at least (LDJac,n), must be
c                set in the subroutine. The derivatives at the current
c                values of y, ys and t.
c
c     LDJac      Integer variable, the leading dimension of the matrix
c                Jac as declared in LIMEX, should not be changed.
c
c     ml         Integer  variable, the  number of  lower diagonals in
c                the Jacobian, should not be changed.
c
c     mu         Integer  variable, the  number of  upper diagonals in
c                the Jacobian, should not be changed.
c
c     Full_or_Band
c                Integer variable, should not be changed.
c                = 0 : Jacobian matrix is a full matrix.
c                = 1 : Jacobian matrix is a band matrix.
c
c     JacInfo    Integer variable, error indicator. If an error in the
c                subroutine  Jacobian  occurs,  set  JacInfo < 0. Then
c                LIMEX aborts with an error  message. If JacInfo => 0,
c                no error is assumed.
c
c     Remarks:   The general form of  this subroutine should look like
c                the  following  code  segment. This  works  for  both
c                matrix representations. The storage scheme for banded
c                matrices is the usual LINPACK/LAPACK format.
c
c                do j = 1, n
c
c                   iUpper = max ( 1, j - mu )
c                   iLower = min ( n, j + ml )
c
c                   if ( Full_or_Band .eq. 0 ) then
c                      k = 0
c                   else
c                      k = mu + 1 - j
c                   end if
c
c                   do i = iUpper, iLower
c
c                      Now define matrix element with indices i, j:
c
c                      Jac(i+k,j) = ...
c
c                   end do
c
c                end do
c
c                If the  numerical  approximation  of the  Jacobian is
c                used, one may supply a dummy routine. If the Jacobian
c                is a  sparse  matrix, it  can be  more  efficient, to
c                compute  the numerical approximation  of the Jacobian
c                also  by this  routine, even  if  for  such  problems
c                LIMEX4.3B should be preferable.
c
c     Other subroutines
c
c     From LIMEX_Dense:
c     -----------------
c
c     Comp_Herm  Computes  right  derivatives  approximations and  the
c                coefficients of the Hermite interpolation.
c
c     Eval_Herm  Evaluates the Hermite interpolation polynomial.
c
c     JacPlot    Generates a PostScript plot of the Jacobian matrix.
c
c     From the BLAS library:
c     ----------------------
c
c     daxpy      Adds a scalar multiple of a  vector to another vector
c                (used only in LIMEX4_3A_Dense.f).
c
c     dcopy      Copies real vectors
c
c     dnrm2      Computes the Euclidean norm of a real vector
c
c     dscal      Scales a real vector by a constant
c
c     From the LAPACK Library:
c     ------------------------
c
c     dgbtrf     Factorizes a banded real matrix
c
c     dgbtrs     Solves a banded real linear system
c
c     dgetrf     Factorizes a general real matrix
c
c     dgetrs     Solves a general real linear system
c
*/
#define MAX_NO_EQNS 4096
#define MAX_ROW_TAB 7


typedef void (*Fcn)(int* n, int* nz,
                    double* t, double* y, double* dy,
                    double* B, int* ir, int* ic,
                    int* info);

typedef void (*Jac)(int* n,
                    double* t, double* y, double* dy,
                    double* J, int* ldJ, int* ml, int* mu,
                    int* full_or_band, int* info);

// original LIMEX API
/*
extern void limex_(
            int*      n,
            Fcn       fcn,
            Jac       jac,
            double*   t0,
            double*   T,
            double*   y0,
            double*   dy0,
            double*   rTol,
            double*   aTol,
            double*   h,
            int*      iOpt,
            double*   rOpt,
            int*      iPos,
            int*      ifail);
*/

// extended LIMD API: iOpt(31):
//                      0 = B (n x n) identity,
//                      1 = B const.,
//                      2 = B may be variable (as in original API)
extern void limd_(
            int*      n,
            Fcn       fcn,
            Jac       jac,
            double*   t0,
            double*   T,
            double*   y0,
            double*   dy0,
            double*   rTol,
            double*   aTol,
            double*   h,
            int*      iOpt,
            double*   rOpt,
            int*      iPos,
            int*      ifail);


// extenteded API (as above) with extrapolation entry if iOpt(32)==1
extern void limdherm_(
            int*      n,
            Fcn       fcn,
            Jac       jac,
            double*   t0,
            double*   T,
            double*   y0,
            double*   dy0,
            double*   rTol,
            double*   aTol,
            double*   h,
            int*      iOpt,
            double*   rOpt,
            int*      iPos,
            int*      ifail,
            int*      kOrder,
            double*   Dense,        // returns the comp_herm_() output array 'Dense' in this special version of LIMEX
            double*   t1,
            double*   t2);


extern void eval_herm_(
            int*      n,
            int*      maxEqns,
            double*   Dense,        // has to be the output array 'Dense' of comp_herm_() !
            int*      kOrder,
            double*   tFac,         // tFac = (tEval - t1) / (t2 - t1)
            double*   yInterp);     // y(.) at interpolation point ( == y(tEval) )

//
/*
extern void comp_herm_(
            int*      n,
            int*      maxEqns,
            double*   Dense,
            int*      kOrder,
            int*      ipt,
            int*      nj,
            double*   Work);
*/

///

// corresponding Hermite interpolation evaluating according to 'Dense' tableau
//
// NOTE: Currently not re-entrant since the 'Dense' array is modified by hermine_().
//       Therefore, it is commented out.
/*
extern void hermine_(
            int*      n,
            int*      kOrder,
            double*   Dense,
            double*   t1,
            double*   t2,
            double*   tEval,
            double*   yEval);
*/
