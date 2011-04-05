// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-03-17 td
// last changed:
//
#ifndef __YE_OLDE_PARKIN_CORE_H
#define __YE_OLDE_PARKIN_CORE_H

// #include "addpkg/dlib/logger.h"
#include "common/Constants.h"
#include "common/PARKINLog.h"

#include "linalg/QRDecomp.h"

#include "UserFunc.h"


namespace PARKIN
{
    struct YeOldeParkinWk
    {
        Real        eps, fc, zscal, xstep, cond;
        Real        fcmin;
        unsigned    iter, itmax;

        YeOldeParkinWk() :
            eps(0.0), fc(0.0), zscal(-1.0), xstep(0.0), cond(0.0),
            fcmin(1.0e-2),
            iter(0), itmax(50)
        { }
    };

    ///
    ///

    class YeOldeParkinCore
    {
        public:
            // c'tor
            YeOldeParkinCore();

            // d'tor
            ~YeOldeParkinCore();

            // initial settings and problem definition
            void setIOpt(IOpt const& iopt);
            void setWk(YeOldeParkinWk const& wk);
            void setProblem(UserFunc* fun);
            int initialise( unsigned m,
                            Vector const& x, Vector const& xscal,
                            Vector const& fobs, Vector const& fscal,
                            Real const rtol,
                            IOpt const& iopt, YeOldeParkinWk const& wk
                          );

            // start iteration
            int run();

            // statistical analyse
            int analyse();

            // get (current iteration) information
            Vector          getSolution();
            YeOldeParkinWk  getWk();
            void            printCounter();

        private:
            // copy c'tor
            YeOldeParkinCore(YeOldeParkinCore const&);
            // assignment
            YeOldeParkinCore const& operator= (YeOldeParkinCore const&);

            // helper routines
            Vector call_FCN(Vector const&, int&);
            Matrix call_JAC(Vector const&, int&);
            int initialise_timer_and_counter();
            int initiation();
            int set_internal_parameters();
            int scale_model_and_measurements();
            int prepare_GN_method();
            int initial_scaling();
            int compute_residual_vector();
            int approximate_negative_jacobian();
            int solve_linear_system();
            int control_steplength();
            int control_steplength_output();
            int predict_damping_factor();
            int optimise_rank();
            int optimise_rank_sub621();
            int compute_next_trial_step();
            int compute_next_sub1();
            int compute_next_sub2(unsigned);
            int reduce_damping_factor();
            int prepare_for_next_iteration_step();
            int exit_solution_entry1();
            int exit_solution_entry2();
            int exit_sub82();
            int exit_sub88();
            int print_final_residual_and_counts();
            int compute_final_simulation_and_sensitivity();

            //
            //          _mprmon =   0      1      2      3      4       5       6
            //  dlib::log_level =  LNONE  LINFO  LVERB  LTALK  LGABBY  LDEBUG  LTRACE
            // logical units for logging; log levels by _mprerr, _mprmon, _mprsol, _mptim
            dlib::logger    _luerr, _lumon, _lusol, _lutim;
            static
            dlib::log_level _loglvl[];

            // data spaces
            IOpt            _iopt;
            YeOldeParkinWk  _wk;

            UserFunc*       _fun;

            //
            unsigned    _n, _m;
            Vector      _x, _xa, _xw, _pivot, _v;
            Vector      _dx, _dxh, _dx1, _dx1a, _diag;
            Vector      _z, _zscal;
            Vector      _f, _fd, _fh, _u;
            Matrix      _A, _AH;
            Real        _fc, _fca, _fc1, _fcmin, _fcmin2;
            Real        _cond, _d1, _sk, _tmin;
            Real        _sumf, _sumfa, _sumx, _sumxa, _sfc1, _sfc2;
            unsigned    _div, _konv, _level;
            unsigned    _irank, _iranka, _irkdec, _irkopt;
            unsigned    _iter, _ic, _itol, _icall, _iscal;
            unsigned    _irkmax, _itmax, _icmax, _divm;
            unsigned    _ifgn, _ijgn;
            Real        _epmach, _epsu, _zscald, _xstepm, _qkap, _skap;
            Real        _tol, _tolf, _tolj, _tolmin, _tolmax, _smax;
            unsigned    _ndecom, _nsolve;

            bool        _lpos, _lscal, _kt;

            QRDecomp    _qrA;

    };

}

/*
      SUBROUTINE PARK11 (NDE,NCEQ,N,NIPAR,M,LDIM,ITPM1,NSW,NRW,NIW,
     @                  TP,YINIT,X,IPARA,Z,IZ,IPZ,PIVOT,A,AH,D,DX,
     @                  DXH,DX1,DX1A,V,XA,XW,HH,F,FD,FH,U,SW,RW,IW,
     @                  KFLAG,IITER,RITER)
C
C
C----------------------------------------------------------------
C  SUBROUTINE PARK11 PERFORMS GAUSS-NEWTON METHOD FOR A DISCRETE
C  QUADRATIC ERROR FUNCTIONAL
C
C
C  REFERENCES:
C=============
C
C  /1/  U. NOWAK, P. DEUFLHARD:
C       TOWARDS PARAMETER IDENTIFICATION FOR LARGE CHEMICAL
C       REACTION SYSTEMS
C       (UNIV. HEIDELBERG, SFB 123: TECHN. REP. ... (1982))
C       IN:
C       P. DEUFLHARD, E. HAIRER (ED.):
C       NUMERICAL TREATMENT OF INVERSE PROBLEMS IN DIFFERENTIAL AND
C       INTEGRAL EQUATIONS
C       BIRKHAEUSER: PROGRESS IN SCIENTIFIC COMPUTING ... (1983)
C
C  /2/  P. DEUFLHARD, V. APOSTOLESCU:
C       A STUDY OF THE GAUSS-NEWTON METHOD FOR THE SOLUTION
C       OF NONLINEAR LEAST SQUARES PROBLEMS
C       IN: SPECIAL TOPICS OF APPLIED MATHEMATICS
C           (ED. BY J.FREHSE, D.PALLASCHKE, U.TROTTENBERG)
C            P. 129 - 150 (1980)
C
C----------------------------------------------------------------
C
C
C DATE OF LATEST CHANGE: FEB. 22, '83
C
C
C**********************************************************************
C                                                                     *
C                                                                     *
C  INPUT PARAMETERS:                                                  *
C-------------------                                                  *
C                                                                     *
C   NDE        NUMBER OF ODE'S                                        *
C   NCEQ       NUMBER OF CHEMICAL EQUATIONS                           *
C              (NUMBER OF KINETIC PARAMETERS)                         *
C   N          NUMBER OF PARAMETERS TO BE IDENTIFIED                  *
C   NQ1H       = N*(N+1)/2                                            *
C   NIPAR      DIMENSION OF INTEGER ARRAY IPARA                       *
C              NIPAR=NQ INDICATES,THAT NO LINEARLY DEPENDENT PARAM.   *
C              ARE DEFINED:  NLD:=0; ELSE: NLD:=NIPAR-3*NQ            *
C   M          NUMBER OF MEASUREMENTS                                 *
C   LDIM       NUMBER OF NONZEROS IN JACOBIAN OF ODE                  *
C   ITPM       NUMBER OF DIFFERENT MEASUREMENT POINTS PLUS 2 FOR      *
C              TSTART AND TFINAL                                      *
C   NRW        DIMENSION OF REAL WORK SPACE FOR TRAJECTORIES          *
C   NSW        DIMENSION OF REAL WORK SPACE FOR SENSITIVITIES         *
C   NIW        DIM. OF INTEGER WORK SPACE                             *
C   TP(ITPM)   HALTING POINTS FOR INTEGRATOR (TSTART,TMEASURE1,       *
C              TMEASURE2,...,TFINAL)                                  *
C   YINIT(NDE) INITIAL VALUES FOR ODE                                 *
C   X(N)       INITIAL ESTIMATE OF PARAMETERS TO BE IDENTIFIED        *
C   IPARA(NIPAR) INDICIES OF P. TO BE IDENTIFIED (INDEPENDENT AND     *
C                LINEARLY DEPENDENT) AND POINTERS                     *
C                ( X(I) = RK(IPARA(I)) , I = 1,2,...,NQ )             *
C   Z(M)       VALUES OF OBSERVED EXPERIMENTAL DATA (MEASUREMENTS)    *
C   IZ(M)      SPECIES INDICES OF MEASUREMENTS                        *
C   IPZ(ITPM)  POINTERS TO Z AND IZ (SEE STORAGE DESCRIPTION BELOW)   *
C   ATA(NQ1H)  REAL WORK ARRAY OF LENGTH N*(N+1)                      *
C   DX,DXH,DX1,DX1A,V,XA,XW,F,FH:                                     *
C              REAL WORK ARRAYS OF LENGTH N                           *
C   RW(NRW)    REAL WORK SPACE FOR TRAJECTORIES                       *
C   SW(NSW)    REAL WORK SPACE FOR SENSITIVITIES                      *
C   IW(NIW)    INTEGER WORK ARRAY                                     *
C                                                                     *
C   KFLAG      PRINT PARAMETER (FOR GAUSS NEWTON ITERATION)           *
C              .EQ.-1: NO PRINT                                       *
C                   0: PRINT SOLUTION ONLY                            *
C                  +1: PRINT OF PROBLEM DESCRIPTION,                  *
C                      INITIAL GUESS OF PARAMETERS,                   *
C                      ITERATIVE VALUES OF LEVEL FUNCTIONS,           *
C                      SOLUTION, VARIANCE OF RESIDUALS AND SOLUTION   *
C                      CONFIDENCE INTERVALS                           *
C                  +2: ADDITIONALLY ALL ITERATES X, VARIANCE-COVAR.-  *
C                      MATRIX AND CORRELATION COEFFICIENTS            *
C                  +3: INITIAL AND FINAL RESIDUALS, SCALING OF        *
C                      MEASUREMENTS ADDITIONALLY                      *
C                  +4: ALL RESIDUALS ADDITIONALLY                     *
C                                                                     *
C                                                                     *
C   IITER      INTEGER ARRAY OF LENGTH 10 FOR ADDITIONAL INFORMATION  *
C                                                                     *
C   IITER(1):  PRINT PARAMETER FOR INTEGRATOR                         *
C              .EQ.0: ERROR MESSAGES ONLY                             *
C              .EQ.1: INTEGRATION MONITOR                             *
C              .EQ.2: ENHANCED INTEGRATION MONITOR                    *
C                                                                     *
C   IITER(2):  SETS PRINT PARAMETER FOR INTEGRATOR ACCORDING TO       *
C              ACTUAL ITERATION NUMBER                                *
C              SET TO IITER(1) IF ITER.LE.IITER(2)                    *
C              SET TO  0  IF ITER.GT.IITER(2)                         *
C                                                                     *
C   IITER(3):  INDICATES TYP OF PARKIN CALL                           *
C              .EQ.0: DO PARAMETER IDENTIFICATION                     *
C              .EQ.1: DO SIMULATION ONLY; COMPUTE SMOOTH TRAJECTORY   *
C                     (WRITE SOLUTION TO DATABASE 'OUT')              *
C              .EQ.3: DO SIMULATION AND SENSITIVITY ANALYSIS (FOR     *
C                     INDICATED PARAMETERS) FOR SMOOTH TRAJECTORY     *
C                     (WRITE SOLUTION TO DATABASE 'OUT' AND (RELATIVE)*
C                      SENSITIVITIES TO DATABASE 'HELP')              *
C                                                                     *
C   IITER(4):  INDICATES TYP OF SIMULATION FOR FINAL PARAMETER        *
C              ITERATES                                               *
C              .EQ.0: NO SIMULATION (-> NO SOLUTION TRAJECTORY ON     *
C                                       DATABASE 'OUT')               *
C              .EQ.1,3: DO SIMULATION WITH FINAL PARAMETER            *
C                       ITERATES AS DESCRIBED FOR IITER(3)            *
C                                                                     *
C   IITER(5):  INDICATES TYP OF DATA OUTPUT FOR PLOT                  *
C              .EQ.0: DATA OUTPUT ACCORDING TO IITER(4)               *
C              .EQ.1: ADDITIONALLY INITIAL TRAJECTORY ON DATABASE OUT *
C              .EQ.2: WRITE SOLUTION OF SMOOTH TRAJECTORY FOR ALL     *
C                     PARAMETER ITERATES TO DATABASE 'HELP'           *
C                     (-> IITER(4).LE.2)                              *
C                                                                     *
C   IITER(6):  MAXIMUM NUMBER OF ITERATIONS (-> ITMAX)                *
C                                                                     *
C   IITER(7):  INDICATES USE OF POSITIVITY OF PARAMETERS              *
C              .EQ.0: ITERATION IN LOGARITHM OF PARAMETERS            *
C              .EQ.1: NO USE OF POSITIVITY CONSTRAINT                 *
C                                                                     *
C   IITER(8):  INDICATES TYPE OF SCALING FOR FUNCTIONAL               *
C              .EQ.0: INTERNAL SCALING WITH DIAGONAL WEIGHTING MATRIX *
C              .EQ.1: USE USER SUPPLIES SUBROUTINE SCALUS             *
C                                                                     *
C                                                                     *
C   RITER      REAL ARRAY OF LENGTH 10 FOR ADDITIONAL INFORMATION     *
C                                                                     *
C   RITER(1): =EPSU :  PRESCIBED RELATIVE PRECISION FOR SOLUTION      *
C                      .EQ.0: APPROPRIATE EPSU INTERNALLY SELECTED    *
C                                                                     *
C   RITER(2): =FC1 :   STARTING VALUE FOR RELAXIATON FACTOR           *
C                      .EQ.0: APPROPRIATE FC1 INTERNALLY SELECTED     *
C                                                                     *
C   RITER(3): =ZSCALD: THRESHOLD FOR (INTERNAL) SCALING OF MEASURE-   *
C                      MENTS (REL.ERROR FOR ALL MEASUREMENTS Z WITH:  *
C                      Z.GT.ZMAX*ZSCALD                               *
C                      .EQ.0: APPROPRIATE ZSCALD INTERNALLY SELECTED  *
C                                                                     *
C   RITER(4): =HMAX:   MAXIMUM PERMITTED STEPSIZE FOR INTEGRATOR FOR  *
C                      FINAL SIMULATION (IITER(3) OR IITER(4)=1/2/3/4)*
C                      .EQ.0: NO STEPSIZE RESTRICTION                 *
C                                                                     *
C                                                                     *
C                                                                     *
C  OUTPUT PARAMETERS                                                  *
C-------------------                                                  *
C                                                                     *
C    KFLAG     ERROR FLAG                                             *
C              .GE.-1: SUCCESSFULL IDENTIFICATION (KFLAG NOT ALTERED) *
C              .EQ.-2: IDENTIFICATION FAILED                          *
C                                                                     *
C                  -1  TERMINATION, SINCE ITERATION DIVERGES          *
C                      INITIAL GUESS TOO BAD OR MODEL TOO FAR         *
C                      FROM BEING COMPATIBLE                          *
C                  -2  TERMINATION AFTER ITMAX ITERATIONS             *
C                      ( AS INDICATED BY INPUT PARAMETER ITMAX )      *
C                  -3  TERMINATION, SINCE RELAXATION                  *
C                      STRATEGY DID NOT SUCCEED                       *
C                   IN CASE OF FAILURE:                               *
C                      - USE BETTER INITIAL GUESS                     *
C                      - OR REFINE MODEL                              *
C                      - OR TURN TO GENERAL OPTIMIZATION ROUTINE      *
C                                                                     *
C                                                                     *
C   IITER(19):  NUMBER OF UNUSED INTEGER WORK SPACE LOCATIONS         *
C               (FOR IITER(19).LT.0 :                                 *
C                  ENLARGE WORK SPACE TO SPEED UP LINEAR ALGEBRA      *
C   IITER(20):  NUMBER OF UNUSED REAL WORK SPACE LOCATIONS            *
C               (FOR IITER(20).LT.0 :                                 *
C                  ENLARGE WORK SPACE TO SPEED UP LINEAR ALGEBRA      *
C                                                                     *
C                                                                     *
C**********************************************************************
C                                                                     *
C                                                                     *
C                                                                     *
C  FURTHER INTERNAL PARAMETERS                                        *
C-----------------------------                                        *
C                                                                     *
C  LPOS       USE POSITIVITY OF RATE CONSTANTS FOR PARAMETER ITERATION*
C                                                                     *
C  FC1        STARTING VALUES OF RELAXIATION FACTOR                   *
C  FC         ACTUAL VALUE OF RELAXATION FACTOR                       *
C  FCMIN      MINIMUM PERMITTED VALUE OF RELAXATION FACTOR            *
C  ITMAX      MAXIMUM PERMITTED NUMBER OF ITERATIONS                  *
C                                                                     *
C  COND       MAXIMUM PERMITTED SUB-CONDITION NUMBER                  *
C             OF JACOBIAN MATRIX                                      *
C  IRANK      STARTING VALUE OF PSEUDO-RANK OF JACOBIAN MATRIX        *
C  IRKMAX     MAXIMUM PERMITTED PSEUDO-RANK OF JACOBIAN MATRIX        *
C  ICMAX      MAXIMUM PERMITTED NUMBER OF RANK-REDUCTION-STEPS        *
C  IDIVM      MAXIMUM PERMITTED NUMBER OF SUCCESSIVE                  *
C             NON-CONTRACTIVE STEPS                                   *
C  TOLF       REQUIRED REL. ACCURACY FOR INTEGRATOR FOR EVALUATION    *
C             OF MODEL FUNCTION F (SOLVE IVP)                         *
C  TOLJ       REQUIRED REL. ACCURACY FOR INTEGRATOR FOR EVALUATION    *
C             OF VARIATIONAL EQUATION W.R.T. PARAMETERS X             *
C                                                                     *
C                                                                     *
C                                                                     *
C**********************************************************************
C                                                                     *
C                                                                     *
C                                                                     *
C  WORKSPACE PARKIN VERSION 1.1:                                      *
C===============================                                      *
C                                                                     *
C  WITH:                                                              *
C                                                                     *
C  N:    NUMBER OF PARAMETERS TO BE ESTIMATED                         *
C  NCEQ: NUMBER OF KINETIC PARAMETER (=NUMB. OF CHEMICAL EQUA.)       *
C  NDE:  NUMBER OF FIRST ORDER DIFFERENTIAL EQUATION (=N. OF SPECIES) *
C  LDIM: NUMBER OF NONZEROS IN JACOBIAN                               *
C  NLU:  NUMBER OF NONZEROS IN DECOMPOSED JACOBIAN                    *
C  JM:   MAX. NUMBER OF ROWS IN EXTRAPOLATION TABLE                   *
C        ( JM:=5 FOR THIS VERSION)                                    *
C                                                                     *
C                                                                     *
C (1) WORKSPACE FOR GN-METHOD                                         *
C----------------------------                                         *
C                                                                     *
C     M*N + N*N + 9*N + 5*M                                           *
C                                                                     *
C                                                                     *
C                                                                     *
C (2) WORKSPACE FOR INTEGRATION                                       *
C------------------------------                                       *
C                                                                     *
C                                                                     *
C  (2.1) FOR TRAJECTORY                                               *
C        RW(NRW): NRW = (JM + 7)*NDE + LDIM + NLU + NDE               *
C                       RULE OF THUMB: NRW = 34 * NDE                 *
C                                                                     *
C  (2.2) FOR SENSITIVITIES                                            *
C        SW(NSW): NSW = (JM + 5)*NCN                                  *
C                                                                     *
C                                                                     *
C  (2.3) FOR POINTER                                                  *
C        IW(NIW): NIW = 14*NDE + LDIM + NLU + NDE                     *
C                       RULE OF THUMB: NIW = 36 * NDE                 *
C                                                                     *
C**********************************************************************
C                                                                     *
C                                                                     *
C                                                                     *
C  STORAGE DESCRIPTION FOR MEASUREMENTS                               *
C--------------------------------------                               *
C                                                                     *
C  MEASUREMENT VALUES ARE STORED LINEARLY IN ARRAY Z(M). FIRST ALL    *
C  MEASUREMENTS OF THE FIRST MEASUREMENT POINT, THEN ALL OF THE       *
C  SECOND MEASUREMENT POINT ETC.. FINALLY ALL MEASUREMENTS OF THE     *
C  LAST MEASUREMENT POINT.                                            *
C  THE SPECIES INDEX OF A VALUE OF ARRAY Z IS STORED IN ARRAY IZ AT   *
C  THE SAME POSITION ( IZ(J)=K MEANS: THE MEASUREMENT Z(J) BELONGS    *
C  TO SPECIES NO. K )                                                 *
C  THE MEASUREMENT POINTS (T-VALUES) ARE STORED IN ARRAY TP, WHERE    *
C  THE FIRST VALUE OF TP IS THE TIME BELONGING TO THE INITIAL VALUES. *
C  SO, TP(2) CONTAINS THE FIRST MEASUREMENT POINT                     *
C  THE ARRAY IPZ CONTAINS POINTERS FROM THE T-VALUES TO THE FIRST     *
C  MEASUREMENT OF ARRAY Z (AND IZ) WHICH BELONGS TO THIS MEASUREMENT  *
C  POINT. SO:                                                         *
C                                                                     *
C  TP(J)  IS THE J-1 MEASUREMENT POINT.                               *
C  Z(IPZ(J))  IS THE FIRST MEASUREMENT OF M. POINT J-1 (=TP(J))       *
C  IZ(IPZ(J))  IS THE SPECIES INDEX OF Z(IPZ(J))                      *
C  IPZ(J+1)-IPZ(J)  IS THE NUMBER OF MEASUREMENTS AT TP(J)            *
C                                                                     *
C                                                                     *
C  STORAGE DESCRIPTION FOR PARAMETERS                                 *
C------------------------------------                                 *
C                                                                     *
C                                                                     *
C                                                                     *
C                                                                     *
C                                                                     *
C**********************************************************************
C
C
*/

#endif // __YE_OLDE_PARKIN_CORE_H
