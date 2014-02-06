// Copyright (C) 2010 - 2012
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2012-10-04 td
// last changed:
//

#include "MultipleShootingGN.h"
#include "odelib/LIMEX_A.h"

using namespace PARKIN;

///

dlib::log_level MultipleShootingGN::_loglvl[] = {
                  dlib::LNONE
                , dlib::LINFO
                , dlib::LVERB
                , dlib::LTALK
                , dlib::LGABBY
                , dlib::LDEBUG
                , dlib::LTRACE
                , dlib::LALL
};

///

//---------------------------------------------------------------------------
MultipleShootingGN::MultipleShootingGN() :
    _luerr("err.MS"), _lumon("mon.MS"), _lusol("sol.MS"), // _lutim("tim.MS"),
    _ode(0),
    _m(0), _n(0),
    _odeSolver( new LIMEX_A() )
{
    IOpt  iopt = IOpt();

    dlib::set_all_logging_levels(dlib::LNONE);
    // _lusol.set_level(dlib::LNONE);
    _luerr.set_logger_header(&print_parkin_logger_header);
    _lumon.set_logger_header(&print_parkin_logger_header);
    _lusol.set_logger_header(&print_parkin_logger_header);

    setIOpt(iopt);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
MultipleShootingGN::~MultipleShootingGN()
{
    delete _odeSolver;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
MultipleShootingGN::setLogStream(std::ostream& out)
{
    dlib::set_all_logging_output_streams(out);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
MultipleShootingGN::setIOpt(IOpt const& iopt)
{
    _luerr.set_level(_loglvl[iopt.mprerr]);
    _lumon.set_level(_loglvl[iopt.mprmon]);
    _lusol.set_level(_loglvl[iopt.mprsol]);

    _iopt = iopt;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
MultipleShootingGN::setProblem(FirstOrderODESystem* ode)
{
    _ode = ode;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Vector
MultipleShootingGN::getNodes()
{
    return _tnodes;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Matrix
MultipleShootingGN::getSolution()
{
    return _X;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Real
MultipleShootingGN::getPeriod()
{
    return _period;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Matrix
MultipleShootingGN::getFloquetMultipliers()
{
    // bool qinit = true;
    // bool qjcrfr = false;
    unsigned    jacgen, low, igh;
    int         iauto, ifail = 0;
    Vector      t1, t2, dE;


    jacgen = _iopt.jacgen;
    iauto = _iopt.iauto;


    // Comp. of the trajectories
    _HH = call_IVPSOL( _tnodes, _period, _X, _XU, ifail );

    ++_kount;

    if ( ifail != 0 )
    {
        // Singular trajectory
        _ierr = 82;
        _FM = Matrix();

        return _FM;
    }

    // Scale unkowns X
    compute_scaling_XW( _xthr );

    if ( iauto > 0 )
    {
        _periodw = _period;
    }

    // iranka = _irank;

    // Jacobian matrix generation
    // --------------------------
    if ( jacgen != 3 )
    {
        _new = 0;
        ifail = 0;

        if ( !(jacgen > 0) )
        {
            // Diff. approx. of Wronskian matrices G(1), ..., G(m1)
            // ----------------------------------------------------
            compute_fd_derivative_G( _X, _period, ifail );

            _kount += _n;
        }
        else
        {
            // Wronskian G(1), ..., G(m1) by integration of variational equation
            // -----------------------------------------------------------------
            solve_var_eq_for_G( _X, _period, ifail );

            _kount += _n;
        }

        // if ( ifail != 0 )
    }
    else
    {
        // Rank-1 updates of Wronskian matrices G(1), ..., G(m1)
        // -----------------------------------------------------
        ++_new;
        compute_rank1_update_G( iauto );
    }

    // Comp. of sensitivity matrix E
    // -----------------------------
    if ( _irank != 0 )
    {
        // Storing row scaling vector
        for (long j = 1; j <= (long)_n; ++j)
        {
            dE(j) = SMALL / _XW(j,1);
        }

        // Scaled matrix product of Wronskian matrices
        _E = multiply_G( dE );

        /// if ( floquet ) exit_with_solution();
    }
    else
    {
        // *** if rank == 0
        _ierr = 81;
        _FM = Matrix();

        return _FM;
    }


    for (long k = 1; k <= (long)_n; ++k)
    {
        Real s = 1.0 / dE(k);

        _E.set_colm(k) = s * _E.colm(k);
    }


    // Using ``standard´´ software for computation of eigenvalues
    //  for computation of multipliers of Wronskian
    balance( _E, low, igh, t2 );
    orthes( low, igh, _E, t2 );
    hqr( low, igh, _E, t1, t2, ifail );

    if ( ifail != 0 )
    {
        _ierr = (unsigned)ifail;
    }

    _FM.set_colm(1) = t1;
    _FM.set_colm(2) = t2;

    return _FM;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
IOpt
MultipleShootingGN::getIOpt()
{
    return _iopt;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
int
MultipleShootingGN::initialise(
                                 Vector const&  tnodes,
                                 Matrix const&  X,
                                 Real const     period,
                                 Real const     rtol,
                                 IOpt const&    iopt
                              )
{
    bool    qsucc = false;
    bool    qrank1 = iopt.qrank1;
    // bool    qfcstr = false;
    int     jacgen = iopt.jacgen;
    int     rscal = iopt.rscal;

    _n = X.nr();
    _m = X.nc();

    _period = period;

    _tolmin = 10.0*_n*EPMACH;
    _tolf = _tolj = rtol;
    _ierr = 0;

    setIOpt(iopt);

    ///

    _tnodes = tnodes;
    _X = X;

    if ( !qsucc )
    {
        printl( _lumon, dlib::LINFO,
                "\n    %s\n    %s\n    %s\n    %s\n\n",
                "****** PERIOD : Multiple Shooting ******",
                "* Method for periodic solution of ODEs *",
                "* via nonlinear least-squares approach *",
                "****************************************"
              );
    }

    check_init();
    if ( _ierr != 0 ) { return _ierr; }


    if ( rscal == 0 )  { _iopt.rscal  = 1; }
    if ( jacgen == 0 ) { _iopt.jacgen = 2; }


    if ( !qsucc )
    {
        printl( _lumon, dlib::LINFO,
                " %s %d\n %s %d\n %s %d\n\n %s %9.2e\n %s %9.2e\n\n",
                "Number of ODE equations                   (N):    ", _n,
                "Number of shooting nodes in unit interval (M):    ", _m,
                "Prescribed relative precision             (TOLF): ", _tolf,
                "                                          (TOLJ): ", _tolj
              );

        std::string shooting;
        if ( _m == 2)
            shooting = "Single";
        else
            shooting = "Multiple";
        printl( _lumon, dlib::LINFO,
                " %s shooting method is being used.\n\n", shooting.c_str()
              );

        std::string jacg, rsmode;
        if ( jacgen == 1 )
            jacg = "a user function";
        else if ( jacgen == 2 )
            jacg = "numerical differentiation (without feedback strategy)";
        else if ( jacgen == 3 )
            jacg = "numerical differentiation (feedback strategy applied)";
        printl( _lumon, dlib::LINFO,
                " The Jacobian is supplied by\n %s.\n", jacg.c_str()
              );
        rsmode = ( _iopt.norowscal == true) ? "inhibited" : "allowed";
        printl( _lumon, dlib::LINFO,
                " Automatic row scaling of the Jacobian is %s.\n\n", rsmode.c_str()
              );
    }

    _nonlin = _iopt.nonlin;

    if ( !qsucc )
    {
        printl( _lumon, dlib::LINFO,
                " Rank-1 updates are %s.\n",
                ((qrank1 == true)? "allowed" : "inhibited")
              );

        std::string nonlintype;
        if ( _nonlin == 1 ) nonlintype = "linear";
        else if ( _nonlin == 2 ) nonlintype = "mildly nonlinear";
        else if ( _nonlin == 3 ) nonlintype = "highly nonlinear";
        else if ( _nonlin == 4 ) nonlintype = "extremely nonlinear";

        printl( _lumon, dlib::LINFO,
                " Problem is specified as being %s.\n", nonlintype.c_str()
              );

        _niter = _ny = _kount = _new = 0;
    }

    return 0;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
MultipleShootingGN::check_init()
{
    _ierr = 0;

    // Checking dimensional parameters N,and M
    if ( (_n <= 0) || (_m <= 1) || (_m != (unsigned)_tnodes.nr()))
    {
        printl( _luerr, dlib::LINFO,
                " Error: %s\n        %s\n        %s%d%s%d%s%d\n",
                "Bad or inconsistent input to dimensional parameters supplied.",
                "Choose N and M positive.",
                "Input is: N = ", _n, " M = ", _m, " #TNODES = ", _tnodes.nr()
              );
        _ierr = 20;
    }

    // Problem type specified by user
    _nonlin = _iopt.nonlin;
    if ( _nonlin == 0 ) _nonlin = 3;
    _iopt.nonlin = _nonlin;

    // Checking of user-precribed RTOL
    if ( (_tolf <= 0.0) || (_tolj <= 0.0) )
    {
        printl( _luerr, dlib::LINFO,
                " Error: Nonpositive RTOL supplied.\n"
              );
        _ierr = 21;
    }

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
int
MultipleShootingGN::fire()
{
    const Real  fcs = 0.7;
    const Real  redh = 1.0e-2;

    Real        /*hstart,*/ eph = 1.0, condh/*, cond1 = 1.0*/;
    Real        sens1 = 0.0/*, rsmall = 0.0*/;
    Real        tol;
    Real        sumxa, /* conva,*/ fch, fcmin2, fcminh;
    Real        th, fcnum, fcdnm;
    unsigned    n1, m1;
    unsigned    jacgen, jacgenv, levl, /*mode,*/ nred;
    unsigned    irkmax, iranka = 0;
    bool        qsucc/*, qinisc*/;      // qsucc: called by successive one-step
    bool        qinit;
    bool        qiter, qjcrfr;      // qjcrfr: Jacobian refresh flag
    bool        /*qgenj,*/ qnext, qredrnk, qrepet;
    bool        qredu, qred;
    int         iauto, ibroy, ifail;
    Real        ddp;
    Matrix      DDX;
    Vector      u, dx1, t1, t2, dE;

    if ( _ode == 0 )
    {
        _ierr = 999;
        printl( _luerr, dlib::LINFO,
                " Error: %s\n        %s\n",
                "No problem (FirstOrderODESystem) set so far.",
                "Please use setProblem(&ode) first.");
        return _ierr;
    }

    fcnum = fcdnm = sens1 = 0.0;
    condh = 1.0;

    qsucc   = false;
    iauto   = _iopt.iauto;
    ibroy   = ( _iopt.qrank1 == true ) ? 1 : 0;
    jacgenv = _iopt.jacgen;
    // mode    = _iopt.mode;

    n1      = ( iauto > 0 ) ? (_n + 1) : _n;
    m1      = ( _m > 0 ) ? (_m - 1) : 0;
    irkmax  = std::min(_n, _m);
    fcmin2  = _fcmin * _fcmin;
    fcminh  = std::sqrt( _fcmin );
    // hstart  = (_tnodes(2) - _tnodes(1)) * _period * redh;

    if ( _fc < _fcmin ) _fc = _fcmin;
    if ( _fc > 1.0 ) _fc = 1.0;

    // Initial preparations
    qjcrfr = qrepet = false;
    qiter = qinit = true;
    ifail = 0;

    // Miscellaneous preparations of first iteration step
    if ( !qsucc )
    {
        _niter = _kount = 0;

        irkmax =
        iranka = _irank;
        _fca   = _fc;
        _conv  = 0.0;

        // Print monitor header
        printl( _lumon, dlib::LVERB,
                "\n\n %s\n %s\n",
                "******************************************************************",
                "   It   Ny     Normf          Normx         Damp.Fct.   New   Rank");

        // Comp. initial residual vector
        _HH = call_IVPSOL( _tnodes, _period, _X, _XU, ifail );

        ++_kount;

        if ( ifail != 0 )
        {
            // Singular trajectory
            _ierr = 82;
            qiter = false;
        }

    }
    else
    {
        // qinisc = false;
    }

    // Main iteration loop
    // ===================

    while ( qiter )
    {
        // Startup of iteration step
        if ( !qjcrfr )
        {
            // Scale unkowns X
            compute_scaling_XW( _xthr );

            if ( iauto > 0 )
            {
                _periodw = _period;
            }

            // qinisc = false;

            if ( !qinit )
            {
                // Save DXQ (and dpq, if neccessary)
                _DXQA = _DXQ;

                if ( iauto > 0 )
                {
                    _dpqa = _dpq;
                }

                // Prelim. pseudo-rank
                iranka = _irank;

                if ( (_irank >= 0) && (_fc > fcminh) )
                {
                    _irank = _n;
                }

                // Est. damping factor (a posteriori)
                // ----------------------------------
                th = _fc - 1.0;
                DDX = _DXQ + th * _DX;

                if ( iauto == 1 )
                {
                    ddp = _dpq + th * _dp;
                }

                fcnum = scalarprod( _DX,_DX, _dp,_dp, iauto );
                fcdnm = scalarprod( DDX,DDX, ddp,ddp, iauto );

                fch = std::sqrt( fcnum / fcdnm ) * _fc * _fc * 0.5;

                // Decisision criterion Jacobian update technique
                //   jacgen = 0:    Numerical differentiation,
                //   jacgen = 1,2:  Integration of Var.Eqn.,
                //   jacgen = 3:    Rank-1 updating
                // ----------------------------------------------
                if ( ibroy != 0 )
                {
                    jacgen = ( ((_fc < _fca) && (_new > 0)) ||
                               (fch < _fc*_sigma) ||
                               (eph*redh > _eps) ||
                               (_irank > iranka) ) ? jacgenv : 3;
                }

                _fca = _fc;

                if ( _nonlin > 1 )
                {
                    _fc = std::min( 1.0, fch );
                }

            } // if ( !qinit )

        } // if ( !qjcrfr )

        qjcrfr = false;
        iranka = _irank;

        // Jacobian matrix generation
        // --------------------------
        if ( jacgen != 3 )
        {
            _new = 0;
            ifail = 0;

            if ( !(jacgen > 0) )
            {
                // Diff. approx. of Wronskian matrices G(1), ..., G(m1)
                // ----------------------------------------------------
                compute_fd_derivative_G( _X, _period, ifail );

                _kount += _n;
            }
            else
            {
                // Wronskian G(1), ..., G(m1) by integration of variational equation
                // -----------------------------------------------------------------
                if ( !qinit )
                {
                    // Adaption of _tolj and _tolf
                    _tolj = std::sqrt( _sumx );
                    if ( _tolj > redh ) _tolj = redh;
                    if ( _tolj < _tolmin ) _tolj = _tolmin;

                    _tolf = tol * 0.1;
                }

                solve_var_eq_for_G( _X, _period, ifail );

                _kount += _n;
            }
        }
        else
        {
            // Rank-1 updates of Wronskian matrices G(1), ..., G(m1)
            // -----------------------------------------------------
            ++_new;
            compute_rank1_update_G( iauto );
        }

        // Comp. of sensitivity matrix E
        // -----------------------------
        if ( _irank != 0 )
        {
            // Storing row scaling vector
            for (long j = 1; j <= (long)_n; ++j)
            {
                dE(j) = SMALL / _XW(j,1);
            }

            // Scaled matrix product of Wronskian matrices
            _E = multiply_G( dE );

            /// if ( floquet ) exit_with_solution();

            // Internal row and column scaling
            for (long k = 1; k <= (long)_n; ++k)
            {
                Real s = _XW(k,1);

                // for (long j = 1; j >= (long)_n; ++j)
                // {
                //     _E(j,k) = -s * _E(j,k);
                // }
                _E.set_colm(k) = -s * _E.colm(k);

                _E(k,k) += SMALL;
            }

            // Extended matrix E
            if ( iauto != 0 )
            {
                for (long j = 1; j <= (long)m1; ++j)
                {
                    long j1 = j+1;
                    Real tj1 = _tnodes(j1) * _period;
                    Real dt = _tnodes(j1) - _tnodes(j);

                    t1 = _XU.colm(j);
                    t2 = call_FCN( tj1, t1, ifail );

                    _FP.set_colm(j) = dt * t2;

                }

                t2 = compute_condensed_RHSp( 1, _FP, t2, dE );

                _E.set_colm(n1) = -_periodw * t2;
            }

            // Monitor for actually applied maximum rank
            if ( irkmax < _irank )
            {
                irkmax = _irank;
            }

            // Save curr. residual matrix _HH (and _FP, if applicable)
            if ( !qrepet )
            {
                _HHA = _HH;

                if ( iauto != 0 )
                {
                    _FPA = _FP;
                }
            }
        }

        qnext = false;
        qredrnk = true;

        // Central part of iteration step

        // Pseudo-rank reduction loop
        // ==========================

        while ( qredrnk )
        {
            // QR Decomposition of sensitivity _E
            // ----------------------------------
            if ( _irank > 0 )
            {
                condh = _cond;
                compute_qrdcmp_E( qrepet, condh, ifail );

                if ( ifail != 0 ) { _ierr = 80; qiter = false; break; }
            }

            // Comp. of condensed right-hand side u
            // ------------------------------------
            if ( !qrepet )
            {
                if ( _irank > 0 )
                {
                    u = compute_condensed_RHSp( 1, _HH, u, dE );
                }
            }

            // LSQ Solution of linear system
            // -----------------------------
            if ( _irank > 0 )
            {
                dx1 = solve_E( qrepet, u, ifail );

                if ( ifail != 0 ) { _ierr = 81; qiter = false; break; }
            }

            if ( !qrepet && (_irank != 0) )
            {
                // _qu.zeros(_n);
                _qu.set_row(1,_irank) = u.row(1,_irank);
            }

            // Descaling of solution dx1
            for (long j = 1; j <= (long)_n; ++j)
            {
                _DX(j,1) = dx1(j) * _XW(j,1);
            }
            if ( iauto > 0 )
            {
                    _dp = _periodw * dx1(n1);
            }

            // Successive computation of _DX(N,2), ..., _DX(N,M)
            // -------------------------------------------------
            compute_rest_of( _DX, _dp, 1, _HH, iauto );

            // ----------------------------------------------
            // Iterative refinement sweep _ny = 1, ..., nymax
            // ----------------------------------------------
            tol = compute_refinement_sweep(
                                            qrepet, iauto, levl,
                                            _DX, _dp, dE, tol, ifail
                                          );

            if ( ifail != 0 ) { _ierr = 82; qiter = false; break; }

            _tolf = tol;

            // Eval. scaled std. level fun.
            // ----------------------------
            _sumf = 0.0;

            for (long j = 1; j <= (long)m1; ++j)
            {
                long j1 = j+1;

                for (long k = 1; k <= (long)_n; ++k)
                {
                    Real tmp = _HH(k,j) / _XW(k,j1);

                    _sumf += (tmp * tmp);
                }
            }

            // ---------------------------------------------------
            // Projection Moore-Penrose pseudo-inverse of Jacobian
            // ---------------------------------------------------
            if ( !((iauto == 0) || (_irank < _n)) )
            {
                // Vector       pivot = _qrE.getPivot();
                unsigned    pvs = _qrE.getPivot()(n1); // _pivot(n1);
                Vector      v = _qrE.getMatH().colm(_n);

                for (long j = 1; j <= (long)_n; ++j)
                {
                    unsigned pj = _qrE.getPivot()(j); // _pivot(j);

                    if (pj != n1)
                    {
                        _XTG(pj, 1) = -v(j) * _XW(pj, 1);
                    }
                    else
                    {
                        _ptg = -v(j) * _periodw;
                    }
                }

                if ( pvs != n1 )
                {
                    _XTG(pvs, 1) = _XW(pvs, 1);
                }
                else
                {
                    _ptg = _periodw;
                }

                _HH.zeros(_n,m1);

                compute_rest_of( _XTG, _ptg, 1, _HH, iauto );
                tol = compute_refinement_sweep(
                                                qrepet, iauto, levl,
                                                _XTG, _ptg, dE, tol, ifail
                                              );

                if ( ifail != 0 ) { _ierr = 83; qiter = false; break; }

                Real s  = scalarprod(  _DX,_XTG,  _dp,_ptg, iauto );
                Real st = scalarprod( _XTG,_XTG, _ptg,_ptg, iauto );

                s /= st;

                _DX = _DX + (-s) * _XTG;
                _dp = _dp + (-s) * _ptg;
            }

            // Eval. scaled natural level fun. and scaled max. error norm _conv
            // ----------------------------------------------------------------
            _sumx = scalarprod( _DX,_DX, _dp,_dp, iauto );
            _conv = 0.0;
            for (long j = 1; j <= (long)m1; ++j)
            {
                for (long k = 1; k <= (long)_n; ++k)
                {
                    Real s = std::fabs( _DX(k,j) ) / _XW(k,j);

                    if ( _conv < s ) { _conv = s; }
                }
            }
            if ( iauto != 0 )
            {
                Real s = std::fabs( _dp ) / _periodw;

                if ( _conv < s ) { _conv = s; }
            }

            // ---------------------------
            // Ordinary G-N correction _DX
            // ---------------------------
            _XA = _X;
            // _DX = _DXQ;
            if ( iauto != 0 )
            {
                _perioda = _period;
                // _dp = _dpq;
            }

            // Eval. of subcondition and sensitivity number
            sumxa = _sumx;
            // conva = _conv;
            // cond1 = 1.0;
            sens1 = 0.0;
            if ( _irank != 0 )
            {
               sens1 = std::fabs( _qrE.getDiag()(1) );
               // cond1 = sens1 / std::fabs( _qrE.getDiag()(_irank) );
            }

            // Est. damping factor (a priori)
            // ------------------------------
            nred = 0;
            qredu = false;

            if ( !((_niter == 0) || ((_nonlin == 1) && (iauto == 0)) ) )
            {
                if ( !(((_new > 0) || (_irank < _n && iranka < _n)) && !qrepet) )
                {
                    // Denom. a priori (Full Rank)
                    DDX = _DX - _DXQA;
                    if ( iauto == 1 )
                    {
                        ddp = _dp - _dpqa;
                    }

                    fcdnm = scalarprod( DDX,DDX, ddp,ddp, iauto );

                    // Denom. a priori (projected)
                    if ( !(_irank < _n) )
                    {
                        if ( iauto != 0 )
                        {
                            Real sum1 = scalarprod( _DXQA,_XTG, _dpqa,_ptg, iauto );
                            Real sum2 = scalarprod(  _XTG,_XTG,  _ptg,_ptg, iauto );
                            Real  del = sum1*sum1/sum2;

                            fcdnm -= del;
                        }

                        // New damp. fact.
                        _fc = (fcdnm > (fcnum*fcmin2)) ?
                                    _fca * std::sqrt(fcnum / fcdnm) :
                                    _fca / _fcmin;
                    }

                } // if ( _new == 0 || qrepet )

                qredu = ( _fc < _fcmin );

                if ( _fc > fcs ) { _fc = 1.0; }

            } // if ( (_niter != 0) && (_nonlin != 1) )

            qrepet = false;

            if ( !qredu )
            {
                log_iteration_vals1( sumxa );

                nred = 0;
                qred = true;
                // Damping factor reduction loop
                // =============================
                while ( qred )
                {
                    // ----------------------------------------------
                    // Prelim. new iterate, exit if prob. only linear
                    // ----------------------------------------------
                    qinit = false;

                    _X = _XA + _fc * _DX;

                    if ( iauto > 0 ) { _period = _perioda + _fc * _dp; }

                    for (long j = 1; j <= (long)_n; ++j)
                    {
                        _X(j, _m) = _X(j, 1);
                    }

                    if ( _nonlin == 1 ) { _ierr = 0; qiter = false; break; }

                    // Plus corresp. new residual
                    // --------------------------
                    _HH = call_IVPSOL( _tnodes, _period, _X, _XU, ifail );

                    ++_kount;

                    if ( (ifail == 1) || (ifail == 2) )
                    {
                        _fc *= 0.5;

                        if ( _fc < _fcmin )
                        {
                            qiter = false; break;
                        }

                        break;
                    }

                    levl = 1;

                    // Comp. of condensed right-hand side u
                    // ------------------------------------
                    if ( !qrepet )
                    {
                        if ( _irank > 0 )
                        {
                            u = compute_condensed_RHSp( 1, _HH, u, dE );
                        }
                    }

                    // LSQ Solution
                    // ------------
                    if ( _irank > 0 )
                    {
                        dx1 = solve_E( qrepet, u, ifail );
                    }

                    // Descaling of solution dx1
                    for (long j = 1; j <= (long)_n; ++j)
                    {
                        _DXQ(j,1) = dx1(j) * _XW(j,1);
                    }
                    if ( iauto > 0 )
                    {
                        _dpq = _periodw * dx1(n1);
                    }

                    // Successive computation of _DXQ(N,2), ..., _DXQ(N,M)
                    // ---------------------------------------------------
                    compute_rest_of( _DXQ, _dpq, 1, _HH, iauto );

                    // -----------------------------------------------
                    // Iterative refinement sweeps _ny = 1, ..., nymax
                    // -----------------------------------------------
                    tol = compute_refinement_sweep(
                                                    qrepet, iauto, levl,
                                                    _DXQ, _dpq, dE, tol, ifail
                                                  );

                    if ( ifail != 0 ) { _ierr = 82; qiter = false; break; }

                    _tolf = tol;

                    // Eval. scaled std. level fun.
                    // ----------------------------
                    _sumf = 0.0;

                    for (long j = 1; j <= (long)m1; ++j)
                    {
                        long j1 = j+1;

                        for (long k = 1; k <= (long)_n; ++k)
                        {
                            Real tmp = _HH(k,j) / _XW(k,j1);

                            _sumf += (tmp * tmp);
                        }
                    }

                    // ---------------------------------------------------
                    // Projection Moore-Penrose pseudo-inverse of Jacobian
                    // ---------------------------------------------------
                    if ( !((iauto == 0) || (_irank < _n)) )
                    {
                        Real s  = scalarprod( _DXQ,_XTG, _dpq,_ptg, iauto);
                        Real st = scalarprod( _XTG,_XTG, _ptg,_ptg, iauto);

                        s /= st;

                        _DXQ = _DXQ - s * _XTG;
                        _dpq = _dpq - s * _ptg;
                    }

                    // Eval. scaled natural level fun. and
                    // scaled max. error norm _conv
                    // -----------------------------------
                    _sumx = scalarprod( _DXQ,_DXQ, _dpq,_dpq, iauto );
                    _conv = 0.0;
                    for (long j = 1; j <= (long)m1; ++j)
                    {
                        for (long k = 1; k <= (long)_n; ++k)
                        {
                            Real s = std::fabs( _DX(k,j) ) / _XW(k,j);

                            if ( _conv < s ) { _conv = s; }
                        }
                    }
                    if ( iauto != 0 )
                    {
                        Real s = std::fabs( _dp ) / _periodw;

                        if ( _conv < s ) { _conv = s; }
                    }

                    // ------------------------------
                    // Simplified G-N correction _DXQ
                    // ------------------------------

                    //
                    // Rank independent convergence test
                    //
                    if ( (_conv <= _eps) && (irkmax == _n) )
                    {
                        _ierr = 0; qiter = false; break;
                    }

                    //
                    // Natural monotonicity test
                    //
                    if ( _sumx > sumxa )
                    {
                        // Output of curr. iterate
                        log_iteration_vals2();

                        // Eval. reduced damping factor
                        // ----------------------------
                        th = std::sqrt( _sumx/ sumxa );
                        th = std::sqrt( 1.0 + 8.0*(th + _fc - 1.0)/_fc ) - 1.0;

                        _fc /= th;

                        ++nred;
                        // Rank reduction if damping is too small
                        qredu = (_fc < _fcmin) || ((_new > 0) && (nred > 1));
                    }
                    else
                    {
                        qnext = true;
                    }

                    qred = !(qnext || qredu);

                } // while ( qred )

                if ( (!qredrnk) || (!qiter) ) break;

            } // if ( !qredu )

            // End of damp. fact. reduction loop
            // ---------------------------------

            if ( qredu )
            {
                // Restore former values for repeating step
                // ----------------------------------------
                qrepet = true;
                levl = 0;

                if ( iauto != 0 )
                {
                    _FP = _FPA;
                    _period = _perioda;
                }

                _X = _XA;
                _XU = _X + _HHA;
                _HH = _HHA;

                printl( _lumon, dlib::LVERB,
                        "  %4d %38s %7.5f     %2d   %4d\n",
                        _niter, "not accepted damping factor:", _fc, _new, _irank
                     );

                if ( _niter == 0 ) { _fc = _fcmin; }

                jacgen = jacgenv;

                if ( _new > 0 )
                {
                    qjcrfr = true;
                    qredu = false;
                    _irank = _n;
                }
                else
                {
                    // Pseudo-rank reduction
                    // ---------------------
                    qrepet = true;

                    for (long j = 1; j <= (long)_irank; ++j)
                    {
                        u(j) = _qu(j);
                    }

                    --_irank;
                }

            } // if ( qredu )

            qredrnk = qredu;

        } // while ( qredrnk )

        if ( !qiter ) break;

        if ( qnext )
        {
            // ----------------------------------
            // Prep. to start next iteration step
            // ----------------------------------
            ++_niter;
            levl = 0;

            // Print curr. values
            log_iteration_vals2();

            // Exit if in single-step mode, or max. _niter reached
            if ( _niter > _nitmax ) { _ierr = 2; qiter = false; break; }

        } // if ( qnext )

    } // while ( qiter )


    ++_niter;


    _X = _X + _DXQ;

    if ( iauto > 0 ) { _period = _period + _dpq; }

    for (long j = 1; j <= (long)_n; ++j)
    {
        _X(j, _m) = _X(j, 1);
    }

    _tolj = tol;


    return 0;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Matrix
MultipleShootingGN::call_IVPSOL(Vector const& tnodes, Real period,
                                Matrix const& XX, Matrix& XU, int& ifail)
{
    Matrix      HH;
    double*    y = new double[_n];

    ifail = 0;
    HH.zeros(_n, _m-1);

    _odeSolver -> setRTol( _tolf );

    for (unsigned j = 1; j < _m; ++j)
    {
        double      tj  = period*tnodes(j);
        double      tj1 = period*tnodes(j+1);
        unsigned    k = 0;

        k = 0;
        while (k < _n)
        {
            y[k++] = XX(k,j);
        }

        dynamic_cast<LIMEX_A*>(_odeSolver) -> resetSuccessiveCallFlag();

        ifail = _odeSolver -> integrate(_n, y, tj, tj1);

        if (ifail == 0)
        {
            k = 0;
            while (k < _n)
            {
                double tmp = y[k++];

                XU(k,j) = tmp;
                HH(k,j) = tmp - XX(k,j+1);
            }
        }
        else
        {
            break;
        }
    }

    delete[] y;

    return HH;
}
//---------------------------------------------------------------------------
Vector
MultipleShootingGN::call_FCN(Real t, Vector const& y, int& ifail)
{
    Vector fcn;
    double* dy = new double[_n];
    double* yy = new double[_n];

    for (unsigned j = 1; j <= _n; ++j)
    {
        yy[j-1] = y(j);
    }

    _ode -> computeDerivatives(t, yy, dy, &ifail);

    fcn.zeros(_n);
    for (unsigned j = 1; j <= _n; ++j)
    {
        fcn(j) = *dy++;
    }

    delete[] yy;
    delete[] dy;
    ifail = 0;

    return fcn;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
MultipleShootingGN::compute_scaling_XW(Real& xthr)
{
    const Real red = 1.0e-2;

    Real xmax;

    _XW.zeros(_n,_m);

    for (unsigned j = 1; j <= _n; ++j)
    {
        _XW(j,1) = std::fabs( _X(j,1) );
    }

    // Arithmetic mean for _XW(_n,2), ..., _XW(_n,_m)
    for (unsigned k = 2; k <= _m; ++k)
    {
        for (unsigned j = 1; j <= _n; ++j)
        {
            _XW(j,k) = 0.5* ( std::fabs( _X(j,k) ) + std::fabs( _XU(j,k-1) ) );
        }
    }

    // Threshold determination
    for (unsigned j = 1; j <= _n; ++j)
    {
        xmax = 0.0;

        for (unsigned k = 1; k <= _m; ++k)
        {
            Real tmp = _XW(j,k);

            if ( xmax < tmp ) xmax = tmp;
        }

        if ( xmax < xthr*red ) xmax = xthr;

        xmax *= red;

        if ( xmax < SMALL ) xmax = 1.0;

        for (unsigned k = 1; k < _m; ++k)
        {
            if ( xmax > _XW(j,k) ) _XW(j,k) = xmax;
        }

        _XW(j,_m) = _XW(j,1);
    }

    xthr = xmax;

    return;
    // return xmax;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Real
MultipleShootingGN::scalarprod(Matrix const& X1,
                               Matrix const& X2,
                               Real p1,
                               Real p2,
                               int iauto)
{
    Real result = 0.0;
    Real rl = _tnodes(_m) - _tnodes(1);

    for (unsigned k = 1; k <= _m; ++k)
    {
        Real s, sum = 0.0;

        if ( k == 1 )
        {
            s = _tnodes(2) - _tnodes(1);
        }
        else if ( k == _m )
        {
            s = _tnodes(_m) - _tnodes(_m-1);
        }
        else
        {
            s = _tnodes(k+1) - _tnodes(k-1);
        }

        s /= rl;

        for (unsigned j = 1; j <= _n; ++j)
        {
            Real s1 = X1(j,k) / _XW(j,k);
            Real s2 = X2(j,k) / _XW(j,k);

            sum += (s1*s2);
        }

        result += (s*sum);
    }

    result *= 0.5;

    if ( iauto == 1 )
    {
        result += ( (p1/_periodw) * (p2/_periodw) );
    }

    return result;
}
//---------------------------------------------------------------------------

//===========================================================================

//---------------------------------------------------------------------------
// Difference approximation of Wronskian matrices _G[0], ..., _G[_m-2]
//---------------------------------------------------------------------------
void
MultipleShootingGN::compute_fd_derivative_G(Matrix const& XX,
                                            Real period,
                                            int& ifail)
{
    double* xjj = new double[_n];

    ifail = 0;

    _odeSolver -> setRTol( _tolf );

    for (unsigned jj = 1; jj < _m; ++jj)
    {
        double      tja = period * _tnodes(jj);
        double      tj1 = period * _tnodes(jj+1);

        for (unsigned j = 1; j <= _n; ++j)
        {
            double      th, s;
            unsigned    k = 0;

            k = 0;
            while (k < _n)
            {
                xjj[k++] = XX(k,jj);
            }

            th = xjj[j-1];

            s = _XW(j,jj) * _reldif;
            s *= ( (th < 0.0) ? -1.0 : 1.0 );

            xjj[j-1] = th + s;

            dynamic_cast<LIMEX_A*>(_odeSolver) -> resetSuccessiveCallFlag();

            ifail = _odeSolver -> integrate(_n, xjj, tja, tj1);

            if (ifail == 0)
            {
                k = 0;
                while (k < _n)
                {
                    double tmp = xjj[k++];

                    _G[jj-1](k,j) = ( tmp - _XU(k,jj) ) / s;
                }
            }
            else
            {
                delete[] xjj;
                return;
            }
        }
    }

    delete[] xjj;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Integration of variational equation
//---------------------------------------------------------------------------
void
MultipleShootingGN::solve_var_eq_for_G( Matrix const& XX,
                                        Real period,
                                        int& ifail)
{
    double                tstart = 0.0;
    double                tend = 1.0;
    ODESolver::Grid     dummy;
    unsigned           nint = _n + _n*_n;
    double*            Y = new double[nint];

    ifail = 0;

    _odeSolver -> setRTol( _tolf );

    dynamic_cast<LIMEX_A*>(_odeSolver) ->
            setODESystem(
                            *_ode,
                            tstart, dummy, tend
                         );

    for (unsigned j = 1; j < _m; ++j)
    {
        unsigned    k;
        double      tj  = period * _tnodes(j);
        double      tj1 = period * _tnodes(j+1);

        k = 0;
        while (k < _n)
        {
            Y[k] = XX(k+1,j);
            ++k;
        }
        while (k < nint)
        {
            // set the remaining _n x _n part to the identity matrix
            Y[k] = ((k+1) % _n == k/_n) ? 1.0 : 0.0;
            ++k;
        }

        dynamic_cast<LIMEX_A*>(_odeSolver) -> resetSuccessiveCallFlag();

        ifail = _odeSolver -> integrate(nint, Y, tj, tj1);

        if (ifail == 0)
        {
            for (unsigned i = 1; i <= _n; ++i)
            {
                k = 0;
                while (k < _n)
                {
                    // access element with index: ( _n*(i-1) + k+1 )
                    double tmp = Y[_n*(i-1) + ++k];

                    _G[j-1](i,k) = tmp;
                    // or is...   _G[j-1](k,i)=tmp;  ...rather correct??!?
                }
            }
        }
        else
        {
            break;
        }
    }

    delete[] Y;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Rank-1 updates of Wronskian matrices _G[0], ..., _G[_m-2]
//---------------------------------------------------------------------------
void
MultipleShootingGN::compute_rank1_update_G(int iauto)
{
    Vector  dxj;
    Real    fch = _fca - 1.0;

    dxj.zeros(_n);

    for (unsigned k = 1; k < _m; ++k)
    {
        Real dnm = 0.0;

        for (unsigned j = 1; j <= _n; ++j)
        {
            Real tmp = _DX(j,k) / _XW(j,k);

            dxj(j) = tmp / _XW(j,k);

            dnm += (tmp * tmp);
        }

        dnm *= _fca;

        if ( dnm == 0.0 ) continue;

        for (unsigned i = 1; i <= _n; ++i)
        {
            Real t = dxj(i) / dnm;

            for (unsigned j = 1; j <= _n; ++j)
            {
                Real s = _G[k-1](j,i);

                if ( s == 0.0 ) continue;

                s += t * ( _HH(j,k) + fch*_HHA(j,k) );

                if ( iauto > 0 )
                {
                    s += t * _fca * _dp * ( _FPA(j,k) - _FP(j,k) );
                }

                _G[k-1](j,i) = s;
            }
        }
    }
    // return;
}
//---------------------------------------------------------------------------

//===========================================================================

//---------------------------------------------------------------------------
Matrix
MultipleShootingGN::multiply_G(Vector const& dE)
{
    Matrix E = dE.diag();

    // E.zeros(_n,_n);

    for (unsigned jj = 1; jj < _m; ++jj)
    {
        unsigned j = _m - jj - 1;

        for (unsigned i = 1; i <= _n; ++i)
        {
            Vector v;

            v.zeros(_n);

            for (unsigned k = 1; k <= _n; ++k)
            {
                Real s = 0.0;

                for (unsigned ell = 1; ell <= _n; ++ell)
                {
                    s += ( E(i,ell) * _G[j](ell,k) );
                }

                v(k) = s;
            }

            E.set_rowm(i) = v.t();

            // for (unsigned k = 1; k <= _n; ++k)
            // {
            //    E(i,k) = t1(k)
            // }
        }
    }

    return E;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Vector
MultipleShootingGN::compute_condensed_RHSp(unsigned jin,
                                           Matrix const& HH,
                                           Vector const& u,
                                           Vector const& dE)
{
    Vector v,w;
    unsigned m1 = _m - 1;

    if ( jin > m1 )
    {
        return u;
    }

    _BG = dE.diag();
    w.zeros(_n);

    for (unsigned j = 1; j <= _n; ++j)
    {
        w(j) = dE(j) * HH(j,m1);
    }

    if ( (m1 == 1) || (jin == m1) )
    {
        return w;
    }

    unsigned m2 = m1 - 1;
    v.zeros(_n);

    for (unsigned jj = jin; jj <= m2; ++jj)
    {
        unsigned j = m2 + jin - jj;

        for (unsigned i = 1; i <= _n; ++i)
        {
            Real s1 = w(i);

            for (unsigned k = 1; k <= _n; ++k)
            {
                Real s2 = 0.0;

                for (unsigned ell = 1; ell <= _n; ++ell)
                {
                    s2 += _BG(i,ell) * _G[j](ell,k);
                }

                v(k) = s2;
            }

            for (unsigned k = 1; k <= _n; ++k)
            {
                s1 += v(k) * HH(k,j);

                _BG(i,k) = v(k);
            }

            w(i) = s1;
        }
    }

    return w;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
MultipleShootingGN::compute_rest_of(Matrix& DX, Real dp,
                                    unsigned jin, Matrix const& HH,
                                    int iauto )
{
    // Vector       v;
    Vector       w = DX.colm(1);
    unsigned    m1 = _m - 1;

    // v.zeros(_n);

    for (unsigned j = 1; j <= m1; ++j)
    {
        unsigned j_1 = j - 1;
        unsigned j1  = j + 1;

        for (unsigned i = 1; i <= _n; ++i)
        {
            Real s = (jin <= j) ? HH(i,j) : 0.0;

            if ( iauto > 0 )
            {
                s += dp * _FP(i,j);
            }

            for (unsigned k = 1; k <= _n; ++k)
            {
                s += w(k) * _G[j_1](i,k);
            }

            // v(i) = s;

            DX(i, j1) = s;
        }

        // w = v;
        w = DX.colm(j1);
    }

    // return w;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Real
MultipleShootingGN::compute_refinement_sweep(bool qrepet, int iauto, int levl,
                                             Matrix& DXQ, Real& dpq,
                                             Vector const& dE, Real tol,
                                             int& ierr)
{
    const Real  redh = 1.0e-2;

    unsigned    jn = 1, jin = _m, ja = 0;
    Real         sigdel = 10.0;
    Real         sigdlh = 0.0;
    Real         th, tolh, corr, eph = _eps;
    Vector       dx1, du, rf;
    Matrix       DHH, DDX;
    Real         ddp = 0.0;

    _ny = 0;
    ierr = 0;

    if ( (_nymax == 0) || _irank < _n)
    {
        // some warning/error message
        return 0.0;
    }

    dx1.zeros(_n);
    du.zeros(_n);
    rf.zeros(_m);
    DHH.zeros(_n,_m);
    DDX.zeros(_n,_m);

    while( jn != _m )
    {
        if (_ny != 0)
        {
            for (unsigned j = jn; j < _m; ++j)
            {
                unsigned j1  = j + 1;
                unsigned j_1 = j - 1;

                for (unsigned i = 1; i <= _n; ++i)
                {
                    Real s = _HH(i,j);

                    if ( iauto == 0 )
                    {
                        s += dpq * _FP(i,j);
                    }

                    for (unsigned k = 1; k <= _n; ++k)
                    {
                        s += _G[j_1](i,k) * DXQ(k,j);
                    }

                    DHH(i,j) = s - DXQ(i,j1);
                }
            }

            if ( _irank > 0 )
            {
                du = compute_condensed_RHSp(jin, DHH, du, dE);
            }

        } // if ( _ny != 0 )

        for (unsigned i = 1; i <= _n; ++i)
        {
            du(i) += dE(i) * ( DXQ(i,1) - DXQ(i,_m) );
        }

        if ( _irank > 0 )
        {
            dx1 = solve_E( qrepet, du, ierr );
        }

        corr = 0.0;

        for (unsigned ell = 1; ell <= _n; ++ell)
        {
            Real s = dx1(ell);

            if ( corr < std::fabs(s) ) { corr = std::fabs(s); }

            s *= _XW(ell,1);

            DDX(ell,1) = s;
            DXQ(ell,1) += s;
        }

        if (iauto != 0)
        {
            Real s = dx1(_n+1);
            if ( corr < std::fabs(s) ) { corr = std::fabs(s); }
            s *= _periodw;
            ddp = s;
            dpq += s;
        }

        if ( corr >= eph ) { /* eph = corr; */ ierr = 3; return tol; }

        rf(1) = corr;

        // Recursive computation of _DDX(_n,2), ..., _DDX(_n,_m)
        compute_rest_of( DDX, ddp, jin, DHH, iauto );

        // Refinement of _DXQ(_n,2), ..., _DXQ(_n,_m)
        for (unsigned j = 2; j <= _m; ++j)
        {
            corr = 0.0;

            for (unsigned i = 1; i <= _n; ++i)
            {
                Real s = DDX(i,j);

                DXQ(i,j) += s;

                s = std::fabs(s) / _XW(i,j);

                if ( corr < s ) { corr = s; }
            }

            rf(j) = corr;
        }

        // Determination of sweep index jn
        ja = jn;

        for (unsigned j = 1; j <= _m; ++j)
        {
            if ( rf(j) > eph ) { break; }
            jn = j;
        }

        ++_ny;

        if ( jn <= ja ) { ierr = 1; return tol; }

        if ( jn != _m )
        {
            jin = jn;

            if ( (_ny > 1) || (levl == 0) ) { continue; }
        }

        // Determination and adaptation of parameters tol and reldif
        if ( (levl == 0) || (_ny > 1) ) { continue; }

        for (unsigned j = 1; j < _m; ++j)
        {
            Real s = 0.0;

            if ( rf(j) != 0.0 )
            {
                s = rf(j+1) / rf(j);
            }
            if ( sigdlh < s )
            {
                sigdlh = s;
            }

            rf(j) = s;
        }

        sigdel = std::max(sigdlh, sigdel);
        th = tol * sigdel;

        if ( th > redh ) { ierr = 2; return tol; }

        if ( th > eph )
        {
            eph = th;
        }

        tolh = _eps / sigdel;

        if (tolh > _tolmin )
        {
            tolh = _tolmin;
        }

        tol = tolh;
        _reldif = std::sqrt( tol / sigdel );

    }

    ierr = 0;

    return tol;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
MultipleShootingGN::compute_qrdcmp_E(
    bool        qrepet,
    Real&        cond,
    int&        ifail
    )
{
    // unsigned irankc, mcon = 0;

//std::cerr << "*** MultipleShootingGN::compute_qrdcmp_E ***" << std::endl;
//std::cerr << "qrepet = " << (int)qrepet << std::endl;
//std::cerr << "cond   = " << cond << std::endl;
//std::cerr << "ifail  = " << ifail << std::endl;
//std::cerr << std::endl;
//std::cerr << "_irank = " << _irank << std::endl;
//std::cerr << std::endl;
//std::cerr << "_E = (pre)" << std::endl;
//std::cerr << _E.t() << std::endl;
//std::cerr << "_qrE = (pre)" << std::endl;
//std::cerr << _qrE.getMat().t() << std::endl;
//std::cerr << std::endl;

    if ( !qrepet )
    {
        _qrE = _E.factorQR(_irank, cond);
        // _qrE = _E.factorQRcon(mcon, _irank, cond);
    }
    else
    {
        _qrE.setNewRank(_irank);
    }

    _irank  = _qrE.getRank();
    cond    = _qrE.getSubCond();
    // irankc  = _qrE.getRankc();
    // condc   = _qrE.getSubCondc();
    // sens    = std::fabs( _qrE.getDiag()(irankc+1) );

    ifail = _qrE.getError().getIerr();

//std::cerr << "_irank  = " << _irank << std::endl;
//std::cerr << "cond    = " << cond << std::endl;
//std::cerr << "ifail   = " << ifail << std::endl;
//std::cerr << std::endl;
//std::cerr << "_E = (post)" << std::endl;
//std::cerr << _E.t() << std::endl;
//std::cerr << "_qrE = (post)" << std::endl;
//std::cerr << _qrE.getMat().t() << std::endl;
//std::cerr << "*** MultipleShootingGN::compute_qrdcmp_E ***" << std::endl;
//std::cerr << std::endl;

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Vector
MultipleShootingGN::solve_E(bool qrepet, Vector& b, int& ifail)
{
    Vector x;

//std::cerr << "*** MultipleShootingGN::solve_E ***" << std::endl;
//std::cerr << "qrepet = " << (int)qrepet << std::endl;
//std::cerr << "ifail  = " << ifail << std::endl;
//std::cerr << std::endl;
//std::cerr << "_qrE = (pre)" << std::endl;
//std::cerr << _qrE.getMat().t() << std::endl;
//std::cerr << "b = (pre)" << std::endl;
//std::cerr << b.t() << std::endl;
//std::cerr << "x = (pre)" << std::endl;
//std::cerr << x.t() << std::endl;
//std::cerr << std::endl;

    if ( !qrepet )
    {
        _qrE.solve(b,x);
    }
    else
    {
        _qrE.solveR(b,x);
    }

//std::cerr << "_qrE = (post)" << std::endl;
//std::cerr << _qrE.getMat().t() << std::endl;
//std::cerr << "b = (post)" << std::endl;
//std::cerr << b.t() << std::endl;
//std::cerr << "x = (post)" << std::endl;
//std::cerr << x.t() << "   (" << length_squared(x) << ")" << std::endl;
//std::cerr << "*** MultipleShootingGN::solve_E ***" << std::endl;
//std::cerr << std::endl;

    ifail = _qrE.getError().getIerr();

    return x;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
MultipleShootingGN::log_iteration_vals1(Real sumx)
{
    int monprio  = _lumon.level().priority;
    int verbprio = dlib::LVERB.priority;
    //

    if ( monprio != verbprio )
    {
        printl( _lumon, dlib::LALL,
                " %s\n",
                "   It   Ny     Normf          Normx         Damp.Fct.   New   Rank"
              );
    }
    else
    {
        printl( _lumon, dlib::LGABBY,
                " %s\n",
                "   It   Ny     Normf          Normx                     New   Rank"
              );
    }

    //

    if ( _niter == 0 )
    {
        printl( _lumon, dlib::LVERB,
                "  %4d  %3d   %14.7e   %14.7e               %2d   %4d\n",
                _niter, _ny, _sumf, sumx, _new, _irank
              );
    }
    else
    {
        printl( _lumon, dlib::LVERB,
                "  %4d  %3d   %14.7e   %14.7e               %2d   %4d\n",
                _niter, _ny, _sumf, sumx, _new, _irank
              );
    }

    //

    if ( (monprio != verbprio) && (_niter != 0) )
    {
        printl( _lumon, dlib::LALL,
                "  %4d  %3d   %14.7e   %14.7e    %7.5f    %2d   %4d\n",
                _niter, _ny, _sumf, sumx, _fc, _new, _irank
              );
    }

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
MultipleShootingGN::log_iteration_vals2()
{
//    printl( _lumon, dlib::LVERB,
//            " %s\n",
//            "   It   Ny     Normf          Normx         Damp.Fct."

    printl( _lumon, dlib::LVERB,
            "  %4d  %3d   %14.7e   %14.7e    %7.5f\n",
            _niter, _ny, _sumf, _sumx, _fc
          );

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
MultipleShootingGN::exchange(unsigned& j, unsigned& k, unsigned& ell,
                             unsigned& m, Matrix& A)
{
    unsigned    n = A.nr();
    double      tmp;

    for (unsigned i = 1; i <= ell; ++i)
    {
        tmp = A(i,j);   A(i,j) = A(i,m);   A(i,m) = tmp;
    }

    for (unsigned i = k; i <= n; ++i)
    {
        tmp = A(j,i);   A(j,i) = A(m,i);   A(m,i) = tmp;
    }

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
/*
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALANCE,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
*/
//---------------------------------------------------------------------------
void
MultipleShootingGN::balance(/* unsigned n = A.nr(), */
                            Matrix& A,
                            unsigned& low, unsigned& high, Vector& dscal)
{
/*
C
C     THIS SUBROUTINE BALANCES A REAL MATRIX AND ISOLATES
C     EIGENVALUES WHENEVER POSSIBLE.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        A CONTAINS THE INPUT MATRIX TO BE BALANCED.
C
C     ON OUTPUT:
C
C        A CONTAINS THE BALANCED MATRIX;
C
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT A(I,J)
C          IS EQUAL TO ZERO IF
C           (1) I IS GREATER THAN J AND
C           (2) J=1,...,LOW-1 OR I=IGH+1,...,N;
C
C        SCALE CONTAINS INFORMATION DETERMINING THE
C           PERMUTATIONS AND SCALING FACTORS USED.
C
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
C                 = D(J,J),      J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
C     THEN 1 TO LOW-1.
C
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.
C
C     THE ALGOL PROCEDURE EXC CONTAINED IN BALANCE APPEARS IN
C     PRBALA  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
C     K,L HAVE BEEN REVERSED.)
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     :::::::::: RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION.
C                RADIX = 16.0D0 FOR LONG FORM ARITHMETIC
C
*/
    unsigned    ell, k;
    double      b2, rdx = std::numeric_limits<double>::radix;
    bool        qconv, qswap = true;


    b2 = rdx*rdx;

    k = 1;
    ell = A.nr();


    while ( qswap )
    {

        for (unsigned jj = 1; jj <= ell; ++jj)
        {
            unsigned j = ell + 1 - jj;

            qswap = true;

            for (unsigned i = 1; i <= ell; ++i)
            {
                if ( i == j ) continue;

                if ( A(j,i) != 0.0 ) { qswap = false; break; }
            }

            if ( qswap )
            {
                // m = ell;

                dscal(ell) = j;

                if ( j != ell ) { exchange(j,k,ell,ell,A); }

                if ( ell == 1 )
                {
                    low = k;
                    high = ell;

                    return;
                }

                --ell;

                break;
            }
        }

    }

    while ( qswap )
    {

        for (unsigned j = k; j <= ell; ++j)
        {
            qswap = true;

            for (unsigned i = k; i <= ell; ++i)
            {
                if ( i == j ) continue;

                if ( A(i,j) != 0.0 ) { qswap = false; break; }
            }

            if ( qswap )
            {
                // m = k;

                dscal(k) = j;

                if ( j != k ) { exchange(j,k,ell,k,A); }

                ++k;

                break;
            }
        }

    }

    for (unsigned i = k; i <= ell; ++i)
    {
        dscal(i) = 1.0;
    }

    qconv = false;

    while ( !qconv )
    {
        double c, r, s, g, f;

        qconv = true;

        for (unsigned i = k; i <= ell; ++i)
        {
            c = r = 0.0;

            for (unsigned j = k; j <= ell; ++j)
            {
                if ( j == i ) continue;

                c += std::fabs( A(j,i) );
                r += std::fabs( A(i,j) );
            }

            ///
            /// guard against zero c or r due to underflow
            ///
            if ( (c == 0.0) || (r == 0.0) ) continue;

            g = r / rdx;
            f = 1.0;
            s = c + r;

            while ( !(c >= g) )
            {
                f *= rdx;
                c *= b2;
            }

            g = r * rdx;

            while ( !(c < g) )
            {
                f /= rdx;
                c /= b2;
            }

            ///
            /// now balance
            ///
            if ( (c + r)/f >= 0.95*s ) continue;

            g = 1.0 / f;
            dscal(i) *= f;
            qconv = false;

            for (unsigned j = k; j <= (unsigned)A.nr(); ++j)
            {
                A(i,j) *= g;
            }

            for (unsigned j = 1; j <= ell; ++j)
            {
                A(j,i) *= f;
            }
        }
    }


    low = k;
    high = ell;


    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
/*
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ORTHES,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
*/
//---------------------------------------------------------------------------
void
MultipleShootingGN::orthes(unsigned low, unsigned high,
                           Matrix& A, Vector& d)
{
/*
C
C     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  PRBALA.  IF  PRBALA  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N;
C
C        A CONTAINS THE INPUT MATRIX.
C
C     ON OUTPUT:
C
C        A CONTAINS THE HESSENBERG MATRIX.  INFORMATION ABOUT
C          THE ORTHOGONAL TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLE UNDER THE
C          HESSENBERG MATRIX;
C
C        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
*/
    unsigned mp, n = A.nr();
    unsigned la = high - 1;
    unsigned kp1 = low + 1;

    if ( la < kp1 ) { return; }


    for (unsigned m = kp1; m <= la; ++m)
    {
        double  f, g, h = 0.0;
        double  scal = 0.0;

        d(m) = 0.0;

        ///
        /// scale column (algol (?) tol then not needed)
        ///
        for (unsigned i = m; i <= high; ++i)
        {
            scal += std::fabs( A(i,m-1) );
        }

        if (scal == 0.0 ) continue;

        mp = m + high;

        for (unsigned ii = m; ii <= high; ++ii)
        {
            unsigned i = mp - ii;

            d(i) = A(i,m-1) / scal;
            h += d(i) * d(i);
        }

        g = - ((d(m) >= 0.0) ? 1.0 : -1.0) * std::sqrt(h);

        h -= d(m) * g;

        d(m) -= g;

        ///
        /// form (I - (U*Ut)/h) * A
        ///
        for (unsigned j = m; j <= n; ++j)
        {
            f = 0.0;

            for (unsigned ii = m; ii <= high; ++ii)
            {
                unsigned i = mp - ii;

                f += d(i) * A(i,j);
            }

            f /= h;

            for (unsigned i = m; i <= high; ++i)
            {
                A(i,j) -= f * d(i);
            }
        }

        ///
        /// form (I - (U*Ut)/h) * A * (I - (U*Ut)/h)
        ///
        for (unsigned i = 1; i <= high; ++i)
        {
            f = 0.0;

            for (unsigned jj = m; jj <= high; ++jj)
            {
                unsigned j = mp - jj;

                f += d(j) * A(i,j);
            }

            f /= h;

            for (unsigned j = m; j <= high; ++j)
            {
                A(i,j) -= f * d(j);
            }
        }

        d(m) *= scal;

        A(m,m-1) = scal * g;
    }

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
/*
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR,
C     NUM. MATH. 14, 219-231(1970) BY MARTIN, PETERS, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 359-371(1971).
C
*/
//---------------------------------------------------------------------------
void
MultipleShootingGN::hqr(unsigned low, unsigned high, Matrix& H,
                        Vector& wr, Vector& wi, int& ifail)
{
/*
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A REAL
C     UPPER HESSENBERG MATRIX BY THE QR METHOD.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  PRBALA.  IF  PRBALA  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N;
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX.  INFORMATION ABOUT
C          THE TRANSFORMATIONS USED IN THE REDUCTION TO HESSENBERG
C          FORM BY  ELMHES  OR  PRORTH, IF PERFORMED, IS STORED
C          IN THE REMAINING TRIANGLE UNDER THE HESSENBERG MATRIX.
C
C     ON OUTPUT:
C
C        H HAS BEEN DESTROYED.  THEREFORE, IT MUST BE SAVED
C          BEFORE CALLING  PRHQR  IF SUBSEQUENT CALCULATION AND
C          BACK TRANSFORMATION OF EIGENVECTORS IS TO BE PERFORMED;
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N;
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
C
*/
    unsigned    en, mp2, m, ell, k = 1;
    unsigned    na, enm2, its, n = H.nr();
    double      machep = std::numeric_limits<double>::epsilon();
    double      p, q, r, s, t, w, x, y, zz, norm = 0.0;
    bool        qlast, qnext = true;

    ifail = 0;

    for (unsigned i = 1; i <= n; ++i)
    {
        for (unsigned j = k; j <= n; ++j)
        {
            norm += std::fabs( H(i,j) );
        }

        k = i;

        if ( (low <= i) && (i <= high) ) continue;

        wr(i) = H(i,i);
        wi(i) = 0.0;
    }

    en = high;
    t = 0.0;

    while ( qnext )
    {
        ///
        /// search for next eigenvalues
        ///
        if ( en < low ) { /* exit */ return; }

        its = 0;
        na = en - 1;
        enm2 = na - 1;

        while ( true )
        {
            for (unsigned ellell = low; ellell <= en; ++ellell)
            {
                ell = en + low - ellell;

                if ( ell == low ) { break; }

                s = std::fabs( H(ell-1,ell-1) ) + std::fabs( H(ell,ell) );

                if ( s == 0.0 ) { s = norm; }
                if ( std::fabs( H(ell,ell-1) ) <= machep * s ) { break; }
            }

            ///
            /// form shift
            ///
            x = H(en,en);

            if ( ell == en )
            {
                /// one root found
                wr(en) = x + t;
                wi(en) = 0.0;

                en = na;

                break;
            }

            y = H(na,na);
            w = H(en,na) * H(na,en);

            if ( ell == na )
            {
                /// two roots found
                p = 0.5*(y - x);
                q = p*p + w;
                zz = std::sqrt(std::fabs( q ));
                x += t;
                if ( q < 0.0 )
                {
                    /// complex pair
                    wr(na) = wr(en) = x + p;
                    wi(na) = zz;
                    wi(en) = -zz;
                }
                else
                {
                    /// real pair
                    zz = p + std::fabs(zz)*((p >= 0.0) ? 1.0 : -1.0);
                    wr(en) = wr(na) = x + zz;
                    if ( zz != 0.0 ) wr(en) = x - w / zz;
                    wi(en) = wi(na) = 0.0;
                }

                en = enm2;

                break;
            }

            if ( its >= 30 )
            {
                /// error: no convergence within 30 steps
                ifail = en;

                return;
            }

            if ( !((its != 10) && (its != 20)) )
            {
                ///
                /// form exceptional shift
                ///
                t += x;

                for (unsigned i = low; i <=en; ++i)
                {
                    H(i,i) -= x;
                }

                s = std::fabs( H(en,na) ) + std::fabs( H(na,enm2) );
                y = x = 0.75 * s;
                w = -0.4375 * s * s;
            }

            ++its;

            ///
            /// look for two consecutive small sub-diagonal elements
            ///
            for (unsigned mm = ell; mm <= enm2; ++mm)
            {
                m = enm2 + ell - mm;

                zz = H(m,m);
                r = x - zz;
                s = y - zz;
                p = H(m,m+1) + (r*s - w) / H(m+1,m);
                q = H(m+1,m+1) - zz - r - s;
                r = H(m+2,m+1);
                s = std::fabs(p) + std::fabs(q) + std::fabs(r);
                p /= s;
                q /= s;
                r /= s;
                if ( m == ell ) { break; }

                if ( std::fabs( H(m,m-1) ) * (std::fabs( q ) + std::fabs( r )) <=
                       machep * std::fabs( p ) *
                       (std::fabs( H(m-1,m-1) ) + std::fabs( zz ) + std::fabs( H(m+1,m+1) ))
                   )
                {
                    break;
                }
            }

            mp2 = m + 2;

            for (unsigned i = mp2; i <= en; ++i)
            {
                H(i,i-2) = 0.0;

                if ( i == mp2 ) continue;

                H(i,i-3) = 0.0;
            }

            ///
            /// double qr step involving rows ell to en
            ///                     and columns m to en
            ///
            for (unsigned k = m; k <= na; ++k)
            {
                qlast = !( k != na );

                if ( !(k == m) )
                {
                    p = H(k,k-1);
                    q = H(k+1,k-1);
                    r = ( !qlast ) ? H(k+2,k-1) : 0.0;
                    x = std::fabs( p ) + std::fabs( q ) + std::fabs( r );
                    if ( x == 0.0 ) { break; }
                    p /= x;
                    q /= x;
                    r /= x;
                }

                s = ( (p >= 0.0) ? 1.0 : -1.0 )*std::sqrt(p*p + q*q + r*r);

                if ( !(k == m) )
                {
                    H(k,k-1) = -s * x;
                }
                else
                {
                    if ( ell != m ) { H(k,k-1) = - H(k,k-1); }
                }

                p += s;
                x = p / s;
                y = q / s;
                zz = r / s;
                q /= p;
                r /= p;

                ///
                /// row modification
                ///
                for (unsigned j = k; j <= en; ++j)
                {
                    p = H(k,j) + q * H(k+1,j);

                    if ( !qlast )
                    {
                        p += r * H(k+2,j);
                    }

                    H(k+2,j) -= p * zz;
                    H(k+1,j) -= p * y;
                    H(k,j) -= p * x;
                }

                unsigned jj = ( en <= (k+3) ) ? en : k+3;

                ///
                /// column modification
                ///
                for (unsigned i = ell; i <= jj; ++i)
                {
                    p = x * H(i,k) + y * H(i,k+1);

                    if ( !qlast )
                    {
                        p += zz * H(i,k+2);
                    }

                    H(i,k+2) -= p * r;
                    H(i,k+1) -= p * q;
                    H(i,k) -= p;
                }

            }

        }

    }

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
