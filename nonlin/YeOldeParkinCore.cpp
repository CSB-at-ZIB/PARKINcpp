// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-03-23 td
// last changed:
//

#include "YeOldeParkinCore.h"

using namespace PARKIN;

///

dlib::log_level YeOldeParkinCore::_loglvl[] = {
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
YeOldeParkinCore::YeOldeParkinCore() :
    _luerr("err.PK"), _lumon("mon.PK"), _lusol("sol.PK"), _lutim("tim.PK"),
    // _iopt(IOpt()), _wk(YeOldeParkinWk()),
    _fun(0),
    _x()
{
    IOpt            iopt = IOpt();
    YeOldeParkinWk  wk = YeOldeParkinWk();

    dlib::set_all_logging_levels(dlib::LNONE);
    // _lusol.set_level(dlib::LNONE);
    _luerr.set_logger_header(&print_parkin_logger_header);
    _lumon.set_logger_header(&print_parkin_logger_header);
    _lusol.set_logger_header(&print_parkin_logger_header);

    setIOpt(iopt);
    setWk(wk);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
YeOldeParkinCore::~YeOldeParkinCore()
{
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
YeOldeParkinCore::setIOpt(IOpt const& iopt)
{
    _luerr.set_level(_loglvl[iopt.mprerr]);
    _lumon.set_level(_loglvl[iopt.mprmon]);
    _lusol.set_level(_loglvl[iopt.mprsol]);

    _iopt = iopt;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
YeOldeParkinCore::setWk(YeOldeParkinWk const& wk)
{
    _wk = wk;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
YeOldeParkinCore::setProblem(UserFunc* fun)
{
    _fun = fun;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Vector
YeOldeParkinCore::getSolution()
{
    return _x;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
YeOldeParkinWk
YeOldeParkinCore::getWk()
{
    return _wk;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
YeOldeParkinCore::printCounter()
{
    printl( _lumon, dlib::LINFO,
            "\n\n\n  %s %6d\n  %s %6d\n  %s %6d\n  %s %6d\n  %s %6d\n\n\n",
            "GN-ITERATIONS:", _iter,
            "GN-F-EVAL.   :", _ifgn,
            "GN-J-EVAL.   :", _ijgn,
            "DECOMP CALLS :", _ndecom,
            "SOLVE CALLS  :", _nsolve
          );
	// ... and a few more counting variables!
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
int
YeOldeParkinCore::initialise(
    unsigned                m,
    Vector const&           x,
    Vector const&           xscal,
    Vector const&           fobs,
    Vector const&           fscal,
    Real const              rtol,
    IOpt const&             iopt,
    YeOldeParkinWk const&   wk
    )
{
    _n = x.nr();
    _m = fobs.nr();

    _xa = _x = x;

    if ( (unsigned) xscal.nr() == _n ) _xw = xscal;
    else _xw.ones(_n);

    _dx.zeros(_n);
    _dxh.zeros(_n);
    _dx1.zeros(_n);
    _dx1a.zeros(_n);
    _v.zeros(_n);
    _pivot.zeros(_n);
    _diag.zeros(_n);

    _z = fobs;

    if ( (unsigned) fscal.nr() == _m ) _zscal = fscal;
    else _zscal.ones(_m);

    _f.zeros(_m);
    _fd.zeros(_m);
    _fh.zeros(_m);
    _u.zeros(_m);

    _epsu = (wk.eps == 0.0) ? rtol : wk.eps;

    setIOpt(iopt);
    setWk(wk);

    return 0;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
int
YeOldeParkinCore::run()
{
    int done;

	// C
	// C  INITIALIZE TIME-COUNTING AND COUNT PARAMETERS
	// C-----------------------------------------------

	initialise_timer_and_counter();

	// C------------
	// C  INITIATION
	// C------------

	initiation();

	// C-------------------------
	// C  SET INTERNAL PARAMETERS
	// C-------------------------

	set_internal_parameters();

	// C---------------------------------------------------
	// C  SCALING OF MODEL FUNCTION VALUES AND MEASUREMENTS
	// C---------------------------------------------------

	scale_model_and_measurements();

	// C
	// C  IDENTIFICATION OF PARAMETERS OR NOT
	// C-------------------------------------
	// C
	//      IF(IITER(3).NE.0) GOTO 8000

	// C-----------------------------------------
	// C  PREPARATIONS TO START GN-METHOD
	// C-----------------------------------------

	prepare_GN_method();

	// C-----------------
	// C  INITIAL SCALING
	// C-----------------

	initial_scaling();

	// C
	// C*****************************************************************
	// C
	// C------------------------------
	// C  START OF MAIN ITERATION LOOP
	// C------------------------------

L1:
	// C
	// C  COMPUTATION OF RESIDUAL VECTOR
	// C

	while ( (done = compute_residual_vector()) == 1 ) ;

	switch ( done )
	{
		case  20: break;
		case  44: goto L44;
		case  93: goto L93;     /// fc == fcmin
		default : return -999;
	}


	_sumfa = _sumf;
	_itol  = 0;

L2:
	// time_init();
	// iscal = 2; itype = 1; infoco = 0;
    _iscal = 2;

	// C-----------------------------------------------
	// C  APPROXIMATION OF THE NEGATIVE JACOBIAN MATRIX
	// C-----------------------------------------------
	// C  (NUMERICAL INTEGRATION OF VARIATIONAL EQUATION)

	while ( approximate_negative_jacobian() == 25 ) ;

	// C-------------------------------
	// C  SOLUTION OF THE LINEAR SYSTEM
	// C-------------------------------

	switch ( solve_linear_system() )
	{
		case   0: break;
		case   8: goto L8;      /// konv > 1
		case  93: goto L93;     /// rank == 0
		default : break;
	}

	// C--------------------
	// C  STEPLENGTH CONTROL
	// C--------------------

	control_steplength();

L61:
	// C-----------------------------
	// C  TRIAL VALUE OF NEXT ITERATE
	// C-----------------------------

	switch ( compute_next_trial_step() )
	{
		case   1: goto L1;
		case  92: goto L92;     /// iter > itmax
		case 621: goto L621;    /// fc < fcmin
		default : return -999;
	}

	// C------------------------------------------------------------------

L44:
	// C  (BEST) LINEAR LEAST SQUARES SOLUTION (SIMPLIFIED GN-CORRECTION)
	// C-----------------------------------------------------------------

	switch ( prepare_for_next_iteration_step() )
	{
		case   1: goto L1;
		case   2: goto L2;      /// (fc != 1) || (tol <= s)
		case   6: break;
		case   8: goto L8;      /// konv > 1
		case  91: goto L91;     /// div > divmax
		case  94: goto L94;     /// reduct > reductmax
		default : return -999;
	}

	// C  REDUCTION OF RELAXATION FACTOR
	// C--------------------------------

	switch( reduce_damping_factor() )
	{
		case  61: goto L61;
		case  93: goto L93;
		case 621: goto L621;
		default : return -999;
	}

	// C------------------------------------------------------------------

L621:
	optimise_rank_sub621();

	do {
        predict_damping_factor();
	}
	while ( optimise_rank() == 51 );

	control_steplength_output();

	goto L61;

	// C
	// C******************************************************************
	// C

L91:
    printl( _luerr, dlib::LINFO,
            "\n\n %s\n %s\n %s\n\n\n",
            "TERMINATION SINCE ITERATION DIVERGES,",
            "INITIAL GUESS TOO BAD OR MODEL TOO FAR",
            "FROM BEING COMPATIBLE"
          );
    exit_solution_entry2();
    return -1;

L92:
    --_iter;
    printl( _luerr, dlib::LINFO,
            "\n\n %s\n %3d  %s\n\n\n",
            "USER-PRESCRIBED TERMINATION AFTER",
            _iter, "ITERATIONS"
          );
    exit_solution_entry2();
    return -2;

L93:
    printl( _luerr, dlib::LINFO,
           "\n\n %s\n %s\n\n\n",
           "TERMINATION SINCE RELAXATION STRATEGY",
           "DID NOT SUCCEED"
          );
    exit_solution_entry2();
    return -3;

L94:
    --_ic;
    printl( _luerr, dlib::LINFO,
            "\n\n %s\n %5d %s\n\n\n",
            "TERMINATION SINCE MORE THAN",
            _ic, "RANK REDUCTIONS PERFORMED"
          );
	exit_solution_entry2();
	return -4;


L8:
	exit_solution_entry1();

	return 0;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
int
YeOldeParkinCore::analyse()
{
    return 0;
}
//---------------------------------------------------------------------------

//===========================================================================

Vector
YeOldeParkinCore::call_FCN(Vector const& x, int& ifail)
{
    if ( _lpos == true )
    {
        long   n = x.nr();
        Vector y(n);

        for(long k = 1; k <= n; ++k) y(k) = std::exp( x(k) );

        return _fun -> fcn( y, ifail );
    }

    return _fun -> fcn( x, ifail );
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Matrix
YeOldeParkinCore::call_JAC(Vector const& x, int& ifail)
{

    if ( _lpos == true )
    {
        long   n = x.nr();
        Vector y(n);

        for(long k = 1; k <= n; ++k) y(k) = std::exp( x(k) );

        Matrix J = _fun -> jac( y, ifail );

        return J * y.diag();
    }

    if ( _iscal == 2 )
    {
        long   n = x.nr();
        Matrix J = _fun -> jac( x, ifail );

        for (long k = 1; k <= n; ++k) J.set_colm(k) = J.colm(k) * _xw(k);

        return J;
    }

    return _fun -> jac( x, ifail );
}

//===========================================================================

// --------------------------------------------------------------------------

int YeOldeParkinCore::initialise_timer_and_counter()
{
    _ifgn = _ijgn = _ndecom = _nsolve = 0;

	return 0;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::initiation()
{
    _lscal  = true;
    _kt     = false;
    _irkopt = 0;

	return 0;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::set_internal_parameters()
{
    // _epsu   = _wk.eps;

    _fcmin  = _wk.fcmin;
    _fc1    = _wk.fc;
    if ( _fc1 < _fcmin ) _fc1 = 0.1;

    _zscald = _wk.zscal;
    if ( _zscald < 0.0 ) _zscald = 1.0e-8;

    _xstepm = _wk.xstep;
    if ( _xstepm < 1.0) _xstepm = 1.1;

    _cond = _wk.cond;
    if ( _cond < 1.0 ) _cond = 1.0e16;
    if ( _cond > 1e16) _cond = 1.0e16;

    _fcmin2 = _fcmin * _fcmin;

    Real const epsmin = 1.0e-5;
    Real const epsmax = 1.0e-2;
    _tolmin = 1.0e-7;
    _tolmax = 1.0e-3;

    _qkap = 1.0e-3;

    if ( (_epsu < epsmin) || (_epsu > epsmax) ) _epsu = epsmax;

    _epmach = EPMACH;
    _tol    = _epsu * 1.0e-1;
    _itmax  = _wk.itmax;
    _divm   = 3;

    _lpos = _iopt.lpos;

	// for (int k = 1; k <= n; ++k) { rk(ipara(k) = _x(k); }

	if ( _iopt.mode != 0 ) return 0; // if ( iiter(3) != 0 ) return 0;
	// if ( kprint >= 1 ) cout << "... headline ..." << endl;
	// if ( kprint < 3 ) return 0;

	printl( _lumon, dlib::LINFO,
            "\n\n              %s\n              %s\n              %s\n\n",
            "****************************************",
            "*****    GAUSS NEWTON ITERATION    *****",
            "****************************************"
          );

    printl( _lumon, dlib::LVERB,
            "\n %s\n %s\n\n\n %s %8.4f\n %s %8.4f\n\n %s %3d\n\n %s %10.2e\n\n %s %10.2e\n",
            "INTERNAL PARAMETER FOR GN-METHOD:",
            "=================================",
            "RELAXATION FACTOR:     START:", _fc1,
            "                     MINIMUM:", _fcmin,
            "MAXIMUM NUMBER OF ITERATIONS:", _itmax,
            "INTEGRATOR TOLERANCE:", _tol,
            "PRESCRIBED RELATIVE ACCURACY FOR SOLUTION:", _epsu
          );

    if ( _qkap > 0.0 )
    {
      printl( _lumon, dlib::LVERB,
              " %s\n",
              "(INTERNALLY ADAPTED TO CONVERGENCY RATE)"
            );
    }

	return 0;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::scale_model_and_measurements()
{
    Real zsmin, zmax = _z(1);

	for (unsigned j = 1; j <= _m; ++j)
	{
	    if ( _z(j) > zmax ) zmax = _z(j);
	}
	zsmin = zmax * _zscald;

	for (unsigned j = 1; j <= _m; ++j)
	{
	    _zscal(j) = _z(j);
	    if ( _zscal(j) < zsmin ) _zscal(j) = zsmin;
	}

	// if ( iiter(8) == 1 ) user_scaling();

	// if ( kprint <= 2 ) return 0;
	if ( _iopt.mode != 0 ) return 0; // if ( iiter(3) != 0 ) return 0;

	// if ( iiter(8) == 0 ) cout << " ... " << zscal << endl;
	// if ( iiter(8) == 1 ) cout << " --- " << endl;

    printl( _lumon, dlib::LVERB,
            "\n %s\n\n",
            "INTERNAL SCALING FOR MEASUREMENTS:"
          );

    char line[256]; *line='\0';
    for (unsigned j = 1; j <= _m; ++j)
    {
        std::sprintf( line, "%s %7d %13.6e", line, j, _zscal(j) );

        if ( j%3 == 0 )
        {
            printl( _lumon, dlib::LVERB, " %s\n", line );
            *line='\0';
        }
    }
    if ( *line != '\0' ) printl( _lumon, dlib::LVERB, " %s\n", line );

	return 0;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::prepare_GN_method()
{
    _irank  = _n;
    _irkmax = _n;
    _icmax  = _itmax;
    _skap   = 0.0;
    _epmach *= 10.0 * _n;

    _fc = _fc1;

    _fca    = _fc;
    _ic     = 0;
    _iranka = _irank;
    _iter   = 0;
    _div    = 0;
    _sfc1   = 0.0;
    _tmin   = 1.0;
    _konv   = 0;
    // _tsave  = _tp(...);

    _tolf  = _tol;
    _tolj  = _tol;
    _itol  = 1;
    _icall = 0;
    _iscal = 0;

	return 0;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::initial_scaling()
{
	for (unsigned k = 1; k <= _n; ++k)
	{
	    Real s = _x(k);
	    if ( _lpos == true ) s = (s > 0.0) ? std::log(s) : -1.0e38;
	    _x(k) = s;
	    if ( s < _epmach ) s = _epmach;
	    if ( _lpos == true ) s = 1.0;
	    _xw(k) = s;
	}

	// if ( kprint < 1 ) return 0;

	if ( _lpos == true )
	{
	    printl( _lumon, dlib::LINFO,
                "\n\n       %s\n       %s\n       %s\n\n\n",
                "------------------------------------",
                "|  ITERATION IN LOGARITHMIC SCALE  |",
                "------------------------------------"
              );
	}

	printl( _lumon, dlib::LINFO,
            "\n\n %s\n\n",
            "INITIAL GUESS OF PARAMETERS:"
          );
    for (unsigned k = 1; k <= _n; ++k)
    {
        printl( _lumon, dlib::LINFO,
                " %8d: %15.6e\n", k, _x(k)
              );
    }
    printl( _lumon, dlib::LINFO, "\n" );

	return 0;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::compute_residual_vector()
{
    int ifail = 0;

	// fjev11();
    _fd = call_FCN( _x, ifail );

    ++_ifgn;

	// if ( kflag < 0 )
	if ( ifail != 0 )
	{
		// if ( kprint >= 1 ) cout << "... message ..." << endl;
		printl( _lumon, dlib::LVERB,
                " %s\n",
                "INTEGRATION FAILED; RESTART WITH FC=FC/4.0"
              );

		if ( _fc == _fcmin ) return 93;

		_fc *= 0.25;

		if ( _fc < _fcmin ) _fc = _fcmin;

		for (unsigned k = 1; k <= _n; ++k)
		{
		    _x(k) = _xa(k) + _fc * _dx(k);
		}

		// if ( nld != 0 )
		// {
		//	    for (int k = 1; k <= n; ++k)
		//	    {
		//		    if ( ld2 <= n ) continue;
		//		    for (int l = ld1; l <= ld2; ++l) { }
		//	    }
		// }

		// for (int k = 1; k <= n; ++k) { }

        _iscal = 0;

		return 1;
	}

    _sumf = 0.0;

	// if ( (_ifgn == 1) && (kprint >= 3) )
	if ( _ifgn == 1 )
	{
	    printl( _lumon, dlib::LTALK,
                "\n %s\n\n               %s        %s      %s  %s\n",
                "RESIDUALS FOR INITIAL GUESS:",
                "YMODEL", "YMEASURE", "RESIDUAL", "SCALED RESIDUAL"
              );
	}

    for (unsigned j = 1; j <= _m; ++j)
    {
        Real s, sh;
        s = sh = _fd(j) - _z(j);

        if ( _lscal == true ) sh /= _zscal(j);
        if ( _ifgn == 1 )
        {
            printl( _lumon, dlib::LTALK,
                    " %8d: %13.4e %13.4e %13.4e %13.4e\n",
                    j, _fd(j), _z(j), s, sh
                  );
        }
        if ( _lscal == true ) s = sh;

        _f(j) = _u(j) = s;
        _sumf += s*s;
    }

	if ( (_iter == 0) || (_itol != 0) ) return 20;


	// C  FIRST MONOTONICITY TEST  (LEVEL-FUNCTION: T(X|I))
	// C---------------------------------------------------

    _kt = false;    // _kt = 1;

	if ( _sumf <= _sumfa ) _level = 1;

	return 44;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::approximate_negative_jacobian()
{
    int ifail = 0;

	// fjev11();
    _A = call_JAC( _x, ifail );

    ++_ijgn;

	// if ( kflag < 0 )
	if ( ifail != 0 )
	{
	    _iscal = 0;

	    printl( _lumon, dlib::LINFO,
                "\n %s\n %s %d\n\n",
                "EVALUATION OF GN-JACOBIAN FAILED;",
                "RESTART WITH ISCAL =", _iscal
              );

		return 25;
	}

//std::cerr << "\n ***** YeOldeParkinCore::approximate_negative_jacobian() *****\n";
//std::cerr << " _A = " << std::endl;
//std::cerr << _A;
//std::cerr << "\n *************************************************************\n";

	for (unsigned j = 1; j <= _m; ++j)
	{
		for (unsigned k = 1; k <= _n; ++k)
		{
		    _A(j,k) = - _A(j,k);
		}
	}

	if ( _lscal == true )
	{
		for (unsigned j = 1; j <= _m; ++j)
		{
		    Real scalh = _zscal(j);

			for (unsigned k = 1; k <= _n; ++k)
			{
			    _A(j,k) /= scalh;
			}
		}
	}

	return 0;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::solve_linear_system()
{
    _smax   = 0.0;
    _irkdec = 0;
    _level  = 0;
    _kt     = false;    // _kt = 0;

	_u = _f;  // for (unsigned j = 1; j <= _m; ++j) _u(j) = _f(j);


	// deccon();  // Householder triangularisation
    _qrA = _A.factorQR( _irank, _cond, 0 );

//std::cerr << "\n ***** YeOldeParkinCore::solve_linear_system() ***** \n";
//std::cerr << " _qrA.getRank() = " << _qrA.getRank() << std::endl;
//std::cerr << " _qrA.getMatH() = " << std::endl;
//std::cerr << _qrA.getMatH();
//std::cerr << "\n *************************************************** \n";

    ++_ndecom;

    _irank = _qrA.getRank();
    // _cond  = _qrA.getSubCond();
    _pivot = _qrA.getPivot();
    _diag  = _qrA.getDiag();
    _AH    = _qrA.getMatH();

	if ( _irank == 0 ) return 93;

    _d1 = std::fabs( _diag(1) );
    _sk = _d1 / std::fabs( _diag(_irank) );

	// solcon();            // (best) linear least squares solution
    _qrA.solve( _u, _v );   // solve QR v = u

    ++_nsolve;

	_fh = _u; // for (unsigned j = 1; j <= _m; ++j) _fh(j) = _u(j);

	_sumxa = _tmin  = 0.0;

	if ( _iter != 0 ) { _sfc1 = _sfc2 = 0.0; }

	for (unsigned k = 1; k <= _n; ++k)
	{
	    Real t = _v(k);
	    Real h = _xw(k);
	    Real s = t*h;
	    _sumxa += t*t;
	    _tmin = std::max( _tmin, std::fabs(t) );
		if ( _iter != 0 )
		{
		    t     -= _dx1a(k)/h;
		    _sfc2 += t*t;
		    t     =  _dx(k)/h;
		    _sfc1 += t*t;
		}
		_xa(k)  = _x(k);
		_dxh(k) = s;
	}

    //
	// Test of Accuracy
	// - naturally scaled L-infinity norm -
	//
	if ( _tmin <= _epsu ) ++_konv;
	if ( _iter > 0 ) _skap = std::sqrt( _sumxa / _sfc1 );
	if ( _konv > 1 ) return 8;

	if ( _iter != 0 )
	{
	    printl( _lumon, dlib::LVERB,
                "\n\n  %s\n",
                "RANK OPTIMISATION:"
              );
	}

	return 0;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::control_steplength()
{
	if ( _iter > 0 )
	{
		do {
			predict_damping_factor();
		} while ( optimise_rank() == 51 );

		return control_steplength_output();
	}

    //

    if ( _kt == true ) _fc = _fcmin;  // if ( _kt < 0 ) ...

    _dx = _dxh; // for (unsigned k = 1; k <= _n; ++k) _dx(k) = _dxh(k);

	return control_steplength_output();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int YeOldeParkinCore::control_steplength_output()
{
    if ( _irank != _n )
    {
        char line[4096]; *line='\0';

        std::sprintf( line, "SELECTED COLUMNS: ");

        for (unsigned k = 1; k <= _irank; ++k)
        {
            std::sprintf( line, "%s %4d", line, (unsigned)_pivot(k));
        }
        printl( _lumon, dlib::LVERB, "    %s\n", line );
    }

    printl( _lumon, dlib::LVERB,
            "\n\n%s\n\n\n",
            "*******************************************************************"
          );

    printl( _lumon, dlib::LINFO,
            "\n\n\n  %s   %s    %s  %s   %s   %s  %s   %s\n\n",
            "IT", "LEVELX", "RELFC", "RANK", "LEVELF",
            "SUBCOND(J)", "SENS.", "KAPPA"
          );
    printl( _lumon, dlib::LINFO,
            " %3d %10.4e %7s %3d  %10.4e  %8.2e %8.2e %8.2e\n",
            _iter, _sumxa, " ", _irank, _sumf,
            _sk, _d1, _skap
          );

    ++_iter;

	return 0;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::predict_damping_factor()
{
    Real const fcs1 = 0.7;

	if ( _irank != _n )
	{
		for (unsigned l = 1; l <= _n; ++l)
		{
		    unsigned k = _pivot(l);
		    _v(l) = _dx1a(k) / _xw(k);
		}

        _AH = _qrA.getMatH();

//std::cerr << "\n ***** YeOldeParkinCore::predict_damping_factor() ***** \n";
//std::cerr << " _AH =  (" << _AH.nr() << " x " << _AH.nc() << ")" << std::endl;
//std::cerr << _AH;
//std::cerr << "\n ****************************************************** \n";

        Real     t = 0.0;
        unsigned rk1 = _irank + 1;

		for (unsigned k = rk1; k <= _n; ++k)
		{
		    Real s = _v(k);

			for (unsigned l = 1; l <= _irank; ++l) s -= _v(l) * _AH(l,k);

			if ( k != rk1 )
			{
			    // unsigned k1 = k - 1;
				for (unsigned l = rk1; l < k; ++l) s -= _v(l) * _AH(l,k);
			}

			s /= _AH(k,k);
			t += s*s;
			_v(k) = s;
		}
		_sfc2 -= t;
	}

    // _fc = _fca / _fcmin;
    // if ( _sfc2 > (_sfc1*_fcmin2) ) _fc = std::sqrt( _sfc1 / _sfc2 ) * _fca;
    _fc = ( _sfc2 > (_sfc1*_fcmin2) ) ?
                std::sqrt( _sfc1 / _sfc2 ) * _fca :
                _fca / _fcmin;

    if ( _fc > fcs1 ) _fc = 1.0;

    _dx = _dxh; // for (unsigned k = 1; k <= _n; ++k) _dx(k) = _dxh(k);

	return 0;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::optimise_rank()
{
    Real t = _fc * _fc * _sumxa;

	if ( _irkdec == 1 ) return 54;

    printl( _lumon, dlib::LVERB,
            " %3s %10.4e  %5.3f  %3d %12s %8.2e -->  %8.2e\n",
            " ", _sumxa, _fc, _irank, " ", _sk, t
          );

	if ( t >= _smax )
	{
	    _irkopt = _irank;

	    if ( t > _smax ) _smax = t;
    }

	return optimise_rank_sub621();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int YeOldeParkinCore::optimise_rank_sub621()
{
    _kt = true;     // _kt = -1;
	_u  = _fh;      // for (int j = 1; j <= _m; ++j) _u(j) = _fh(j);

    --_irank;

	int entry = 0;
	if ( (_irank  > 0) && (_fc < 1.0) ) {
	    entry = 1;
	}
	else {
	    _irank = _irkopt;
	    _irkdec = 1;
	    if ( _irkopt == _n ) entry = 2;
	}

	switch ( entry )
	{
		default:
		case 0 :
		case 1 : _qrA.setNewRank( _irank );      // deccon();
                 ++_ndecom;
                 _diag  = _qrA.getDiag();
                 _pivot = _qrA.getPivot();
                 // _irank = _qrA.getRank();
                 // _cond  = _qrA.getSubCond();
                 _AH    = _qrA.getMatH();

//std::cerr << "\n ***** YeOldeParkinCore::optimise_rank_sub621() ***** \n";
//std::cerr << " entry          = " << entry << std::endl;
//std::cerr << " _irank         = " << _irank << std::endl;
//std::cerr << " _qrA.getRank() = " << _qrA.getRank() << std::endl;
//std::cerr << " _qrA.getMatH() = " << std::endl;
//std::cerr << _qrA.getMatH() << std::endl;
//std::cerr << "\n **************************************************** \n";

		case 2 : _d1 = std::fabs( _diag(1) );
                 _sk = _d1 / std::fabs( _diag(_irank) );
                 _qrA.solveR( _u, _v );                   // solcon();
                 ++_nsolve;
                 break;
	}

    _sumxa = _tmin = 0.0;

    if ( _iter > 0 ) _sfc2 = 0.0;

	for (unsigned k = 1; k <= _n; ++k)
	{
	    Real t = _v(k);
	    Real h = _xw(k);
	    Real s = t*h;

	    _sumxa += t*t;
	    _tmin  = std::max( _tmin, std::fabs(t) );

	    if ( _iter != 0 )
	    {
	        t     -= _dx1a(k) / h;
	        _sfc2 += t*t;
	    }

	    _dxh(k) = s;
	}

	return 51;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::compute_next_trial_step()
{
    int done;

    _x = _xa + _fc * _dx;
	// for (unsigned k = 1; k <= _n; ++k) _x(k) = _xa(k) + _fc*_dx(k);

	while ( (done=compute_next_sub1()) == 0 );

	return done;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int YeOldeParkinCore::compute_next_sub1()
{
	// if ( nld != 0 )
	// {
	//      for (int k = 1; k <= n; ++k)
	//      {
	//		    if ( ld2 <= n ) continue;
	//		    if ( apply_trafo && (fabs(s) > 50.0) )
	//			    return compute_next_sub2(k);
    //
	//		    for (int l = ld1; l <= ld2; ++l) { }
	//	    }
	// }

	for (unsigned k = 1; k <= _n; ++k)
	{
	    Real s = _x(k);

		if ( _lpos == true )
		{
			if ( std::fabs(s) > 50.0 )
				return compute_next_sub2(k);
		}
	}

    _sumfa = _sumf;

    _iscal = (_fc == 1.0) ? 2 : 1;

	return (_iter <= _itmax) ? 1 : 92;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int YeOldeParkinCore::compute_next_sub2(unsigned kk)
{
    if ( _x(kk) > 0.0) _fc = ( 49.0 - _xa(kk)) / _dx(kk);
    if ( _x(kk) < 0.0) _fc = (-49.0 - _xa(kk)) / _dx(kk);

    printl( _lumon, dlib::LINFO,
            "\n    %s\n    %s %6.3f\n\n",
            "TRIAL VALUE OUT OF RANGE;",
            "RELFC REDUCED TO:", _fc
          );

	if ( _fc < _fcmin ) return 621;

    _x = _xa + _fc * _dx;
	// for (unsigned k = 1; k <= _n; ++k) _x(k) = _xa(k) + _fc * _dx(k);

	return 0;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

int YeOldeParkinCore::reduce_damping_factor()
{
    printl( _lumon, dlib::LINFO,
            " %3s %10.4e  %5.3f %5s %10.4e\n",
            " ", _sumx, _fc, " ", _sumf
          );

    Real t = (std::sqrt(_sumx/_sumxa) + _fc - 1.0) / _fc;

    _fc /= (std::sqrt(8.0*t + 1.0) - 1.0);

	if ( _fc >= _fcmin ) return 61;

	if ( _irank <= 1 )
	{
		if ( _fca <= _fcmin ) return 93;

		_fc = _fcmin;
	}

	// C  REDUCTION OF RANK
	// C-------------------

    --_iter;

	return 621;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::prepare_for_next_iteration_step()
{
    if ( _kt == true )      // if ( _kt < 0 )
        _qrA.solveR( _u, _v );
    else
        _qrA.solve( _u, _v );
	// solcon();  // (best) linear least squares solution
    ++_nsolve;

	_fh = _u; // for (unsigned j = 1; j <= _m; ++j) _fh(j) = _u(j);

    _sumx = 0.0;

	for (unsigned k = 1; k <= _n; ++k)
	{
	    Real t = _v(k);
	    Real h = _xw(k);
	    Real s = t*h;

	    _sumx += t*t;
	    _dx1(k) = s;
	}

	// C  SECOND MONOTONICITY TEST
	// C  IN TERMS OF NATURAL LEVEL-FUNCTION
	// C------------------------------------

    if ( _sumx <= _sumxa ) ++_level;

    if ( (_sumx <= (_sumxa*_qkap)) && (_fc == 1.0) ) ++_konv;
	if ( _konv > 1 ) return 8;

	if ( _level <= 0 ) return 6;

	// C----------------------------------------------------
	// C  PREPARATIONS TO START THE FOLLOWING ITERATION STEP
	// C----------------------------------------------------

	// C  TEST ON LOCAL CONTRACTION OF MAPPING

    ++_div;
    if ( (_sumxa < _sfc1) || (_irank != _n) ) _div = 0;
	if ( _div > _divm ) return 91;

    //

    printl( _lumon, dlib::LINFO,
            " %3s %10.4e  %5.3f %5s %10.4e\n",
            " ", _sumx, _fc, " ", _sumf
          );

    printl( _lumon, dlib::LVERB,
            "\n  %s\n\n",
            "PARAMETER ITERATES:"
          );
    for (unsigned k = 1; k <= _n; ++k)
    {
        printl( _lumon, dlib::LVERB,
                " %8d:  %15.6e\n", k, _x(k)
              );
    }

    // if ( (nld > 0) && (kprint >= 3) )
    // {
    //     cout << param << endl;
    //     cout << rk << endl;
    // }

    printl( _lumon, dlib::LGABBY,
            "\n\n  %s\n\n               %s        %s      %s  %s\n",
            "ACTUAL RESIDUAL VECTOR:",
            "YMODEL", "YMEASURE", "RESIDUAL", "SCALED RESIDUAL"
          );

    for (unsigned j = 1; j <= _m; ++j)
    {
        Real s = _fd(j) - _z(j);
        Real sh = s;
        if ( _lscal == true ) sh /= _zscal(j);

        printl( _lumon, dlib::LGABBY,
                " %8d: %13.4e %13.4e %13.4e %13.4e\n",
                j, _fd(j), _z(j), s, sh
              );
    }

    _fca = _fc;

    if ( (_irank < _n) && (_iranka >= _irank) ) ++_ic;
    if ( _ic > _icmax ) return 94;

    _iranka = _irank;

    _irank = std::min(_irkmax, _n);

	// C  RESCALING

	_dx1a = _dx1;

	for (unsigned k = 1; k <= _n; ++k)
	{
	    // _dx1a(k) = _dx1(k);
	    Real t = 0.5*( std::fabs(_x(k)) + std::fabs(_xa(k)) );
		if ( t < std::fabs(_dx(k)) ) continue;
		if ( t < SMALL ) t = 1.0;
		if ( _lpos == true ) t = 1.0;
		_xw(k) = t;
	}

    Real sh, s = std::sqrt( _sumxa / _n );

    _itol = 0;

	if ( (_fc != 1.0) || (_tol <= s) ) return 2;

    if ( s < _tolmin ) s = _tolmin;
    if ( s > _tolmax ) s = _tolmax;
    if ( s <= 1.0e-2) sh = 1.0e-2;
    if ( s <= 1.0e-3) sh = 1.0e-3;
    if ( s <= 1.0e-4) sh = 1.0e-4;
    if ( s <= 1.0e-5) sh = 1.0e-5;
    if ( s <= 1.0e-6) sh = 1.0e-6;
    if ( s <= 1.0e-7) sh = 1.0e-7;
	if ( sh == _tol ) return 2;
	if ( sh != _tol ) return 2;

    _tol = _tolf = _tolj = sh;
    _itol = _iscal = 1;

	// C  WRITE MESSAGE THAT TOLERANCE HAS BEEN CHANGED

	printl( _lumon, dlib::LINFO,
            "\n    %s %10.2e\n\n",
            "INTEGRATOR TOLERANCE REDUCED TO:", _tol
          );

	return 1;
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::exit_solution_entry1()
{
	if ( _iopt.mprmon < 2 ) return exit_sub88();
	if ( _iopt.mprmon < 1 ) return print_final_residual_and_counts();

    printl( _lumon, dlib::LINFO,
            " %3s %10.4e  %5.3f %5s %10.4e\n",
            " ", _sumx, _fc, " ", _sumf
          );

    printl( _lumon, dlib::LINFO,
        "\n\n%s\n%s\n\n\n",
        "*******************************************************************",
        "*******************************************************************"
      );

    printl( _lumon, dlib::LINFO,
           "\n %s %3d  %s\n %4d  %s\n %4d  %s\n\n\n",
           "SOLUTION OBTAINED AFTER", _iter,
           "ITERATIONS BY",
           _ifgn, "FUNCTION EVALUATIONS AND",
           _ijgn, "JACOBIAN EVAL."
          );

    printl( _lumon, dlib::LINFO,
           "\n    %s\n    %s\n\n",
           "SOLUTION:",
           "========="
          );

	return exit_sub82();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int YeOldeParkinCore::exit_solution_entry2()
{
    printl( _lumon, dlib::LINFO,
            "\n    %s\n\n",
            "FINAL PARAMETER ITERATES:"
          );

	return exit_sub82();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int YeOldeParkinCore::exit_sub82()
{
	for (unsigned k = 1; k <= _n; ++k)
	{
	    printl( _lumon, dlib::LINFO,
                " %8d: %15.6e\n", k, _x(k)
              );
	}
	// if (nld > 0) {
    //     cout << param << endl;
    //     cout << rk << endl;
	// }

	printl( _lumon, dlib::LINFO,
            "\n\n  %s   %s %12.3e\n",
            "FINALLY ACHIVED ACCURACY :",
            "  EPS =", _tmin
          );

    char line[256]; *line='\0';

    std::sprintf( line, "%s %10.6f\n",
                 "KAPPA =", _skap
                );
    if ( _sfc1 >  1.0e-35 ) std::sprintf(line, "%s\n", "NOT AVAILABLE");
    if ( _skap <= 0.0 ) std::sprintf(line, "%s\n", "NOT AVAILABLE");

	printl( _lumon, dlib::LINFO,
            "\n  %s   %s\n",
            "INCOMPATIBILITY FACTOR   :", line
          );

	return exit_sub88();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int YeOldeParkinCore::exit_sub88()
{
    printl( _lumon, dlib::LINFO,
            "\n\n    %s\n\n  %s          %s\n\n",
            "FINAL PARAMETER ESTIMATES:",
            "J.TH PARAM", "VALUE"
          );
    for (unsigned k = 1; k <= _n; ++k)
    {
        printl( _lumon, dlib::LINFO,
               "  %6d     %19.10e\n", k, _x(k)
              );
    }
    printl( _lumon, dlib::LINFO,
            "\n\n"
          );

	return print_final_residual_and_counts();
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

int YeOldeParkinCore::print_final_residual_and_counts()
{
	if ( _iopt.mprmon < 3 )
		return compute_final_simulation_and_sensitivity();

    printl( _lumon, dlib::LINFO,
            "\n\n %s\n\n               %s        %s      %s  %s\n",
            "FINAL RESIDUALS:",
            "YMODEL", "YMEASURE", "RESIDUAL", "SCALED RESIDUAL"
          );

	for (unsigned j = 1; j <= _m; ++j)
	{
	    Real s = _fd(j) - _z(j);
	    Real sh = s;

	    if ( _lscal == true ) sh /= _zscal(j);

	    printl( _lumon, dlib::LINFO,
                " %8d: %13.4e %13.4e %13.4e %13.4e\n",
                j, _fd(j), _z(j), s, sh
              );
	}

	// cout << "... final counting statistics ..." << endl;
	// cout << gnt << ftime << xjtime << dtime << stime << otime << endl;

    printl( _lumon, dlib::LINFO,
            "\n\n %s\n %s\n\n %s %6d\n %s %6d\n %s %6d\n %s %6d\n\n",
            "FINAL COUNTS:",
            "=============",
            "GN-F-EVAL.:   ", _ifgn,
            "GN-J-EVAL.:   ", _ijgn,
            "DECOMP CALLS: ", _ndecom,
            "SOLVE CALLS:  ", _nsolve
          );
	// ... and a few more counting variables!

	return compute_final_simulation_and_sensitivity();
}

// --------------------------------------------------------------------------

int YeOldeParkinCore::compute_final_simulation_and_sensitivity()
{
    return 0;

    /*
	if ( iiter(4) == 0 ) return 0;

	// sens00();
	// fjev11();

	if ( kprint <= 0 ) return 0;
	if ( typ != 1 ) return 0;

	// C  FINAL RESIDUALS
	// C-----------------
	for (int j = 1; j <= m; ++j) { }

	// C  SCALE OBSERVATION MATRIX
	// C--------------------------
	for (int k = 1; k <= n; ++k) {
		for (int j = 1; j <= m; ++j) { }
	}

	// C  DECOMPOSITION
	// C---------------
	// deccon();

	if ( rank != n ) return 0;

	//
	//stat();
	//

	return 0;
	*/
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


