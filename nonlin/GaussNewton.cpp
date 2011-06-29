// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-10-13 td
// last changed:
//

#include "GaussNewton.h"

using namespace PARKIN;

///

dlib::log_level GaussNewton::_loglvl[] = {
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
GaussNewton::GaussNewton() :
    _luerr("err.GN"), _lumon("mon.GN"), _lusol("sol.GN"), _lutim("tim.GN"),
    // _iopt(IOpt()), _wk(GaussNewtonWk()),
    _fun(0),
    _n(0), _m(0), _mfit(0), _mcon(0),
    _x(), _xscal(), _fi(),
    _rtol(0.0), _tolmin(0.0),
    _nitmax(50), _nonlin(0), _irank(0), _irankc(0), _ierr(0),
    _maxmn(0),
    _AA(), _A(),
    _dx(), _dxq(), _xa(), _xw(), _dxqa(), _delxq(), _eta(),
    _fmodel(), _f(), _fa(), _fw(), _qu(), _rq(), _rqkeep(),
    _fc(0.0), _fcmin(0.0), _sigma(0.0), _fca(0.0), _cond(0.0),
    _conv(0.0), _sumx(0.0), _sumxs(0.0), _dlevf(0.0),
    _niter(0), _ncorr(0), _nfcn(0), _njac(0), _nfcnj(0), _nrejr1(0), _new(0),
    _qbdamp(false),
    _vcv(), _xl(), _xr(),
    _qrA()
{
    IOpt            iopt = IOpt();
    GaussNewtonWk   wk = GaussNewtonWk();

    dlib::set_all_logging_levels(dlib::LNONE);
    // _lusol.set_level(dlib::LNONE);
    _luerr.set_logger_header(&print_parkin_logger_header);
    _lumon.set_logger_header(&print_parkin_logger_header);
    _lusol.set_logger_header(&print_parkin_logger_header);

    _AA.zeros(1,1);
    _A.zeros(1,1);

    setIOpt(iopt);
    setWk(wk);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
GaussNewton::~GaussNewton()
{
    // nothing to do?!
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::setIOpt(IOpt const& iopt)
{
    _luerr.set_level(_loglvl[iopt.mprerr]);
    _lumon.set_level(_loglvl[iopt.mprmon]);
    _lusol.set_level(_loglvl[iopt.mprsol]);

    _iopt = iopt;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::setWk(GaussNewtonWk const& wk)
{
    _wk = wk;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::setProblem(UserFunc* fun)
{
    _fun = fun;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Vector
GaussNewton::getSolution()
{
    if ( _iopt.lpos == true )
    {
        long   n = _x.nr();
        Vector y(n);

        for (long j = 1; j <= n; ++j) y(j) = std::exp( _x(j) );

        return y;
    }

    return _x;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
std::vector<Vector>
GaussNewton::getSolutionIter()
{
    if ( _iopt.lpos == true )
    {
        std::vector<Vector> yiter( _xiter.size() );

        for (unsigned k = 0; k < _niter; ++k)
        {
            Vector y = _xiter[k];
            for (long j = 1; j <= y.nr(); ++j) y(j) = std::exp( y(j) );
            yiter[k] = y;
        }

        return yiter;
    }

    return _xiter;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
GaussNewtonWk
GaussNewton::getWk()
{
//        Real    fcbnd, ajdel, ajmin, etadif, etaini;
//        Real    dlevf, sumx, prec, skap;
//        Real    fcstart, fcmin, sigma, cond;
//        long    niter, nitmax, irank;
//        int     ifail;

    _wk.dlevf  = _dlevf;
    _wk.sumx   = _sumx;
    _wk.cond   = _cond;
    _wk.niter  = _niter;
    _wk.nitmax = _nitmax;
    _wk.irank  = _irank;
    _wk.ncorr  = _ncorr;
    _wk.nrejr1 = _nrejr1;
    _wk.njac   = _njac;
    _wk.nfcn   = _nfcn;
    _wk.nfcnj  = _nfcnj;

    if (_sigma2 > 0.0)
    {
        _wk.sigma2 = _sigma2;
        _wk.xl     = _xl;
        _wk.xr     = _xr;
        _wk.vcv    = _vcv;
    }

    return _wk;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::printCounter()
{
    if ( _ierr == 0 )
    {
        printl(_lumon, dlib::LINFO,
               "\n\n         %s%s%s"
               "\n         %s%7i%s"
               "\n         %s%7i%s"
               "\n         %s%7i%s"
               "\n         %s%7i%s"
               "\n         %s%7i%s"
               "\n         %s%7i%s"
               "\n         %s\n\n",
               "******  Counter * ", "Gauss-Newton", " ******",
               "***  Gauss-Newton iter.: ",  _niter, "  ***",
               "***  Corrector steps   : ",  _ncorr, "  ***",
               "***  Rejected rk-1 st. : ", _nrejr1, "  ***",
               "***  Jacobian eval.    : ",   _njac, "  ***",
               "***  Function eval.    : ",   _nfcn, "  ***",
               "***  ...  for Jacobian : ",  _nfcnj, "  ***",
               "*************************************"
        );
    }
    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
int
GaussNewton::initialise(
    unsigned                m,
    Vector const&           x,
    Vector const&           xScal,
    Vector const&           fObs,
    Vector const&           fScal,
    Real const              rtol,
    IOpt const&             iopt,
    GaussNewtonWk const&    wk
    )
{
    bool        qsucc  = wk.qsucc;
    bool        qrank1 = iopt.qrank1;
    bool        qfcstr = (wk.fcstart > 0.0);
    int         jacgen = iopt.jacgen;
    int         rscal = iopt.rscal;
    unsigned    minmn;

    _m = m;

    if ( iopt.lpos == true )
    {
        long n = x.nr();

        _x.zeros(n);
        _xscal.ones(n);
        for (long j = 1; j <= n; ++j) _x(j) = std::log( x(j) );
    }
    else
    {
        _x = x;
        _xscal = xScal;
    }

    _fi = fObs;

    _n = x.nr();
    _mfit = fObs.nr();

    _tolmin = 10.0*_n*EPMACH;
    _rtol = rtol;
    _ierr = 0;

    setIOpt(iopt);
    setWk(wk);

    if ( !qsucc )
        printl( _lumon, dlib::LINFO,
                "\n    %s\n    %s\n    %s\n    %s\n\n",
                "********* PARKIN : Gauß-Newton *********",
                "* Method for the solution of nonlinear *",
                "* least-squares problems               *",
                "****************************************"
              );

    check_init();
    if ( _ierr != 0 ) return _ierr;

    if ( !qsucc )
    {
        _xiter.clear();
        _sumxall.clear();
        _dlevfall.clear();
        _sumxqall.clear();
        _tolall.clear();
        _fcall.clear();
    }

    _mcon = _m - _mfit;
    _maxmn = std::max(_m, _n);
    // minmn = std::min(_m,_n);

    if ( rscal == 0 )  _iopt.rscal = 1;
    if ( jacgen == 0 ) _iopt.jacgen = 2;


    if ( !qsucc )
    {
        printl( _lumon, dlib::LINFO,
                " %s  %d\n %s  %d\n %s  %d\n\n %s %9.2e\n\n",
                "Number of parameters to be estimated        (N): ", _n,
                "Number of data to fit, e.g. observations (MFIT): ", _mfit,
                "Number of equality constraints           (MCON): ", _mcon,
                "Prescribed relative precision            (XTOL): ", _rtol
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
    if ( _iopt.boundeddamp == 0 )
        _qbdamp = (_nonlin == 4);
    else if ( _iopt.boundeddamp == 1 )
        _qbdamp = true;
    else if ( _iopt.boundeddamp == 2 )
        _qbdamp = false;

    if ( _qbdamp && (_wk.fcbnd < 1.0) ) _wk.fcbnd = 10.0;

    if ( qrank1 && (_m > _n) )
    {
        printl( _luerr, dlib::LVERB,
                "\n\n Warning: Broyden steps are disallowed for\n"
                    "          the overdetermined system.\n"
              );

        qrank1 = _iopt.qrank1 = false;
    }

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
        if ( _qbdamp )
            printl( _lumon, dlib::LINFO,
                    " Bounded damping strategy is active;"
                    " factor is %e\n\n", _wk.fcbnd
                  );
        else
            printl(_lumon, dlib::LINFO,
                    " Bounded damping strategy is off.\n\n"
                  );
    }

    // Maximum permitted number of iteration steps
    _nitmax = _wk.nitmax;
    if ( _nitmax <= 0 ) _nitmax = 50;
    _wk.nitmax = _nitmax;

    if ( !qsucc )
        printl(_lumon, dlib::LINFO,
                " %s  %d\n",
                "Maximum permitted number of iteration steps: ", _nitmax
              );

    // Initial damping factor for highly nonlinear problems
    if ( !qfcstr )
        _wk.fcstart = (_nonlin == 4) ? 1.0e-4 : 1.0e-2;

    // Minimal permitted dampnig factor
    if ( _wk.fcmin <= 0.0 )
        _wk.fcmin = (_nonlin == 4) ? 1.0e-8 : 1.0e-2;
    _fcmin = _wk.fcmin;

    // Broyden-update decision parameter SIGMA
    if ( _wk.sigma < 1.0 ) _wk.sigma = 2.0;
    if ( !qrank1 ) _wk.sigma = 10.0 / _fcmin;
    _sigma = _wk.sigma;

    // Starting value of damping factor (FCMIN <= FC <= 1.0)
    if ( (_nonlin <= 2) && !qfcstr )
        _fc  = 1.0;         // for linear or mildly nonlinear problems
    else
        _fc = _wk.fcstart;  // for highly or extremely nonlinear problems

    _wk.fcstart = _fc;

    // Initial rank
    _irank = _wk.irank;
    minmn = std::min(_m,_n);
    if ( (_irank <= 0) || _irank > minmn ) _wk.irank = minmn;

    // Maximum permitted subcondition number of Jacobian A
    _cond = wk.cond;
    if ( _cond < 1.0 ) _cond = 1.0 / EPMACH;
    _wk.cond = _cond;

    if ( !qsucc )
    {
        printl( _lumon, dlib::LVERB,
                " %s\n %s %9.2e\n %s %9.2e\n %s %9.2e\n %s  %d\n %s %9.2e\n\n",
                "Internal parameters:",
                "Starting value for damping factor (FCSTART): ", _wk.fcstart,
                "Minimum allowed damping factor      (FCMIN): ", _fcmin,
                "Rank-1 updated decision parameter   (SIGMA): ", _sigma,
                "Initial Jacobian pseudo-rank        (IRANK): ", _wk.irank,
                "Maximum permitted subcondition       (COND): ", _cond
              );

        _wk.fw.zeros(_m);
        _fw.zeros(_m);
        for (long j = 1; j <= (long)_mfit; ++j)
            if ( ((Real)SMALL <= fScal(j)) && (fScal(j) <= (Real)GREAT) )
            {
                _wk.fw(_mcon+j) = _fw(_mcon+j) = 1.0 / fScal(j);
            }
            else
            {
                _wk.fw(_mcon+j) = _fw(_mcon+j) = 1.0;
                printl( _luerr, dlib::LVERB,
                        " %s fscal(%d) = %10.2e %s\n",
                        "Warning: Bad scaling value", j, fScal(j),
                        "replaced by 1.0"
                      );
            }

        _niter = _ncorr = _nfcn = _njac = _nfcnj = _nrejr1 = _new = 0;
    }

    return 0;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::check_init()
{
    _ierr = 0;

    // Checking dimensional parameters N,M, and MFIT
    if ( (_n <= 0) || (_m <= 0) || (_mfit < 0) || (_mfit > _m) )
    {
        printl( _luerr, dlib::LINFO,
                " Error: %s\n        %s\n        %s%d%s%d%s%d\n",
                "Bad or inconsistent input to dimensional parameters supplied.",
                "Choose N and M positive, and MFIT <= M.",
                "Input is: N = ", _n, " M = ", _m, " MFIT = ", _mfit
              );
        _ierr = 20;
    }

    // Problem type specified by user
    _nonlin = _iopt.nonlin;
    if ( _nonlin == 0 ) _nonlin = 3;
    _iopt.nonlin = _nonlin;

    // Checking of user-precribed RTOL
    std::string action, bound;
    if ( _rtol <= 0.0 )
    {
        printl( _luerr, dlib::LINFO,
                " Error: Nonpositive RTOL supplied.\n"
              );
        _ierr = 21;
    }
    else
    {
        Real tolmax = 1.0e-1;
        action = bound = "";
        if ( _rtol < _tolmin )
        {
            _rtol = _tolmin;
            action = "increased";
            bound = "smallest";
        }
        if ( _rtol > tolmax )
        {
            _rtol = tolmax;
            action = "decreased";
            bound = "largest";
        }
        if ( !action.empty() )
            printl( _luerr, dlib::LVERB,
                    " Warning: User prescribed RTOL %s"
                    " to reasonable %s value RTOL = %11.2e\n",
                    action.c_str(), bound.c_str(), _rtol
                  );
    }

    // Test user prescribed accuracy and scaling on proper values
    if ( _n <= 0 ) return;
    //
    Real defscl = (_nonlin >= 3) ? _rtol : 1.0;
    //
    for (long k = 1; k <= (long)_n; ++k)
    {
        Real xscl;
        xscl = _xscal(k);
        action = bound = "";

        if ( xscl < 0.0 )
        {
            printl( _luerr, dlib::LINFO,
                    " Error: Negative value in xscal(%d) supplied.\n", k
                  );
            _ierr = 22;
        }

        if ( xscl == 0.0 ) xscl = defscl;

        if ( (xscl > 0.0) && (xscl < SMALL) )
        {
            xscl = (Real)SMALL;
            action = "increased";
            bound = "too small";
        }

        if ( xscl > GREAT )
        {
            xscl = (Real)GREAT;
            action = "decreased";
            bound = "too big";
        }

        if ( !action.empty() )
            printl( _luerr, dlib::LVERB,
                    " Warning: xscal(%d) = %9.2e %s, %s to %9.2e\n",
                    k, _xscal(k), bound.c_str(), action.c_str(), xscl
                  );

        _xscal(k) = xscl;
    }

    return;
}
//---------------------------------------------------------------------------



//---------------------------------------------------------------------------
int
GaussNewton::computeSensitivity()
{
    Real        ajdel, ajmin;
    Real        epdiff = std::sqrt(10.0 * EPMACH);
    Real        etadif, etaini, etamax, etamin;
    unsigned    jacgen, iscal;
    bool        qsucc, qscale, qinisc = true;
    int         ifail = 0;

    if ( _fun == 0 )
    {
        _ierr = 999;
        printl( _luerr, dlib::LINFO,
                " Error: %s\n        %s\n",
                "No problem (fcn/jac pair) set so far.",
                "Please use setProblem(&fun) first.");
        return _ierr;
    }

    if ( (_n <= 0) || (_m <= 0) )
    {
        _ierr = 998;
        printl( _luerr, dlib::LINFO,
                " Error: %s\n        %s\n        %s\n",
                "More input is needed:",
                "Vector x where to compute, and scaling values.",
                "Please use initialise(...) first.");
        return _ierr;

    }
    qsucc  = _wk.qsucc;
    qscale = (_iopt.norowscal != true);
    jacgen = _iopt.jacgen;
    iscal  = _iopt.iscal;

    _ierr = 0;

    if ( !qsucc )  // computeSensitivity() only available before run() is called
    {
        ajdel = ajmin = etadif = etaini = etamax = etamin = 0.0;

        // Numerical differentation related initialisation
        if ( jacgen == 2 )
        {
            ajdel = _wk.ajdel;
            if ( ajdel <= SMALL ) ajdel = std::sqrt(10.0 * EPMACH);
            ajmin = _wk.ajmin;
        }
        else if ( jacgen == 3 )
        {
            etadif = _wk.etadif;
            if ( etadif <= SMALL ) etadif = 1.0e-6;
            etaini = _wk.etaini;
            if ( etaini <= SMALL ) etaini = 1.0e-6;
            etamax = std::sqrt(epdiff);
            etamin = epdiff*etamax;
        }

        if ( jacgen == 3 )
        {
            Vector v;
            v.ones(_n); _eta.zeros(_n);
            _eta = etaini * v;
        }

        _AA.zeros(_m,_n);
        _A.zeros(_m,_n);

        _xa = _x;

        compute_scaling_xw(iscal, qinisc);

        //
        if ( jacgen == 1 )
        {
            _fmodel.zeros(_m);

            //_timon(2);
            _AA = call_JAC( _x, ifail );    // _fun->jac(_x, ifail);
            //_timoff(2);
        }
        else
        {
            //_timon(1);
            _fmodel = call_FCN( _x, ifail );
            //_timoff(1);

            if (ifail != 0) { _ierr = 182; return _ierr; }

            _lumon.set_level( dlib::LDEBUG );

            //_timon(2);
            if ( jacgen == 3 ) compute_jcf_AA(etamin, etamax, etadif, ifail);
            if ( jacgen == 2 ) compute_jac_AA(ajdel, ajmin, ifail);
            //_timoff(2);

            _lumon.set_level( _loglvl[_iopt.mprmon] );
        }

        if ( (jacgen == 1) && (ifail <  0) ) { _ierr = 183; return _ierr; }
        if ( (jacgen != 1) && (ifail != 0) ) { _ierr = 182; return _ierr; }

        // Copy Jacobian to work matrix _A(_m2,_n)
        // _A = -1.0 * _AA;   //  * _xw.diag();    // _A(1:_m2, 1:_n) = _AA(1:_m2, 1:_n)

        //_A.scale_columns(_xw);
        for (long k = 1; k <= (long)_n; ++k)
        {
            for (long j = 1; j <= (long)_m; ++j)
            {
                _A(j,k) = _AA(j,k) * _xw(k);
            }
        }

//std::cerr << "*** GaussNewton::computeSensitivity ***" << std::endl;
//std::cerr << " _fmodel = " << std::endl;
//std::cerr << _fmodel.t() << std::endl;
//std::cerr << " _xw = " << std::endl;
//std::cerr << _xw.t() << std::endl;
//std::cerr << " _AA = " << std::endl;
//std::cerr << _AA.t() << std::endl;
//std::cerr << " _A = " << std::endl;
//std::cerr << _A.t() << std::endl;
//std::cerr << " _fw = " << std::endl;
//std::cerr << _fw.t() << std::endl;
//std::cerr << "*** GaussNewton::computeSensitivity ***" << std::endl;
//std::cerr << std::endl;

        // Row scaling of _A(_m,_n)
        if ( qscale ) compute_row_scaling_A(); // else _fw.ones(_m);

        _qrA = _A.factorQRcon(_mcon, _irank, _cond);
    }

    return _ierr;
}
//---------------------------------------------------------------------------
QRconDecomp
GaussNewton::getSensitivity()
{
    return _qrA;
}
//---------------------------------------------------------------------------
Matrix
GaussNewton::getSensitivityMatrix()
{
    return _A;
}
//---------------------------------------------------------------------------



//---------------------------------------------------------------------------
int
GaussNewton::run()
{
    Real        ajdel, ajmin, cond1, condco = 0.0;
    Real        conva, dmue, dxnrm;
    Real        epdiff, etadif, etaini, etamax, etamin;
    Real        fcdnm, fck2, fcmin2, fcminh, fcbnd;
    Real        fcnump, fcnmp2, fcnumk, fch, fcbh, fcredu;
    Real        rsmall, sens1 = 0.0;
    Real        skap, sumxa, sumxk, th, dlevxa;
    Real        prec, dxanrm;
    unsigned    iscal, mode, nred, minmn, iterm, jacgen;
    unsigned    rscal, iranka, mfout, modefi, ifccnt, mconh;
    bool        qsucc, qgenj, qinisc, qjcrfr, qlinit, qredrnk;
    bool        qiter, qrepet, qnext, qredu, qscale, qrank1, qred;
    int         ifail;
    Vector      t1, t2, t3;

    if ( _fun == 0 )
    {
        _ierr = 999;
        printl( _luerr, dlib::LINFO,
                " Error: %s\n        %s\n",
                "No problem (fcn/jac pair) set so far.",
                "Please use setProblem(&fun) first.");
        return _ierr;
    }

    ajdel = ajmin = etadif = etaini = etamax = etamin = 0.0;
    fck2 = sumxk = sumxa = fcnump = fch = 0.0;
    mfout = iranka = ifccnt = 0;

    // Control flags and control integers
    qsucc  = _wk.qsucc;                                     // _wk[Wk::QSUCC];
    qscale = (_iopt.norowscal != true);                     // _iopt[IOpt::NOROWSCAL] != true;
    qrank1 = _iopt.qrank1;                                  // _iopt[IOpt::QRANK1];
    iscal  = _iopt.iscal;                                   // _iopt[IOpt::ISCAL];
    rscal  = _iopt.rscal;
    mode   = _iopt.mode;                                    // _iopt[IOpt::MODE];
    iterm  = _iopt.iterm;                                   // _iopt[IOpt::ITERM];
    jacgen = _iopt.jacgen;                                  // _iopt[IOpt::JACGEN];

    // _mcon = _m - _mfit;

    // Derivated internal parameters
    minmn  = std::min(_m, _n);
    fcmin2 = _fcmin * _fcmin;
    fcminh = std::sqrt(_fcmin);
    rsmall = std::sqrt(10.0 * _rtol);

    if ( _fc < _fcmin ) _fc = _fcmin;
    if ( _fc > 1.0 ) _fc = 1.0;

    // Initial preparations
    qjcrfr = qlinit = qrepet = false;
    qiter = true;
    ifail = 0;

    // _d.zeros(_n);                                           // = dlib::zeros_matrix<Real>(_n,1);
    // _pivot.zeros(_n);                                       // = dlib::zeros_matrix<Real>(_n,1);

    fcbnd = ( _qbdamp ) ? _wk.fcbnd : 0.0;                  // fcbnd = ( _qbdamp ) ? _wk[Wk::FCBND] : 0.0;

    // Numerical differentation related initialisation
    if ( jacgen == 2 )
    {
        ajdel = _wk.ajdel;                                  // _wk[Wk::AJDEL];
        if ( ajdel <= SMALL ) ajdel = std::sqrt(10.0 * EPMACH);
        ajmin = _wk.ajmin;                                  // _wk[Wk::AJMIN];
    }
    else if ( jacgen == 3 )
    {
        etadif = _wk.etadif;                                // _wk[Wk::ETADIF];
        if ( etadif <= SMALL ) etadif = 1.0e-6;
        etaini = _wk.etaini;                                // _wk[Wk::ETAINI];
        if ( etaini <= SMALL ) etaini = 1.0e-6;
        epdiff = std::sqrt(10.0 * EPMACH);
        etamax = std::sqrt(epdiff);
        etamin = epdiff*etamax;
    }

    // Miscellaneous preparations of first iteration step
    if ( !qsucc )
    {
        _sumxall.resize(_nitmax+1);
        _dlevfall.resize(_nitmax+1);
        _sumxqall.resize(_nitmax+1);
        _tolall.resize(_nitmax+1);
        _fcall.resize(_nitmax+1);

        _niter    = _ncorr = _nrejr1 = _nfcn = _nfcnj = ifccnt = 0;
        _wk.niter = _niter;                                 // _wk[Wk::NITER] = _niter;
        _xiter.clear(); _xiter.resize(_nitmax+2);
        _xiter[0] = _x;
        _rq.zeros(_m);                                      // = dlib::zeros_matrix<Real>(_m,1);
        _rqkeep.zeros(_m);                                  // = dlib::zeros_matrix<Real>(_m,1);
        _AA.zeros(_m,_n);
        _A.zeros(_m,_n);
        // _fmodel.zeros(_m);
        iranka = _irank;
        qgenj  = qinisc = true;
        fck2   = _fca   = _fc;
        sumxk  = _conv  = 0.0;
        if ( jacgen == 3 )                                  // dlib::ones_matrix<Real>(_n,1);
        {
            Vector v;
            v.ones(_n); _eta.zeros(_n);
            _eta = etaini * v;
        }
        _xa    = _x;
        // _xw.ones(_n);

        // Print monitor header
        printl( _lumon, dlib::LVERB,
                "\n\n\n %s\n\n %s\n",
                "*******************************************************************",
                "   It      Normf               Normx       Damp.Fct.   New     Rank"
              );

        // Startup step
        // Computation of residual vector
        //_timon(1);
        _fmodel = call_FCN( _x, ifail );        // _fun->fcn(_x, ifail);
        //_timoff(1);

        ++_nfcn;

        if ( ifail != 0 )
        {
            _ierr = 82;
            qiter = false;
        }
        else
        {
            // _f = dlib::zeros_matrix<Real>(_m,1);
            // dlib::set_rowm( _f, range(_mcon+1, _mcon+_mfit) ) =
            //         dlib::rowm( _fmodel, range(_mcon+1, _mcon+_mfit) ) - dlib::rowm( _fi, range(1, _mfit) );
            // dlib::set_rowm( _f, range(1, _mcon) ) =
            //         dlib::rowm( _fmodel, range(1, _mcon) );
            _f.zeros(_m);
            _f.set_row( _mcon+1, _mcon+_mfit ) = _fmodel.row( _mcon+1, _mcon+_mfit ) - _fi.row( 1, _mfit );
            if ( _mcon > 0 ) _f.set_row( 1, _mcon ) = _fmodel.row( 1, _mcon );
        }

        ///

        if ( qscale ) compute_scaling_fw(rscal);

        ///
    }
    else
    {
        qinisc = false;
    }

    // Main iteration loop
    // ===================

    while ( qiter )
    {
        // Startup of iteration step
        if ( !qjcrfr )
        {
            compute_scaling_xw(iscal, qinisc);

            qinisc = false;

            if ( _niter != 0 )
            {
                // Prelininary pseudo-rank
                iranka = _irank;
                _dxqa = _dxq;

                // Computation of linear residuum
                if ( (_m > _n) || (iranka != _m) )
                {
                    // _rq = dlib::zeros_matrix<Real>(_m,1);
                    // dlib::set_rowm( _rq, range(_mcon+1,_m) ) =
                    //    dlib::subm( _AA, range(_mcon+1,_m) , range(1,_n) ) * dlib::rowm( _dxqa, range(1,_n) ) +
                    //    dlib::rowm(  _f, range(_mcon+1,_m) );
                    _rq.zeros(_m);
                    _rq.set_row( _mcon+1, _m ) =
                            _AA.subm( _mcon+1, _m, 1, _n ) * _dxqa.row( 1, _n ) + _f.row( _mcon+1, _m );
                    _rqkeep = _rq;
                    mfout = _mfit;
                }
                else
                {
                    mfout = 0;
                }

                // A posteriori estimate of damping factor
                sumxk  =  // Denominator for kappa (SKAP) estimate
                // fcnump = dlib::length_squared(
                //            dlib::pointwise_multiply(_dx, dlib::reciprocal(_xw)) );
                fcnump = length_squared( _dx.dotdiv(_xw) );
                th = _fc - 1.0;
                // fcdnm = dlib::length_squared(
                //            dlib::pointwise_multiply((_dxqa + th*_dx), dlib::reciprocal(_xw)) );
                fcdnm = length_squared( (_dxqa + th*_dx).dotdiv(_xw) );
                fch = _fc * _fc * 0.5 * std::sqrt(fcnump / fcdnm);

                // Computation of the numerator of the damping factor predictor
                // fcnmp2 = dlib::length_squared(
                //                     dlib::pointwise_multiply( _dxqa, dlib::reciprocal(_xw)) );
                fcnmp2 = length_squared( _dxqa.dotdiv(_xw) );
                fcnump *= fcnmp2;

                // Decision critrion for Jacobian updating technique
                qgenj = ( (_fc < _fca) && (_new > 0) ) ||
                        (fch < _fc*_sigma) || (iranka != _n) || (_m > _n) || !qrank1;
                fck2 = _fca = _fc;
                _irank = minmn;
                if ( _nonlin > 1 ) _fc = std::min(1.0, fch);
            } // if ( _niter != 0 )

        } // if ( !qjcrfr )

        qjcrfr = false;

        // Jacobian matrix ( stored to _AA(_m2,_n) )
        // Jacobian computation by JAC -or- diff.approx (qgenj==true) -or- rank-1 update (qgenj==false)
        if ( qgenj )
        {
            _new = 0;
            if ( jacgen == 1 )
            {
                //_timon(2);
                _AA = call_JAC( _x, ifail );    // _fun->jac(_x, ifail);
                //_timoff(2);
            }
            else
            {
                //_timon(2);
                if ( jacgen == 3 ) compute_jcf_AA(etamin, etamax, etadif, ifail);
                if ( jacgen == 2 ) compute_jac_AA(ajdel, ajmin, ifail);
                //_timoff(2);
            }

            ++_njac;

           if ( (jacgen == 1) && (ifail <  0) ) { _ierr = 83; break; }
           if ( (jacgen != 1) && (ifail != 0) ) { _ierr = 82; break; }
        }
        else
        {
            ++_new;
            compute_rank1_up_AA();
        }

        // Copy Jacobian to work matrix _A(_m2,_n)
        _A = -1.0 * _AA;    // _A(1:_m2, 1:_n) = _AA(1:_m2, 1:_n)
        //_A.scale_columns(_xw);
        for (long k = 1; k <= (long)_n; ++k)
            for (long j = 1; j <= (long)_m; ++j)
            {
                _A(j,k) = _A(j,k) * _xw(k);
            }
        // Row scaling of _A(_m,_n)
        if ( qscale ) compute_row_scaling_A(); else _fw.ones(_m);  // = dlib::ones_matrix<Real>(_m,1);
        // Save and scale _f(_m)
        _fa = _f;
         t3 = _f.dotmult(_fw);  // dlib::pointwise_multiply(_f, _fw);
        // Scaling of linear residuum _rq(_m)
        if ( (_m > _n) || (iranka != _m) )
            // dlib::set_rowm( _rq, range(1,_m) ) =
            //     dlib::pointwise_multiply( dlib::rowm( _rq, range(1,_m) ), dlib::rowm( _fw, range(1,_m) ) );
            _rq.set_row(1,_m) = _rq.row(1,_m).dotmult(_fw.row(1,_m));

        qnext = false;
        qredrnk = true;

        // Central part of iteration step

        // Pseudo-rank reduction loop
        // ==========================

        while ( qredrnk )
        {
            // Decomposition of _A
            cond1 = _cond;
            mconh = _mcon;
            //_timon(3);
            compute_qrdcmp_A(qrepet, mconh, cond1, condco, sens1, ifail);
            //_timoff(3);

            qlinit = true;
            if ( ifail != 0 ) { _ierr = 80; qiter = false; break; }
            // save irankc, sens1, condco from qrdcmp_A ?!??

            // Solution of linear system
            //_timon(4);
            t2 = solve_A(qrepet, t3, ifail);
            //_timoff(4);

            if ( ifail != 0 ) { _ierr = 81; qiter = false; break; }

            //
            if ( (!qrepet) && (_irank != 0) )
            {
                // dlib::set_rowm( _qu, range(1,_m) ) = dlib::rowm( t3, range(1,_m) );
                _qu.zeros(_m);
                _qu.set_row(1,_m) = t3.row(1,_m);
            }

            // Evaluation of scaled natural level function SUMX
            //  scaled max err norm CONV, (scaled) std level function DLEVF,
            //  and computation of ordinary Gauss-Newton correction _dx(_n)
            _wk.sumx  = _sumx;
            _wk.dlevf = _dlevf;
            _dx    = compute_levels_and_dx(t2);
            _xa    = _x;
            sumxa  = _sumx;
            dlevxa = std::sqrt(sumxa / _n);
            conva  = _conv;
            dxanrm = wnorm(_dx);
            _sumxall[_niter]  = dlevxa;
            _dlevfall[_niter] = _dlevf;

            // A priori estimate of damping factor FC
            qredu = false;

            if ( (_niter != 0) && (_nonlin != 1) )
            {
                if ( (_new == 0) || qrepet )
                {
                    // Computation of denominator of a priori estimate
                    if ( (_m > _n) || (iranka != _m) )
                    {
                        //_timon(4);
                        _delxq = solve_A(qrepet, _rq, ifail);
                        //_timoff(4);

                        if ( ifail != 0 ) { _ierr = 81; qiter = false; break; }

                        // fcdnm = dlib::length_squared(
                        //             dlib::pointwise_multiply((_dx-_dxq), dlib::reciprocal(_xw)) - _delxq );
                        fcdnm = length_squared( (_dx - _dxq).dotdiv(_xw) - _delxq );
                    }
                    else
                    {
                        // fcdnm = dlib::length_squared(
                        //            dlib::pointwise_multiply((_dx-_dxq), dlib::reciprocal(_xw)) );
                        fcdnm = length_squared( (_dx - _dxq).dotdiv(_xw) );
                    }
                    if ( _irank != _n )
                    {
                        // Rank-deficient case
                        // t1 = dlib::pointwise_multiply( _dxqa, dlib::reciprocal(_xw) );
                        t1     = _dxqa.dotdiv(_xw);
                        fcdnm -= _qrA.projectIntoSubspace(t1);
                    }
                    fcdnm *= _sumx;

                    // New damping factor
                    if ( fcdnm > fcnump*fcmin2 )
                    {
                        dmue = _fca * std::sqrt(fcnump / fcdnm);
                        _fc = std::min(dmue, 1.0);
                        if ( dmue > 10.0 ) ++ifccnt;
                    }
                    else
                    {
                        _fc = 1.0;
                        dmue = -1.0;
                        if ( _fca >= 1.0 ) ++ifccnt;
                    }

                    printl( _lumon, dlib::LDEBUG,
                            " %s\n"
                            " %s%10u             %s%18.10e\n"
                            " %s%18.10e     %s%18.10e\n"
                            " %s%18.10e     %s%18.10e\n"
                            " %s%18.10e     %s%18.10e\n"
                            " %s\n",
                            "+++ a priori estimate +++",
                            " ifccnt = ", ifccnt,
                            " fc     = ", _fc,
                            " fca    = ", _fca,
                            " dmue   = ", dmue,
                            " fcnump = ", fcnump,
                            " fcdnm  = ", fcdnm,
                            " sumx   = ", _sumx,
                            " conv   = ", _conv,
                            "+++++++++++++++++++++++++"
                          );

                    if ( _qbdamp )
                    {
                        fcbh = _fca * fcbnd;
                        if ( _fc > fcbh )
                        {
                            _fc = fcbh;
                            printl( _lumon, dlib::LGABBY,
                                    " %s\n",
                                    "*** incr. rest. act. (a priori) ***"
                                  );
                        }
                        fcbh = _fca / fcbnd;
                        if ( _fc < fcbh )
                        {
                            _fc = fcbh;
                            printl( _lumon, dlib::LGABBY,
                                    " %s\n",
                                    "*** decr. rest. act. (a priori) ***"
                                  );
                        }
                    }

                } // if ( _new == 0 || qrepet )

                qredu = (_fc < _fcmin);

            } // if ( _niter != 0 && _nonlin != 1 )

            qrepet = false;

            if ( !qredu )
            {
                // Save natural level for later
                fcnumk = _sumx;

                //_timon(5);
                log_iteration_vals1(dlevxa); // dlib::LVERB
                //_timoff(5);

                nred = 0;
                qred = true;
                // Damping factor reduction loop
                // =============================
                while ( qred )
                {
                    // Preliminary new iterate
                    _x = _xa + _fc*_dx;
                    _fcall[_niter] = _fc;

                    if ( _nonlin == 1 ) { _ierr = 0; qiter = false; break; }

                    // Computation of residual vector
                    //_timon(1);
                    _fmodel = call_FCN( _x, ifail );     // _fun->fcn(_x, ifail);
                    //_timoff(1);

                    ++_nfcn;

                    if ( ifail < 0 ) { _ierr = 82; qiter = false; break; }
                    //
                    if ( (ifail == 1) || (ifail == 2) )
                    {
                        fcredu = ( ifail == 1 ) ? 0.5 : _fmodel(1);

                        if ( (fcredu <= 0.0) || (fcredu >= 1.0) ) { _ierr = 83; qiter = false; break; }

                        printl( _lumon, dlib::LVERB,
                                " %8s%2d%41s%5.3f%4s%2d        %4d\n",
                                "  ", _niter,
                                "fcn could not be evaluated", _fc,
                                "  ", _new, _irank
                              );
                        fch = _fc;
                        _fc *= fcredu;
                        if ( fch > _fcmin ) _fc = std::max(_fc, _fcmin);
                        if ( _qbdamp )
                        {
                            fcbh = fch / fcbnd;
                            if ( _fc < fcbh )
                            {
                                _fc = fcbh;
                                printl( _lumon, dlib::LGABBY,
                                        " %s\n",
                                        "*** decr. rest. act. (fcn redu.) ***"
                                      );
                            }
                        }
                        if ( _fc < _fcmin ) { _ierr = 3; qiter = false; break; }

                        break;
                    }

                    // dlib::set_rowm( _f, range(_mcon+1,_mcon+_mfit) ) =
                    //         dlib::rowm( _fmodel, range(_mcon+1,_mcon+_mfit) ) - dlib::rowm( _fi, range(1,_mfit) );
                    // dlib::set_rowm( _f, range(1,_mcon) ) =
                    //         dlib::rowm( _fmodel, range(1,_mcon) );
                    // dlib::set_rowm( t3, range(1,_m) ) =
                    //         dlib::pointwise_multiply( dlib::rowm(_f, range(1,_m) ), dlib::rowm( _fw, range(1,_m) ) );
                    _f.set_row( _mcon+1, _mcon+_mfit ) = _fmodel.row( _mcon+1, _mcon+_mfit ) - _fi.row( 1, _mfit );
                    if ( _mcon > 0 ) _f.set_row( 1, _mcon ) = _fmodel.row( 1, _mcon );
                    t3.set_row( 1, _m ) = _f.row( 1, _m ).dotmult( _fw.row(1,_m) );

                    // Solution of linear (_mfit,_n)-least squares problem
                    //_timon(4);
                    t2 = solve_A(qrepet, t3, ifail);
                    //_timoff(4);

                    if ( ifail != 0 ) { _ierr = 81; qiter = false; break; }

                    // Evaluation of SUMX, CONV, DLEVF, and reduced Gauss-Newton correction _dxq
                    _dxq = compute_levels_and_dx(t2);
                    _sumxqall[_niter] = std::sqrt(_sumx / _n);
                    dxnrm = wnorm(_dxq);

                    // Convergence test
                    _tolall[_niter] = dxnrm;

                    if ( iterm == 0 )
                    {
                        if ( ((dxnrm <= _rtol) && (ifccnt >= 3)) ||
                             (dxanrm <= _rtol)
                           ) { _ierr = 0; qiter = false; break; }
                    }
                    else if ( iterm == 1 )
                    {
                        if ( (dxnrm <= _rtol) && (ifccnt >= 3)
                           ) { _ierr = 0; qiter = false; break; }
                    }
                    else if ( iterm == 2 )
                    {
                        if ( (dxanrm <= _rtol)
                           ) { _ierr = 0; qiter = false; break; }
                    }

                    _fca = _fc;

                    // Natural monotonicity test
                    if ( _sumx > sumxa )
                    {
                        // Output of iterate
                        //_timon(5);
                        log_iteration_vals2( std::sqrt(_sumx/_n) , _niter , "*" ); // dlib::LTALK
                        //_timoff(5);

                        // Evaluation of reduced damping factor
                        th = _fca - 1.0;
                        // fcdnm = dlib::length_squared(
                        //            dlib::pointwise_multiply((_dxq+th*_dx), dlib::reciprocal(_xw)) );
                        fcdnm = length_squared( (_dxq + th*_dx).dotdiv(_xw) );
                        _fc = _fca * _fca * 0.5 * std::sqrt(fcnumk / fcdnm);
                        if ( _qbdamp )
                        {
                            fcbh = _fca / fcbnd;
                            if ( _fc < fcbh )
                            {
                                _fc = fcbh;
                                printl( _lumon, dlib::LGABBY,
                                        " %s\n",
                                        "*** decr. rest. act. (a post) ***"
                                      );
                            }
                        }

                        ++_ncorr;
                        ++nred;
                        ifccnt = 0;

                        // Rank reduction if damping factor is too small
                        qredu = (_fc < _fcmin) || ((_new > 0) && (nred > 1));
                    }
                    else
                    {
                        qnext = true;
                    }

                    qred = !(qnext || qredu);

                } // while ( qred )

                if ( (!qredrnk) || (!qiter) ) break;

                // End of damping factor reduction loop
                // ====================================
            } // if ( !qredu )

            if ( qredu )
            {
                // Restore former values for repeating step
                ++_nrejr1;
                _x = _xa;
                // dlib::set_rowm( _f, range(1,_m) ) = dlib::rowm( _fa, range(1,_m) );
                _f.set_row( 1, _m ) = _fa.row( 1, _m );
                _dxq = _dxqa;

                printl( _lumon, dlib::LVERB,
                        "    %2d %40s%5.3f     %2d     %4d\n",
                        _niter, "not accepted damping factor: ", _fc, _new, _irank
                     );

                ifccnt = 0;
                _fca = fck2;
                if ( _niter == 0 ) _fc = _fcmin;
                if ( _new > 0 )
                {
                    qgenj = qjcrfr = true;
                    qredu = false;
                    _fc = fch;
                    _irank = minmn;
                    // dlib::set_rowm( _rq, range(1,_m) ) = dlib::rowm( _rqkeep, range(1,_m) );
                    _rq.set_row( 1, _m ) = _rqkeep.row( 1, _m );
                }
                else
                {
                    // Pseudo-rank reduction
                    qrepet = true;
                    t3.set_row( 1, _m )  = _qu.row( 1, _m );
                    --_irank;
                    if ( _irank == 0 ) { _ierr = 3; qiter = false; break; }
                }
            } // if ( qredu )

            qredrnk = qredu;

        } // while ( qredrnk )

        if ( !qiter ) break;
        // End of pseudo-rank reduction loop
        // =================================

        if ( qnext )
        {
            // Preparations to start the next iteration step
            //_timon(5);
            log_iteration_vals2( std::sqrt(_sumx/_n) , _niter+1 , "*" ); // dlib::LTALK
            //_timoff(5);

            // print natural level of current iteration,
            // and (it) if in single-step mode
            _sumxs = _sumx;
            _sumx = sumxa;

            //_timon(5);
            if ( _niter != 0 )
                log_solout(mfout, _xa, _rqkeep.row(_mcon+1,_m), 2, dlib::LVERB);
            else if ( _niter == 0 )
                log_solout(0, _xa, _rq, 1, dlib::LINFO);
            //_timoff(5);

            ++_niter;
            _wk.niter = _niter;                             // _wk[Wk::NITER] = _niter;
            _xiter[_niter] = _x;

            if ( _niter >= _nitmax ) { _ierr = 2; qiter = false; break; }

            _fca = _fc;

            if ( mode == 1 ) { _ierr = -1; _wk.qsucc = true; /*_wk[Wk::QSUCC] = true;*/ return _ierr; }

        } // if ( qnext )

    } // while ( qiter )

    // End of main iteration loop
    // ==========================

    // Exits
    // -----

    // Solution exit
    //
    if ( _ierr == 0 )
    {
        if ( _nonlin != 1 )
        {
            _x = _x + _dxq;
            _xiter[_niter+1] = _x;  // one more update to store!!
            if ( _irank < minmn ) _ierr = 1;

            //_timon(5);
            log_iteration_vals2( std::sqrt(_sumx /_n) , _niter+1 , "."); // dlib::LVERB
            //_timoff(5);

            if ( _ierr == 0 )
            {
                printl( _lumon, dlib::LINFO,
                        "\n\n\n %s\n %s %3d %s\n",
                        "Solution of _nonlinear_ least squares problem",
                        "obtained within", _niter+1, " iteration steps."
                      );
            }
        }
        else
        {
            printl( _lumon, dlib::LINFO,
                    "\n\n\n %s\n\n %s\n",
                    "Least squares solution of _linear_ system.",
                    "No estimate avaible for achieved relative accuracy."
                  );
        }

    } // if ( _ierr == 0 )

    // Fail exit messages
    //
    if ( _ierr == 1 )
    {
        printl( _luerr, dlib::LINFO,
                " %s\n",
                "Iteration terminates at stationary point."
              );
    }
                                                            // prec = _wk[Wk::PREC] = _wk[Wk::SKAP] = -1.0;
    prec = _wk.prec = _wk.skap = -1.0;

    if ( (_ierr == 0 || _ierr == 1) && _nonlin != 1 )
    {
        if ( sumxk < _tolmin )
            skap = -1.0;
        else if ( sumxa  < _tolmin )
            skap = std::sqrt(_tolmin);
        else
            skap = std::sqrt(sumxa / sumxk);

        if ( skap >= 0.0 )
        {
            printl( _lumon, dlib::LINFO,
                      "\n Incompatibility factor kappa %10.3e"
                      "\n (sumxa =%10.3e , sumxk =%10.3e , tolmin =%10.3e)\n", skap, sumxa, sumxk, _tolmin
                  );
        }
        else
        {
            printl( _lumon, dlib::LINFO,
                    "\n Incompatibility factor kappa not available\n"
                  );
        }

        if ( _ierr == 0 && skap < 1.0 )
        {
            prec = (skap >= 0.0) ? std::max(std::sqrt(sumxa/_n)*skap/(1.0 - skap), EPMACH) : EPMACH;

            printl( _lumon, dlib::LINFO,
                    "\n Achieved relative accuracy %10.3e\n", prec
                  );
        }

        _wk.prec = prec;                                    // _wk[Wk::PREC] = prec;
        _wk.skap = skap;                                    // _wk[Wk::SKAP] = skap;

    } // if ( (_ierr == 0 || _ierr == 1) && _nonlin != 1 )

    _rtol = _wk.prec;                                       // _wk[Wk::PREC];

    //
    if ( _ierr == 2 )
    {
        printl( _luerr, dlib::LINFO,
                " %s %3d %s\n",
                "Iteration terminates after NITMAX =", _nitmax, "steps."
              );
    }
    //
    if ( _ierr == 3 )
    {
        printl( _luerr, dlib::LINFO,
                " %s\n",
                "Gauss-Newton method fails to converge."
              );
    }
    //
    if ( _ierr == 80 )
    {
        printl( _luerr, dlib::LINFO,
                " %s %5d %s\n",
                "Error", ifail,
                "given by linear solver (decomp)."
              );
    }
    //
    if ( _ierr == 81 )
    {
        printl( _luerr, dlib::LINFO,
                " %s %5d %s\n",
                "Error", ifail,
                "given by linear solver (solve)."
              );
    }
    //
    if ( _ierr == 82 )
    {
        printl( _luerr, dlib::LINFO,
                " %s %5d %s\n",
                "Error", ifail,
                "given by user function FCN."
              );
    }
    //
    if ( _ierr == 83 )
    {
        printl( _luerr, dlib::LINFO,
                " %s %5d %s\n",
                "Error", ifail,
                "given by user function JAC."
              );
    }

    if ( _ierr >= 80 && _ierr <= 83 ) _wk.ifail = ifail;    // _wk[Wk::IFAIL] = ifail;

    if ( (_ierr == 82 || _ierr == 83) && _niter <= 1 )
    {
        printl( _luerr, dlib::LINFO,
                " %s\n",
                "Try to find a better initial guess."
              );
    }

    // Common exit
    //
    if ( _mcon > 0 )
    {
        printl( _lumon, dlib::LINFO,
                "\n\n   Subcondition (    1,%4d) of constrained part %10.3e\n",
                _irankc, condco
              );
        printl( _lumon, dlib::LINFO,
                "\n\n   Subcondition ( %4d,%4d) of least squares part %10.3e\n",
                _irankc+1, _irank, cond1
              );
    }
    else
    {
        printl( _lumon, dlib::LINFO,
                "\n\n   Subcondition (    1,%4d) of least squares part %10.3e\n",
                _irank, cond1
              );
    }

    printl( _lumon, dlib::LINFO,
            "\n\n   Sensitivity (lsq) %10.3e\n\n",
            sens1
          );

    _sumxs = _sumx;
    _sumx = sumxa;

    //_timon(5);
    if ( _niter != 0 )
        log_solout(mfout, _xa, _rqkeep.row(_mcon+1,_m), 2, dlib::LVERB);
    else if ( _niter == 0 )
        log_solout(0, _xa, _rq, 1, dlib::LINFO);
    //_timoff(5);

    ++_niter;
    _wk.niter = _niter;                                     // _wk[Wk::NITER] = _niter;

    modefi = (_ierr == 0) ? 3 : 4;

    //_timon(5);
    log_solout(0, _x, _rq, modefi, dlib::LINFO);
    //_timoff(5);

    _xscal = _xw;

    return _ierr;
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
Vector
GaussNewton::call_FCN(Vector const& x, int& ifail)
{

    if ( _iopt.lpos == true )
    {
        long   n = x.nr();
        Vector y(n);

        for(long j = 1; j <= n; ++j) y(j) = std::exp( x(j) );

        return _fun->fcn( y, ifail );
    }

    return _fun->fcn( x, ifail );
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Matrix
GaussNewton::call_JAC(Vector const& x, int& ifail)
{

    if ( _iopt.lpos == true )
    {
        long   n = x.nr();
        Vector y(n);

        for(long j = 1; j <= n; ++j) y(j) = std::exp( x(j) );

        Matrix J = _fun->jac( y, ifail );

        return J * y.diag();
    }

    return _fun->jac( x, ifail );
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
void
GaussNewton::compute_scaling_xw(unsigned iscal, bool /*qinisc*/)
{
    if ( iscal == 1 )
    {
        _xw = _xscal;
    }
    else
    {
        _xw.zeros(_n); // = dlib::zeros_matrix<Real>(_n,1);

        for (long j = 1; j <= _xw.nr(); ++j)
        {
            _xw(j) = std::max( _xscal(j),
                               std::max( 0.5*(std::fabs(_x(j)) + std::fabs(_xa(j))), SMALL)
                             );
        }
    }

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::compute_scaling_fw(unsigned rscal)
{
    if ( rscal == 3 )
    {
        for (long j = 1; j <= _fw.nr(); ++j)
        {
            Real r = std::fabs( _f(j) );

            _fw(j) = ( (SMALL <= r) && (r <= GREAT)) ?
                            1.0 / _f(j) :
                            SMALL;         // 1.0 ???;
        }
    }
    else if ( rscal == 2 )
    {
        for (long j = 1; j <= _fw.nr(); ++j)
        {
            Real r = std::fabs( _fw(j) / _f(j) );

            _fw(j) = ( (SMALL <= r) && (r <= GREAT) ) ?
                            _fw(j) / _f(j) :
                            _fw(j);
        }
    }

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::compute_row_scaling_A()
{
    Real    s1;
    long    mcon = _mcon;

    for (long k = 1; k <= mcon; ++k)
    {
        // s1 = dlib::max( dlib::abs( dlib::rowm( _A, k ) ));
        s1 = _A.rowm(k).norm("inf");

        if ( s1 > 0.0 )
        {
            s1     = 1.0 / s1;
            _fw(k) = s1;
            // dlib::set_rowm( _A, k ) = dlib::rowm( _A, k ) * s1;
            _A.set_rowm(k) = _A.rowm(k) * s1;
        }
        else
        {
            _fw(k) = 1.0;
        }
    }

    for (long k = mcon+1; k <= _A.nr(); ++k)
    {
        s1 = _fw(k);

        // dlib::set_rowm( _A, k ) = dlib::rowm( _A, k ) * s1;
        _A.set_rowm(k) = _A.rowm(k) * s1;
    }

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::compute_qrdcmp_A(
    bool        qrepet,
    unsigned&   mcon,
    Real&       cond,
    Real&       condc,
    Real&       sens,
    int&        ifail
    )
{

//std::cerr << "*** GaussNewton::compute_qrdcmp_A ***" << std::endl;
//std::cerr << "qrepet = " << (int)qrepet << std::endl;
//std::cerr << "mcon   = " << mcon << std::endl;
//std::cerr << "cond   = " << cond << std::endl;
//std::cerr << "condc  = " << condc << std::endl;
//std::cerr << "sens   = " << sens << std::endl;
//std::cerr << "ifail  = " << ifail << std::endl;
//std::cerr << std::endl;
//std::cerr << "_irank = " << _irank << std::endl;
//std::cerr << std::endl;
//std::cerr << "_A = (pre)" << std::endl;
//std::cerr << _A.t() << std::endl;
//std::cerr << "_qrA = (pre)" << std::endl;
//std::cerr << _qrA.getMat().t() << std::endl;
//std::cerr << std::endl;

    if ( !qrepet )
        _qrA = _A.factorQRcon(mcon, _irank, cond);
    else
        _qrA.setNewRank(_irank);

    _irank  = _qrA.getRank();
    _irankc = _qrA.getRankc();
    cond    = _qrA.getSubCond();
    condc   = _qrA.getSubCondc();
    sens    = std::fabs( _qrA.getDiag()(_irankc+1) );

    ifail = _qrA.getError().getIerr();

//std::cerr << "_irank  = " << _irank << std::endl;
//std::cerr << "_irankc = " << _irankc << std::endl;
//std::cerr << "cond    = " << cond << std::endl;
//std::cerr << "condc   = " << condc << std::endl;
//std::cerr << "sens    = " << sens << std::endl;
//std::cerr << "ifail   = " << ifail << std::endl;
//std::cerr << std::endl;
//std::cerr << "_A = (post)" << std::endl;
//std::cerr << _A.t() << std::endl;
//std::cerr << "_qrA = (post)" << std::endl;
//std::cerr << _qrA.getMat().t() << std::endl;
//std::cerr << "*** GaussNewton::compute_qrdcmp_A ***" << std::endl;
//std::cerr << std::endl;

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Vector
GaussNewton::solve_A(bool qrepet, Vector& b, int& ifail)
{
    Vector x;

//std::cerr << "*** GaussNewton::solve_A ***" << std::endl;
//std::cerr << "qrepet = " << (int)qrepet << std::endl;
//std::cerr << "ifail  = " << ifail << std::endl;
//std::cerr << std::endl;
//std::cerr << "_qrA = (pre)" << std::endl;
//std::cerr << _qrA.getMat().t() << std::endl;
//std::cerr << "b = (pre)" << std::endl;
//std::cerr << b.t() << std::endl;
//std::cerr << "x = (pre)" << std::endl;
//std::cerr << x.t() << std::endl;
//std::cerr << std::endl;

    if ( !qrepet )  _qrA.solve(b,x);
    else            _qrA.solveR(b,x);

//std::cerr << "_qrA = (post)" << std::endl;
//std::cerr << _qrA.getMat().t() << std::endl;
//std::cerr << "b = (post)" << std::endl;
//std::cerr << b.t() << std::endl;
//std::cerr << "x = (post)" << std::endl;
//std::cerr << x.t() << "   (" << length_squared(x) << ")" << std::endl;
//std::cerr << "*** GaussNewton::solve_A ***" << std::endl;
//std::cerr << std::endl;

    ifail = _qrA.getError().getIerr();

    return x;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Vector
GaussNewton::compute_levels_and_dx(Vector const& dx1)
{
    Vector dxq;

    // 1 Descaling of solution dx1 (stored to dxq)
    dxq = dx1.dotmult(_xw);

    // 2 Evalutation of scaled natural level function SUMX
    //   and scaled maximum error norm CONV
    _conv = dx1.norm("inf");
    _sumx = length_squared(dx1);

    // 3 Evaluation of (scaled) standard level function DLEVF
    _dlevf = std::sqrt( length_squared(_f)/_m );

    return dxq;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::compute_jac_AA(Real ajdelta, Real ajmin, int& ifail)
{
    Vector  fu;
    Real    w, u;
    long    m = _fmodel.nr();
    long    n = _x.nr();
    int     su;

    ifail = 0;
    _AA.zeros(m,n);

    for (long k = 1; k <= n; ++k)
    {
        printl( _lumon, dlib::LDEBUG,
                " %s: %s%04d\n",
                "jac", "computing difference quotient of parameter #", k);

        w = _x(k);

        su = ( w < 0.0 ) ? -1 : 1;
        u  = std::max( std::max(std::fabs(w), ajmin), _xw(k) );
        u *= ajdelta * su;

        _x(k) = w + u;

        fu = call_FCN( _x, ifail );      // _fun->fcn( _x, ifail );
        ++_nfcnj;

        _x(k) = w;

        if ( ifail != 0 ) break;

        _AA.set_colm(k) = (1.0/u) * Matrix(fu - _fmodel);
    }

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::compute_jcf_AA(
    Real etamin,
    Real etamax,
    Real etadif,
    int& ifail
)
{
    const Real small2 = 0.1;

    Vector  fu;
    Real    w, u, sumd, hg, fhj;
    bool    is, qexit, qfine;
    int     su;

    long m = _fmodel.nr();
    long n = _x.nr();

    ifail = su = 0;
    w = u = sumd = hg = fhj = 0.0;
    _AA.zeros(m,n);

    for (long k = 1; k <= n; ++k)
    {
        printl( _lumon, dlib::LDEBUG,
                " %s: %s%04d\n",
                "jcf", "computing difference quotient of parameter #", k);

        is = qexit = qfine = false;
        while ( !qfine )
        {
            w  = _x(k);
            su = ( w < 0.0 ) ? -1 : 1;

            u = _eta(k) * _xw(k) * su;

            _x(k) = w + u;

            fu = call_FCN( _x, ifail );      // _fun->fcn( _x, ifail );
            ++_nfcnj;

            _x(k) = w;

            if ( ifail != 0 ) { qexit = true; break; }

            sumd = 0.0;
            for (long j = 1; j <= m; ++j)
            {
                hg = std::max( std::fabs(_fmodel(j)), std::fabs(fu(j)) );
                fhj = fu(j) - _fmodel(j);
                if ( hg != 0.0 ) sumd += (fhj/hg)*(fhj/hg);
                _AA(j,k) = fhj / u;
            }
            sumd = std::sqrt( sumd / m );
            qfine = true;
            if ( (sumd != 0.0) && (is == false) )
            {
                _eta(k) = std::min(
                            etamax,
                            std::max( etamin, _eta(k)*std::sqrt(etadif/sumd) )
                          );
                is = true;
                qfine = ( (_conv < small2) || (sumd >= etamin) );
            }
        }

        if ( qexit ) break;
    }

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::compute_rank1_up_AA()
{
    Vector  dxj, dxf;
    Real    dnm, s1;

    long m = _AA.nr();
    long n = _AA.nc();

    dxj  = _dx.dotdiv(_xw);
    dnm  = length_squared( dxj.row(1,n) );
    dxj  = dxj.dotdiv(_xw);
    dnm *= _fca;
    if ( dnm != 0.0 )
    {
        s1  = _fca - 1.0;
        dxf = _f.row(1,m) + _fa.row(1,m)*s1;
        for (long k = 1; k <= n; ++k)
        {
            s1 = dxj(k) / dnm;
            _AA.set_colm(k) = Matrix( _AA.colm(k) + s1 * dxf.row(1,m) );
        }
    }

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Real
GaussNewton::wnorm(Vector const& z) const
{
    return std::sqrt( length_squared( z.dotdiv(_xw) ) / z.nr() );
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::log_iteration_vals1(Real dlevx)
{
    int monprio  = _lumon.level().priority;
    int verbprio = dlib::LVERB.priority;
    //

    if ( monprio != verbprio )
    {
        printl( _lumon, dlib::LALL,
                " %s\n",
                "   It      Normf               Normx       Damp.Fct.   New     Rank"
              );
    }
    else
    {
        printl( _lumon, dlib::LGABBY,
                " %s\n",
                "   It      Normf               Normx                   New     Rank"
              );
    }

    //

    if ( _niter == 0 )
    {
        printl( _lumon, dlib::LVERB,
                "  %4d     %14.7e      %10.3e                %2d     %4d\n",
                _niter, _dlevf, dlevx, _new, _irank
              );
    }
    else
    {
        printl( _lumon, dlib::LVERB,
                "  %4d     %14.7e      %10.3e                %2d     %4d\n",
                _niter, _dlevf, dlevx, _new, _irank
              );
    }

    //

    if ( (monprio != verbprio) && (_niter != 0) )
    {
        printl( _lumon, dlib::LALL,
                "  %4d     %14.7e      %10.3e      %5.3f     %2d     %4d\n",
                _niter, _dlevf, dlevx, _fc, _new, _irank
              );
    }

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::log_iteration_vals2(Real dlevx, unsigned niter, std::string mark)
{
//    printl( _lumon, dlib::LVERB,
//            " %s\n",
//            "   It      Normf           Normx           Damp.Fct."
//          );

    printl( _lumon, dlib::LVERB,
            "  %4d     %14.7e    %1s %10.3e      %5.3f\n",
            niter, _dlevf, mark.c_str(), dlevx, _fc
          );

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::log_solout(
    unsigned                mfit,
    Vector const&           x,
    Vector const&           rq,
    int                     mode,
    dlib::log_level const   loglvl
    )
{
    char line[256];
    long n = _x.nr();
    bool qnorm = true;

    if ( qnorm )
    {
        if ( mode == 1 )
        {
            printl( _lusol, dlib::LINFO,
                    " %s\n %s%d\n\n %s\n",
                    "Start data:",
                    "N = ", n,
                    "Format: iteration number, (x(j), j=1,...,N), Normf, Normx"
                  );
            printl( _lusol, dlib::LINFO,
                    " %s\n",
                    "Initial data:"
                  );
        }
        else if ( mode == 3 )
        {
            printl( _lusol, dlib::LINFO,
                    " %s\n",
                    "Solution data:"
                  );
        }
        else if ( mode == 4 )
        {
            printl( _lusol, dlib::LINFO,
                    " %s\n",
                    "Final data:"
                  );
        }

        printl( _lusol, dlib::LINFO,
                " %5d\n",
                _wk.niter
              );

        *line = '\0';
        for (long l1 = 1, l2 = 0; l1 <= n; ++l1)
        {
            sprintf( line, "%s%18.10e ", line, x(l1) );
            ++l2;
            if ( l2 == 3 )
            {
                printl( _lusol, dlib::LINFO, " %s\n", line);
                *line = '\0';
                l2 = 0;
            }
        }
        if ( *line != '\0' ) printl( _lusol, dlib::LINFO, " %s\n", line );
        printl( _lusol, dlib::LINFO,
                " %18.10e %18.10e\n",
                // _dlevf, sqrt(_sumx / n)
                _wk.dlevf, std::sqrt(_wk.sumx/n)
              );

        if ( mode == 2 )
        {
            if ( mfit != 0 )
            {
                printl( _lusol, loglvl,
                        "\n   %s\n",
                        "Residuum for current iteration:"
                      );
                *line = '\0';
                for (long l1 = 1, l2 = 0; l1 <= (long)mfit; ++l1)
                {
                    sprintf( line, "%s%18.10e ", line, rq(l1) );
                    ++l2;
                    if ( l2 == 3 )
                    {
                        printl( _lusol, loglvl, "   %s\n", line);
                        *line = '\0';
                        l2 = 0;
                    }
                }
                if ( *line != '\0' ) printl( _lusol, loglvl, "   %s\n", line );
            }
        }

        if ( mode == 1 )
            printl( _lusol, loglvl,
                    " %s\n",
                    "Intermediate data:"
                  );
        else if ( mode >= 3 )
            printl( _lusol, loglvl,
                    " %s\n",
                    "End data:"
                  );
    }

    return;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
int
GaussNewton::analyse()
{
    Real        epdiff = std::sqrt(10.0 * EPMACH);
    Real        ajdel, ajmin, etamin, etamax, etadif;
    Real        cond1, condco = 0.0, sens1 = 0.0;
    unsigned    mconh, jacgen = _iopt.jacgen;
    int         ifail;
    bool        qrepet;

    _iopt.qstat = true;
    ajdel = ajmin = etamin = etamax = etadif = 0.0;

   // Numerical differentation related initialisation
    if ( jacgen == 2 )
    {
        ajdel = _wk.ajdel;                                  // _wk[Wk::AJDEL];
        if ( ajdel <= SMALL ) ajdel = std::sqrt(10.0 * EPMACH);
        ajmin = _wk.ajmin;                                  // _wk[Wk::AJMIN];
    }
    else if ( jacgen == 3 )
    {
        etadif = _wk.etadif;                                // _wk[Wk::ETADIF];
        if ( etadif <= SMALL ) etadif = 1.0e-6;
        // epdiff = sqrt(10.0 * EPMACH);
        etamax = std::sqrt(epdiff);
        etamin = epdiff*etamax;
    }

    //_timon(1);
    _fmodel = call_FCN( _x, ifail );         // _fun->fcn(_x, ifail);
    //_timoff(1);

    ++_nfcn;

    if ( ifail != 0 )
    {
        _ierr = 182;
        printl( _luerr, dlib::LERROR,
                "\n %s\n %s\n",
                "Computation of the statistical analysis stopped",
                "since user function FCN failed."
              );
        return _ierr;
    }

    _f.zeros(_m);
    _f.set_row( _mcon+1, _mcon+_mfit ) = _fmodel.row( _mcon+1, _mcon+_mfit ) - _fi.row( 1, _mfit );
    if ( _mcon > 0 ) _f.set_row( 1, _mcon ) = _fmodel.row( 1, _mcon );

    if ( jacgen == 1 )
    {
        // _timon(2);
        _AA = call_JAC( _x, ifail );         // _fun->jac(_x, ifail);
        //_timoff(2);
    }
    else
    {
        //_timon(2);
        if ( jacgen == 3 ) compute_jcf_AA(etamin, etamax, etadif, ifail);
        if ( jacgen == 2 ) compute_jac_AA(ajdel, ajmin, ifail);
        //_timoff(2);
    }

    ++_njac;

    if ( (jacgen == 1) && (ifail <  0) )
    {
       _ierr = 183;
        printl( _luerr, dlib::LERROR,
                "\n %s\n %s\n",
                "Computation of the statistical analysis stopped",
                "since user function JAC failed."
              );
       return _ierr;
    }
    if ( (jacgen != 1) && (ifail != 0) )
    {
        _ierr = 182;
        printl( _luerr, dlib::LERROR,
                "\n %s\n %s\n",
                "Computation of the statistical analysis stopped",
                "since user function FCN failed."
              );
        return _ierr;
    }

    _A = _AA;

    // Decomposition of _A
    cond1  = _cond;
    mconh  = _mcon;
    qrepet = false;
    _irank = std::min(_m,_n);
    //_timon(3);
    compute_qrdcmp_A(qrepet, mconh, cond1, condco, sens1, ifail);
    //_timoff(3);

    if ( ifail != 0 )
    {
        _ierr = 180;
        printl( _luerr, dlib::LERROR,
                "\n %s\n %s\n",
                "Computation of the statistical analysis stopped",
                "since linear system (decomp) failed."
              );
        return _ierr;
    }

    if ( (_mfit - (_n-_mcon)) > 0 )
    {
        compute_statistics( _fmodel.row(_mcon+1,_m), _fi );
    }
    else
    {
        printl( _lumon, dlib::LINFO,
                "\n\n %s\n\n",
                "Statistical analysis only available for overdetermined systems."
              );
    }

    return _ierr;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
GaussNewton::compute_statistics(
    Vector const& ymodel,
    Vector const& y
    )
{
    Matrix rinv;
    Vector res, v;
    Real   sum1;
    long   m = _A.nr();
    long   n = _A.nc();

    _ierr = 0;

    _xl.zeros(n);
    _xr.zeros(n);
    _vcv.zeros(n,n);
    rinv.zeros(n,n);
    v.zeros(n);

    printl( _lumon, dlib::LINFO,
            "\n\n\n     %s\n\n\n",
            "Assuimng the classical linear model:"
          );

    res = ymodel.row(1,_mfit) - y.row(1,_mfit);
    sum1 = res.norm("l2");
    sum1 = (sum1*sum1) / (_mfit - (_irank-_mcon));

    _sigma2 = sum1;

    printl( _lumon, dlib::LINFO,
            " %s\n %s\n\n",
            "Best unbiased estimate of variance and std.dev. of residuals:",
            "-------------------------------------------------------------"
          );
    printl( _lumon, dlib::LINFO,
            " sigma2 =%10.3e    sigma =%10.3e\n",
            _sigma2, std::sqrt(_sigma2)
          );

    // Computation of covariance matrix
    if ( n != 1 )
    {
        // Apply pseudo-inverse to Cholesky decomposition of
        // variance covariance matrix of data error.
        // We assume that this matrix is sigma*I
        m = _mfit + _mcon;

        for (long j = 1; j <= n; ++j)
        {
            _xl.zero();
            _xl(j) = (j <= (long)_mcon) ? 0.0 : 1.0;

            _qrA.solveR(_xl, _xr);

            // rinv.set_colm( j ) = _xr;
            for (long k = 1; k <= n; ++k) rinv(k,j) = _xr(k);
        }

        // Product rinv * rinv^t
        _vcv = _sigma2 * ( rinv * rinv.t() );
    }
    else
    {
        Real d = (_qrA.getDiag())(1);
        _vcv(1,1) = _sigma2 / ( d*d );
    }

    // Pretty printout of the covariance matrix
    printl( _lumon, dlib::LINFO,
            "\n\n\n %s\n %s\n\n",
            "Covariance matrix of parameters",
            "-------------------------------"
          );

    {
        char line[256];
        long jh, nh1, nh2 = 0;

        while ( nh2 != n )
        {
            nh1 = nh2 + 1;
            nh2 = nh1 + 6;
            if ( nh2 > n ) nh2 = n;

            *line = '\0';
            for (long l = nh1; l <= nh2; ++l)
            {
                sprintf(line, "%s%10ld", line, l);
            }
            printl( _lumon, dlib::LINFO, " %s\n", line);

            for (long j = nh1; j <= n; ++j)
            {
                jh = j;
                if ( jh > nh2 ) jh = nh2;

                sprintf(line, " %3ld", j);
                for (long k = nh1; k <= jh; ++k)
                {
                    sprintf( line, "%s%10.2e", line, _vcv(j,k) );
                }
                printl( _lumon, dlib::LINFO, " %s\n", line);
            }
        }
    }


    // Computation and pretty printout of the correlation coefficients
    for (long j = 1; j <= n; ++j)
    {
        v(j) = std::sqrt( _vcv(j,j) );
        rinv(j,j) = (v(j) != 0.0) ? 1.0 : 0.0;
    }

    if ( n != 1 )
    {
        // nm1 = n - 1;
        for (long j = 1; j < n; ++j)
        {
            if ( v(j) == 0.0 )
            {
                for (long l = j+1; l <= n; ++l) rinv(l,j) = 0.0;
            }
            else
            {
                for (long k = j+1; k <= n; ++k)
                {
                    rinv(k,j) = (v(k) != 0.0) ? _vcv(k,j) / ( v(k)*v(j) ) : 0.0;
                }
            }
        }

        printl( _lumon, dlib::LINFO,
                "\n\n\n %s\n %s\n\n",
                "Correlation coefficients",
                "------------------------"
              );

        {
            char line[256];
            long jh, nh1, nh2 = 0;

            while ( nh2 != n )
            {
                nh1 = nh2 + 1;
                nh2 = nh1 + 9;
                if ( nh2 > n ) nh2 = n;

                *line = '\0';
                sprintf(line, "   ");
                for (long l = nh1; l <= nh2; ++l)
                {
                    sprintf(line, "%s%7ld", line, l);
                }
                printl( _lumon, dlib::LINFO, " %s\n", line);

                for (long j = nh1; j <= n; ++j)
                {
                    jh = j;
                    if ( jh > nh2 ) jh = nh2;

                    sprintf(line, " %3ld", j);
                    for (long k = nh1; k <= jh; ++k)
                    {
                        sprintf( line, "%s%7.2f", line, rinv(j,k) );
                    }
                    printl( _lumon, dlib::LINFO, " %s\n", line);
                }
            }
        }
    }

    // Standard error in parameters
    printl( _lumon, dlib::LINFO,
            "\n\n\n %s\n %s\n %s\n\n",
            "Standard deviation of parameters",
            "--------------------------------",
            "  No.  Estimate           sigma(X)"
          );

    for (long j = 1; j <= n; ++j)
    {
        _xr(j) = (_x(j) != 0.0)
                    ? std::fabs( 100.0*v(j)/_x(j) )
                    : std::fabs( 100.0*v(j) );
        printl( _lumon, dlib::LINFO,
                " %4ld  %10.3e   +/-  %10.3e    =%8.2f %%\n",
                j, _x(j), v(j), _xr(j)
              );
    }

    // Associated confidence intervals
    printl( _lumon, dlib::LINFO,
            "\n\n\n %s\n %s\n",
            "Independent confidence intervals",
            "--------------------------------"
          );
    {
        // FISH15: Array containing upper 5% values of fisher(1,l)-distribution
        // (L=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
        //  26,27,28,29,30,32,34,36,38,40,42,44,46,48,50,60,70,80,90,100,200,300)
        Real const fish15[] = {
              0.00,
            161.00, 18.51, 10.13, 7.71, 6.61, 5.99, 5.59,
              5.32,  5.12,  4.96, 4.84, 4.75, 4.67, 4.60, 4.54, 4.49,
              4.45,  4.41,  4.38, 4.35, 4.32, 4.30, 4.28, 4.26, 4.24,
              4.23,  4.21,  4.20, 4.18, 4.17, 4.15, 4.13, 4.11, 4.10,
              4.08,  4.07,  4.06, 4.05, 4.04, 4.03, 4.00, 3.98, 3.96,
              3.95,  3.94,  3.89, 3.88
        };

        long l, lh;
        Real fnma, h;

        l = _mfit - (n - _mcon);
        if ( l <= 30 )          lh = l;
        else if ( l <=  50 )    lh = 15 + long(0.5*l);
        else if ( l <= 100 )    lh = 35 + long(0.1*l);
        else if ( l <  300 )    lh = 45 + long(0.05*l);
        else                    lh = 47;

        fnma = fish15[lh];
        h = std::sqrt( fnma );

        printl( _lumon, dlib::LINFO,
                "   (on 95%%-prob.lev. using F-dist. F(alpha,1,m-n) =%6.2f)\n\n",
                fnma
              );

        for (long j = 1; j <= n; ++j)
        {
            _xr(j) = h*v(j);
            v(j) = 100.0 * _xr(j) / _x(j);
            _xl(j) = _x(j) - _xr(j);
            _xr(j) = _x(j) + _xr(j);

            printl( _lumon, dlib::LINFO,
                    " %4ld   ( %10.3e , %10.3e )\n",
                    j, _xl(j), _xr(j)
                  );
        }
    }

    return;
}
//---------------------------------------------------------------------------

