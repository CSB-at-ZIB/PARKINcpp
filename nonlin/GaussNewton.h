// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-10-19 td
// last changed:
//
#ifndef __GAUSS_NEWTON_H
#define __GAUSS_NEWTON_H

#include <vector>

// #include "addpkg/dlib/logger.h"

#include "common/Types.h"
#include "common/Constants.h"
#include "common/PARKINLog.h"

#include "linalg/QRconDecomp.h"

#include "UserFunc.h"

namespace PARKIN
{

    struct Info
    {
        // Real     ;
        // int     ;
        // bool    ;
    };

    //

    struct GaussNewtonWk
    {
        Vector  fw;
        Real    fcbnd, ajdel, ajmin, etadif, etaini;
        Real    dlevf, sumx, prec, skap;
        Real    fcstart, fcmin, sigma, cond;
        long    niter, nitmax, irank;
        long    ncorr, nrejr1, njac, nfcn, nfcnj;
        int     ifail;
        bool    qsucc;

        GaussNewtonWk() :
            fw(),
            fcbnd(0.0), ajdel(0.0), ajmin(0.0), etadif(0.0), etaini(0.0),
            dlevf(0.0), sumx(0.0), prec(0.0), skap(0.0),
            fcstart(0.0), fcmin(0.0), sigma(0.0), cond(0.0),
            niter(0), nitmax(0), irank(0),
            ncorr(0), nrejr1(0), njac(0), nfcn(0), nfcnj(0),
            ifail(0), qsucc(false)
        { }
    };

    ///
    ///

    class GaussNewton
    {
        public:
            // c'tor
            GaussNewton();
            //GaussNewton(IOpt const& iopt, GaussNewtonWk const& wk);

            // d'tor
            ~GaussNewton();

            // initial settings and problem definitions
            void setIOpt(IOpt const& iopt);
            void setWk(GaussNewtonWk const& wk);
            void setProblem(UserFunc* prob);
            int initialise( unsigned m,
                            Vector const& x, Vector const& xscal,
                            Vector const& fobs, Vector const& fscal,
                            Real const rtol,
                            IOpt const& iopt, GaussNewtonWk const& wk
                          );

            // start iteration
            int run();

            // statistical analysis
            int analyse();

            // get current iteration information
            Vector          getSolution();
            //Info            getInfo();
            GaussNewtonWk   getWk();
            void            printCounter();

            std::vector<Vector> getSolutionIter();

        private:
            // copy c'tor
            GaussNewton(GaussNewton const&);
            // assignmet
            GaussNewton const& operator= (GaussNewton const&);

            // helper routines
            void check_init();
            Vector call_FCN(Vector const&, int&);   // call needed as trafo (e.g. iopt.lpos) entry point
            Matrix call_JAC(Vector const&, int&);   // call needed as trafo (e.g. iopt.lpos) entry point
            void compute_scaling_xw(unsigned, bool);
            void compute_jcf_AA(Real, Real, Real, int&);
            void compute_jac_AA(Real, Real, int&);
            void compute_rank1_up_AA();
            void compute_row_scaling_A();
            void compute_qrdcmp_A(bool, unsigned&, Real&, Real&, Real&, int&);
            Vector solve_A(bool, Vector&, int&);
            Vector compute_levels_and_dx(Vector const&);
            Real wnorm(Vector const&) const;
//            Real compute_norm_projection(Vector const&);
            void log_iteration_vals1(Real);
            void log_iteration_vals2(Real, unsigned, std::string);
            void log_solout(unsigned, Vector const&, Vector const&, int, dlib::log_level const);
            void compute_statistics(Vector const&, Vector const&);

            //          _mprmon =   0      1      2      3      4       5       6
            //  dlib::log_level =  LNONE  LINFO  LVERB  LTALK  LGABBY  LDEBUG  LTRACE
            // logical units for logging; log levels by _mprerr, _mprmon, _mprsol, _mptim
            dlib::logger    _luerr, _lumon, _lusol, _lutim;
            static
            dlib::log_level _loglvl[];

            // data spaces
            IOpt            _iopt;
            GaussNewtonWk   _wk;
            Info            _info;

            // user function with methods _fun->fcn() and _fun->jac()
            UserFunc*       _fun;

            // persistent variables storing all important iteration vectors
            std::vector<Vector> _xiter;
            std::vector<Real>   _sumxall, _dlevfall, _sumxqall;
            std::vector<Real>   _tolall, _fcall;

            // working variables for one run
            unsigned        _n, _m, _mfit, _mcon;
            Vector          _x, _xscal, _fi;
            Real            _rtol, _tolmin;
            unsigned        _nitmax, _nonlin, _irank, _irankc, _ierr;
            unsigned        _maxmn;
            Matrix          _AA, _A;
            Vector          _dx, _dxq, _xa, _xw, _dxqa, _delxq, _eta;
            Vector          _fmodel, _f, _fa, _fw, _qu, _rq, _rqkeep;
            Real            _fc, _fcmin, _sigma, _fca, _cond;
            Real            _conv, _sumx, _sumxs, _dlevf;
            unsigned        _niter, _ncorr, _nfcn, _njac, _nfcnj, _nrejr1, _new;
            bool            _qbdamp;
            Matrix          _vcv;
            Vector          _xl, _xr;

            QRconDecomp     _qrA;
    };

}
#endif // __GAUSS_NEWTON_H
