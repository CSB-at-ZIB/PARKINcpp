// Copyright (C) 2010 - 2012
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2012-10-02 td
// last changed:
//
#ifndef __MULTIPLE_SHOOTING_GN_H
#define __MULTIPLE_SHOOTING_GN_H

#include <vector>

#include "common/Types.h"
#include "common/Constants.h"
#include "common/PARKINLog.h"

#include "linalg/QRconDecomp.h"

#include "odelib/FirstOrderODESystem.h"
#include "odelib/ODESolver.h"
//#include "odelib/LIMEX_A.h"

#include "UserFunc.h"  // for decl Iopt()

namespace PARKIN
{

    class MultipleShootingGN
    {
        public:
            // c'tor
            MultipleShootingGN();

            // d'tor
            ~MultipleShootingGN();


            void setLogStream(std::ostream& out);
            // initial settings and problem definitions
            void setIOpt(IOpt const& iopt);
            void setProblem(FirstOrderODESystem* ode);
            int initialise( Vector const& tnodes, Matrix const& X,
                             Real const period,
                             Real const rtol, IOpt const& iopt
                          );

            // start Gauss-Newton iteration
            int fire();

            // get results (after iteration)
            Vector  getNodes();
            Matrix  getSolution();

            Real    getPeriod();
            Matrix  getFloquetMultipliers();

            //
            IOpt    getIOpt();

        private:
            // copy c'tor
            MultipleShootingGN(MultipleShootingGN const&);
            // assignment
            MultipleShootingGN const& operator= (MultipleShootingGN const&);

            // helper routines
            void check_init();
            Matrix call_IVPSOL(Vector const&, Real, Matrix const&, Matrix&, int&);
            Vector call_FCN(Real, Vector const&, int&);
            Real scalarprod(Matrix const&, Matrix const&, Real, Real, int);
            void compute_scaling_XW(Real&);
            void compute_fd_derivative_G(Matrix const&, Real, int&);
            void solve_var_eq_for_G(Matrix const&, Real, int&);
            void compute_rank1_update_G(int);
            Matrix multiply_G(Vector const&);
            Vector compute_condensed_RHSp(unsigned, Matrix const&,
                                          Vector const&, Vector const&);
            void compute_rest_of(Matrix&, Real, unsigned, Matrix const&, int);
            Real compute_refinement_sweep(bool,int,int,
                                          Matrix&,Real&,Vector const&,
                                          Real,int&
                                         );
            void compute_qrdcmp_E(bool,Real&,int&);
            Vector solve_E(bool,Vector&,int&);
            void log_iteration_vals1(Real);
            void log_iteration_vals2();

            //
            void exchange(unsigned&,unsigned&,unsigned&,unsigned&,Matrix&);
            //
            void balance(Matrix&,unsigned&,unsigned&,Vector&);
            void orthes(unsigned,unsigned,Matrix&,Vector&);
            void hqr(unsigned,unsigned,Matrix&,Vector&,Vector&,int&);

            //          _mprmon =   0      1      2      3      4       5       6
            //  dlib::log_level =  LNONE  LINFO  LVERB  LTALK  LGABBY  LDEBUG  LTRACE
            // logical units for logging; log levels by _mprerr, _mprmon, _mprsol
            dlib::logger            _luerr, _lumon, _lusol;
            static
            dlib::log_level         _loglvl[];

            // data spaces
            IOpt                    _iopt;

            // user ODE system with (compulsory!) methods
            //                  _ode->computeDerivatives() and
            //                  _ode->computeJacobian()
            FirstOrderODESystem*    _ode;

            // working variables for one salvo (i.e. run)
            unsigned                _m, _n;
            Vector                  _tnodes;
            Real                    _period, _dp, _dpq, _dpqa;
            Real                    _perioda, _periodw, _ptg;
            Vector                  _qu;
            Matrix                  _X, _XU, _DX, _DXQ, _DXQA, _XA, _XW;
            Matrix                  _HH, _HHA, _FP, _FPA, _XTG, _BG;
            Matrix                  _E, _FM;

            unsigned                _niter, _ny, _ierr;
            unsigned                _kount, _new, _nitmax, _nymax;
            unsigned                _irank, _nonlin;
            Real                    _fc, _fcmin, _fca, _cond, _sigma;
            Real                    _sumx, _sumf, _conv;
            Real                    _eps, _tolf, _tolj, _tolmin, _xthr;
            Real                    _reldif;

            Matrix*                 _G;
            QRDecomp                _qrE;
            // QRconDecomp     _qrE;
            ODESolver*              _odeSolver;
    };

}

#endif // __MULTIPLE_SHOOTING_GN_H
