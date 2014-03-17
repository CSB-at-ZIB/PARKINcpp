// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-02-04 td
// last changed:
//
#ifndef __LIMEX_A_H
#define __LIMEX_A_H

#include "FirstOrderODESystem.h"
#include "ODESolver.h"
/// #include "ODETrajectory.h"
/// #include "LIMEXTrajectory.h"
/// #include "CubicHermiteTrajectory.h"

///

extern "C"
{
    #include "addpkg/LIMEX4_3A/LIMEX4_3A.h"
}

///

namespace PARKIN
{
    class LIMEX_A : public ODESolver
    {
        public:
            LIMEX_A();
            // LIMEX_A(ODESolver const& other) :
            //        ODESolver(other),
            //        _trajectory(0)
            // { _solverId = ODE_SOLVER_LIMEX_A; }
            virtual ~LIMEX_A();

            virtual LIMEX_A* clone() { return new LIMEX_A(*this); }

            virtual void setODESystem(
                                        FirstOrderODESystem& ode,
                                        double              t0,
                                        Grid const&         y0,
                                        double              tEnd,
                                        int                  bandwidth = 0
                                       );

            virtual void setODESystem(
                                        FirstOrderODESystem& ode,
                                        double              t0,
                                        Grid const&         y0,
                                        Grid const&         refGrid,
                                        double              tEnd,
                                        int                 bandwidth = 0
                                       );


            virtual int integrate();
            virtual int integrate(unsigned n, double* yIni,
                                  double tLeft, double tRight);

            virtual int integrateSensitivitySystem(unsigned nDAE);
            virtual int integrateSensitivitySystem(unsigned nDAE,
                                                    unsigned n, double* yIni,
                                                    double tLeft, double tRight);

            virtual Grid&       getAdaptiveGridPoints();
            virtual Trajectory& getAdaptiveSolution();

            virtual Grid&           getSolutionGridPoints()     { return _solPoints; }
            virtual Trajectory&     getSolutionTrajectory()     { return _solution; }
            virtual Grid&           getDataGridPoints()         { return _datPoints; }
            virtual Trajectory&     getDataTrajectory()         { return _data; }
            virtual ODETrajectory*  getRawTrajectory()  { return _trajectory->clone(); }

            virtual std::string     getErrorMessage(int rc);

            ///

            int integrateWithoutInterpolation();
            int integrateWithoutInterpolation(unsigned n, double* yIni,
                                              double tLeft, double tRight);

            void setODESystem(
                                FcnLimex     fcn,
                                JacLimex     jac,
                                double      t0,
                                Grid const& y0,
                                Grid const& refGrid,
                                double      tEnd,
                                int         bandwidth = 0
                             );

            int             getSuccessiveCallFlag()         { return _iOpt[15]; }
            void            resetSuccessiveCallFlag()       { _iOpt[15] = -1; }


        private:
            int             getLimexInterpolationFlag()     { return _iOpt[31]; }
            void            setLimexInterpolationFlag(int j=1);

        private:
            void initOpt();
            void computeAndSaveLimexTrajectory(double* from, double to, double* y = 0);

            int integrateWithLIMEXInterp();
            int integrateWithLIMEXInterp(unsigned n, double* yIni, double tLeft, double tRight);
            int integrateWithCubicInterp();
            int integrateWithCubicInterp(unsigned n, double* yIni, double tLeft, double tRight);

            int integrateSensLIMEX(unsigned nnDAE);
            int integrateSensLIMEX(unsigned nnDAE, unsigned n, double* yIni, double tLeft, double tRight);
            int integrateSensCubic(unsigned nnDAE);
            int integrateSensCubic(unsigned nnDAE, unsigned n, double* yIni, double tLeft, double tRight);

            int       _n;         // size of differential algebraic system
            FcnLimex  _fcn;       // external function computing f(t,y) and B
            JacLimex  _jac;       // external function computing residual f - B * y'
            double   _t0;        // starting point of integration
            double   _T;         // end point of integration
            double*  _y0;        // initial values of solution at t0
            double*  _dy0;       // initial derivative of solution at t0; if unknown set to zero
            double   _h;         // initial stepsize guess
            int      _iOpt[32];  // input integer option field
            double   _rOpt[5];   // input real value option field (integration control parameter)
            int*     _iPos;      // field of length n: check/garantee positive solution component, if corresponding entry is set to 1
            int      _ifail[3];  // return error indication(s)
            int      _kOrder;            // output of current integration order (in single step mode)
                                        // output of Hermite interpolation tableau (max. size in FORTRAN size definition file...)
            double   _Dense[MAX_NO_EQNS*(2+MAX_ROW_TAB*(MAX_ROW_TAB+1)/2)];
            double   _t1;                // output start of current subinterval in single step mode
            double   _t2;                // output end of current subinterval in single step mode

            Grid        _solPoints;
            Trajectory  _solution;
            Grid        _datPoints;
            Trajectory  _data;

            /// LIMEXTrajectory             _trajectory;
            /// CubicHermiteTrajectory      _trajectory;
            ODETrajectory*              _trajectory;
    };

    ///

    class LIMEXWrapper
    {
        public:
            // LIMEX callback interface:
            //    'static' qualifier needed as calling parameter
            //     of FORTRAN subroutine limex_( ... ) !
            static void xfcn(  int* n, int* nz,
                                double* t, double* y, double* dy,
                                double* B, int* ir, int* ic,
                                int* info);
            static void xjac(  int* n,
                                double* t, double* y, double* dy,
                                double* J, int* ldJ,
                                int* ml, int* mu, int* full_or_band,
                                int* info);

            static void setODE(FirstOrderODESystem& ode)
            {
                _ode = &ode;
            }

        private:
            static FirstOrderODESystem* _ode;
    };

}
#endif // __LIMEX_A_H
