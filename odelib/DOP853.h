// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-12-07 td
// last changed:
//
#ifndef __DOP853_H
#define __DOP853_H

#include <cstdio>
#include <vector>
// #include <sstream>

#include "ODESolver.h"

extern "C"
{
    #include "addpkg/Ode/dop853b.h"
}


namespace PARKIN
{

    class DOP853 : public ODESolver
    {
        public:
            DOP853();
            virtual ~DOP853();

            virtual int         integrate();
            virtual int         integrate(unsigned n, double* yIni,
                                          double xLeft, double Right);
            virtual Grid&       getAdaptiveGridPoints();
            virtual Trajectory& getAdaptiveSolution();

            /// void    fcn (unsigned n, double x, double* y, double* f, double* cd);
            /// void solout (long nr, double xold, double x, double* y, unsigned n, int* irtrn);
            /// {
            ///     double contd8(unsigned j, double s);    // for (0 <= j < n)  if  (xold < s < x)
            ///
            ///     ...
            /// }

            void setODESystem(
                                FcnEqDiff       fcn,
                                double          x0,
                                Grid const&     y0,
                                Grid const&     refGrid,
                                double          xend
                              ); // , SolTrait solout);

            Grid&        getSolutionGridPoints() { return _solPoints; }
            Trajectory&  getSolutionTrajectory() { return _solution; }
            Grid&        getDataGridPoints()     { return _datPoints; }
            Trajectory&  getDataTrajectory()     { return _data; }

        private:
            unsigned    _n;             // dimension of the system
            FcnEqDiff   _fcn;           // function computing the values of f(x,y)
            double      _x0;            // initial x-value
            double*     _y0;            // initial values for y
            double*     _y;             // values for y (work space for integrator)
            double      _xend;          // final x-value (the difference (xend - x) may be positive or negative)
            int         _itol;          // switch for _rtol/_atol: scalar (itol==0), vector (itol==1)

            SolTrait    _solout;        // call-back function providing the numerical solution during integration
            int         _iout;          // switch for controlling the calls of _solout
            FILE*       _fileout;       // message stream

            double      _uround;        // rounding unit
            double      _safe;          // safety factor
            double      _fac1, _fac2;   // parameters for step size selection
            double      _beta;          // for stabilised step size control
            double      _hmax;          // maximal step size
            double      _h;             // initial step size
            long        _nmax;          // maximal number of allowed steps
            int         _meth;          // switch for the choice of the coefficients
            long        _nstiff;        // test for stiffness

            unsigned    _nrdens;        // number of components for which dense output is required
            unsigned*   _icont;         // indices of components for which dense output is required, >= nrdens
            unsigned    _licont;        // declared length of _icont
            double*     _cd;            // client data

            Grid        _solPoints;     // computed, adaptive grid by the integrator
            Trajectory  _solution;      // computed solution trajectories

            Grid        _datPoints;     // the measurement grid points (given!), sorted by setODESystem()
            Trajectory  _data;          // the (interpolated) simulated data
    };

    ///

    class DOP853Wrapper
    {
            typedef ODESolver::GridIterConst GridIterConst;

        public:
            static void solout(long nr, double xold, double x, double* y, unsigned n, int* irtrn);

            static void setObj(DOP853& obj) { _obj = &obj; }

        private:
            static DOP853* _obj;
    };


}
#endif // __DOP853_H
