// Copyright (C) 2010 - 2013
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2013-01-22 td
// last changed:
//
#ifndef __METAN_A_H
#define __METAN_A_H

#include "FirstOrderODESystem.h"
#include "ODESolver.h"

///

extern "C"
{
    #include "addpkg/METAN1_A/METAN1_A.h"
}

///


namespace PARKIN
{

    class METAN_A : public ODESolver
    {
        public:
            METAN_A();
            virtual ~METAN_A();

            virtual METAN_A* clone() { return new METAN_A(*this); }

            virtual void setODESystem(
                                        FirstOrderODESystem& ode,
                                        double              t0,
                                        Grid const&         y0,
                                        double              tEnd,
                                        int                 bandwidth = 0
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

            virtual Grid&           getAdaptiveGridPoints();
            virtual Trajectory&     getAdaptiveSolution();

            virtual Grid&           getSolutionGridPoints()   { return _solPoints; }
            virtual Trajectory&     getSolutionTrajectory()   { return _solution; }
            virtual Grid&           getDataGridPoints()       { return _datPoints; }
            virtual Trajectory&     getDataTrajectory()       { return _data; }
            virtual ODETrajectory*  getRawTrajectory()        { return _trajectory; }

            //

            void setODESystem(
                                FcnMetan        fcn,
                                double         t0,
                                Grid const&    y0,
                                Grid const&    refGrid,
                                double         tEnd,
                                int             bandwidth = 0
                              );

        private:
            int             _n;
            FcnMetan         _fcn;
            double          _t0;
            double*         _y0;
            double*         _y;
            double          _tEnd;
            double          _hMax;
            double          _h;
            int             _kFlag;

            SoutMetan       _sout;

            Grid            _solPoints;
            Trajectory      _solution;
            Grid            _datPoints;
            Trajectory      _data;

            ODETrajectory*  _trajectory;
    };

    ///

    class METANWrapper
    {
        public:
            static void xfcn(int* n, double* t, double* y, double* dy, int* ifail);

            static void xsout(int* n, double* tOld, double* t, double* y, double* dy);

            static void setODE(FirstOrderODESystem& ode) { _ode = &ode; }
            static void setObj(METAN_A& obj) { _obj = &obj; }

        private:
            static FirstOrderODESystem*     _ode;
            static METAN_A*                 _obj;
    };


}

#endif // __METAN_A_H
