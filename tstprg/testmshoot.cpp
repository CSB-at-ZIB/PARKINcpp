#include <iostream>
#include "testcases.h"
#include "parkin.h"


// --------------------------------------------------------------------------

namespace PARKIN
{

    class TestODE : public FirstOrderODESystem
    {
        public:
            TestODE() : _dim(3), _sigma(16.0), _R(153.083), _b(4.0) { }
            virtual ~TestODE() { }

            virtual void computeDerivatives( Real const   t,
                                              Real*          y,
                                              Real*          dy,
                                              int*          info);

            virtual void computeJacobian( Real const      t,
                                            Real*           y,
                                            Real*           dy,
                                            Real*           J,
                                            int*            ldJ,
                                            int*            full_or_band,
                                            int*            info
                                          );

            virtual void computeMassMatrix( Real const     /* t */,
                                              Real*          /* y */,
                                              Real*          B,
                                              int*           ir,
                                              int*           ic
                                            )
            {
                for (int j = 0; j < _dim; ++j)
                {
                    B[j] = 1.0;
                    ir[j] = j+1;
                    ic[j] = j+1;
                }
            }

            virtual int getSystemDimension()
            {
                return _dim;
            }

            virtual int getMassMatrixNz()
            {
                return _dim;
            }

        private:
            int     _dim;
            Real     _sigma;
            Real     _R;
            Real     _b;
    };

    // ----------------------------------------------------------------------

    void TestODE::computeDerivatives(Real const    t,
                                      Real*         y,
                                      Real*         dy,
                                      int*          info)
    {
        dy[0] = _sigma * (y[1] - y[0]);
        dy[1] = y[0] * (_R - y[2]) - y[1];
        dy[2] = y[0] * y[1] - _b + y[2];

        *info = 0;
    }

    // ----------------------------------------------------------------------

    void TestODE::computeJacobian( Real const   t,
                                    Real*        y,
                                    Real*        dy,
                                    Real*        J,
                                    int*        ldJ,
                                    int*        full_or_band,
                                    int*        info)
    {
        *info = -987;

        if (*full_or_band)
        {
            J[0] = - _sigma;
            J[1] =   _sigma;
            J[2] =   0.0;

            J[3] =  _R - y[2];
            J[4] = -1.0;
            J[5] = -y[0];

            J[6] =   y[1];
            J[7] =   y[0];
            J[8] = - _b;

            *info = 0;
        }
    }

}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

using namespace PARKIN;

// --------------------------------------------------------------------------

int testmshoot()
{
    int                n = 3, m = 11;
    Real                period = 0.95, rtol = 1.0e-4;
    MultipleShootingGN  mshoot;
    LIMEX_A             solver;
    TestODE             ode;
    IOpt                iopt = mshoot.getIOpt();
    Vector              tnodes(m);
    Matrix              X(n,m), FM(n,2);
    ODESolver::Grid     z(n);
    double             t0, t1, yIni[3];
    int                 ifail = 0;

    iopt.iauto     = 1;
    // iopt.mode      = 1;         // single step
    iopt.jacgen    = 3; // 1;         // use Jacobian in FCN
    iopt.qrank1    = false;     // no Broyden-rank1 updates
    iopt.nonlin    = 3;         // highly nonlinear problem
    iopt.norowscal = false;     // automatic row scaling ON
    iopt.mprmon    = 2;

    tnodes(1) = 0.0;

    for (int i = 2; i <= m; ++i)
    {
        tnodes(i) = tnodes(i-1) + 1.0/(m-1.0);
    }

    std::cout << "\n tnodes = " << tnodes.t() << std::endl;

    yIni[0] = z[0] = X(1,1) = 0.0;
    yIni[1] = z[1] = X(2,1) = -28.0;
    yIni[2] = z[2] = X(3,1) = 140.0;


    solver.setInterpolationFlag(1);
    solver.setODESystem(ode,0.0,z,10.0);


    for (int i = 2; i <= m; ++i)
    {
        t0 = tnodes(i-1)*period;
        t1 = tnodes(i)*period;

        solver.integrate(n, yIni, t0, t1);

        X(1,i) = yIni[0];
        X(2,i) = yIni[1];
        X(3,i) = yIni[2];
    }


    std::cout << "\n\nPeriod guess:\n" << period << std::endl;
    std::cout << "\nInitial guess: X = " << X.t() << std::endl;

    mshoot.setProblem(&ode);
    mshoot.initialise(tnodes, X, period, rtol, iopt);

    ifail = mshoot.fire();

    std::cout << "mshoot:  rc = " << ifail << std::endl;
    std::cout << "mshoot:   p = " << mshoot.getPeriod() << std::endl;

    FM =  mshoot.getFloquetMultipliers();

    std::cout << "mshoot: flq = " << FM.t() << std::endl;

    return ifail;
}
