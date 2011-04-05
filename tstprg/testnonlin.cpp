#include <iostream>
#include "testcases.h"
#include "parkin.h"


// --------------------------------------------------------------------------

namespace PARKIN
{

    class FCN : public UserFunc
    {
        public:
            virtual ~FCN() { }
            virtual Vector fcn(Vector const& x, unsigned& ifail);
            virtual Matrix jac(Vector const& x, unsigned& ifail);
            Vector  init(Vector& x);
        private:
            Vector _observation;
    };

    // ----------------------------------------------------------------------

    Vector FCN::fcn(Vector const& x, unsigned& ifail)
    {
        Vector  r;
        Real    h;

        r.zeros(12);

        r(1) = x(1) - x(2)*x(2) - x(3)*x(3) + 100.0;

        for (long k = 2; k <= 12; ++k)
        {
            h = ( 2.0*(k - 2.0) - x(3) ) / x(2);
            r(k) = x(1)*exp( -0.5*h*h );
        }

        ifail = 0;
        return r;
    }

    // ----------------------------------------------------------------------

    Matrix FCN::jac(Vector const& x, unsigned& ifail)
    {
        Matrix  J;
        Real    h, exph, hexph;

        J.zeros(12,3);

        J(1,1) = 1.0;
        J(1,2) = -2.0*x(2);
        J(1,3) = -2.0*x(3);

        for (long k = 2; k <= 12; ++k)
        {
            h     = ( 2.0*(k - 2.0) - x(3) ) / x(2);
            exph  = exp( -0.5*h*h );
            hexph = x(1)*exph*h / x(2);

            J(k,1) = exph;
            J(k,2) = hexph*h;
            J(k,3) = hexph;
        }

        ifail = 0;
        return J;
    }

    // ----------------------------------------------------------------------

    Vector FCN::init(Vector& x0)
    {
        Vector  xsol, fit;
        Real    h;

        x0.zeros(3);
        xsol.zeros(3);

        x0(1)   = 1.0;   x0(2)   = 2.0;  x0(3)   = 5.0;
        xsol(1) = 100.0; xsol(2) = 10.0; xsol(3) = 10.0;

        fit.zeros(11);

        for (long k = 1; k <= 11; ++k)
        {
            h      = ( 2.0*(k - 1.0) - xsol(3) ) / xsol(2);
            fit(k) = xsol(1)*exp( -0.5*h*h ) + sin(10.0*k)*1.0e-2;
        }

        return fit;
    }
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

using namespace PARKIN;

// --------------------------------------------------------------------------

int testnonlin()
{
    Real            rtol = 1.0e-4;
    GaussNewton     gn;
    FCN             problem;
    IOpt            iopt;
    GaussNewtonWk   wk;
    Vector          x, xscal, fobs, fscal;
    unsigned        ifail, m;

    iopt.mode      = 1;         // single step
    iopt.jacgen    = 1;         // use Jacobian in FCN
    iopt.qrank1    = false;     // no Broyden-rank1 updates
    iopt.nonlin    = 3;         // highly nonlinear problem
    iopt.norowscal = false;     // automatic row scaling ON
    iopt.mprmon    = 2;

    fobs = problem.init(x);
    std::cout << "\n\nSimulated experimental data:\n" << fobs << std::endl;
    xscal.zeros(x.nr());
    fscal.zeros(fobs.nr());
    m = problem.fcn(x,ifail).nr();

    gn.setProblem(&problem);
    gn.initialise(m,x,xscal,fobs,fscal,rtol,iopt,wk);

    int ierr = -1;
    while ( ierr == -1 )
    {
        ierr = gn.run();
        std::cout << "One Gauß-Newton step." << std::endl;
    }

    gn.analyse();

    // std::cout << problem.jac(x, ifail) << std::endl;

    return 0;
}
