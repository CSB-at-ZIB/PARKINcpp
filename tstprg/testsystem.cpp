// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-12-08 td
// last changed:
//

#include <vector>

#include "linalg/Vector.h"
#include "system/Expression.h"
#include "system/BioSystem.h"
#include "system/BioPAR.h"

#include "nonlin/GaussNewton.h"


using namespace PARKIN;


int testsystem()
{
    /// define and solve forward problem first

    Real                    tstart = 2.0;
    Real                    tend   = 8.0;
    Expression::Param       var;
    BioSystem::Species      species;
    BioSystem::Parameter    param;
    BioSystem               biosys(tstart, tend);

    std::map< std::string, Expression > exprmap;

    species.push_back("x1");
    species.push_back("x2");
    param.push_back("a");
    param.push_back("b");

    exprmap["x1"] = Expression(TIMES, "a", "x2");
    exprmap["x2"] = Expression(TIMES, Expression(MINUS, "b"), "x1");

    // biosys.setSpecies(species);
    biosys.setParameters(param);
    biosys.setODESystem(exprmap);

    biosys.setInitialValue("x1",  5.0);
    biosys.setInitialValue("x2", -3.0);

    Vector meastp(9);
    for (long j=1; j <= meastp.nr(); ++j) meastp(j) = tstart + j*(tend-tstart)/9.0;
    biosys.setMeasurementTimePoints( meastp );

    //
    //  take the following lines in to change the default values
    //  ATOL = 10*EPMACH, and RTOL = 1.0e-12
    //
    //biosys.setSolverATol(EPMACH);
    //biosys.setSolverRTol(1.0e-9);

    Vector syndata, synscal;
    for (unsigned j = 1; j <= 3; ++j)
    {
        var["a"] = 2.0/j;
        var["b"] = j;


        syndata = biosys.computeModel(var, "adaptive");
        // syndata = biosys.computeModel(var, "adaptive");


        Vector tp( biosys.getOdeTrajectoryTimePoints() );
        Vector x1( biosys.getOdeTrajectory(0) );
        Vector x2( biosys.getOdeTrajectory(1) );

        std::cout << " t = " << std::endl;
        std::cout << tp.t() << std::endl;

        std::cout << " sol = " << std::endl;
        std::cout << x1.t() << std::endl;
        std::cout << x2.t() << std::endl;
        std::cout << "=====================" << std::endl;
    }

    /// and, subsequently, prepare and solve for the inverse problem

    Real   rtol = 1.0e-5;
    Vector p, pscal;
    // Vector syndata, synscal;

    // syndata = biosys.getMeasurements();

    std::cout << std::endl;
    std::cout << "Simulated  experimental  data  :  nr() = " << syndata.nr() << std::endl;
    std::cout << syndata << std::endl;

    p.zeros(2);
    p(1) = 1.42; // 2.0/p(2);  // initial guess: very far away from 'true' solution [2/3 3]
    p(2) = 0.8;
    pscal.zeros(2);
    synscal.zeros( syndata.nr() );

    IOpt           iopt;
    GaussNewtonWk  wk;
    GaussNewton    gn;
    BioPAR         prob( &biosys ); // , param );

    iopt.mode      = 0;   // 0:normal run, 1:single step
    iopt.jacgen    = 1;   // 1:user supplied Jacobian, 2:num.diff., 3:num.diff.(with feedback)
    iopt.qrank1    = false;   // allow Broyden rank-1 updates if __true__
    iopt.nonlin    = 4;   // 1:linear, 2:mildly nonlin., 3:highly nonlin., 4:extremely nonlin.
    iopt.norowscal = false;   // allow for automatic row scaling of Jacobian if __false__
    iopt.mprmon    = 2;
    iopt.mprerr    = 1;

    gn.setProblem( &prob );
    gn.initialise( syndata.nr(), p, pscal, syndata, synscal, rtol, iopt, wk );

    gn.run();
    gn.analyse();
    gn.printCounter();

    std::cout << std::endl;
    std::cout << "Inv.Prob. solution = " << std::endl;
    std::cout << gn.getSolution().t() << std::endl;

    return 0;
}
