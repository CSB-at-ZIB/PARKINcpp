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


int testsystem_aux()
{
    /// define and solve forward problem first

    Real                    tstart = 1.0;
    Real                    tend   = 4.0;
    Expression::Param       var;

    BioSystem::Species      species;
    BioSystem::Parameter    param;
    BioSystem::Parameter    par1;
    // BioSystem::Auxiliary    aux;
    BioSystem               biosys(tstart, tend);

    BioSystem::ExpressionMap aux;
    BioSystem::ExpressionMap emap;

    species.push_back("s1");
    species.push_back("s2");

    param.push_back("compartment");
    param.push_back("k1");

    aux["react1"] = Expression( TIMES, "compartment", Expression(TIMES, "k1", "s1") );
    // exprlist.push_back( Expression(MINUS, Expression(TIMES, "compartment", Expression(TIMES, "k1", "s1") ) ) );
    // exprlist.push_back( Expression(TIMES, "compartment", Expression(TIMES, "k1", "s1") ) );
    emap["s1"] = Expression( MINUS, aux["react1"] );
    emap["s2"] = aux["react1"] ;


    // biosys.setSpecies(species);
    biosys.setParameters(param);
    //biosys.setAuxiliary(aux);

    //biosys.setAuxExpressions(auxilist);
    biosys.setODESystem(emap);

    biosys.setInitialValue("s1", 1.5e-4);
    biosys.setInitialValue("s2", 0.0);

    Vector meastp(10);
    for (long j=1; j <= meastp.nr(); ++j) meastp(j) = tstart + j*(tend-tstart)/10.0;
    biosys.setMeasurementTimePoints( meastp );

//    Vector breaktp(4);
//    for (long j=1; j <= breaktp.nr(); ++j) breaktp(j) = tstart + (j-1)*(tend-tstart)/3.0;
//    biosys.setBreakpoints( breaktp );

    //
    //  take the following lines in to change the default values
    //  ATOL = 10*EPMACH, and RTOL = 1.0e-12
    //
    biosys.setSolverRTol(1.0e-5);
    biosys.setSolverATol(1.0e-8);

    Vector syndata, synscal;
//    for (unsigned j = 1; j <= 3; ++j)
//    {
    var["compartment"] = 1.0;
    var["k1"] = 1.5;


    biosys.computeModel(var, "init");
    // syndata = biosys.computeModel(var, "init");


    Vector tp( biosys.getOdeTrajectoryTimePoints() );
    Vector x1( biosys.getOdeTrajectory(0) );
    Vector x2( biosys.getOdeTrajectory(1) );

    std::cout << " t = " << std::endl;
    std::cout << tp.t() << std::endl;

    std::cout << " sol = " << std::endl;
    std::cout << x1.t() << std::endl;
    std::cout << x2.t() << std::endl;
    std::cout << "=====================" << std::endl;
//    }

    /// and, subsequently, prepare and solve for the inverse problem

    BioSystem::MeasurementList  measlist = biosys.getMeasurementList();
    BioSystem                   invBiosys( tstart, tend );

    // invBiosys.setSpecies(species);
    invBiosys.setParameters(param);

    invBiosys.setODESystem(emap);
    invBiosys.setInitialValue( species[0], biosys.getInitialValue(species[0]) );
    invBiosys.setInitialValue( species[1], biosys.getInitialValue(species[1]) );

    invBiosys.setParamValue( param[0], var[param[0]] );
    invBiosys.setParamValue( param[1], var[param[1]] );

    invBiosys.setMeasurementList( meastp, measlist );

    //
    // From this point on only the new object "invBiosys" is referenced
    //

    invBiosys.setSolverRTol(1.0e-6);
    invBiosys.setSolverATol(1.0e-8);

    Real   rtol = 1.0e-5;
    Vector p, pscal;
    // Vector syndata, synscal;


    syndata = invBiosys.getMeasurements();
    synscal.zeros( syndata.nr() );

    std::cout << std::endl;
    std::cout << "Simulated  experimental  data  :  nr() = " << syndata.nr() << std::endl;
    std::cout << syndata << std::endl;

    p.zeros(2);
    p(1) = 0.8; // 2.0/p(2);  // initial guess: very far away from 'true' solution [2/3 3]
    p(2) = 2.0;
    pscal.zeros(2);

    par1.push_back( "compartment" );
    par1.push_back( "k1" );

    IOpt           iopt;
    GaussNewtonWk  wk;
    GaussNewton    gn;
    BioPAR         prob( &invBiosys, par1 );

    iopt.mode      = 0;   // 0:normal run, 1:single step
    iopt.jacgen    = 3;   // 1:user supplied Jacobian, 2:num.diff., 3:num.diff.(with feedback)
    iopt.qrank1    = false;   // allow Broyden rank-1 updates if __true__
    iopt.nonlin    = 4;   // 1:linear, 2:mildly nonlin., 3:highly nonlin., 4:extremely nonlin.
    iopt.norowscal = false;   // allow for automatic row scaling of Jacobian if __false__
    iopt.mprmon    = 2;
    iopt.mprerr    = 1;


    // wk.nitmax = 2;
    wk.cond = 1.0 / 1.0e-3;

    gn.setProblem( &prob );
    gn.initialise( syndata.nr(), p, pscal, syndata, synscal, rtol, iopt, wk );

    gn.run();
    // gn.analyse();
    gn.printCounter();

    std::cout << std::endl;
    std::cout << "Inv.Prob. solution = " << std::endl;
    std::cout << gn.getSolution().t() << std::endl;

    Expression::Param final = invBiosys.getSysPar();

    std::cout << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << "Final parameter value(s)" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << std::setw(20) << "Parameter Name" << "   " << "Value" << std::endl;
    std::cout << std::setw(20) << "--------------" << "   " << "-----" << std::endl;
    for (unsigned j = 0; j < param.size(); ++j)
        std::cout << std::setw(20) << param[j] << "   " << final[param[j]] << std::endl;

    return 0;
}
