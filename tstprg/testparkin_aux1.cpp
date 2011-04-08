// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-12-08 td
// last changed:
//

#include <cmath>
#include <vector>

#include "addpkg/dlib/time_this.h"

#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "system/Expression.h"
#include "system/BioSystem.h"
#include "system/BioPAR.h"

#include "nonlin/YeOldeParkinCore.h"


using namespace PARKIN;


int testparkin_aux()
{
    srand(0);

    /// define and solve forward problem first
    ///
    /// taken from old PARK11 example files / cases
    /// pyridine: chemin.pyex (chemin.xml)

    Real                        tstart = 0.0;
    Real                        tend   = 6.0;
    Expression::Param           var, par;

    BioSystem::Species          species;
    BioSystem::Parameter        param;
    BioSystem::Parameter        par1;
    BioSystem                   biosys(tstart, tend);

    BioRHS::ExpressionMap    aux;
    BioRHS::ExpressionMap    emap;

    //

    species.push_back("A");
    species.push_back("B");
    species.push_back("C");
    species.push_back("D");
    species.push_back("E");
    species.push_back("F");
    species.push_back("G");

    //

    param.push_back("p1");      var["p1"]  =  1.0;
    param.push_back("p2");      var["p2"]  =  0.494;
    param.push_back("p3");      var["p3"]  = 29.0;
    param.push_back("p4");      var["p4"]  =  9.0;
    param.push_back("p5");      var["p5"]  =  0.053;
    param.push_back("p6");      var["p6"]  =  2.0;
    param.push_back("p7");      var["p7"]  =  0.06;
    param.push_back("p8");      var["p8"]  =  5.0;
    param.push_back("p9");      var["p9"]  =  0.02;
    param.push_back("p10");     var["p10"] =  0.5;
    param.push_back("p11");     var["p11"] =  1.19;

    //

    // Reaction / rule:
    //
    //      r1  :  p1 * A
    //      r2  :  p2 * B
    //      r3  :  p3 * B * C
    //      r4  :  p4 * C * C
    //      r5  :  p5 * D
    //      r6  :  p6 * C
    //      r7  :  p7 * D
    //      r8  :  p8 * E
    //      r9  :  p9 * B
    //      r10 :  p10 * D * F
    //      r11 :  p11 * E * F

    aux["r1"] = Expression(TIMES, "p1", "A");

    aux["r2"] = Expression(TIMES, "p2", "B");

    aux["r3"] = Expression(TIMES, "p3", Expression(TIMES, "B", "C"));

    aux["r4"] = Expression(TIMES, "p4", Expression(TIMES, "C", "C"));

    aux["r5"] = Expression(TIMES, "p5", "D");

    aux["r6"] = Expression(TIMES, "p6", "C");

    aux["r7"] = Expression(TIMES, "p7", "D");

    aux["r8"] = Expression(TIMES, "p8", "E");

    aux["r9"] = Expression(TIMES, "p9", "B");

    aux["r10"] = Expression(TIMES, "p10", Expression(TIMES, "D", "F"));

    aux["r11"] = Expression(TIMES, "p11", Expression(TIMES, "E", "F"));

    //

    // ODEs:
    //
    //      C' = ( 1*r1) + (-1*r2) + (-1*r3) + (-1*r8) + ( 1*r9) + ( 1*r10)
    //         = (vi) + (- C * k1 * X * pow(K5 + C,-1)) + (- C * kd) +
    //           (- a1 * C * Y) + (a2 * Z) + (alpha * d1 * Z)
    //
    //      X' = ( 1*r6) + (-1*r7)
    //         = (V3p * M * (1 - X) * pow(K3 - X + 1, -1)) + (- V4 * X * pow(K4 + X, -1))
    //
    //      M' = ( 1*r4) + (-1*r5)
    //         = (C * V1p * (1 - M) * pow(K6 + C, -1) * pow(K1 - M + 1, -1)) +
    //           (- M * V2 * pow(K2 + X, -1))
    //
    //      Y' = (-1*r8) + ( 1*r9) + ( 1*r11) + ( 1*r12) + (-1*r13)
    //         = (- a1 * C * Y) + (a2 * Z) + (alpha *kd * Z) + (vs) + (- d1 * Y)
    //
    //      Z' = ( 1*r8) + (-1*r9) + (-1*r10) + (-1*r11)
    //           (a1 * C * Y) + (- a2 * Z) + (- alpha * d1 * Z) + (- alpha * kd * Z)

    emap["C"] = Expression(
                            PLUS,
                            Expression(TIMES, 1.0, aux["r1"]),
                            Expression(
                                PLUS,
                                Expression(MINUS, aux["r2"]),
                                Expression(
                                    PLUS,
                                    Expression(MINUS, aux["r3"]),
                                    Expression(
                                        PLUS,
                                        Expression(MINUS, aux["r8"]),
                                        Expression(
                                            PLUS,
                                            Expression(TIMES, 1.0, aux["r9"]),
                                            Expression(TIMES, 1.0, aux["r10"])
                                        )
                                    )
                                )
                            )
                        );

    emap["X"] = Expression(
                            PLUS,
                            Expression(TIMES, 1.0, aux["r6"]),
                            Expression(MINUS, aux["r7"])
                        );

    emap["M"] = Expression(
                            PLUS,
                            Expression(TIMES, 1.0, aux["r4"]),
                            Expression(MINUS, aux["r5"])
                        );

    emap["Y"] = Expression(
                            PLUS,
                            Expression(MINUS, aux["r8"]),
                            Expression(
                                PLUS,
                                Expression(TIMES, 1.0, aux["r9"]),
                                Expression(
                                    PLUS,
                                    Expression(TIMES, 1.0, aux["r11"]),
                                    Expression(
                                        PLUS,
                                        Expression(TIMES, 1.0, aux["r12"]),
                                        Expression(MINUS, aux["r13"])
                                    )
                                )
                            )
                        );

    emap["Z"] = Expression(
                            PLUS,
                            Expression(TIMES, 1.0, aux["r8"]),
                            Expression(
                                PLUS,
                                Expression(MINUS, aux["r9"]),
                                Expression(
                                    PLUS,
                                    Expression(MINUS, aux["r10"]),
                                    Expression(MINUS, aux["r11"])
                                )
                            )
                        );


    biosys.setODESystem(emap);
    // biosys.setSpecies(species);  // has now a complete different effect:
                                    // resets to an identity 'emap' mapping
    biosys.setParameters(param);


    biosys.setInitialValue("C", 0.0);
    biosys.setInitialValue("X", 0.0);
    biosys.setInitialValue("M", 0.0);
    biosys.setInitialValue("Y", 1.0);
    biosys.setInitialValue("Z", 1.0);

    Vector meastp(58);
    for (long j=1; j <= meastp.nr(); ++j) meastp(j) = tstart + j*(tend-tstart)/59.0;
    biosys.setMeasurementTimePoints( meastp );

    // Breakpoints / Event Management:
    //
    //  Subdivision of integration interval [t0 T] :  [ t0 = b1, b2, b3, ..., bn-1, bn = T ]
    //
    Vector breaktp;
//    breaktp.zeros(4);
//    for (long j=1; j <= breaktp.nr(); ++j) breaktp(j) = tstart + (j-1)*(tend-tstart)/3.0;
//    biosys.setBreakpoints( breaktp );

    // and succinctly some event handling
    //
/*
    emap.clear();
    emap = biosys.getEvent(1);

    emap[species[3]] = Expression(0.5);

    biosys.setEvent(1, emap);


    std::cout << std::endl;
    for (long j = 0; j < breaktp.nr(); ++j)
    {
        BioRHS::ExpressionMap icmap = biosys.getEvent(j);
        BioRHS::EMapIterConst cBeg = icmap.begin();
        BioRHS::EMapIterConst cEnd = icmap.end();

        std::cout << "Event " << std::setw(2) << j+1 <<
                     "   (at time point t = " << breaktp(j+1) << ")" << std::endl;
        std::cout << "--------" << std::endl;
        for (BioRHS::EMapIterConst it = cBeg; it != cEnd; ++it)
        {
            std::cout << it->first << " :  " << it->second << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
*/

    emap.clear();
    emap = biosys.getODEExpr();

    for (unsigned j = 0; j < emap.size(); ++j)
    {
        std::cout << std::endl;
        std::cout << species[j] << " :  "  << emap[species[j]] << std::endl;
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "--------------------------" << std::endl;
    std::cout << "Current parameter value(s)" << std::endl;
    std::cout << "--------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << std::setw(20) << "Parameter Name" << "   " << "Value" << std::endl;
    std::cout << std::setw(20) << "--------------" << "   " << "-----" << std::endl;
    for (unsigned j = 0; j < param.size(); ++j)
        std::cout << std::setw(20) << param[j] << "   " << var[param[j]] << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    //
    //  take the following lines in to change the default values
    //  ATOL = 10*EPMACH, and RTOL = 1.0e-12
    //
    biosys.setSolverRTol(1.0e-5);
    biosys.setSolverATol(1.0e-8);


    Vector syndata, synscal, vref, utmp;

TIME_THIS_TO( std::cerr << " *** Call: biosys.computeModel() *** " << std::endl;

//for (long n=0; n < 100; ++n)
    vref = biosys.computeModel(var, "init");
    // syndata = biosys.computeModel(var, "init");

, std::cerr );
//std::cerr << std::endl;
//std::cerr << " *** Retn: biosys.computeModel() *** " << std::endl;



    Matrix Jac;
    long   n = species.size();
    long   T = meastp.nr();
    long   j = 0;

    utmp.zeros(3);

    par["K3"] = var["K3"];      utmp(1) = ( par["K3"] );
    par["K4"] = var["K4"];      utmp(2) = ( par["K4"] );
    par["K5"] = var["K5"];      utmp(3) = ( par["K5"] );

    Jac.zeros( n*T, 7 );


TIME_THIS_TO( std::cerr << " *** Call: biosys.computeJacobian() *** " << std::endl;

    Jac.set_colm( 1, 3 ) =
        biosys.computeJacobian( par ) * utmp.diag();  // ... * exp( u=log(par) ).diag()

, std::cerr )


    for (Expression::ParamIterConst it = par.begin(); it != par.end(); ++it)
    {
        Real   h = sqrt(1e-10), ptmp;
        Vector vtmp;

        ptmp = par[it->first];
        par[it->first] = exp( log(ptmp) + h );
        vtmp = biosys.computeModel(par);
        par[it->first] = ptmp;

        Jac.set_colm( 5+j++ ) = (1.0/h)*(vtmp - vref);
    }

    std::cout << " Jacobian Jac (" << Jac.nr() << " x " << Jac.nc() << ") = " << std::endl;
    std::cout << Jac << std::endl;



// exit(42);



    std::cout << std::endl;
    std::cout << " -=-=-=-=-=-=-=- " << std::endl;

    emap.clear();
    emap = biosys.getVarExpr();

    for (unsigned jj = 1, j = 0; j < species.size(); ++j)
    {
        std::string eqn = species[j];

        for (unsigned k = 0; k < species.size(); ++k)
        {
            std::string s = eqn + "_" + species[k];
            std::cout << s << " :  ";
            std::cout << "( " << jj++ << " / " << emap.size() << " )" << std::endl;
            std::cout << emap[s];
            std::cout << std::endl;
        }

        for (unsigned k = 0; k < param.size(); ++k)
        {
            std::string s = eqn + "_" + param[k];
            std::cout << s << " :  " << std::endl;
            std::cout << "( " << jj++ << " / " << emap.size() << " )" << std::endl;
            std::cout << emap[s];
            std::cout << std::endl;
        }
    }

    emap.clear();
    emap = biosys.getODEExpr();

    std::cout << std::endl;
    std::cout << emap[species[species.size()-1]].df( param[param.size()-1] );
    std::cout << std::endl;

    std::cout << " -=-=-=-=-=-=-=- " << std::endl;


    Vector tp( biosys.getOdeTrajectoryTimePoints() );
    Vector x1( biosys.getOdeTrajectory(0) );
    Vector x2( biosys.getOdeTrajectory(1) );
    Vector x3( biosys.getOdeTrajectory(2) );
    Vector x4( biosys.getOdeTrajectory(3) );
    Vector x5( biosys.getOdeTrajectory(4) );

    std::cout << "=====================" << std::endl;
    std::cout << " t = " << std::endl;
    std::cout << tp.t();
    std::cout << "=====================" << std::endl;
    std::cout << " C [a.u] = " << std::endl;
    std::cout << x1.t();
    std::cout << "---------------------" << std::endl;
    std::cout << " X [a.u] = " << std::endl;
    std::cout << x2.t();
    std::cout << "---------------------" << std::endl;
    std::cout << " M [a.u] = " << std::endl;
    std::cout << x3.t();
    std::cout << "---------------------" << std::endl;
    std::cout << " Y [a.u] = " << std::endl;
    std::cout << x4.t();
    std::cout << "---------------------" << std::endl;
    std::cout << " Z [a.u] = " << std::endl;
    std::cout << x5.t();
    std::cout << "=====================" << std::endl;

    ///
    /// and, subsequently, prepare and solve for the inverse problem
    ///

    Real   rtol       = 1.0e-5;
    Real   solverRTol = 1.0e-5;
    Real   solverATol = 1.0e-8;
    Vector p, pscal;
    // Vector syndata, synscal;

    BioSystem::MeasurementList  measlist = biosys.getMeasurementList();
    BioSystem                   invBiosys( tstart, tend );

    invBiosys.setSolverRTol( solverRTol );
    invBiosys.setSolverATol( solverATol );

    invBiosys.setODESystem(emap);
    // invBiosys.setSpecies(species);       // see above comment!
    invBiosys.setParameters(param);

    for (unsigned j = 0; j < species.size(); ++j)
        invBiosys.setInitialValue( species[j], biosys.getInitialValue(species[j]) );

    for (unsigned j = 0; j < param.size(); ++j)
        invBiosys.setParamValue( param[j], biosys.getParamValue(param[j]) );

    //
    // Some kind of bootstrapping ...
    //

    BioSystem::MeasurementList  newMeas;
    Vector                      newtp;
    Real                        sigma = 0.0; // 0.015;   // add 1.5% white noise ...
    unsigned                    sampleSize = 1;

    newMeas.clear();
    newtp.zeros( sampleSize*meastp.nr() );

    for (unsigned tp = 0; tp < (unsigned) meastp.nr(); ++tp )
    {
        for (unsigned k = 0; k < sampleSize; ++k) {
            newtp( sampleSize*tp + k + 1 ) = meastp(tp+1);
            newMeas.push_back( MeasurementPoint() );
        }

        for (unsigned j = 0; j < species.size(); ++j)
        {
            double val = measlist[tp][species[j]].first;

            for (unsigned k = 0; k < sampleSize; ++k)
            {
                Real newval = val + sigma * randn();

                newMeas[sampleSize*tp+k][species[j]] =
                    std::make_pair<Real,Real>( newval, 1.0 );
            }
        }
    }

    //
    // From this point on only the new object "invBiosys" is referenced
    //

    invBiosys.setMeasurementList( newtp, newMeas );

    meastp = newtp; // invBiosys.getMeasurementTimePoints();
    measlist = invBiosys.getMeasurementList();

    syndata = invBiosys.getMeasurements();
    synscal.ones( syndata.nr() );
    synscal = sigma * synscal;

    std::cout << std::endl;
    std::cout << "Simulated  experimental  data  :  nr() = " << syndata.nr() << std::endl;
    // std::cout << syndata << std::endl;

    std::ofstream out;

    out.open("testsystem_aux.csv");

    if ( !out.is_open() )
    {
        std::cerr << "### ERROR: Could not open data file 'testsystem_aux.csv'.\n";
        return -3;
    }

    std::cout << "\n";
    std::cout << "Timepoint [s]";
    out << "Timepoint [s]";
    for (unsigned j = 0; j < species.size(); ++j)
    {
        std::cout << "\t";
        std::cout << species[j] << " [a.u.]";
        out << "\t";
        out << species[j] << " [a.u.]";
    }
    std::cout << "\n";
    out << "\n";
    for (unsigned tp = 0; tp < (unsigned)meastp.nr(); ++tp)
    {
        std::cout << std::fixed << std::setw(10) << meastp(tp+1);
        out << std::fixed << std::setw(10) << meastp(tp+1);
        for (unsigned j = 0; j < species.size(); ++j)
        {
            std::cout << "\t";
            std::cout << std::scientific << measlist[tp][species[j]].first;
            out << "\t";
            out << std::scientific << measlist[tp][species[j]].first;
        }
        std::cout << "\n";
        out << "\n";
    }

    out.close();

//    invBiosys.computeModel(var);
//
//    std::cout << "####### invBiosys.computeModel() results #######" << std::endl;
//    for (unsigned j = 0; j < species.size(); ++j)
//    {
//        std::cout << "  trajector for \"" << species[j] << "\" = " << std::endl;
//        std::cout << invBiosys.getSimTrajectoryPoints(species[j]).t() << std::endl;
//    }
//    std::cout << "################################################" << std::endl;


    //
    // Initial guess for GaussNewton
    //

    p.zeros(3);
    p(1) = 0.5;         // true: 0.2
    p(2) = 0.5;         // true: 0.1
    p(3) = 0.05;        // true: 0.02
    pscal.zeros(3);

    var.clear();

    par1.push_back( "K3" );         var[ "K3" ] = p(1);
    par1.push_back( "K4" );         var[ "K4" ] = p(2);
    par1.push_back( "K5" );         var[ "K5" ] = p(3);

    IOpt                iopt;
    YeOldeParkinWk      wk;
    YeOldeParkinCore    parkin;
    BioPAR              prob( &invBiosys, par1 );

    iopt.mode      = 0;     // 0:normal run, 1:single step
    // iopt.jacgen    = 3;     // 1:user supplied Jacobian, 2:num.diff., 3:num.diff.(with feedback)
    // iopt.qrank1    = false;     // allow Broyden rank-1 updates if __true__
    // iopt.nonlin    = 4;     // 1:linear, 2:mildly nonlin., 3:highly nonlin., 4:extremely nonlin.
    // iopt.norowscal = false;     // allow for automatic row scaling of Jacobian if __false__
    iopt.lpos      = false;      // force solution vector to be positive (all components > 0.0)
                            //          _mprmon =   0      1      2      3      4       5       6
                            //  dlib::log_level =  LNONE  LINFO  LVERB  LTALK  LGABBY  LDEBUG  LTRACE
    iopt.mprmon    = 3;
    iopt.mprerr    = 1;


    // if ( iopt.jacgen > 1 )
        wk.cond = 1.0 / sqrt(solverRTol);
    // wk.itmax = 150;


    parkin.setProblem( &prob );
    parkin.initialise( syndata.nr(), p, pscal, syndata, synscal, rtol, iopt, wk );

    parkin.run();
    // parkin.analyse();
    // parkin.printCounter();


    Expression::Param final = invBiosys.getSysPar();

    std::cout.unsetf( std::ios_base::floatfield );
    std::cout << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << "Final parameter value(s)" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << std::setw(20) << std::right << "Parameter Name" << "      " <<
                 std::setw(15) << std::left << "Value" <<
                 std::setw(20) << std::left << "Initial" <<
                 std::endl;
    std::cout << std::setw(20) << std::right << "--------------" << "      " <<
                 std::setw(15) << std::left << "-----" <<
                 std::setw(20) << std::left << "-------" <<
                 std::endl;
    for (unsigned j = 0; j < param.size(); ++j)
        std::cout << std::setw(20) << std::right << param[j] << "      " <<
                     std::setw(15) << std::left << final[param[j]] <<
                     std::setw(20) << std::left << var[param[j]] <<
                     std::endl;

    return 0;
}
