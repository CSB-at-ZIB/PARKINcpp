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

    param.push_back("p1");      var["p1"]  =  0.1806040385e+1;      //  1.0;
    param.push_back("p2");      var["p2"]  =  0.8935523699e+0;      //  0.494;
    param.push_back("p3");      var["p3"]  =  0.2939611589e+2;      // 29.0;
    param.push_back("p4");      var["p4"]  =  0.9217335023e+1;      //  9.0;
    param.push_back("p5");      var["p5"]  =  0.5812799976e-1;      //  0.053;
    param.push_back("p6");      var["p6"]  =  0.2428563221e+1;      //  2.0;
    param.push_back("p7");      var["p7"]  =  0.6410427954e-1;      //  0.06;
    param.push_back("p8");      var["p8"]  =  0.5560335128e+1;      //  5.0;
    param.push_back("p9");      var["p9"]  =  0.2009778604e-1;      //  0.02;
    param.push_back("p10");     var["p10"] =  0.5772851642e+0;      //  0.5;
    param.push_back("p11");     var["p11"] =  0.2150523182e+1;      //  1.19;

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
    //      A' = -r1 + r9
    //         = - (p1 * A) + (p9 * B)
    //
    //      B' = r1 - r2 - r3 + r7 - r9 + r10
    //         = (p1 * A) - (p2 * B) - (p3 * B * C) + (p7 * D) - (p9 * B) + (p10 * D * F)
    //
    //      C' = r2 - r3 - (2*r4) - r6 + r8 + r10 + (2*r11)
    //         = (p2 * B) - (p3 * B * C) - 2*(p4 * C * C) - (p6 * C) + (p8 * E) + (p10 * D * F) + 2*(p11 * E * F)
    //
    //      D' = r3 - r5 - r7 - r10
    //         = (p3 * B * C) - (p5 * D) - (p7 * D) - (p10 * D * F)
    //
    //      E' = r4 + r5 - r8 - r11
    //         = (p4 * C * C) + (p5 * D) - (p8 * E) - (p11 * E * F)
    //
    //      F' = r3 + r4 + r6 - r10 - r11
    //         = (p3 * B * C) + (p4 * C * C) + (p6 * C) - (p10 * D * F) - (p11 * E * F)
    //
    //      G' = r6 + r7 + r8
    //         = (p6 * C) + (p7 * D) + (p8 * E)

    emap["A"] = Expression(
                            PLUS,
                            Expression(MINUS, aux["r1"]),
                            aux["r9"]
                          );

    emap["B"] = Expression(
                            PLUS,
                            Expression( PLUS,
                                        Expression( PLUS,
                                                    aux["r1"],
                                                    Expression(MINUS, aux["r2"])
                                                  ),
                                        Expression( MINUS, aux["r3"] )
                                      ),
                            Expression( PLUS,
                                        Expression( PLUS,
                                                    aux["r7"],
                                                    Expression(MINUS, aux["r9"])
                                                  ),
                                        aux["r10"]
                                      )
                          );

    emap["C"] = Expression(
                            PLUS,
                            Expression( PLUS,
                                        Expression( PLUS,
                                                    aux["r2"],
                                                    Expression(MINUS, aux["r3"])
                                                  ),
                                        Expression( PLUS,
                                                    Expression(TIMES, -2.0, aux["r4"]),
                                                    Expression(MINUS, aux["r6"])
                                                  )
                                      ),
                            Expression( PLUS,
                                        Expression( PLUS.
                                                    aux["r8"],
                                                    aux["r10"]
                                                  ),
                                        Expression( TIMES, 2.0, aux["r11"] )
                                      )
                          );

    emap["D"] = Expression(
                            PLUS,
                            Expression( PLUS,
                                        aux["r3"],
                                        Expression(MIUNS, aux["r5"])
                                      ),
                            Expression( PLUS,
                                        Expression(MINUS, aux["r7"]),
                                        Expression(MINUS, aux["r10"])
                                      )
                          );

    emap["E"] = Expression(
                            PLUS,
                            Expression( PLUS,
                                        aux["r4"],
                                        aux["r5"]
                                      ),
                            Expression( PLUS,
                                        Expression(MINUS, aux["r8"]),
                                        Expression(MINUS, aux["r11"])
                                      )
                          );

    emap["F"] = Expression(
                            PLUS,
                            Expression( PLUS,
                                        Expression( PLUS,
                                                    aux["r3"],
                                                    aux["r4"]
                                                  ),
                                        aux["r6"]
                                      ),
                            Expression( PLUS,
                                        Expression(MINUS, aux["r10"]),
                                        Expression(MINUS, aux["r11"])
                                      )
                          );


    emap["G"] = Expression(
                            PLUS,
                            Expression( PLUS,
                                        aux["r6"],
                                        aux["r7"]),
                            aux["r8"]
                          );


    biosys.setODESystem(emap);
    // biosys.setSpecies(species);  // has now a complete different effect:
                                    // resets to an identity 'emap' mapping
    biosys.setParameters(param);


    biosys.setInitialValue("A", 1.0);
    biosys.setInitialValue("B", 0.0);
    biosys.setInitialValue("C", 0.0);
    biosys.setInitialValue("D", 0.0);
    biosys.setInitialValue("E", 0.0);
    biosys.setInitialValue("F", 0.0);
    biosys.setInitialValue("G", 0.0);

    Vector meastp(58);
    for (long j=1; j <= meastp.nr(); ++j) meastp(j) = tstart + j*(tend-tstart)/59.0;
    biosys.setMeasurementTimePoints( meastp );

    // Breakpoints / Event Management:
    //
    //  Subdivision of integration interval [t0 T] :  [ t0 = b1, b2, b3, ..., bn-1, bn = T ]
    //
//    Vector breaktp;
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
    vref = biosys.computeModel(var, "adaptive");
    // syndata = biosys.computeModel(var, "adaptive");

, std::cerr );
//std::cerr << std::endl;
//std::cerr << " *** Retn: biosys.computeModel() *** " << std::endl;



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
    Vector x6( biosys.getOdeTrajectory(5) );
    Vector x7( biosys.getOdeTrajectory(6) );

    std::cout << "=====================" << std::endl;
    std::cout << " t = " << std::endl;
    std::cout << tp.t();
    std::cout << "=====================" << std::endl;
    std::cout << " A [a.u] = " << std::endl;
    std::cout << x1.t();
    std::cout << "---------------------" << std::endl;
    std::cout << " B [a.u] = " << std::endl;
    std::cout << x2.t();
    std::cout << "---------------------" << std::endl;
    std::cout << " C [a.u] = " << std::endl;
    std::cout << x3.t();
    std::cout << "---------------------" << std::endl;
    std::cout << " D [a.u] = " << std::endl;
    std::cout << x4.t();
    std::cout << "---------------------" << std::endl;
    std::cout << " E [a.u] = " << std::endl;
    std::cout << x5.t();
    std::cout << "---------------------" << std::endl;
    std::cout << " F [a.u] = " << std::endl;
    std::cout << x6.t();
    std::cout << "---------------------" << std::endl;
    std::cout << " G [a.u] = " << std::endl;
    std::cout << x7.t();
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

    out.open("testparkin_aux.csv");

    if ( !out.is_open() )
    {
        std::cerr << "### ERROR: Could not open data file 'testparkin_aux.csv'.\n";
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

    p.zeros(11);
    p(1)  =  0.1000e+1;        // true:
    p(2)  =  0.4940e+0;        // true:
    p(3)  =  0.2900e+2;        // true:
    p(4)  =  0.9000e+1;        // true:
    p(5)  =  0.5300e-1;        // true:
    p(6)  =  0.2000e+1;        // true:
    p(7)  =  0.6000e-1;        // true:
    p(8)  =  0.5000e+1;        // true:
    p(9)  =  0.2000e-1;        // true:
    p(10) =  0.5000e+0;        // true:
    p(11) =  0.1190e+1;        // true:

    pscal.zeros(11);

    var.clear();

    par1.push_back( "p1" );         var[ "p1" ]  = p(1);
    par1.push_back( "p2" );         var[ "p2" ]  = p(2);
    par1.push_back( "p3" );         var[ "p3" ]  = p(3);
    par1.push_back( "p4" );         var[ "p4" ]  = p(4);
    par1.push_back( "p5" );         var[ "p5" ]  = p(5);
    par1.push_back( "p6" );         var[ "p6" ]  = p(6);
    par1.push_back( "p7" );         var[ "p7" ]  = p(7);
    par1.push_back( "p8" );         var[ "p8" ]  = p(8);
    par1.push_back( "p9" );         var[ "p9" ]  = p(9);
    par1.push_back( "p10" );        var[ "p10" ] = p(10);
    par1.push_back( "p11" );        var[ "p11" ] = p(11);

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
