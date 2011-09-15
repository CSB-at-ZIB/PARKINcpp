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
    ///
    /// taken from SBML example files / cases
    /// cyclin: BIOMD0000000008.xml

    Real                        tstart = 1.0;
    Real                        tend   = 100.0;
    Expression::Param           var;

    BioSystem::Species          species;
    BioSystem::Parameter        param;
    BioSystem::Parameter        par1;
    BioSystem                   biosys(tstart, tend);

    BioSystem::ExpressionMap    aux;
    BioSystem::ExpressionMap    emap;

    //

    species.push_back("C");
    species.push_back("X");
    species.push_back("M");
    species.push_back("Y");
    species.push_back("Z");

    //

    // param.push_back("V1");      var["V1"]  = 0.0;
    param.push_back("V1p");     var["V1p"] = 0.75;
    param.push_back("K1");      var["K1"]  = 0.1;

    param.push_back("V2");      var["V2"]  = 0.25;
    param.push_back("K2");      var["K2"]  = 0.1;

    // param.push_back("V3");      var["V3"]  = 0.0;
    param.push_back("V3p");     var["V3p"] = 0.3;
    param.push_back("K3");      var["K3"]  = 0.2;

    param.push_back("V4");      var["V4"]  = 0.1;
    param.push_back("K4");      var["K4"]  = 0.1;

    param.push_back("K5");      var["K5"]  = 0.02;
    param.push_back("K6");      var["K6"]  = 0.3;

    param.push_back("vi");      var["vi"]  = 0.1;
    param.push_back("k1");      var["k1"]  = 0.5;
    param.push_back("kd");      var["kd"]  = 0.02;
    param.push_back("a1");      var["a1"]  = 0.05;
    param.push_back("a2");      var["a2"]  = 0.05;
    param.push_back("d1");      var["d1"]  = 0.05;
    param.push_back("vs");      var["vs"]  = 0.2;
    param.push_back("alpha");   var["alpha"] = 0.1;

    //

    // Assignment:
    //
    //      V1 := C * V1p * pow( C + K6 , -1 )
    //      V3 := V3p * M

    aux["V1"] = Expression(
                    TIMES,
                    "C",
                    Expression(
                        TIMES,
                        "V1p",
                        Expression(
                            POWER,
                            Expression(PLUS, "C", "K6"),
                            -1.0
                        )
                    )
                );

    aux["V3"] = Expression(TIMES, "V3p", "M");

    //

    // Reaction / rule:
    //
    //      r1  :  vi
    //      r2  :  C * k1 * X * pow( C + K5, -1 )
    //      r3  :  C * kd
    //      r4  :  (1 + -M) * V1 * pow( K1 + -M + 1, -1 )
    //      r5  :  M * V2 * pow( K2 + M, -1 )
    //      r6  :  V3 * (1 + -X) * pow( K3 + -X + 1, -1 )
    //      r7  :  V4 * X * pow( K4 + X, -1 )
    //      r8  :  a1 * C * Y
    //      r9  :  a2 * Z
    //      r10 :  alpha * d1 * Z
    //      r11 :  alpha * kd * Z
    //      r12 :  vs
    //      r13 :  d1 * Y

    aux["r1"] = Expression("vi");

    aux["r2"] = Expression(
                    TIMES,
                    "C",
                    Expression(
                        TIMES,
                        "k1",
                        Expression(
                            TIMES,
                            "X",
                            Expression(
                                POWER,
                                Expression(PLUS, "C", "K5"),
                                -1.0
                            )
                        )
                    )
                );

    aux["r3"] = Expression(TIMES, "C", "kd");

    aux["r4"] = Expression(
                    TIMES,
                    Expression(PLUS, 1.0, Expression(MINUS, "M")),
                    Expression(
                        TIMES,
                        aux["V1"],
                        Expression(
                            POWER,
                            Expression(PLUS, "K1", Expression(PLUS, Expression(MINUS, "M"), 1.0)),
                            -1.0
                        )
                    )
                );

    aux["r5"] = Expression(
                    TIMES,
                    "M",
                    Expression(
                        TIMES,
                        "V2",
                        Expression(
                            POWER,
                            Expression(PLUS, "K2", "M"),
                            -1.0
                        )
                    )
                );

    aux["r6"] = Expression(
                    TIMES,
                    aux["V3"],
                    Expression(
                        TIMES,
                        Expression(PLUS, 1.0, Expression(MINUS, "X")),
                        Expression(
                            POWER,
                            Expression(PLUS, "K3", Expression(PLUS, Expression(MINUS, "X"), 1.0)),
                            -1.0
                        )
                    )
                );

    aux["r7"] = Expression(
                    TIMES,
                    "V4",
                    Expression(
                        TIMES,
                        "X",
                        Expression(
                            POWER,
                            Expression(PLUS, "K4", "X"),
                            -1.0
                        )
                    )
                );

    aux["r8"] = Expression(TIMES, "a1", Expression(TIMES, "C", "Y"));

    aux["r9"] = Expression(TIMES, "a2", "Z");

    aux["r10"] = Expression(TIMES, "alpha", Expression(TIMES, "d1", "Z"));

    aux["r11"] = Expression(TIMES, "alpha", Expression(TIMES, "kd", "Z"));

    aux["r12"] = Expression("vs");

    aux["r13"] = Expression(TIMES, "d1", "Y");

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
    // biosys.setSpecies(species);
    biosys.setParameters(param);
    //biosys.setAuxiliary(aux);
    //biosys.setAuxExpressions(auxilist);

    biosys.setInitialValue("C", 0.0);
    biosys.setInitialValue("X", 0.0);
    biosys.setInitialValue("M", 0.0);
    biosys.setInitialValue("Y", 1.0);
    biosys.setInitialValue("Z", 1.0);

    Vector meastp(50);
    for (long j=1; j <= meastp.nr(); ++j) meastp(j) = tstart + j*(tend-tstart)/50.0;
    biosys.setMeasurementTimePoints( meastp );

    // Breakpoints:
    //
    //  Subdivision of integration interval [t0 T] :  [ t0 = b1, b2, b3, ..., bn-1, bn = T ]
    //
    Vector breaktp(4);
    for (long j=1; j <= breaktp.nr(); ++j) breaktp(j) = tstart + (j-1)*(tend-tstart)/3.0;
    biosys.setBreakpoints( breaktp );


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


    Vector syndata, synscal;


    biosys.computeModel(var, "adaptive");
    // syndata = biosys.computeModel(var, "adaptive");

    std::cout << std::endl;
    std::cout << " -=-=-=-=-=-=-=- " << std::endl;

    emap.clear();
    emap = biosys.getVarExpr();

    for (unsigned jj = 1, j = 0; j < species.size(); ++j)
    {
        std::string eqn = species[j];

        for (unsigned k = 0; k < species.size(); ++k)
        {
            std::string s = eqn + " / " + species[k];
            std::cout << s << " :  ";
            std::cout << "( " << jj++ << " / " << emap.size() << " )" << std::endl;
            std::cout << emap[s];
            std::cout << std::endl;
        }

        for (unsigned k = 0; k < param.size(); ++k)
        {
            std::string s = eqn + " / " + param[k];
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
    // invBiosys.setSpecies(species);
    invBiosys.setParameters(param);

    for (unsigned j = 0; j < species.size(); ++j)
        invBiosys.setInitialValue( species[j], biosys.getInitialValue(species[j]) );

    for (unsigned j = 0; j < param.size(); ++j)
        invBiosys.setParamValue( param[j], biosys.getParamValue(param[j]) );

    //
    // From this point on only the new object "invBiosys" is referenced
    //

    invBiosys.setMeasurementList( meastp, measlist );


    syndata = invBiosys.getMeasurements();
    synscal.zeros( syndata.nr() );

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


    // initial guess for GaussNewton

    p.zeros(3);
    p(1) = 0.5;         // true: 0.2
    p(2) = 0.5;         // true: 0.1
    p(3) = 0.05;        // true: 0.02
    pscal.zeros(3);

    var.clear();

    par1.push_back( "K3" );         var["K3"] = p(1);
    par1.push_back( "K4" );         var["K4"] = p(2);
    par1.push_back( "K5" );         var["K5"] = p(3);

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


    // if ( iopt.jacgen > 1 )
        wk.cond = 1.0 / sqrt(solverRTol);
    // wk.nitmax = 2;


    gn.setProblem( &prob );
    gn.initialise( syndata.nr(), p, pscal, syndata, synscal, rtol, iopt, wk );

    gn.run();
    // gn.analyse();
    gn.printCounter();

    std::vector<Vector> piter = gn.getSolutionIter();

    std::cout << std::endl;
    std::cout << "Inv.Prob. solution iteration" << std::endl;
    std::cout << "----------------------------" << std::endl;
    for (unsigned j = 0; j < piter.size(); ++j)
        std::cout << "it = " << j << "\n" << piter[j].t();
    std::cout << std::endl;

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
