// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-12-08 td
// last changed:
//

#include <cmath>
// #include <cstdlib> // for rand(), RAND_MAX
#include <vector>

#include "addpkg/dlib/time_this.h"

#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "system/Expression.h"
#include "system/BioSystem.h"
#include "system/BioProcessor.h"

// #include "system/BioPAR.h"
//
// #include "nonlin/GaussNewton.h"


using namespace PARKIN;

/*
double randu()  // return a random number ( uniformly distributed in [0,1[ )
{
    return ( std::rand() / (1.0 + RAND_MAX) );
}

double randn()
{
    double u = 0.0, v = 0.0;
    while ( u == 0.0 ) u = randu();
    while ( v == 0.0 ) v = randu();
    return std::sqrt(-2*std::log(u)) * std::cos(2*M_PI*v);
    // or: std::sqrt(-2*std::log(u)) * std::sin(2*M_PI*v);
}
*/

int testsystem_aux()
{
    srand(0);

    /// define and solve forward problem first
    ///
    /// taken from SBML example files / cases
    /// cyclin: BIOMD0000000008.xml

    Real                        tstart = 1.0;
    Real                        tend   = 100.0;
    Expression::Param           par, var, varThres;

    BioSystem::Species          species;
    BioSystem::Parameter        param;
    BioSystem::Parameter        par1;
    BioSystem                   biosys(tstart, tend);

    BioRHS::ExpressionMap    aux;
    BioRHS::ExpressionMap    emap;

    //

    species.push_back("C");
    species.push_back("X");
    species.push_back("M");
    species.push_back("Y");
    species.push_back("Z");

    //

    // param.push_back("V1");      var["V1"]  = 0.0;
    param.push_back("glo_V1p");     var["glo_V1p"] = 0.75;
    param.push_back("r04_K1");      var["r04_K1"]  = 0.1;

    param.push_back("r05_V2");      var["r05_V2"]  = 0.25;
    param.push_back("r05_K2");      var["r05_K2"]  = 0.1;

    // param.push_back("V3");      var["V3"]  = 0.0;
    param.push_back("glo_V3p");     var["glo_V3p"] = 0.3;
    param.push_back("r06_K3");      var["r06_K3"]  = 0.2;

    param.push_back("r07_V4");      var["r07_V4"]  = 0.1;
    param.push_back("r07_K4");      var["r07_K4"]  = 0.1;

    param.push_back("r02_K5");      var["r02_K5"]  = 0.02;
    param.push_back("glo_K6");      var["glo_K6"]  = 0.3;

    param.push_back("r01_vi");      var["r01_vi"]  = 0.1;
    param.push_back("r02_k1");      var["r02_k1"]  = 0.5;
    param.push_back("r03_kd");      var["r03_kd"]  = 0.02;
    param.push_back("r11_kd");      var["r11_kd"]  = 0.02;
    param.push_back("r08_a1");      var["r08_a1"]  = 0.05;
    param.push_back("r09_a2");      var["r09_a2"]  = 0.05;
    param.push_back("r10_d1");      var["r10_d1"]  = 0.05;
    param.push_back("r13_d1");      var["r13_d1"]  = 0.05;
    param.push_back("r12_vs");      var["r12_vs"]  = 0.2;
    param.push_back("r10_alpha");   var["r10_alpha"] = 0.1;
    param.push_back("r11_alpha");   var["r11_alpha"] = 0.1;

    //

    // Assignment:
    //
    //      V1 := C * V1p * pow( C + K6 , -1 )
    //      V3 := V3p * M

    aux["glo_V1"] = Expression(
                    TIMES,
                    "C",
                    Expression(
                        TIMES,
                        "glo_V1p",
                        Expression(
                            POWER,
                            Expression(PLUS, "C", "glo_K6"),
                            -1.0
                        )
                    )
                );

    aux["glo_V3"] = Expression(TIMES, "glo_V3p", "M");

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

    aux["r1"] = Expression("r01_vi");

    aux["r2"] = Expression(
                    TIMES,
                    "C",
                    Expression(
                        TIMES,
                        "r02_k1",
                        Expression(
                            TIMES,
                            "X",
                            Expression(
                                POWER,
                                Expression(PLUS, "C", "r02_K5"),
                                -1.0
                            )
                        )
                    )
                );

    aux["r3"] = Expression(TIMES, "C", "r03_kd");

    aux["r4"] = Expression(
                    TIMES,
                    Expression(PLUS, 1.0, Expression(MINUS, "M")),
                    Expression(
                        TIMES,
                        aux["glo_V1"],
                        Expression(
                            POWER,
                            Expression(PLUS, "r04_K1", Expression(PLUS, Expression(MINUS, "M"), 1.0)),
                            -1.0
                        )
                    )
                );

    aux["r5"] = Expression(
                    TIMES,
                    "M",
                    Expression(
                        TIMES,
                        "r05_V2",
                        Expression(
                            POWER,
                            Expression(PLUS, "r05_K2", "M"),
                            -1.0
                        )
                    )
                );

    aux["r6"] = Expression(
                    TIMES,
                    aux["glo_V3"],
                    Expression(
                        TIMES,
                        Expression(PLUS, 1.0, Expression(MINUS, "X")),
                        Expression(
                            POWER,
                            Expression(PLUS, "r06_K3", Expression(PLUS, Expression(MINUS, "X"), 1.0)),
                            -1.0
                        )
                    )
                );

    aux["r7"] = Expression(
                    TIMES,
                    "r07_V4",
                    Expression(
                        TIMES,
                        "X",
                        Expression(
                            POWER,
                            Expression(PLUS, "r07_K4", "X"),
                            -1.0
                        )
                    )
                );

    aux["r8"] = Expression(TIMES, "r08_a1", Expression(TIMES, "C", "Y"));

    aux["r9"] = Expression(TIMES, "r09_a2", "Z");

    aux["r10"] = Expression(TIMES, "r10_alpha", Expression(TIMES, "r10_d1", "Z"));

    aux["r11"] = Expression(TIMES, "r11_alpha", Expression(TIMES, "r11_kd", "Z"));

    aux["r12"] = Expression("r12_vs");

    aux["r13"] = Expression(TIMES, "r13_d1", "Y");

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

    Vector meastp(20);
    for (long j=1; j <= meastp.nr(); ++j) meastp(j) = tstart + j*(tend-tstart)/20.0;
//    biosys.setMeasurementTimePoints( meastp );

    // Breakpoints / Event Management:
    //
    //  Subdivision of integration interval [t0 T] :  [ t0 = b1, b2, b3, ..., bn-1, bn = T ]
    //
    Vector breaktp;

/*
    breaktp.zeros(4);
    for (long j = 1; j <= breaktp.nr(); ++j) breaktp(j) = tstart + (j-1)*(tend-tstart)/3.0;

TIME_THIS_TO( std::cerr << " *** Call: biosys.setBreakpoints() *** " << std::endl;

    biosys.setBreakpoints( breaktp );

, std::cerr);

    // and succinctly some event handling
    //

    for (long j = 1; j < breaktp.nr()-1; ++j)
    {
        // emap.clear();
        emap = biosys.getEvent(j);

        emap[species[3]] = Expression(0.75);

        biosys.setEvent(j, emap);
    }

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
    {
        biosys.setParamValue(param[j], var[param[j]]);

        std::cout << std::setw(20) << param[j] <<
                             "   " << biosys.getParamValue(param[j]) <<
                                      std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;

    //
    //  take the following lines in to change the default values
    //  ATOL = 10*EPMACH, and RTOL = 1.0e-12
    //
    biosys.setSolverRTol(1.0e-8);
    biosys.setSolverATol(1.0e-8);


    //-----------------------------------------------------------------------------
// #define SENSI
#ifdef SENSI
    double                      rtolS = 1.0e-5;
    long                        nS    = species.size();
    long                        TS    = meastp.nr();
    Vector                      pS, pscalS;
    Vector                      dummyMeasS, dummyScalS;

    dummyMeasS.zeros(nS*TS);
    dummyScalS.ones(nS*TS);

/*
    BioSystem::MeasurementList  mlistS;

    for (long tp = 0; tp < TS; ++tp)
    {
        mlistS.push_back( MeasurementPoint() );
        for (long j = 0; j < (long)species.size(); ++j)
        {
            mlistS[tp][species[j]] = std::make_pair<double,double>(0.0,0.0);
        }
    }

    biosys.setMeasurementList( meastp, mlistS );
*/

    //
    // Compute sensitivities at initial guess of GaussNewton
    //

    pS.zeros(3);
    pS(1) = 0.2;         // true: 0.2
    pS(2) = 0.1;         // true: 0.1
    pS(3) = 0.02;        // true: 0.02
    pscalS.zeros(3);
    pscalS(1) = pS(1);
    pscalS(2) = pS(2);
    pscalS(3) = pS(3);

    par1.clear();

    par1.push_back( "r06_K3" );         // var[ "r06_K3" ] = pS(1);
    par1.push_back( "r07_K4" );         // var[ "r07_K4" ] = pS(2);
    par1.push_back( "r02_K5" );         // var[ "r02_K5" ] = pS(3);

    IOpt           ioptS;
    GaussNewtonWk  wkS;
    GaussNewton    gnS;
    BioPAR         probS( &biosys, par1 );

    ioptS.mode      = 0;     // 0:normal run, 1:single step
    ioptS.jacgen    = 1;     // 1:user supplied Jacobian, 2:num.diff., 3:num.diff.(with feedback)
    ioptS.qrank1    = false;     // allow Broyden rank-1 updates if __true__
    ioptS.nonlin    = 4;     // 1:linear, 2:mildly nonlin., 3:highly nonlin., 4:extremely nonlin.
    ioptS.rscal     = 1;     // 1:use unchanged fscal, 2:recompute/modify fscal, 3:use automatic scaling only
    ioptS.norowscal = false;     // allow for automatic row scaling of Jacobian if __false__
    ioptS.lpos      = false;      // force solution vector to be positive (all components > 0.0)
                            //          _mprmon =   0      1      2      3      4       5       6
                            //  dlib::log_level =  LNONE  LINFO  LVERB  LTALK  LGABBY  LDEBUG  LTRACE
    ioptS.mprmon    = 2;
    ioptS.mprerr    = 1;


    // if ( iopt.jacgen > 1 )
        // wkS.cond = 1.0 / sqrt(solverRTol);
    // wkS.nitmax = 15;


    gnS.setProblem( &probS );
    gnS.initialise( dummyMeasS.nr(), pS, pscalS, dummyMeasS, dummyScalS, rtolS, ioptS, wkS );


TIME_THIS_TO( std::cerr << " *** Call: gn.computeSensitivity() *** " << std::endl;
    std::cerr << "rc = " << gnS.computeSensitivity();
, std::cerr)

    std::cout << std::endl;
    std::cout << "Sensitivity QR diagonal" << std::endl;
    std::cout << "-----------------------" << std::endl;
    std::cout << " qrA.getDiag() = " << std::endl;
    std::cout << gnS.getSensitivity().getDiag().t() << std::endl;

    Matrix reducedA;
    Matrix A = gnS.getSensitivityMatrix();

    std::cout << "  A = (" << A.nr() << " x " << A.nc() << ") " << std::endl;
    std::cout << A << std::endl;

    reducedA.zeros(nS, A.nc());
    for (unsigned k = 0; k < (unsigned)A.nc(); ++k)
    {
        Vector v = A.colm( k + 1 );

        for (unsigned j = 0; j < (unsigned)nS; ++j)
        {
            // long   off = j + 1;
            double s = 0.0;

            for (unsigned tp = 0; tp < (unsigned)TS; ++tp)
            {
                double val = v( tp*nS + j + 1 ); // off);
                s += std::pow(val, 2.0);
            }

            reducedA(j+1,k+1) = std::sqrt(s);
        }
    }

    BioSystem::Species   cSpec  = probS.getSpecies();
    BioSystem::Parameter cParam = probS.getCurrentParameter();

    std::cout << " rA = (" << reducedA.nr() << " x " << reducedA.nc() << ") " << std::endl;
    std::cout << std::setw(10) << "\t";
    for (unsigned k = 0; k < cParam.size(); ++k)
    {
        std::cout << std::setw(10) << cParam[k] << " ";
    }
    std::cout << std::endl << std::endl;
    for (unsigned j = 0; j < cSpec.size(); ++j)
    {
        std::cout << std::setw(10) << cSpec[j] << "\t";
        for (unsigned k = 0; k < cParam.size(); ++k)
        {
            std::cout << std::setw(10) << reducedA(j+1,k+1) << " ";
        }
        std::cout << std::endl << std::endl;
    }

#endif
    //-----------------------------------------------------------------------------


    Vector syndata, synscal, vref, utmp;

/*
TIME_THIS_TO( std::cerr << " *** Call: biosys.computeModel() *** " << std::endl;

//for (long n=0; n < 100; ++n)
    vref = biosys.computeModel(var, "adaptive");
    // syndata = biosys.computeModel(var, "adaptive");

, std::cerr );
*/

//std::cerr << " vref = " << std::endl;
//std::cerr << vref;
//std::cerr << std::endl;
//std::cerr << " *** Retn: biosys.computeModel() *** " << std::endl;
//
//exit(-999);

    Matrix Jac, tmat;
//    long   n = species.size();
//    long   T = biosys.getOdeTrajectoryTimePoints().nr(); // meastp.nr();
    long   j = 0;

    utmp.zeros(3);

    par["r06_K3"] = var["r06_K3"];      utmp(1) = ( par["r06_K3"] );
    par["r07_K4"] = var["r07_K4"];      utmp(2) = ( par["r07_K4"] );
    par["r02_K5"] = var["r02_K5"];      utmp(3) = ( par["r02_K5"] );

//    Jac.zeros( n*T, 7 );


TIME_THIS_TO( std::cerr << " *** Call: biosys.computeJacobian() *** " << std::endl;

    // Jac.set_colm( 1, 3 ) =
    tmat =
        biosys.computeJacobian( par, "adaptive" ); // * utmp.diag();  // ... * exp( u=log(par) ).diag()

, std::cerr )

    Jac.zeros( tmat.nr(), 7 );
    Jac.set_colm( 1 , 3 ) = tmat;

//exit(-41);

    vref = biosys.computeModel(par);

    for (Expression::ParamIterConst it = par.begin(); it != par.end(); ++it)
    {
        Real   h = std::sqrt(1e-10), psav;
        Vector vpert;

        psav = par[it->first];
        par[it->first] = psav + h; // std::exp( std::log(psav) + h );
        vpert = biosys.computeModel(par);
        par[it->first] = psav;

        Jac.set_colm( 5+j++ ) = (1.0/h)*(vpert - vref);
    }

    std::cout << " Jacobian Jac (" << Jac.nr() << " x " << Jac.nc() << ") = " << std::endl;
    std::cout << Jac << std::endl;



// exit(-42);

    biosys.setMeasurementTimePoints( meastp );
    biosys.computeModel(par, "init");


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

    Real   xtol       = 1.0e-5;
    Real   solverRTol = 1.0e-8;
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
    Vector                      newtp, savtp;
    Real                        sigma = 0.015;   // add 1.5% white noise ...
    unsigned                    sampleSize = 3;

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
                Real newval = val + sigma*randn();

                newMeas[sampleSize*tp+k][species[j]] =
                    std::make_pair<Real,Real>( newval, 1.0 );
            }
        }
    }

    //
    // From this point on only the new object "invBiosys" is referenced
    //

    invBiosys.setMeasurementList( newtp, newMeas );

    savtp = meastp;
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


    ///

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

    par1.clear();
    var.clear();
    varThres.clear();

    par1.push_back( "r06_K3" );     var[ "r06_K3" ] = p(1);     varThres[ "r06_K3" ] = EPMACH;
    par1.push_back( "r07_K4" );     var[ "r07_K4" ] = p(2);     varThres[ "r07_K4" ] = EPMACH;
    par1.push_back( "r02_K5" );     var[ "r02_K5" ] = p(3);     varThres[ "r02_K5" ] = EPMACH;


/*
TIME_THIS_TO( std::cerr << " *** Call: invBiosys.computeJacobian() *** " << std::endl;

    tmat =
        invBiosys.computeJacobian( var, "adaptive" );

, std::cerr )
*/


    IOpt           iopt;
//    GaussNewtonWk  wk;
//    GaussNewton    gn;
//    BioPAR         prob( &invBiosys, par1 );
    BioProcessor   proc( &invBiosys, "nlscon" );

    iopt.mode      = 0;     // 0:normal run, 1:single step
    iopt.jacgen    = 3;     // 1:user supplied Jacobian, 2:num.diff., 3:num.diff.(with feedback)
    iopt.qrank1    = false;     // allow Broyden rank-1 updates if __true__
    iopt.nonlin    = 4;     // 1:linear, 2:mildly nonlin., 3:highly nonlin., 4:extremely nonlin.
    iopt.rscal     = 1;     // 1:use unchanged fscal, 2:recompute/modify fscal, 3:use automatic scaling only
    iopt.norowscal = false;     // allow for automatic row scaling of Jacobian if __false__
    iopt.lpos      = false;     // force solution vector to be positive (all components > 0.0)
                            //          _mprmon =   0      1      2      3      4       5       6
                            //  dlib::log_level =  LNONE  LINFO  LVERB  LTALK  LGABBY  LDEBUG  LTRACE
    iopt.itmax     = 45;

    iopt.mprmon    = 2;
    iopt.mprerr    = 1;


    proc.setIOpt( iopt );

    proc.setCurrentParamValues( var );
    proc.setCurrentParamThres( varThres );

    Expression::Param speThres( proc.getCurrentSpeciesThres() );
    for (Expression::ParamIterConst it = speThres.begin();
                                    it != speThres.end(); ++it)
    {
        speThres[it->first] = EPMACH;
    }
    proc.setCurrentSpeciesThres( speThres );

/*
TIME_THIS_TO( std::cerr <<  " *** call: proc.computeSensitivityTrajectories() *** " << std::endl;
    std::cerr << "sensTraj.size() = " << proc.computeSensitivityTrajectories().size() << std::endl;
, std::cerr )
*/

TIME_THIS_TO( std::cerr <<  " *** call: proc.prepareDetailedSensitivities() *** " << std::endl;
    std::cerr << "rc = " << proc.prepareDetailedSensitivities( savtp ) << std::endl;
, std::cerr )


    std::vector<Matrix> sensMatList = proc.getSensitivityMatrices();

    for (unsigned j = 0; j < sensMatList.size(); ++j)
    {
        std::cout << " Sensitivity for timepoint #" << j+1 << " (t = " << savtp(j+1) << ") :" << std::endl;
        std::cout << " mat (" << sensMatList[j].nr() << " x " << sensMatList[j].nc() << ") = " << std::endl;
        std::cout << sensMatList[j] << std::endl;
    }


// exit(-43);

    //
    // !!! !!! !!!
    //

    invBiosys.setMeasurementList( measlist );

    //
    // !!! !!! !!!
    //

    proc.identifyParameters(xtol);

    // Expression::Param final = proc.getIdentificationResults();

/*
    GaussNewtonWk  wk;
    GaussNewton    gn;
    BioPAR         prob( &invBiosys, par1 );

    // if ( iopt.jacgen > 1 )
        // wk.cond = 1.0 / sqrt(solverRTol);
    // wk.nitmax = 15;

    gn.setProblem( &prob );
    gn.initialise( syndata.nr(), p, pscal, syndata, synscal, xtol, iopt, wk );


TIME_THIS_TO( std::cerr << " *** Call: gn.computeSensitivity() *** " << std::endl;
    std::cerr << "rc = " << gn.computeSensitivity();
, std::cerr)

    std::cout << std::endl;
    std::cout << "Sensitivity QR diagonal" << std::endl;
    std::cout << "-----------------------" << std::endl;
    std::cout << " qrA.getDiag() = " << std::endl;
    std::cout << gn.getSensitivity().getDiag().t() << std::endl;


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
*/
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
