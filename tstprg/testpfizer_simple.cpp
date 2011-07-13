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
#include "system/BioPAR.h"

#include "nonlin/GaussNewton.h"


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

int testpfizer_simple()
{
    srand(0);
    std::cout << std::endl;
    std::cout << "**********************************************************************" << std::endl;
    std::cout << "*                                                                    *" << std::endl;
    std::cout << "* Simple test system: Cmplx-Rcptr (w/ central/peripheral compartmnt) *" << std::endl;
    std::cout << "*                                                                    *" << std::endl;
    std::cout << "**********************************************************************" << std::endl;
    std::cout << std::endl;

    /// define and solve forward problem first

    Real                        tstart = -10.0;
    Real                        tend   = 100.0;
    Expression::Param           var, par;

    BioSystem::Species          species;
    BioSystem::Parameter        param;
    BioSystem::Parameter        par1;
    BioSystem                   biosys(tstart, tend);

    BioRHS::ExpressionMap    aux;
    BioRHS::ExpressionMap    emap;

    //

    // species.push_back("Gut");
    species.push_back("V_c");
    species.push_back("V_p");
    species.push_back("Receptor");
    species.push_back("Complex");
    // species.push_back("src_Receptor");
    // species.push_back("snk_Receptor");
    // species.push_back("snk_Complex");

    //

    param.push_back("dose");     var[ "dose"  ] = 35.0; // 47.297;
    // param.push_back("k_a");      var[ "k_a"   ] = 0.0;
    param.push_back("k_out");    var[ "k_out" ] = 0.21984;
    param.push_back("k_cp");     var[ "k_cp"  ] = 0.1229;
    param.push_back("k_pc");     var[ "k_pc"  ] = 0.1642;
    param.push_back("k_on");     var[ "k_on"  ] = 1.0e7;
    param.push_back("k_off");    var[ "k_off" ] = 1.0e-5;
    param.push_back("k_syn");    var[ "k_syn" ] = 2.028;
    param.push_back("k_deg");    var[ "k_deg" ] = 0.3149;

    //

    // Reaction / rule:
    //
    //      re4:  r1  :  k_pc * V_p
    //      re5:  r2  :  k_cp * V_c
    //      re6:  r3  :  k_out * V_c
    //      re7:  r4  :  k_syn
    //      re8:  r5  :  k_deg * Receptor
    //      re11: r6  :  k_on * V_c * Receptor
    //      re12: r7  :  k_off * Complex

    aux["r1"] = Expression(TIMES, "k_pc", "V_p");
    aux["r2"] = Expression(TIMES, "k_cp", "V_c");
    aux["r3"] = Expression(TIMES, "k_out", "V_c");
    aux["r4"] = Expression("k_syn");
    aux["r5"] = Expression(TIMES, "k_deg", "Receptor");
    aux["r6"] = Expression(TIMES, "k_on",
                           Expression(TIMES, "V_c", "Receptor"));
    aux["r7"] = Expression(TIMES, "k_off", "Complex");

    //

    // ODEs:
    //
    //      (V_c)' = r1 - r2 - r3 - r6 + r7
    //             = re4 - re5 -re6 - re11 + re12
    //
    //      (V_p)' = -r1 + r2
    //             = -re4 + re5
    //
    //      (Rpr)' = r4 - r5 - r6 + r7
    //             = re7 - re8 - re11 + re12
    //
    //      (Cmx)' = r6 - r7
    //             = re11 - re12

    emap["V_c"] = Expression(
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
                                    Expression(MINUS, aux["r6"]),
                                    aux["r7"]
                                )
                            )
                        )
                  );

    emap["V_p"] = Expression(
                        PLUS,
                        Expression(TIMES, -1.0, aux["r1"]),
                        aux["r2"]
                  );

    emap["Receptor"] = Expression(
                            PLUS,
                            Expression(TIMES, 1.0, aux["r4"]),
                            Expression(
                                PLUS,
                                Expression(MINUS, aux["r5"]),
                                Expression(
                                    PLUS,
                                    Expression(MINUS, aux["r6"]),
                                    aux["r7"]
                                )
                            )
                       );

    emap["Complex"] = Expression(
                            PLUS,
                            Expression(TIMES, 1.0, aux["r6"]),
                            Expression(MINUS, aux["r7"])
                        );


    biosys.setODESystem(emap);
    // biosys.setSpecies(species);  // has now a complete different effect:
                                    // resets to an identity 'emap' mapping
    biosys.setParameters(param);


    biosys.setInitialValue("V_p",       0.0);
    biosys.setInitialValue("V_c",       0.0);
    biosys.setInitialValue("Receptor", 10.0);
    biosys.setInitialValue("Complex",   0.0);

    Vector meastp(200);
    for (long j=1; j <= meastp.nr(); ++j) meastp(j) = tstart + j*(tend-tstart)/200.0;
    biosys.setMeasurementTimePoints( meastp );

    // Breakpoints / Event Management:
    //
    //  Subdivision of integration interval [t0 T] :  [ t0 = b1, b2, b3, ..., bn-1, bn = T ]
    //
    Vector breaktp;

    breaktp.zeros(3);
    // for (long j = 1; j <= breaktp.nr(); ++j) breaktp(j) = tstart + (j-1)*(tend-tstart)/2.0;
    breaktp(1) = tstart;
    breaktp(2) = 0.0;
    breaktp(3) = tend;

TIME_THIS_TO( std::cerr << " *** Call: biosys.setBreakpoints() *** " << std::endl;

    biosys.setBreakpoints( breaktp );

, std::cerr);

    // and succinctly some event handling
    //

    {
        emap = biosys.getEvent(1);

        emap["V_c"] = Expression("dose");

        biosys.setEvent(1, emap);
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

    ///

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
    Real rTol = 1.0e-7;
    Real aTol = 1.0e-9;
    Real xTol = 1.0e-3;

    biosys.setSolverRTol(rTol);
    biosys.setSolverATol(aTol);


    //-----------------------------------------------------------------------------
#define SENSI
#ifdef SENSI
    double                      xTolS = xTol;
    long                        nS    = species.size();
    long                        qS    = param.size();
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
    par1.clear();
    par1 = param;

    pS.zeros( qS );
    pscalS.zeros( qS );

    for (long j = 1; j <= qS; ++j )
    {
        pS(j) = pscalS(j) = var[ par1[j-1] ];
    }

    IOpt           ioptS;
    GaussNewtonWk  wkS;
    GaussNewton    gnS;
    BioPAR         probS( &biosys, par1 );

    ioptS.mode      = 0;     // 0:normal run, 1:single step
    ioptS.jacgen    = 3;     // 1:user supplied Jacobian, 2:num.diff., 3:num.diff.(with feedback)
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
        wkS.cond = 1.0 / /*std::sqrt*/(rTol);
    // wkS.nitmax = 15;


    gnS.setProblem( &probS );
    gnS.initialise( dummyMeasS.nr(), pS, pscalS, dummyMeasS, dummyScalS, xTolS, ioptS, wkS );


TIME_THIS_TO( std::cerr << " *** Call: gn.computeSensitivity() *** " << std::endl;
    std::cerr << "rc = " << gnS.computeSensitivity();
, std::cerr)


    Matrix reducedA;
    Matrix A = gnS.getSensitivityMatrix();

    // std::cout << "  A = (" << A.nr() << " x " << A.nc() << ") " << std::endl;
    // std::cout << A << std::endl;

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

            reducedA(j+1,k+1) = std::sqrt(s / TS);
        }
    }

    Vector pivot = gnS.getSensitivity().getPivot();
    Vector sensDiag = gnS.getSensitivity().getDiag();
    Vector subCond = sensDiag;

    BioSystem::Species   cSpec  = probS.getSpecies();
    BioSystem::Parameter cParam = probS.getCurrentParameter();

    for (long j = 1; j <= sensDiag.nr(); ++j)
    {
        subCond(j) = (sensDiag(j) != 0.0) ? std::fabs( sensDiag(1) ) / std::fabs( sensDiag(j) ) : 0.0;
    }

    std::cout << std::endl;
    std::cout << "Sensitivity QR diagonal" << std::endl;
    std::cout << "-----------------------" << std::endl;
    std::cout << " qrA.getDiag() = " << std::endl;
    std::cout << sensDiag.t() << std::endl;
    std::cout << " corresp. Sub-Condition = " << std::endl;
    std::cout << " " << subCond.t(); // << std::endl;


    for (unsigned j = 0; j < cParam.size(); ++j)
    {
        std::cout << std::setw(8) << cParam[pivot(j+1)-1];
    }

    std::cout << std::endl << std::endl;
    std::cout << std::endl;


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

TIME_THIS_TO( std::cerr << " *** Call: biosys.computeModel() *** " << std::endl;

//for (long n=0; n < 100; ++n)
    vref = biosys.computeModel(var, "init");
    // syndata = biosys.computeModel(var, "init");

, std::cerr );

//std::cerr << " vref = " << std::endl;
//std::cerr << vref;
//std::cerr << std::endl;
//std::cerr << " *** Retn: biosys.computeModel() *** " << std::endl;
//
//exit(-999);

/*
    Vector tp( biosys.getOdeTrajectoryTimePoints() );
    Vector x1( biosys.getOdeTrajectory(0) );
    Vector x2( biosys.getOdeTrajectory(1) );
    Vector x3( biosys.getOdeTrajectory(2) );
    Vector x4( biosys.getOdeTrajectory(3) );

    std::cout << "=====================" << std::endl;
    std::cout << " t = " << std::endl;
    std::cout << tp.t();
    std::cout << "=====================" << std::endl;
    std::cout << " " << species[0] << " [a.u] = " << std::endl;
    std::cout << x1.t();
    std::cout << "---------------------" << std::endl;
    std::cout << " " << species[1] << " [a.u] = " << std::endl;
    std::cout << x2.t();
    std::cout << "---------------------" << std::endl;
    std::cout << " " << species[2] << " [a.u] = " << std::endl;
    std::cout << x3.t();
    std::cout << "---------------------" << std::endl;
    std::cout << " " << species[3] << " [a.u] = " << std::endl;
    std::cout << x4.t();
    std::cout << "=====================" << std::endl;
*/
    ///
    /// and, subsequently, prepare and solve for the inverse problem
    ///

    Real   reconXTol  = xTol;
    Real   solverRTol = rTol;
    Real   solverATol = aTol;
    Vector p, pscal;
    // Vector syndata, synscal;

    BioSystem::MeasurementList  measlist = biosys.getMeasurementList();
    BioSystem                   invBiosys( tstart, tend );

    invBiosys.setSolverRTol( solverRTol );
    invBiosys.setSolverATol( solverATol );

    invBiosys.setODESystem(emap);
    // invBiosys.setSpecies(species);       // see above comment!
    invBiosys.setParameters(param);
    invBiosys.setBreakpoints(breaktp);

    for (unsigned j = 0; j < species.size(); ++j)
        invBiosys.setInitialValue( species[j], biosys.getInitialValue(species[j]) );

    for (unsigned j = 0; j < param.size(); ++j)
        invBiosys.setParamValue( param[j], biosys.getParamValue(param[j]) );

    for (unsigned j = 0; j < (unsigned)breaktp.nr(); ++j)
        invBiosys.setEvent( j, biosys.getEvent(j) );

    //
    // Some kind of bootstrapping ...
    //

    BioSystem::MeasurementList  newMeas;
    Vector                      newtp;
    Real                        sigma = 0.025;   // add 2.5% white noise ...
    unsigned                    sampleSize = 1;
    unsigned                    nSpecies = 2;

    newMeas.clear();
    newtp.zeros( sampleSize*meastp.nr() );

    for (unsigned tp = 0; tp < (unsigned) meastp.nr(); ++tp )
    {
        for (unsigned k = 0; k < sampleSize; ++k) {
            newtp( sampleSize*tp + k + 1 ) = meastp(tp+1);
            newMeas.push_back( MeasurementPoint() );
        }

        for (unsigned j = 0; j < nSpecies /* species.size() */; ++j)
        {
            double val = measlist[tp][species[j]].first;

            for (unsigned k = 0; k < sampleSize; ++k)
            {
                Real newval = val + sigma*randn();

                if (newval >= 0.0)
                {
                    newMeas[sampleSize*tp+k][species[j]] =
                        std::make_pair<Real,Real>( newval, 1.0 );
                }
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

    out.open("testpfizer_simple.csv");

    if ( !out.is_open() )
    {
        std::cerr << "### ERROR: Could not open data file 'testpfizer_simple.csv'.\n";
        return -3;
    }

    std::cout << "\n";
    std::cout << "Timepoint [s]";
    out << "Timepoint [s]";
    for (unsigned j = 0; j < nSpecies /* species.size() */; ++j)
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
        for (unsigned j = 0; j < nSpecies /* species.size() */; ++j)
        {
            if ( measlist[tp].count(species[j]) > 0 )
            {
                std::cout << "\t";
                std::cout << std::scientific << measlist[tp][species[j]].first;
                out << "\t";
                out << std::scientific << measlist[tp][species[j]].first;
            }
            else
            {
                std::cout << "\t";
                std::cout << std::setw(13) << " ";
                out << "\t";
                out << std::setw(13) << " ";
            }
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


    std::ifstream in;

    in.open("iv7mg_b.csv");
    // in.open("testpfizer_simple.csv");

    if ( !in.is_open() )
    {
        std::cerr << "### ERROR: Could not open measurement data file 'iv7mg_b.csv'.\n";
        return -4;
    }

    ODESolver::Grid    mTimepoint;
    BioSystem::Species mSpecies;
    std::string        line;
    long               tp = 0;

    mSpecies.clear();
    newMeas.clear();

    while ( std::getline(in, line) )
    {
        std::stringstream ss(line);
        std::string       id, unit, field;

        if ( mSpecies.empty() )
        {
            ss >> id >> unit;  // we do not check for id=='Timepoint' here...
            std::cout << id << " " << unit << "\t";

            while ( ss.good() )
            {
                ss >> id >> unit;
                std::cout << id << " " << unit << "\t";
                mSpecies.push_back(id);
            }
            std::cout << std::endl;

            continue;
        }
        else
        {
            long j = -1;

            while (std::getline(ss, field, '\t'))
            {
                std::stringstream fs(field);
                Real              measVal;

                if ( !(fs >> measVal) )
                {
                    std::cout << ">   . / .    <" << "\t";

                    if (j < 0)
                    {
                        break;
                    }
                    else
                    {
                        ++j;
                        continue;
                    }
                }

                std::cout << ">" << measVal << "<" << "\t";

                if (j < 0)
                {
                    mTimepoint.push_back(measVal);
                    newMeas.push_back( MeasurementPoint() );
                    tp = newMeas.size()-1;
                }
                else
                {
                    newMeas[tp][mSpecies[j]] = std::make_pair<Real,Real>(measVal, 0.0);
                }

                ++j;
            }

            if ( newMeas[tp].empty() )
            {
                newMeas.pop_back();
                mTimepoint.pop_back();
            }

            std::cout << std::endl;
        }
    }

    in.close();

    std::cout << "Read " << newMeas.size() << " timepoint(s) with measurement data." << std::endl;

    invBiosys.setMeasurementList( Vector(mTimepoint), newMeas );

    syndata = invBiosys.getMeasurements();
    synscal.ones( syndata.nr() );
    synscal = sigma * synscal;

    std::cout << std::endl;
    std::cout << "Measured  experimental  data  :  nr() = " << syndata.nr() << std::endl;
    // std::cout << syndata << std::endl;


    //
    // Initial guess for GaussNewton
    //

    long q = param.size();

    par1.clear();
    p.zeros( q );
    pscal.zeros( q );

    for (long j = 1; j <= q; ++j)
    {
        par1.push_back( param[j-1] );
        p(j) = pscal(j) = var[ par1.back() ];
    }

    p(5) = pscal(5) = 1.0;      var[ par1[4] ] = p(5);
    p(6) = pscal(6) = 1.0;      var[ par1[5] ] = p(6);


    IOpt           iopt;
    GaussNewtonWk  wk;
    GaussNewton    gn;
    BioPAR         prob( &invBiosys, par1 );

    iopt.mode      = 0;     // 0:normal run, 1:single step
    iopt.jacgen    = 3;     // 1:user supplied Jacobian, 2:num.diff., 3:num.diff.(with feedback)
    iopt.qrank1    = false;     // allow Broyden rank-1 updates if __true__
    iopt.nonlin    = 3;     // 1:linear, 2:mildly nonlin., 3:highly nonlin., 4:extremely nonlin.
    iopt.rscal     = 1;     // 1:use unchanged fscal, 2:recompute/modify fscal, 3:use automatic scaling only
    iopt.norowscal = false;     // allow for automatic row scaling of Jacobian if __false__
    iopt.lpos      = true;      // force solution vector to be positive (all components > 0.0)
                            //          _mprmon =   0      1      2      3      4       5       6
                            //  dlib::log_level =  LNONE  LINFO  LVERB  LTALK  LGABBY  LDEBUG  LTRACE
    iopt.mprmon    = 2;
    iopt.mprerr    = 1;


    // if ( iopt.jacgen > 1 )
        wk.cond = 1.0 / /*std::sqrt*/(solverRTol);
    // wk.nitmax = 15;


    gn.setProblem( &prob );
/*
    gn.initialise( syndata.nr(), p, pscal, syndata, synscal, reconXTol, iopt, wk );

TIME_THIS_TO( std::cerr << " *** Call: gn.computeSensitivity() *** " << std::endl;
    std::cerr << "rc = " << gn.computeSensitivity();
, std::cerr)

    std::cout << std::endl;
    std::cout << "Sensitivity QR diagonal" << std::endl;
    std::cout << "-----------------------" << std::endl;
    std::cout << " qrA.getDiag() = " << std::endl;
    std::cout << gn.getSensitivity().getDiag().t() << std::endl;
*/

    gn.initialise( syndata.nr(), p, pscal, syndata, synscal, reconXTol, iopt, wk );
    gn.run();
    // gn.analyse();
    gn.printCounter();

    std::vector<Vector> piter = gn.getSolutionIter();
    Vector psol = gn.getSolution();

    std::cout << std::endl;
    std::cout << "Inv.Prob. solution iteration" << std::endl;
    std::cout << "----------------------------" << std::endl;
    for (unsigned j = 0; j < piter.size(); ++j)
        std::cout << "it = " << j << "\n" << piter[j].t();
    std::cout << std::endl;


    Expression::Param final = invBiosys.getSysPar();
    /*
    Expression::Param final = var;
    for (long j = 1; j <= psol.nr(); ++j)
    {
        final[ par1[j-1] ] = psol(j);
    }
    */

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
