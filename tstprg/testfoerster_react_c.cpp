// Copyright (C) 2012
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2012-03-05 td
// last changed:
//

#include <cmath>
// #include <cstdlib> // for rand(), RAND_MAX
#include <vector>

#include "addpkg/dlib/time_this.h"
#include "addpkg/gnuplot_cpp/gnuplot_i.hpp"  // for graphical gnuplot output on the fly

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

///

struct Experiment
{
    Real                tsync;
    Vector              breaktp;
    Expression::Param   y0;
    std::string         meas_fname;
};

///

int testfoerster_react_c()
{
    srand(0);
    std::cout << std::endl;
    std::cout << "**********************************************************************" << std::endl;
    std::cout << "*                                                                    *" << std::endl;
    std::cout << "* Immunoassay system: Reaction Scheme C (modified Scheme B w/red. P) *" << std::endl;
    std::cout << "*                                                                    *" << std::endl;
    std::cout << "**********************************************************************" << std::endl;
    std::cout << std::endl;

    const Real                    u = 1.6605389;  // * 10^(-27) [kg]    atomic unit mass
    const Real                  N_A = 6.0221413;  // * 10^(+23) [1/mol] Avogadro constant

    /// define and solve forward problem first

    Real                        taccu   = 0.0;
    Real                        ti      = 90.0;
    Real                        t_ini[] = { 10.0, 10.0, 19.0 };
    Real                        tf      = 9 * 60.0;

    Real                        P0[]    = { 5.0, 10.0, 20.0 };
    Real                        A0      = 75.0 / (1.5 * (0.375 + 0.75) * u * N_A);
    Real                        B0      = 10.0;
    Real                        Ps0     = 10.0;

    std::string                 measf[] = {
                                            "Foerster2_faser501_5ng.dat",
                                            "Foerster2_faser501_10ng.dat",
                                            "Foerster2_faser501_20ng.dat"
                                          };

    const unsigned              nexp = 3; /// 2;
    Experiment                  foerster[nexp];

    taccu = 0.0;

    for (unsigned j = 0; j < nexp; ++j)
    {
        Real ttmp = ti + t_ini[j];

        foerster[j].tsync = taccu + ttmp;

        foerster[j].y0.clear();
        foerster[j].y0["A"]   =  A0;
        foerster[j].y0["P"]   =  (1000.0/314.0) * P0[j];
        foerster[j].y0["AP"]  =  0.0;
        foerster[j].y0["APs"] =  0.0;
        foerster[j].y0["B"]   =  3.0 * B0;
        foerster[j].y0["AB"]  =  0.0;
        foerster[j].y0["As"]  =  0.0;
        foerster[j].y0["Ps"]  =  3.0 * Ps0;
        foerster[j].y0["S"]   =  0.0;
        foerster[j].y0["Sw"]  =  0.0;

        foerster[j].breaktp.zeros(3);
        foerster[j].breaktp(1) = 0.0;
        foerster[j].breaktp(2) = ti;
        foerster[j].breaktp(3) = ttmp + tf;

        foerster[j].meas_fname = measf[j];

        taccu += foerster[j].breaktp(3);
    }

    //

    Real                        tstart = foerster[0].breaktp(1);
    Real                        tend   = taccu;
    Expression::Param           var, par, parThres;

    BioSystem::Species          species;
    BioSystem::Parameter        param;
    BioSystem::Parameter        par1;
    BioSystem                   biosys(tstart, tend);

    BioRHS::ExpressionMap       aux;
    BioRHS::ExpressionMap       emap;
    BioRHS::ExprTypeMap         tmap;

    //

    species.push_back("A");     // mobiler Antikoerper-Konzentration
    species.push_back("P");     // Antigen(P4)-Konzentration
    species.push_back("AP");    // mobiler Antigen-Antikoerper-Komplex
    species.push_back("APs");   // AG-Derivat-AK-Komplex an Saeule (immobil)
    species.push_back("B");     // Verunreinigungskonzentration an Saeule
    species.push_back("AB");    // Konkurenz-Komplex an Saeule
    species.push_back("As");    // bindungsfaehiger AK in Loesung
    species.push_back("Ps");    // P4-Derivat-Konzentration an Saeule (wichtig, falls nicht im Ueberschuss!)
    species.push_back("S");     // Messsignal
    species.push_back("Sw");    // Schalter fuer Inkubationsstop

    //

    param.push_back("k_1");     var[ "k_1" ] = 1.0e-2;
    param.push_back("k_2");     var[ "k_2" ] = 1.0e-2;
    param.push_back("k_3");     var[ "k_3" ] = 1.0e-2;
    param.push_back("k_4");     var[ "k_4" ] = 1.0e-2;
    param.push_back("k_5");     var[ "k_5" ] = 1.0e-2;
    param.push_back("k_6");     var[ "k_6" ] = 1.0e-2;
    param.push_back("k_7");     var[ "k_7" ] = 1.0e-2;
    param.push_back("k_8");     var[ "k_8" ] = 1.0e-2;
    param.push_back("b");       var[  "b"  ] = 5.0; // 10.0;
    param.push_back("slp");     var[ "slp" ] = 1.0e+3;
    /*
    param.push_back("A0");      var[  "A0" ] = 75.0 / (1.5 * (0.375 + 0.75) * u * N_A);
    param.push_back("P0");      var[  "P0" ] = (1000.0 / 314.0) *
                                                // 5 ;
                                                // 10;
                                                 20;
                                                // 40;
    */

    // std::string fname = "Foerster2_faser501_5ng.dat";
    // std::string fname = "Foerster2_faser501_10ng.dat";
    // std::string fname = "Foerster2_faser501_20ng.dat";

    //

    // Reaction / rule:
    //
    //      r1  :  k_1 * A * P         A + P      --->   AP
    //      r2  :  k_2 * AP            AP         --->   A + P
    //      r3  :  k_3 * As (* Ps)     As (+ Ps)  --->   APs
    //      r4  :  k_4 * APs           APs        --->   As (+ Ps)
    //      r5  :  k_5 * As (* B)      As (+ B)   --->   AB
    //      r6  :  k_6 * AB            AB         --->   As (+ B)
    //      r7  :  k_7 * A             A          --->   As
    //      r8  :  k_8 * As            As         --->   A

    aux["r1"] = Expression(TIMES, "k_1", Expression(TIMES, "A", "P"));
    aux["r2"] = Expression(TIMES, "k_2", "AP");
    aux["r3"] = Expression(TIMES, "Sw", Expression(TIMES, "k_3", Expression(TIMES, "As", "Ps")));
    aux["r4"] = Expression(TIMES, "Sw", Expression(TIMES, "k_4", "APs"));
    aux["r5"] = Expression(TIMES, "Sw", Expression(TIMES, "k_5", Expression(TIMES, "As", "B")));
    aux["r6"] = Expression(TIMES, "Sw", Expression(TIMES, "k_6", "AB"));
    aux["r7"] = Expression(TIMES, "Sw", Expression(TIMES, "k_7", "A"));
    aux["r8"] = Expression(TIMES, "Sw", Expression(TIMES, "k_8", "As"));

    //

    // ODEs:
    //
    //  if t > -ti:
    //
    //        (A)' = - r1 + r2
    //
    //        (P)' = - r1 + r2
    //
    //       (AP)' = + r1 - r2
    //
    // only if  t >= -t_ini (via Sw value):
    //
    //        (A)' = - r1 + r2 - r7 + r8
    //
    //        (P)' = - r1 + r2
    //
    //       (AP)' = + r1 - r2
    //
    //      (APs)' = + r3 - r4
    //
    //        (B)' = - r5 + r6
    //
    //       (AB)' = + r5 - r6
    //
    //       (As)' = - r3 + r4 - r5 + r6 + r7 - r8
    //
    //       (Ps)' = - r3 + r4

    //

    emap["A"]   = Expression(
                        PLUS,
                        Expression(MINUS, aux["r2"], aux["r1"]),
                        Expression(MINUS, aux["r8"], aux["r7"])
                  );

    emap["P"]   = Expression(
                        MINUS, aux["r2"], aux["r1"]
                        // PLUS,
                        // Expression(TIMES, -1.0, aux["r1"]),
                        // aux["r2"]
                  );

    emap["AP"]  = Expression(
                        MINUS, aux["r1"], aux["r2"]
                  );

    emap["APs"] = Expression(
                        MINUS, aux["r3"], aux["r4"]
                  );

    emap["B"]   = Expression(
                        MINUS, aux["r6"], aux["r5"]
                  );

    emap["AB"]  = Expression(
                        MINUS, aux["r5"], aux["r6"]
                  );

    emap["As"]  = Expression(
                        PLUS,
                        Expression(MINUS, aux["r4"], aux["r3"]),
                        Expression(
                                PLUS,
                                Expression(MINUS, aux["r6"], aux["r5"]),
                                Expression(MINUS, aux["r7"], aux["r8"])
                        )
                  );

    emap["Ps"]  = Expression(
                        MINUS, aux["r4"], aux["r3"]
                  );

    emap["S"]   = Expression(
                        MINUS,
                        "S",
                        Expression(PLUS, "b", Expression(TIMES, "slp", "APs"))
                  );

    emap["Sw"]  = Expression( 0.0 );


    tmap["S"] = 2;



    biosys.setODESystem(emap);
    biosys.setODETypes(tmap);
    // biosys.setSpecies(species);  // has now a complete different effect:
                                    // resets to an identity 'emap' mapping
    biosys.setParameters(param);


    Expression::ParamIterConst y0Beg = foerster[0].y0.begin();
    Expression::ParamIterConst y0End = foerster[0].y0.end();

    for (Expression::ParamIterConst it = y0Beg; it != y0End; ++it)
    {
        biosys.setInitialValue(it->first, it->second);
    }

/*
    Vector meastp(200);
    for (long j=1; j <= meastp.nr(); ++j) meastp(j) = tstart + j*(tend-tstart)/200.0;
    biosys.setMeasurementTimePoints( meastp );
*/

    std::vector<Real> allbp;

    allbp.clear();
    taccu = 0.0;

    for (unsigned j = 0; j < nexp; ++j)
    {
        long nbp = foerster[j].breaktp.nr();

        for (long k = 1; k < nbp; ++k)
        {
            allbp.push_back( taccu + foerster[j].breaktp(k) );
        }

        taccu += foerster[j].breaktp(nbp);
    }

    allbp.push_back( taccu );


    // Breakpoints / Event Management:
    //
    //  Subdivision of integration interval [t0 T] :  [ t0 = b1, b2, b3, ..., bn-1, bn = T ]
    //

//TIME_THIS_TO( std::cerr << " *** Call: biosys.setBreakpoints() *** " << std::endl;

    Vector breaktp(allbp);

    biosys.setBreakpoints( breaktp );

//, std::cerr);

    // and succinctly some event handling
    //

    long jj = 0;

    for (unsigned j = 0; j < nexp; ++j)
    {
        emap = biosys.getEvent(jj+1);

        emap["Sw"] = Expression( 1.0 );

        if ( j == 0 )
        {
            biosys.setEvent(jj+1, emap);
        }
        else
        {
            y0Beg = foerster[j].y0.begin();
            y0End = foerster[j].y0.end();

            aux.clear();

            for (Expression::ParamIterConst it = y0Beg; it != y0End; ++it)
            {
                aux[it->first] = Expression( it->second );
            }

            biosys.setEvent(jj,   aux);
            biosys.setEvent(jj+1, emap);
        }

        jj += (foerster[j].breaktp.nr() - 1);
    }

    //
    //
    //

    std::cout << std::endl;
    for (long j = 0; j < breaktp.nr(); ++j)
    {
        BioRHS::ExpressionMap icmap = biosys.getEvent(j);
        BioRHS::EMapIterConst cBeg = icmap.begin();
        BioRHS::EMapIterConst cEnd = icmap.end();

        std::cout << "Event " << std::setw(2) << j+1 <<
                     "   (at time t = " << breaktp(j+1) << ")" << std::endl;
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
    tmap = biosys.getODETypes();

    for (unsigned j = 0; j < emap.size(); ++j)
    {
        std::cout << std::endl;
        std::cout << ((tmap[species[j]] == 2) ? " 0 * " : "     ");
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
    Real rTol = 1.0e-4;
    Real aTol = 1.0e-8;
    Real xTol = 1.0e-4;

    biosys.setSolverRTol(rTol);
    biosys.setSolverATol(aTol);


    //-----------------------------------------------------------------------------

    Vector syndata, synscal, vref, utmp;

TIME_THIS_TO( std::cerr << " *** Call: biosys.computeModel() *** " << std::endl;

//for (long n=0; n < 100; ++n)
    vref = biosys.computeModel(var, "adaptive");
    // syndata = biosys.computeModel(var, "adaptive");

, std::cerr );

    Vector meastp = biosys.getOdeTrajectoryTimePoints();

std::cerr << " meastp = " << std::endl;
std::cerr << meastp.t();
std::cerr << std::endl;

std::cerr << " vref = " << std::endl;
std::cerr << vref;
std::cerr << std::endl;

std::cerr << " *** Retn: biosys.computeModel() *** " << std::endl;

    //
    //
    //

    Gnuplot g1("lines");

//    //
//    // Styles
//    //
//    g1.reset_plot();
//    std::cout << std::endl << std::endl;
//    std::cout << "*** showing styles" << std::endl;
//
//    std::cout << "sine in points" << std::endl;
//    g1.set_pointsize(0.8).set_style("linespoints");
//    g1.plot_equation("sin(x)","lpoints");
//
//    std::cout << "sine in impulses" << std::endl;
//    g1.set_style("impulses");
//    g1.plot_equation("sin(x)","impulses");
//
//    std::cout << "sine in steps" << std::endl;
//    //g1.set_style("steps");
//    g1.set_style("steps").plot_equation("sin(x)","steps");
//    //g1.replot();

    //

    BioSystem::Species    spe   = biosys.getSpecies();
    ODESolver::Trajectory tra   = biosys.getOdeTrajectory();
    Vector                tratp = biosys.getOdeTrajectoryTimePoints();
    std::vector<Real>     gnutp;

    gnutp.clear();
    for (long j = 0; j < tratp.nr(); ++j)
    {
        gnutp.push_back( tratp(j+1) );
    }

    g1.reset_plot();
    g1.set_grid();
    for (unsigned j = 0; j < spe.size(); ++j)
    {
        if ( spe[j] != "S" )
        {
            g1.plot_xy( gnutp, tra[j], spe[j] );
        }
    }

    //g1.replot();

    //

/*
    std::cout << std::endl;
    std::cout << "Press ENTER to close Gnuplot window." << std::endl;

    while ( std::cin.get() != '\n' )
        ;

exit(-999);
*/

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

    // BioSystem::MeasurementList  measlist = biosys.getMeasurementList();
    BioSystem                   invBiosys( tstart, tend );

    invBiosys.setSolverRTol( solverRTol );
    invBiosys.setSolverATol( solverATol );

    invBiosys.setODESystem(emap);
    invBiosys.setODETypes(tmap);
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
/*
    Real                        sigma = 0.025;   // add 2.5% white noise ...
    unsigned                    sampleSize = 3;
    unsigned                    nSpecies = 2;

    newMeas.clear();
    newtp.zeros( sampleSize*meastp.nr() );

    for (unsigned tp = 0; tp < (unsigned) meastp.nr(); ++tp )
    {
        for (unsigned k = 0; k < sampleSize; ++k) {
            newtp( sampleSize*tp + k + 1 ) = meastp(tp+1);
            newMeas.push_back( MeasurementPoint() );
        }

        // for (unsigned j = 0; j < species.size(); ++j)
        for (unsigned j = 0; j < nSpecies; ++j)
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
*/

/*
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
    // for (unsigned j = 0; j < species.size(); ++j)
    for (unsigned j = 0; j < nSpecies; ++j)
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
        // std::cout << std::fixed << std::setw(10) << meastp(tp+1);
        // out << std::fixed << std::setw(10) << meastp(tp+1);
        std::cout << std::scientific << std::setprecision(3) << std::setw(11) << meastp(tp+1);
        out << std::scientific << std::setprecision(3) << std::setw(11) << meastp(tp+1);
        // for (unsigned j = 0; j < species.size(); ++j)
        for (unsigned j = 0; j < nSpecies; ++j)
        {
            if ( measlist[tp].count(species[j]) > 0 )
            {
                std::cout << "\t";
                std::cout << std::scientific << std::setprecision(6) << measlist[tp][species[j]].first;
                out << "\t";
                out << std::scientific << std::setprecision(6) << measlist[tp][species[j]].first;
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
*/

//    invBiosys.computeModel(var);
//
//    std::cout << "####### invBiosys.computeModel() results #######" << std::endl;
//    for (unsigned j = 0; j < species.size(); ++j)
//    {
//        std::cout << "  trajectory for \"" << species[j] << "\" = " << std::endl;
//        std::cout << invBiosys.getSimTrajectoryPoints(species[j]).t() << std::endl;
//    }
//    std::cout << "################################################" << std::endl;

    ODESolver::Grid    mRawMeas;
    ODESolver::Grid    mTimepoint;
    BioSystem::Species mSpecies;
    long               currtp = 0;

    newMeas.clear();

    for (unsigned jexp = 0; jexp < nexp; ++jexp)
    {
        Real          toff  = foerster[jexp].tsync;
        std::string   fname = foerster[jexp].meas_fname;
        std::ifstream in;

        in.open(fname.c_str());
        // in.open("Foerster2_faser501_5ng.dat");
        // in.open("testfoerster_simple.csv");

        if ( !in.is_open() )
        {
            std::cerr << "### ERROR: Could not open measurement data file '" << fname << "'.\n";
            return -4;
        }

        std::string line;

        mSpecies.clear();
        mSpecies.push_back("S");

        while ( std::getline(in, line) )
        {
            std::stringstream ss(line);
            std::string       id, unit, field;

            /*
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
            */
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
                        mTimepoint.push_back(toff + measVal);
                        newMeas.push_back( MeasurementPoint() );
                        currtp = newMeas.size()-1;
                    }
                    else
                    {
                        mRawMeas.push_back(measVal);
                        newMeas[currtp][mSpecies[j]] = std::make_pair<Real,Real>(measVal, 0.0);
                    }

                    ++j;
                }

                if ( newMeas[currtp].empty() )
                {
                    newMeas.pop_back();
                    mTimepoint.pop_back();
                    mRawMeas.pop_back();
                }

                std::cout << std::endl;
            }
        }

        in.close();
    }

    std::cout << "Read " << newMeas.size() << " timepoint(s) with measurement data." << std::endl;

    invBiosys.setMeasurementList( /*Vector*/(mTimepoint), newMeas );

    //
    //
    //

    Gnuplot g2("linespoints");

    for (unsigned j = 0; j < spe.size(); ++j)
    {
        if ( spe[j] == "S" )
        {
            g2.set_grid().set_pointsize(1.5);
            g2.plot_xy( gnutp, tra[j], spe[j] );
        }
    }

    g2.set_style("points").set_pointsize(2.0);
    g2.plot_xy( mTimepoint, mRawMeas, "Measurements" );

    //

    /*
    std::cout << "Press ENTER to close Gnuplot window." << std::endl;

    while ( std::cin.get() != '\n' )
        ;
    */

    //
    // Initial guess for GaussNewton
    //

    long q = 9; // 8; // param.size();

    par.clear();
    par1.clear();
    p.zeros( q );
    pscal.zeros( q );

    for (long j = 1; j <= q; ++j)
    {
        par1.push_back( param[j-1] );
        p(j) = pscal(j) = var[ par1.back() ];
             par[ par1.back() ] = p(j);
        parThres[ par1.back() ] = std::fabs( p(j) ) * 1.0e-3; // 1.0e-2; // EPMACH;
    }

    // p(5) = pscal(5) = 1.0;      var[ par1[4] ] = p(5);
    // p(6) = pscal(6) = 1.0;      var[ par1[5] ] = p(6);

    BioProcessor proc( &invBiosys, "nlscon" );
    IOpt         iopt = proc.getIOpt();

//    GaussNewtonWk  wk;
//    GaussNewton    gn;
//    BioPAR         prob( &invBiosys, par1 );

    iopt.mode        = 0;     // 0:normal run, 1:single step
    iopt.jacgen      = 3;     // 1:user supplied Jacobian, 2:num.diff., 3:num.diff.(with feedback)
    iopt.qrank1      = false;     // allow Broyden rank-1 updates if __true__
    iopt.nonlin      = 3;     // 1:linear, 2:mildly nonlin., 3:highly nonlin., 4:extremely nonlin.
    iopt.rscal       = 1;     // 1:use unchanged fscal, 2:recompute/modify fscal, 3:use automatic scaling only
    iopt.norowscal   = false;     // allow for automatic row scaling of Jacobian if __false__
    iopt.transf      = 1;
    // iopt.lpos      = true;      // force solution vector to be positive (all components > 0.0)
                            //          _mprmon =   0      1      2      3      4       5       6
                            //  dlib::log_level =  LNONE  LINFO  LVERB  LTALK  LGABBY  LDEBUG  LTRACE
    iopt.mprmon    = 2;
    iopt.mprerr    = 1;

    iopt.itmax = 75;


    // if ( iopt.jacgen > 1 )
    //    wk.cond = 1.0 / /*std::sqrt*/(solverRTol);
    // wk.nitmax = 15;


    proc.setIOpt( iopt );

    proc.setCurrentParamValues( par );
    proc.setCurrentParamThres( parThres );

    Expression::Param speThres( proc.getCurrentSpeciesThres() );
    for (Expression::ParamIterConst it = speThres.begin();
                                    it != speThres.end(); ++it)
    {
        speThres[it->first] = 0.1; // solverRTol * 10.0; // EPMACH;
    }
    speThres["S"] = 1000.0;

    proc.setCurrentSpeciesThres( speThres );



    Expression::Param finalParam;
    int               rc = proc.identifyParameters( reconXTol );

    // if ( 0 == rc )
    {
        finalParam = proc.getIdentificationResults();
    }
    /*
    Expression::Param finalParam = var;
    for (long j = 1; j <= psol.nr(); ++j)
    {
        finalParam[ par1[j-1] ] = psol(j);
    }
    */

    std::cout.unsetf( std::ios_base::floatfield );
    std::cout << std::scientific << std::setprecision(4);
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
    {
        std::string pname = param[j];

        std::cout << std::setw(20) << std::right << pname << "      " <<
                     std::setw(15) << std::left << finalParam[pname] <<
                     std::setw(20) << std::left << var[pname] <<
                     std::endl;

        if ( finalParam[pname] == 0.0 )
        {
            finalParam[pname] = var[pname];
        }
    }

    //

    BioProcessor::TrajectoryMap traMap;

    proc.setCurrentParamValues(finalParam);

    traMap = proc.computeModel();
    tratp  = proc.getAdaptiveTimepoints();

    gnutp.clear();

    for (long j = 0; j < tratp.nr(); ++j)
    {
        gnutp.push_back( tratp(j+1) );
    }


    for (unsigned j = 0; j < spe.size(); ++j)
    {
        if ( spe[j] == "S" )
        {
            g2.set_style("linespoints").set_pointsize(1.5);
            g2.plot_xy( gnutp, traMap[spe[j]], spe[j] );
        }
    }

    //

    std::cout << std::endl;
    std::cout << "Press ENTER to close Gnuplot window." << std::endl;

    while ( std::cin.get() != '\n' )
        ;

    return rc;
}
