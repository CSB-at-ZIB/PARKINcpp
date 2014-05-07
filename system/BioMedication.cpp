// Copyright (C) 2010 - 2014
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2014-01-17 td
// last changed:
//
#include <cmath>

#include "BioMedication.h"

using namespace PARKIN;

//----------------------------------------------------------------------------
BioMedication::BioMedication() :
    _adm()
{
}
//----------------------------------------------------------------------------
BioMedication::BioMedication(MedicationList const& medList) :
    _adm()
{
    _adm = medList;
}
//----------------------------------------------------------------------------
BioMedication::~BioMedication()
{
}
//----------------------------------------------------------------------------
BioMedication::BioMedication(BioMedication const& m)
{
    if ( this != &m )
    {
        _adm = m._adm;
    }
}
//----------------------------------------------------------------------------
BioMedication const&
BioMedication::operator= (BioMedication const& m)
{
    if ( this != &m )
    {
        _adm = m._adm;
    }

    return *this;
}
//----------------------------------------------------------------------------
void
//Vector
BioMedication::computeMedication(double* y, long n, double* dy)
{
    // Vector df;

    if ( _adm.size() <= 0 )
    {
        return;
    }

    if ( n == 0 )
    {
        return;
        // return df;
    }

    double    t = y[0];            // current time in the simulation
    int    nAdm = _adm.size();     // number of medication administrations

    for (int j = 0; j < nAdm; ++j)
    {
        double     sumDosing = 0.0;
        DosingList dosing    = _adm[j].dosing;

        long   nDose    = dosing.size();
        long   indexAdm = _adm[j].index;
        double clRate   = _adm[j].clRate;
        double beta     = _adm[j].beta;

        if ( indexAdm <= 0 )
        {
            continue;
        }

        dy[indexAdm] = 0.0;

        for (long k = 0; k < nDose; ++k)
        {
            double tAccu = t - dosing[k].first;  // first: dosing time point,  second: actual dose

            if ( tAccu > 0.0 )
            {
                // Medication release according to Gamma distribution
                // (with shape parameter alpha := 2), and
                // using the superposition principle for repeated administration
                //
                sumDosing += dosing[k].second * (beta * beta * tAccu * exp( -beta * tAccu ) );
            }
        }

        dy[indexAdm] = sumDosing - clRate * y[indexAdm];
    }

    // return df;
}
//----------------------------------------------------------------------------
