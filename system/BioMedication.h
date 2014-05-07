// Copyright (C) 2010 - 2014
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2014-01-17 td
// last changed:
//
#ifndef __BIO_MEDICATION_H
#define __BIO_MEDICATION_H

#include <utility>  // for `std::pair'
#include <string>
#include <vector>

#include "common/Types.h"
//#include "linalg/Matrix.h"
//#include "linalg/Vector.h"

namespace PARKIN
{
    typedef std::vector< std::pair<Real,Real> > DosingList;

    struct Medication
    {
        long         index;
        std::string  adm_comp;
        Real         beta;
        Real         clRate;
        DosingList   dosing;
    };

    ///

    class BioMedication
    {
        public:

            typedef std::vector< Medication >   MedicationList;

            /// c'tor
            BioMedication();
            BioMedication(MedicationList const& medList);
            /// d'tor
            ~BioMedication();

            BioMedication(BioMedication const& other);
            BioMedication const& operator= (BioMedication const& other);

            //

            MedicationList const&
            getMedicationList() const
            {
                return _adm;
            }

            //

            void
            setMedicationList(MedicationList const& medList)
            {
                _adm = medList;
            }

            void
            addMedication(Medication const& med)
            {
                _adm.push_back(med);
            }

            //

            // Vector
            void
            computeMedication(double* y, long n, double* dy);

        private:

            MedicationList _adm;

    };

}

#endif // __BIO_MEDICATION_H
