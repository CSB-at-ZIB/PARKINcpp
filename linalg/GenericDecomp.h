// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 29.03.2010 td
// last changed:
//

#ifndef __GENERIC_DECOMP_H
#define __GENERIC_DECOMP_H

#include <cmath>

#include "common/Types.h"
#include "Matrix.h"
#include "Vector.h"

namespace PARKIN
{
    // static ErrType dcmpErr(0);

    class GenericDecomp
    {
        public:
            // Contructor
            GenericDecomp() : _err(new ErrType(0)), _valid(false) { }

            // Destructor
            virtual ~GenericDecomp() { delete _err; _err = 0; }

            // Initiate decomposition of matrix
            virtual ErrType const& decompose(Matrix const& A) = 0;

            // Solve for right-hand side vector b
            virtual ErrType const& solve(Matrix const& A, Vector const& b, Vector& x) = 0;
            virtual ErrType const& solve(Vector& b, Vector& x) const = 0;

            //
            virtual void getFirstFactor(Matrix& M1) const = 0;
            virtual void getSecondFactor(Matrix& M2) const = 0;

            //
            bool isValid() const { return _valid; }
            ErrType const& getError() const { return *_err; }

            // Copy constructor
            GenericDecomp(GenericDecomp const& mat)
            {
                if (this != &mat)
                {
                    int ierr = mat.getError().getIerr();
                    std::string msg = mat.getError().getMsg();

                    _err = new ErrType(ierr,msg);
                    _valid = mat._valid;
                }
            }

            // Assignment opaerator
            GenericDecomp const& operator= (GenericDecomp const& mat)
            {
                if (_err == 0) _err = new ErrType(0);
                if (this != &mat)
                {
                    int ierr = mat.getError().getIerr();
                    std::string msg = mat.getError().getMsg();

                    _err->setIerr(ierr);
                    _err->setMsg(msg);
                    _valid = mat._valid;
                }
                return *this;
            }

        private:
//            // Copy constructor
//            GenericDecomp(GenericDecomp const& mat) : _err(mat._err), _valid(mat.isValid()) { }
//            // Assignment opaerator
//            GenericDecomp const& operator= (GenericDecomp const&) { return *this; }

        protected:
            ErrType*    _err;
            bool        _valid;
    };

}

#endif // __GENERIC_DECOMP_H

