// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-02-21 td
// last changed:
//

//#include <sstream>
#include "BioRHS.h"

using namespace PARKIN;

//----------------------------------------------------------------------------
BioRHS::BioRHS() :
    _rhs_type(),
    _rhs(), _drhs(), _species(), _parameter(),
    _num_table_entries(2), _data_table(new double*[2])
    /// , _Fz(0), _Fp(0)
{
    _data_table[0] = _data_table[1] = 0;
}
////----------------------------------------------------------------------------
//BioRHS::BioRHS(long n) :
//    _rhs(), _drhs()
//{
//    _species.clear();
//
//    for (long j = 1; j <= n; ++j)
//    {
//        std::ostringstream s;
//        s << "_" << j << "_";
//        _species.push_back( s.str() );
//    }
//
//    _parameter.clear();
//}
//----------------------------------------------------------------------------
BioRHS::BioRHS(ExpressionMap const& emap) :
    _rhs_type(),
    _rhs(emap), _drhs(), _species(), _parameter(),
    _num_table_entries(2), _data_table(new double*[2])
    /// _Fz(0), _Fp(0)
{
    EMapIterConst eBeg = emap.begin();
    EMapIterConst eEnd = emap.end();

    _species.clear();
    _rhs_type.clear();

    for (EMapIterConst it = eBeg; it != eEnd; ++it)
    {
        _species.push_back( it->first );
        _rhs_type[it->first] = 1;
    }

    StrIterConst sBeg = _species.begin();
    StrIterConst sEnd = _species.end();
    for (EMapIterConst it = eBeg; it != eEnd; ++it)
    {
        // species offsets start at 1 (0 is for time!)
        long       a = 0;
        Expression expr = it->second;

        expr.off("odeTime", a++);

        for (StrIterConst itSpe = sBeg; itSpe != sEnd; ++itSpe)
            expr.off(*itSpe, a++);
    }

    _parameter.clear();

    _data_table[0] = _data_table[1] = 0;

    /// reserveJacobianMemory();
    computeDerivativeExpression();
}
//----------------------------------------------------------------------------
BioRHS::BioRHS(Species const& species) :
    _rhs_type(),
    _rhs(), _drhs(), _species(), _parameter(),
    _num_table_entries(2), _data_table(new double*[2])
    /// , _Fz(0), _Fp(0)
{
    StrIterConst sBeg = species.begin();
    StrIterConst sEnd = species.end();

    _rhs_type.clear();
    _rhs.clear();

    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        _rhs_type[*it] = 1;
        _rhs[*it] = Expression(*it);
    }


    EMapIterConst rBeg = _rhs.begin();
    EMapIterConst rEnd = _rhs.end();

    _species.clear();

    for (EMapIterConst it = rBeg; it != rEnd; ++it)
    {
        _species.push_back( it->first );
    }

    StrIterConst spBeg = _species.begin();
    StrIterConst spEnd = _species.end();
    for (EMapIterConst it = rBeg; it != rEnd; ++it)
    {
        // species offsets start at 1 (0 is for time!)
        long       a = 0;
        Expression expr = it->second;

        expr.off("odeTime", a++);

        for (StrIterConst itSpe = spBeg; itSpe != spEnd; ++itSpe)
            expr.off(*itSpe, a++);
    }

    _parameter.clear();

    _data_table[0] = _data_table[1] = 0;

    /// reserveJacobianMemory();
    computeDerivativeExpression();
}
//----------------------------------------------------------------------------
BioRHS::BioRHS(ExpressionMap const& emap,
               Parameter const& param) :
    _rhs_type(),
    _rhs(emap), _drhs(), _species(), _parameter(param),
    _num_table_entries(2), _data_table(new double*[2])
    /// , _Fz(0), _Fp(0)
{
    EMapIterConst eBeg = emap.begin();
    EMapIterConst eEnd = emap.end();

    _species.clear();
    _rhs_type.clear();

    for (EMapIterConst it = eBeg; it != eEnd; ++it)
    {
        _species.push_back( it->first );
        _rhs_type[it->first] = 1;
    }

    StrIterConst sBeg = _species.begin();
    StrIterConst sEnd = _species.end();
    StrIterConst vBeg = _parameter.begin();
    StrIterConst vEnd = _parameter.end();
    for (EMapIterConst it = eBeg; it != eEnd; ++it)
    {
        // species offsets start at 1 (0 is for time!)
        long       a = 0;
        long       b = 0;
        Expression expr = it->second;

        expr.off("odeTime", a++);

        for (StrIterConst itSpe = sBeg; itSpe != sEnd; ++itSpe)
            expr.off(*itSpe, a++);
        for (StrIterConst itPar = vBeg; itPar != vEnd; ++itPar)
            expr.off(*itPar, --b);
    }

    _data_table[0] = _data_table[1] = 0;

    /// reserveJacobianMemory();
    computeDerivativeExpression();
}
//----------------------------------------------------------------------------
BioRHS::~BioRHS()
{
    // for (long j = 1; j < _num_table_entries; ++j) delete[] _data_table[j];
    delete[] _data_table;

    /// delete[] _Fz;
    /// delete[] _Fp;
}
//----------------------------------------------------------------------------
BioRHS::BioRHS(BioRHS const& r)
{
    if ( this != &r )
    {
        _rhs_type  = r._rhs_type;
        _rhs       = r._rhs;
        _drhs      = r._drhs;
        _species   = r._species;
        _parameter = r._parameter;

        _num_table_entries = r._num_table_entries;
        _data_table = new double*[_num_table_entries];
        for (long j = 0; j < _num_table_entries; ++j)
            _data_table[j] = r._data_table[j];

        /// _Fz = _Fp = 0;
        /// reserveJacobianMemory();
        // computeDerivativeExpression();
    }
}
//----------------------------------------------------------------------------
BioRHS const&
BioRHS::operator= (BioRHS const& r)
{
    if ( this != &r )
    {
        _rhs_type  = r._rhs_type;
        _rhs       = r._rhs;
        _drhs      = r._drhs;
        _species   = r._species;
        _parameter = r._parameter;

        // for (long j = 1; j < _num_table_entries; ++j) delete[] _data_table[j];
        delete[] _data_table;

        _num_table_entries = r._num_table_entries;
        _data_table = new double*[_num_table_entries];
        for (long j = 0; j < _num_table_entries; ++j)
            _data_table[j] = r._data_table[j];

        /// _Fz = _Fp = 0;
        /// reserveJacobianMemory();
        // computeDerivativeExpression();
    }

    return *this;
}
//----------------------------------------------------------------------------
/*
void
BioRHS::reserveJacobianMemory()
{
    long n = _species.size() + 1L;
    long q = _parameter.size();

    delete[] _Fz;
    delete[] _Fp;

    _Fz = new double[n*n];

    if ( q > 0L )
    {
        _Fp = new double[n*q];
    }
}
*/
//----------------------------------------------------------------------------
void
BioRHS::computeDerivativeExpression() // Expression::Param const& var)
{
    EMapIterConst   rBeg = _rhs.begin();
    EMapIterConst   rEnd = _rhs.end();
    StrIterConst    sBeg = _species.begin();
    StrIterConst    sEnd = _species.end();
    StrIterConst    pBeg = _parameter.begin();
    StrIterConst    pEnd = _parameter.end();

    _drhs.clear();

    for (EMapIterConst itRHS = rBeg; itRHS != rEnd; ++itRHS)
    {
        Expression expr = itRHS->second;

        for (StrIterConst itSpe = sBeg; itSpe != sEnd; ++itSpe)
        {
            // species offsets start at 1 (and not at 0 !)
            long        a = 0;
            long        b = 0;
            std::string s = (itRHS->first) + " / " + *itSpe;

            _drhs[s] = expr.df( *itSpe );

            _drhs[s].off("odeTime", a++);

            for (StrIterConst it = sBeg; it != sEnd; ++it)
                _drhs[s].off(*it,a++);
            for (StrIterConst it = pBeg; it != pEnd; ++it)
                _drhs[s].off(*it,--b);
        }

        for (StrIterConst itPar = pBeg; itPar != pEnd; ++itPar)
        // for (Expression::ParamIterConst itPar = vBeg; itPar != vEnd; ++itPar)
        {
            // species offsets start at 1 (and not at 0 !)
            long        a = 0;
            long        b = 0;
            std::string s = (itRHS->first) + " / " + *itPar; // itPar->first;

            _drhs[s] = expr.df( *itPar );

            _drhs[s].off("odeTime", a++);

            for (StrIterConst it = sBeg; it != sEnd; ++it)
                _drhs[s].off(*it,a++);
            for (StrIterConst it = pBeg; it != pEnd; ++it)
                _drhs[s].off(*it,--b);
        }
    }
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void
BioRHS::setRHSType(ExprTypeMap const& tmap)
{
    TMapIterConst tBeg = tmap.begin();
    TMapIterConst tEnd = tmap.end();

    for (TMapIterConst it = tBeg; it != tEnd; ++it)
    {
        _rhs_type[it->first] = it->second;
    }

    return;
}
//----------------------------------------------------------------------------
/*
BioRHS::ExprTypeMap const&
BioRHS::getRHSType() const
{
    return _rhs_type;
}
*/
//----------------------------------------------------------------------------
void
BioRHS::setRHS(ExpressionMap const& emap)
{
    _rhs = emap;

    _species.clear();
    _rhs_type.clear();

    EMapIterConst rBeg = _rhs.begin();
    EMapIterConst rEnd = _rhs.end();

    for (EMapIterConst it = rBeg; it != rEnd; ++it)
    {
        _species.push_back( it->first );
        _rhs_type[ it->first ] = 1;
    }

    StrIterConst sBeg = _species.begin();
    StrIterConst sEnd = _species.end();
    StrIterConst pBeg = _parameter.begin();
    StrIterConst pEnd = _parameter.end();

    for (EMapIterConst it = rBeg; it != rEnd; ++it)
    {
        // species offsets start at 1 (0 is for time!)
        long       a = 0;
        long       b = 0;
        Expression expr = it->second;

        expr.off("odeTime", a++);

        for (StrIterConst itSpe = sBeg; itSpe != sEnd; ++itSpe)
            expr.off(*itSpe, a++);
        for (StrIterConst itPar = pBeg; itPar != pEnd; ++itPar)
            expr.off(*itPar, --b);
    }

    // _parameter.clear();

    /// reserveJacobianMemory();
    computeDerivativeExpression();
}
//----------------------------------------------------------------------------
void
BioRHS::resetSpecies(BioRHS::Species const& species)
{
    StrIterConst sBeg = species.begin();
    StrIterConst sEnd = species.end();

    _rhs.clear();
    _rhs_type.clear();

    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        _rhs_type[*it] = 1;
        _rhs[*it] = Expression(*it);
    }


    EMapIterConst rBeg = _rhs.begin();
    EMapIterConst rEnd = _rhs.end();

    _species.clear();

    for (EMapIterConst it = rBeg; it != rEnd; ++it)
    {
        _species.push_back( it->first );
    }

    StrIterConst spBeg = _species.begin();
    StrIterConst spEnd = _species.end();
    for (EMapIterConst it = rBeg; it != rEnd; ++it)
    {
        // species offsets start at 1 (0 is for time!)
        long       a = 0;
        Expression expr = it->second;

        expr.off("odeTime", a++);

        for (StrIterConst itSpe = spBeg; itSpe != spEnd; ++itSpe)
            expr.off(*itSpe, a++);
    }

    // _parameter.clear();

    /// reserveJacobianMemory();
    computeDerivativeExpression();
}
//----------------------------------------------------------------------------
BioRHS::Species const&
BioRHS::getSpecies() const
{
    return _species;
}
//----------------------------------------------------------------------------
void
BioRHS::setParameters(Expression::Param const& param)
{
    Expression::ParamIterConst pBeg = param.begin();
    Expression::ParamIterConst pEnd = param.end();

    _parameter.clear();

    for (Expression::ParamIterConst it = pBeg; it != pEnd; ++it)
        _parameter.push_back( it->first );

    EMapIterConst rBeg = _rhs.begin();
    EMapIterConst rEnd = _rhs.end();
    StrIterConst  vBeg = _parameter.begin();
    StrIterConst  vEnd = _parameter.end();
    for (EMapIterConst it = rBeg; it != rEnd; ++it)
    {
        long       b = 0;
        Expression expr = it->second;

        for (StrIterConst itPar = vBeg; itPar != vEnd; ++itPar)
            expr.off(*itPar, --b);
    }

    /// reserveJacobianMemory();
    computeDerivativeExpression();
}
//----------------------------------------------------------------------------
void
BioRHS::setParameters(Parameter const& param)
{
    // _parameter.clear();
    _parameter = param;

    EMapIterConst rBeg = _rhs.begin();
    EMapIterConst rEnd = _rhs.end();
    StrIterConst  vBeg = _parameter.begin();
    StrIterConst  vEnd = _parameter.end();
    for (EMapIterConst it = rBeg; it != rEnd; ++it)
    {
        long       b = 0;
        Expression expr = it->second;

        for (StrIterConst itPar = vBeg; itPar != vEnd; ++itPar)
            expr.off(*itPar, --b);
    }

    /// reserveJacobianMemory();
    computeDerivativeExpression();
}
//----------------------------------------------------------------------------
BioRHS::Parameter const&
BioRHS::getParameters() const
{
    return _parameter;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void
BioRHS::b(double* y, int* nz, double* B, int* ir, int* ic)
{
    int            k   = 0;
    int            idx = 0;
    StrIterConst   sBeg = _species.begin();
    StrIterConst   sEnd = _species.end();

    // first component is for time variable!
     B[k] = 1.0;
    ir[k] = ic[k] = ++idx;
    ++k;

    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        ++idx;

        switch ( _rhs_type[*it] )
        {
            /*
            case 3 :
                --idx;
            */
            case 2 :
                break;

            case 1 :
            default:
                 B[k] = 1.0;
                ir[k] = ic[k] = idx;

                ++k;

                break;
        };
    }

    *nz = k;

    return;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
Vector
BioRHS::f(Expression::Param& par,
          long n, double* dy)
{
    /// long                rules = 0;
    Vector              dz; // ( _species.size() );
    StrIterConst        sBeg = _species.begin();
    StrIterConst        sEnd = _species.end();
    /// Expression::Param   aux = par;

/*
    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        if ( _rhs_type[*it] == 3 )
        {
            ++rules;
            par[*it] = _rhs[*it].eval(aux);
        }
    }
*/
    if ( n == 0 )
    {
        // species index start at 2 (1 is for time!)
        long j = 0;
        dz.zeros( _species.size() + 1 ); // - rules );

        dz(++j) = 1.0;

        for (StrIterConst it = sBeg; it != sEnd; ++it)
        {
            /*
            if ( _rhs_type[*it] == 3 )
            {
                dz(++j) = 0.0;
            }
            else
            */
            {
                dz(++j) = _rhs[*it].eval(par);
            }
        }

        return dz;
    }

    ///

    *(dy++) = 1.0;

    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        /*
        if ( _rhs_type[*it] == 3 )
        {
            *(dy++) = 0.0;
        }
        else
        */
        {
            *(dy++) = _rhs[*it].eval(par);
        }
    }

    return dz;
}
//----------------------------------------------------------------------------
Vector
BioRHS::f(double* y,
          long n, double* dy)
{
    /// long           k, rules = 0;
    Vector         dz; // ( _species.size() );
    StrIterConst   sBeg = _species.begin();
    StrIterConst   sEnd = _species.end();
    /// double         z[_species.size()];

    _data_table[0] = y;

/*
    k = 0;
    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        if ( _rhs_type[*it] == 3 )
        {
            ++rules;
            z[k] = _rhs[*it].eval( _data_table );
        }

        ++k;
    }

    k = 0;
    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        if ( _rhs_type[*it] == 3 )
        {
            _data_table[0][k] = z[k];
        }

        ++k;
    }
*/

    if ( n == 0 )
    {
        // species index start at 2 (1 is for time!)
        long   j = 0;
        dz.zeros( _species.size() + 1 ) ; // - rules );

        dz(++j) = 1.0;  // dummy equation y' = 1 (time!)

        for (StrIterConst it = sBeg; it != sEnd; ++it)
        {
            /*
            if ( _rhs_type[*it] == 3 )
            {
                dz(++j) = 0.0;
            }
            else
            */
            {
                dz(++j) = _rhs[*it].eval( _data_table );
            }
        }

        return dz;
    }

    ///

    *dy++ = 1.0;  // empty dummy equation y' = 1 (time!)

    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        /*
        if ( _rhs_type[*it] == 3 )
        {
            *(dy++) = 0.0;
        }
        else
        */
        {
            *(dy++) = _rhs[*it].eval( _data_table );
        }
    }

    return dz;
}
//----------------------------------------------------------------------------
Matrix
BioRHS::Jf(Expression::Param& par,
           long n, double* J, long ldJ)
{
    Matrix Fz; // (n,n);

    if ( n == 0 )
    {
        // long
        n = _species.size() + 1L;
        Fz.zeros( n, n );

        Fz(1,1) = 0.0;

        for (long j = 2; j <= n; ++j)
        {
            for (long k = 2; k <= n; ++k)
            {
                std::string s = _species[j-2] + " / " + _species[k-2];

                Fz(j,k) = _drhs[s].eval(par);
            }
        }

        return Fz;
    }


    ///


    // J[ 0 ] = 0.0;

    for (long j = 0; j < n; ++j)
    {
        J[ j ] = 0.0;
    }


    for (long k = 1; k < n; ++k)    // swapping inner and outer loop for speed
    {
        long kk = ldJ*k;

        J[ kk ] = 0.0;

        for (long j = 1; j < n; ++j)  // j still runs over the rows in mind (!)
        {
            std::string s = _species[j-1] + " / " + _species[k-1];

            J[ kk + j ] = _drhs[s].eval(par);
        }
    }

    return Fz;
}
//----------------------------------------------------------------------------
Matrix
BioRHS::df(Expression::Param const& var, Matrix const& Zp,
           Expression::Param& par)
{
    long                n = Zp.nr();
    long                q = Zp.nc();
    Matrix              FzWork; // ( n, n );
    Matrix              FpWork; // ( n, q );

    if ( ( (n-1) > (long)_species.size()   ) ||
         (     q > (long)_parameter.size() )
       )
    {
       return Matrix();
    }

    FzWork.zeros(n,n);
    FpWork.zeros(n,q);

    FzWork(1,1) = 0.0;

    for (long j = 2; j <= n; ++j)
    {
        std::string s = _species[j-2] + " / ";

        for (long k = 2; k <= n; ++k)
        {
            std::string t = s + _species[k-2];

            FzWork(j,k) = _drhs[t].eval(par);
        }

        Expression::ParamIterConst it = var.begin();

        for (long k = 1; k <= q; ++k)
        {
            // std::string t = s + _parameter[k-1];
            std::string t = s + ((it++)->first);

            FpWork(j,k) = _drhs[t].eval(par);
        }
    }

    // Matrix dZ = FzWork * Zp + FpWork;

    return  FzWork * Zp  +  FpWork;  //  dZ;
}
//----------------------------------------------------------------------------
Matrix
BioRHS::df(Expression::Param const& var, double* Zp, double* y,
           long n, long q, double* dZ)
{
    Matrix DummyMat;
    // Matrix              Fz( n, n );
    // Matrix              Fp( n, q );

    if ( ( (n-1) > (long) _species.size()) ||
         (     q > (long) var.size())
       )
    {
       return DummyMat;
    }


    double* FzSave = new double[n*n];
    //obsolete: double* ZpSave = new double[n*q];
    double* FpSave = new double[n*q];

    double* FzWork    = FzSave; // _Fz;
    //obsolete: double* Zptmp = ZpSave;
    double* FpWork    = FpSave; // _Fp;


    _data_table[0] = y;


    *FzWork++  = 0.0;

    for (long k = 1; k < n; ++k)
    {
        *FzWork++ = 0.0;
    }

    for (long k = 0; k < q; ++k)
    {
        *FpWork++ = 0.0;

        //obsolete: *Zptmp++ = *Zp++;
    }


    for (long j = 1; j < n; ++j)
    {
        std::string s = _species[j-1] + " / ";

        *FzWork++ = 0.0;

        for (long k = 1; k < n; ++k)
        {
            std::string t = s + _species[k-1];

            *FzWork++ = _drhs[t].eval( _data_table );
        }


        Expression::ParamIterConst it = var.begin();

        for (long k = 0; k < q; ++k)
        {
            std::string t = s + ((it++)->first);

            *FpWork++ = _drhs[t].eval( _data_table );

            //obsolete: *Zptmp++ = *Zp++;
        }
    }


    FzWork = FzSave; // _Fz; // FzSave;
    //obsolete: Zp = ZpSave;
    FpWork = FpSave; // _Fp; // FpSave;


    for (long k = 0; k < q; ++k)
    {
        double*   Fpk = FpWork++;
        double* Fztmp = FzWork;

        for (long j = 0; j < n; ++j)
        {
            double* Zptmp = Zp;

            *dZ = *Fpk;

            for (long l = 0; l < n; ++l)
            {
                *dZ += (*Fztmp++) * (*Zptmp++);
            }

            ++dZ;
            Fpk += q;
        }

        Zp += n;
    }


    delete[] FpSave;
    //obsolete: delete[] ZpSave;
    delete[] FzSave;


    return  DummyMat;
}
//----------------------------------------------------------------------------
