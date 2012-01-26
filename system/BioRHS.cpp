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
    _rhs(), _drhs(), _species(), _parameter(),
    _num_table_entries(2), _data_table(new double*[2])
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
    _rhs(emap), _drhs(), _species(), _parameter(),
    _num_table_entries(2), _data_table(new double*[2])
{
    EMapIterConst eBeg = emap.begin();
    EMapIterConst eEnd = emap.end();

    _species.clear();

    for (EMapIterConst it = eBeg; it != eEnd; ++it)
    {
        _species.push_back( it->first );
    }

    StrIterConst sBeg = _species.begin();
    StrIterConst sEnd = _species.end();
    for (EMapIterConst it = eBeg; it != eEnd; ++it)
    {
        long       a = 0;
        Expression expr = it->second;

        for (StrIterConst itSpe = sBeg; itSpe != sEnd; ++itSpe)
            expr.off(*itSpe, a++);
    }

    _parameter.clear();

    _data_table[0] = _data_table[1] = 0;

    computeDerivativeExpression();
}
//----------------------------------------------------------------------------
BioRHS::BioRHS(Species const& species) :
    _rhs(), _drhs(), _species(), _parameter(),
    _num_table_entries(2), _data_table(new double*[2])
{
    StrIterConst sBeg = species.begin();
    StrIterConst sEnd = species.end();

    _rhs.clear();

    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
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
        long       a = 0;
        Expression expr = it->second;

        for (StrIterConst itSpe = spBeg; itSpe != spEnd; ++itSpe)
            expr.off(*itSpe, a++);
    }

    _parameter.clear();

    _data_table[0] = _data_table[1] = 0;

    computeDerivativeExpression();
}
//----------------------------------------------------------------------------
BioRHS::BioRHS(ExpressionMap const& emap,
               Parameter const& param) :
    _rhs(emap), _drhs(), _species(), _parameter(param),
    _num_table_entries(2), _data_table(new double*[2])
{
    EMapIterConst eBeg = emap.begin();
    EMapIterConst eEnd = emap.end();

    _species.clear();

    for (EMapIterConst it = eBeg; it != eEnd; ++it)
    {
        _species.push_back( it->first );
    }

    StrIterConst sBeg = _species.begin();
    StrIterConst sEnd = _species.end();
    StrIterConst vBeg = _parameter.begin();
    StrIterConst vEnd = _parameter.end();
    for (EMapIterConst it = eBeg; it != eEnd; ++it)
    {
        long       a = 0;
        long       b = 0;
        Expression expr = it->second;

        for (StrIterConst itSpe = sBeg; itSpe != sEnd; ++itSpe)
            expr.off(*itSpe, a++);
        for (StrIterConst itPar = vBeg; itPar != vEnd; ++itPar)
            expr.off(*itPar, --b);
    }

    _data_table[0] = _data_table[1] = 0;

    computeDerivativeExpression();
}
//----------------------------------------------------------------------------
BioRHS::~BioRHS()
{
    // for (long j = 1; j < _num_table_entries; ++j) delete[] _data_table[j];
    delete[] _data_table;
}
//----------------------------------------------------------------------------
BioRHS::BioRHS(BioRHS const& r)
{
    if ( this != &r )
    {
        _rhs       = r._rhs;
        _drhs      = r._drhs;
        _species   = r._species;
        _parameter = r._parameter;

        _num_table_entries = r._num_table_entries;
        _data_table = new double*[_num_table_entries];
        for (long j = 0; j < _num_table_entries; ++j)
            _data_table[j] = r._data_table[j];

        // computeDerivativeExpression();
    }
}
//----------------------------------------------------------------------------
BioRHS const&
BioRHS::operator= (BioRHS const& r)
{
    if ( this != &r )
    {
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

        // computeDerivativeExpression();
    }

    return *this;
}
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
            long        a = 0;
            long        b = 0;
            std::string s = (itRHS->first) + " / " + *itSpe;

            _drhs[s] = expr.df( *itSpe );

            for (StrIterConst it = sBeg; it != sEnd; ++it)
                _drhs[s].off(*it,a++);
            for (StrIterConst it = pBeg; it != pEnd; ++it)
                _drhs[s].off(*it,--b);
        }

        for (StrIterConst itPar = pBeg; itPar != pEnd; ++itPar)
        // for (Expression::ParamIterConst itPar = vBeg; itPar != vEnd; ++itPar)
        {
            long        a = 0;
            long        b = 0;
            std::string s = (itRHS->first) + " / " + *itPar; // itPar->first;

            _drhs[s] = expr.df( *itPar );

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
BioRHS::setRHS(ExpressionMap const& emap)
{
    _rhs = emap;

    _species.clear();

    EMapIterConst rBeg = _rhs.begin();
    EMapIterConst rEnd = _rhs.end();

    for (EMapIterConst it = rBeg; it != rEnd; ++it)
    {
        _species.push_back( it->first );
    }

    StrIterConst sBeg = _species.begin();
    StrIterConst sEnd = _species.end();
    StrIterConst pBeg = _parameter.begin();
    StrIterConst pEnd = _parameter.end();

    for (EMapIterConst it = rBeg; it != rEnd; ++it)
    {
        long       a = 0;
        long       b = 0;
        Expression expr = it->second;

        for (StrIterConst itSpe = sBeg; itSpe != sEnd; ++itSpe)
            expr.off(*itSpe, a++);
        for (StrIterConst itPar = pBeg; itPar != pEnd; ++itPar)
            expr.off(*itPar, --b);
    }

    // _parameter.clear();

    computeDerivativeExpression();
}
//----------------------------------------------------------------------------
void
BioRHS::resetSpecies(BioRHS::Species const& species)
{
    StrIterConst sBeg = species.begin();
    StrIterConst sEnd = species.end();

    _rhs.clear();

    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
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
        long       a = 0;
        Expression expr = it->second;

        for (StrIterConst itSpe = spBeg; itSpe != spEnd; ++itSpe)
            expr.off(*itSpe, a++);
    }

    // _parameter.clear();

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
Vector
BioRHS::f(Expression::Param& par,
          long n, double* dy)
{
    StrIterConst   sBeg = _species.begin();
    StrIterConst   sEnd = _species.end();

    if ( n == 0 )
    {
        long   j = 0;
        Vector dz( _species.size() );

        for (StrIterConst it = sBeg; it != sEnd; ++it)
        {
            dz(++j) = _rhs[*it].eval(par);
        }

        return dz;
    }

    ///

    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        *(dy++) = _rhs[*it].eval(par);
    }

    return Vector();
}
//----------------------------------------------------------------------------
Vector
BioRHS::f(double* y,
          long n, double* dy)
{
    StrIterConst   sBeg = _species.begin();
    StrIterConst   sEnd = _species.end();

    _data_table[0] = y;

    if ( n == 0 )
    {
        long   j = 0;
        Vector dz( _species.size() );

        for (StrIterConst it = sBeg; it != sEnd; ++it)
        {
            dz(++j) = _rhs[*it].eval( _data_table );
        }

        return dz;
    }

    ///

    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        *dy++ = _rhs[*it].eval( _data_table );
    }

    return Vector();
}
//----------------------------------------------------------------------------
Matrix
BioRHS::Jf(Expression::Param& par,
           long n, double* J, long ldJ)
{
    if ( n == 0 )
    {
        // long
        n = _species.size();
        Matrix  Fz( n, n );

        for (long j = 1; j <= n; ++j)
        {
            for (long k = 1; k <= n; ++k)
            {
                std::string s = _species[j-1] + " / " + _species[k-1];

                Fz(j,k) = _drhs[s].eval(par);
            }
        }

        return Fz;
    }

    ///

    for (long j = 0; j < n; ++j)
    {
        for (long k = 0; k < n; ++k)
        {
            std::string s = _species[j] + " / " + _species[k];

            J[ ldJ*k + j ] = _drhs[s].eval(par);
        }
    }

    return Matrix();
}
//----------------------------------------------------------------------------
Matrix
BioRHS::df(Expression::Param const& var, Matrix const& Zp,
           Expression::Param& par)
{
    long                n = Zp.nr();
    long                q = Zp.nc();
    Matrix              Fz( n, n );
    Matrix              Fp( n, q );

    if ( (n > (long)_species.size()) ||
         (q > (long)_parameter.size())
       )
    {
       return Matrix();
    }

    for (long j = 1; j <= n; ++j)
    {
        std::string s = _species[j-1] + " / ";

        for (long k = 1; k <= n; ++k)
        {
            std::string t = s + _species[k-1];

            Fz(j,k) = _drhs[t].eval(par);
        }

        Expression::ParamIterConst it = var.begin();

        for (long k = 1; k <= q; ++k)
        {
            // std::string t = s + _parameter[k-1];
            std::string t = s + ((it++)->first);

            Fp(j,k) = _drhs[t].eval(par);
        }
    }

    // Matrix dZ = Fz * Zp + Fp;

    return  Fz * Zp  +  Fp;  //  dZ;
}
//----------------------------------------------------------------------------
Matrix
BioRHS::df(Expression::Param const& var, double* Zp, double* y,
           long n, long q, double* dZ)
{
    // Matrix              Fz( n, n );
    // Matrix              Fp( n, q );

    if ( (n > (long) _species.size()) ||
         (q > (long) var.size())
       )
    {
       return Matrix();
    }


    double* FzSave = new double[n*n];
    double* ZpSave = new double[n*q];
    double* FpSave = new double[n*q];

    double* Fz    = FzSave;
    double* Zptmp = ZpSave;
    double* Fp    = FpSave;


    _data_table[0] = y;


    for (long j = 0; j < n; ++j)
    {
        std::string s = _species[j] + " / ";

        for (long k = 0; k < n; ++k)
        {
            std::string t = s + _species[k];

            *Fz++ = _drhs[t].eval( _data_table );
        }


        Expression::ParamIterConst it = var.begin();

        for (long k = 0; k < q; ++k)
        {
            std::string t = s + ((it++)->first);

            *Fp++ = _drhs[t].eval( _data_table );

            *Zptmp++ = *Zp++;
        }
    }

    Fz = FzSave;
    Zp = ZpSave;
    Fp = FpSave;

    for (long j = 0; j < n; ++j)
    {
        for (long k = 0; k < q; ++k)
        {
            double* Fztmp = Fz;
            double* Zptmp = Zp + k;

            *dZ = *Fp++;

            for (long l = 0; l < n; ++l)
            {
                *dZ += (*Fztmp++) * (*Zptmp);
                Zptmp += q;
            }

            dZ++;
        }

        Fz += n;
    }


    delete[] FpSave;
    delete[] ZpSave;
    delete[] FzSave;


    return  Matrix();
}
//----------------------------------------------------------------------------
