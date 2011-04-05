// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-12-09 td
// last changed:
//
#include <limits>
#include <cmath>

#include "ExprNode.h"

using namespace PARKIN;

///

//---------------------------------------------------------------------------
Real
Literal::eval(Param& /*x*/) const
{
    return _val;
}
//---------------------------------------------------------------------------
Real
Literal::eval(double* /*y*/[]) const
{
    return _val;
}
//---------------------------------------------------------------------------
Expression
Literal::df(std::string const& /*n*/) const
{
    return Expression(0.0);
}
//---------------------------------------------------------------------------
void
Literal::prt(std::ostream& s) const
{
    s << _val;
}
//---------------------------------------------------------------------------
bool
Literal::eq(Real v) const
{
    return (_val == v);
}
//---------------------------------------------------------------------------
void
Literal::off(std::string const& /*n*/, long const /*a*/)
{
}
//---------------------------------------------------------------------------

///

//---------------------------------------------------------------------------
Real
Variable::eval(Param& x) const
{
    return x[_name];
}
//---------------------------------------------------------------------------
Real
Variable::eval(double* y[]) const
{
    const long zero = 0;
    const long one = 1;
    if ( _off < zero ) return *(y[1] - _off - one);
    return *(y[0] + _off);
}
//---------------------------------------------------------------------------
Expression
Variable::df(std::string const& n) const
{
    if (_name.compare(n) == 0) return Expression(1.0);
    return Expression(0.0);
}
//---------------------------------------------------------------------------
void
Variable::prt(std::ostream& s) const
{
    s << _name;
}
bool
Variable::eq(Real v) const
{
    return false;
}
//---------------------------------------------------------------------------
void
Variable::off(std::string const& n, long const a)
{
    if (_name.compare(n) == 0) _off = a;
}
//---------------------------------------------------------------------------

///

//---------------------------------------------------------------------------
//std::map<std::string, long> UnaryExpr::_opMap;
std::map<ExprNodeType, std::string> UnaryExpr::_opMap;
//---------------------------------------------------------------------------
Real
UnaryExpr::eval(Param& x) const
{
    double tmp;

    tmp = _expr.eval(x);

    return opEval(tmp);
}
//---------------------------------------------------------------------------
Real
UnaryExpr::eval(double* y[]) const
{
    Real tmp;

    tmp = _expr.eval( y );

    return opEval(tmp);
}
//---------------------------------------------------------------------------
Real
UnaryExpr::opEval(Real const& v) const
{
    switch(_op)
    {
        case MINUS   : return -1.0 * v;                     // _expr.eval(x);
        case ABS     : return std::fabs( v );               // _expr.eval(x) );
        case ARCCOS  : return std::acos( v );               // _expr.eval(x) );
        case ARCOSH  : return ::acosh( v );                 // _expr.eval(x) );
        case ARCSIN  : return std::asin( v );               // _expr.eval(x) );
        case ARSINH  : return ::asinh( v );                 // _expr.eval(x) );
        case ARCTAN  : return std::atan( v );               // _expr.eval(x) );
        case ARTANH  : return ::atanh( v );                 // _expr.eval(x) );
        case CEIL    : return std::ceil( v );               // _expr.eval(x) );
        case COS     : return std::cos( v );                // _expr.eval(x) );
        case COSH    : return std::cosh( v );               // _expr.eval(x) );
        case EXP     : return std::exp( v );                // _expr.eval(x) );
        // case DELTA   : return delta( v );                   // _expr.eval(x) );
        case FLOOR   : return std::floor( v );              // _expr.eval(x) );
        case LN      : return std::log( v );                // _expr.eval(x) );
        case LOG     : return std::log10( v );              // _expr.eval(x) );
        case SIGN    : return sign( v );                    // _expr.eval(x) );
        case SIN     : return std::sin( v );                // _expr.eval(x) );
        case SINH    : return std::sinh( v );               // _expr.eval(x) );
        case TAN     : return std::tan( v );                // _expr.eval(x) );
        case TANH    : return std::tanh( v );               // _expr.eval(x) );
        default      : return std::numeric_limits<Real>::signaling_NaN();
    }
}
//---------------------------------------------------------------------------
Expression
UnaryExpr::df(std::string const& n) const
{
    Expression e = _expr.df(n);

    if ( e.eq(0.0) ) return Expression(0.0);

    //switch (_opMap[_op])
    switch (_op)
    {
        case MINUS :    if ( e.eq(1.0) ) return Expression(-1.0);
                        return Expression( MINUS,  e );

        case ABS   :    if ( e.eq(1.0) ) return Expression( SIGN, _expr );
                        return Expression( TIMES,  Expression(SIGN, _expr), e );

        case EXP   :    if ( e.eq(1.0) ) return Expression( EXP, _expr );
                        return Expression( TIMES,  Expression(EXP, _expr), e );

        case LN    : return Expression( DIVIDE, e, _expr);

        case LOG   : return Expression(
                                DIVIDE,
                                Expression(DIVIDE, e, _expr),
                                std::log(10.0)
                            );

        // case SIGN  : return Expression( TIMES,  Expression(DELTA, _expr), e );

        case SIN   :    if ( e.eq(1.0) ) return Expression( COS, _expr );
                        return Expression( TIMES,  Expression(COS, _expr), e );

        case COS   :    if ( e.eq(1.0) ) return Expression( MINUS,  Expression(SIN, _expr) );
                        return Expression( MINUS,  Expression(TIMES, Expression(SIN, _expr), e) );

        case TAN   : return Expression(
                                DIVIDE,
                                e,
                                Expression(
                                    TIMES,
                                    Expression(COS, _expr),
                                    Expression(COS, _expr)
                                )
                            );

        default    : return Expression(0.0);
    }
}
//---------------------------------------------------------------------------
void
UnaryExpr::prt(std::ostream& s) const
{
    s << "(" << _opMap[_op] << " " << _expr << ")";
}
//---------------------------------------------------------------------------
bool
UnaryExpr::eq(Real v) const
{
    return false;
}
//---------------------------------------------------------------------------
void
UnaryExpr::off(std::string const& n, long const a)
{
    _expr.off(n,a);
}
//---------------------------------------------------------------------------

///

//---------------------------------------------------------------------------
Real
BinaryExpr::eval(Param& x) const
{
    switch(_op)
    {
        case '+' : return _lexpr.eval(x) + _rexpr.eval(x);
        case '-' : return _lexpr.eval(x) - _rexpr.eval(x);
        case '*' :  if ( _lexpr.eq(0.0) || _rexpr.eq(0.0) ) return 0.0;
                    return _lexpr.eval(x) * _rexpr.eval(x);
        case '/' : return _lexpr.eval(x) / _rexpr.eval(x);
        case '^' : return std::pow( _lexpr.eval(x), _rexpr.eval(x) );
        default  : return std::numeric_limits<Real>::signaling_NaN();
    }
}
//---------------------------------------------------------------------------
Real
BinaryExpr::eval(double* y[]) const
{
    switch(_op)
    {
        case '+' : return _lexpr.eval(y) + _rexpr.eval(y);
        case '-' : return _lexpr.eval(y) - _rexpr.eval(y);
        case '*' :  if ( _lexpr.eq(0.0) || _rexpr.eq(0.0) ) return 0.0;
                    return _lexpr.eval(y) * _rexpr.eval(y);
        case '/' : return _lexpr.eval(y) / _rexpr.eval(y);
        case '^' : return std::pow( _lexpr.eval(y), _rexpr.eval(y) );
        default  : return std::numeric_limits<Real>::signaling_NaN();
    }
}
//---------------------------------------------------------------------------
Expression
BinaryExpr::df(std::string const& n) const
{
    Expression le = _lexpr.df(n);
    Expression re = _rexpr.df(n);

    switch (_op)
    {
        /// (f + g)' = f' + g'
        case '+' :  if ( le.eq(0.0) )
                    {
                        if ( re.eq(0.0) ) return Expression(0.0);
                        return re;
                    }
                    if ( re.eq(0.0) )
                    {
                        if ( le.eq(0.0) ) return Expression(0.0);
                        return le;
                    }
                    return Expression( PLUS, le, re );

        /// (f - g)' = f' - g'
        case '-' :  if ( le.eq(0.0) )
                    {
                        if ( re.eq(0.0) ) return Expression(0.0);
                        return Expression(MINUS, re);
                    }
                    if ( re.eq(0.0) )
                    {
                        if ( le.eq(0.0) ) return Expression(0.0);
                        return le;
                    }
                    return Expression( MINUS, le, re );

        /// (f*g)' = f'*g + f*g'
        case '*' :  if ( le.eq(0.0) )
                    {
                        if ( re.eq(0.0) ) return Expression(0.0);
                        if ( re.eq(1.0) ) return _lexpr;
                        return Expression( TIMES, _lexpr, re);
                    }
                    if ( le.eq(1.0) )
                    {
                        if ( re.eq(0.0) ) return _rexpr;
                        if ( re.eq(1.0) ) return Expression( PLUS, _rexpr, _lexpr );
                        return Expression( PLUS, _rexpr, Expression(TIMES, _lexpr, re) );
                    }
                    if ( re.eq(0.0) )
                    {
                        if ( le.eq(0.0) ) return Expression(0.0);
                        if ( le.eq(1.0) ) return _rexpr;
                        return Expression( TIMES, le, _rexpr );
                    }
                    if ( re.eq(1.0) )
                    {
                        if ( le.eq(0.0) ) return _lexpr;
                        if ( le.eq(1.0) ) return Expression( PLUS, _rexpr, _lexpr );
                        return Expression( PLUS, Expression(TIMES, le, _rexpr), _lexpr );
                    }
                    return Expression(
                                PLUS,
                                Expression( TIMES, le,     _rexpr ),
                                Expression( TIMES, _lexpr, re     )
                           );

        /// (f/g)' = (f'*g - f*g')/(g^2)
        case '/' :  if ( le.eq(0.0) )
                    {
                        if ( re.eq(0.0) ) return Expression(0.0);
                        if ( re.eq(1.0) ) return Expression(
                                                        DIVIDE,
                                                        Expression(MINUS, _lexpr),
                                                        Expression(TIMES, _rexpr, _rexpr)
                                                 );
                        return Expression(
                                    DIVIDE,
                                    Expression( MINUS, Expression( TIMES, _lexpr, re ) ),
                                    Expression( TIMES, _rexpr, _rexpr )
                               );
                    }
                    if ( le.eq(1.0) )
                    {
                        if ( re.eq(0.0) ) return Expression( DIVIDE, 1.0, _rexpr );
                        if ( re.eq(1.0) ) return Expression(
                                                    DIVIDE,
                                                    Expression(MINUS, _rexpr, _lexpr),
                                                    Expression(TIMES, _rexpr, _rexpr)
                                                 );
                        return Expression(
                                    DIVIDE,
                                    Expression(
                                        MINUS,
                                        _rexpr,
                                        Expression( TIMES, _lexpr, re )
                                    ),
                                    Expression( TIMES, _rexpr, _rexpr )
                               );
                    }
                    if ( re.eq(0.0) )
                    {
                        if ( le.eq(0.0) ) return Expression(0.0);
                        if ( le.eq(1.0) ) return Expression( DIVIDE, 1.0, _rexpr );
                        return Expression( DIVIDE, le, _rexpr );
                    }
                    if ( re.eq(1.0) )
                    {
                        if ( le.eq(0.0) ) return Expression(
                                                    DIVIDE,
                                                    Expression(MINUS, _lexpr),
                                                    Expression(TIMES, _rexpr, _rexpr)
                                                 );
                        if ( le.eq(1.0) ) return Expression(
                                                    DIVIDE,
                                                    Expression(MINUS, _rexpr, _lexpr),
                                                    Expression(TIMES, _rexpr, _rexpr)
                                                 );
                        return Expression(
                                    DIVIDE,
                                    Expression(
                                        MINUS,
                                        Expression( TIMES, le, _rexpr ),
                                        _lexpr
                                    ),
                                    Expression( TIMES, _rexpr, _rexpr )
                               );
                    }
                    return Expression(
                                DIVIDE,
                                Expression(
                                    MINUS,
                                    Expression( TIMES, le,     _rexpr ),
                                    Expression( TIMES, _lexpr, re     )
                                ),
                                Expression( TIMES, _rexpr, _rexpr )
                           );

        /// (f^g)' = [exp(g*log f)]' = (f^g) * g' * (log f) + f^(g - 1) * g * f'
        case '^' :  if ( le.eq(0.0) )
                    {
                        if ( re.eq(0.0) ) return Expression(0.0);
                        if ( re.eq(1.0) ) return Expression(
                                                        TIMES,
                                                        Expression(POWER, _lexpr, _rexpr),
                                                        Expression(LN, _lexpr)
                                                 );
                        return Expression(
                                    TIMES,
                                    re,
                                    Expression(
                                        TIMES,
                                        Expression(POWER, _lexpr, _rexpr),
                                        Expression(LN, _lexpr)
                                    )
                               );
                    }
                    if ( le.eq(1.0) )
                    {
                        if ( re.eq(0.0) ) return Expression(
                                                    TIMES,
                                                    _rexpr,
                                                    Expression( POWER, _lexpr, Expression(MINUS, _rexpr, 1.0) )
                                                 );
                        if ( re.eq(1.0) ) return Expression(
                                                    PLUS,
                                                    Expression(
                                                        TIMES,
                                                        _rexpr,
                                                        Expression( POWER, _lexpr, Expression(MINUS, _rexpr, 1.0) )
                                                    ),
                                                    Expression(
                                                        TIMES,
                                                        Expression(POWER, _lexpr, _rexpr),
                                                        Expression(LN, _lexpr)
                                                    )
                                                 );
                        return Expression(
                                    PLUS,
                                    Expression(
                                        TIMES,
                                        _rexpr,
                                        Expression( POWER, _lexpr, Expression(MINUS, _rexpr, 1.0) )
                                    ),
                                    Expression(
                                        TIMES,
                                        re,
                                        Expression(
                                            TIMES,
                                            Expression(POWER, _lexpr, _rexpr),
                                            Expression(LN, _lexpr)
                                        )
                                    )
                               );
                    }
                    if ( re.eq(0.0) )
                    {
                        if ( le.eq(0.0) ) return Expression(0.0);
                        if ( le.eq(1.0) ) return Expression(
                                                        TIMES,
                                                        _rexpr,
                                                        Expression(
                                                            POWER,
                                                            _lexpr,
                                                            Expression(MINUS, _rexpr, 1.0)
                                                        )
                                                 );
                        return Expression(
                                    TIMES,
                                    le,
                                    Expression(
                                        TIMES,
                                        _rexpr,
                                        Expression( POWER, _lexpr, Expression(MINUS, _rexpr, 1.0) )
                                    )
                               );
                    }
                    if ( re.eq(1.0) )
                    {
                        if ( le.eq(0.0) ) return Expression(
                                                    TIMES,
                                                    Expression(POWER, _lexpr, _rexpr),
                                                    Expression(LN, _lexpr)
                                                 );
                        if ( le.eq(1.0) ) return Expression(
                                                    PLUS,
                                                    Expression(
                                                        TIMES,
                                                        _rexpr,
                                                        Expression( POWER, _lexpr, Expression(MINUS, _rexpr, 1.0) )
                                                    ),
                                                    Expression(
                                                        TIMES,
                                                        Expression(POWER, _lexpr, _rexpr),
                                                        Expression(LN, _lexpr)
                                                    )
                                                 );
                        return Expression(
                                    PLUS,
                                    Expression(
                                        TIMES,
                                        le,
                                        Expression(
                                            TIMES,
                                            _rexpr,
                                            Expression( POWER, _lexpr, Expression(MINUS, _rexpr, 1.0) )
                                        )
                                    ),
                                    Expression(
                                        TIMES,
                                        Expression(POWER, _lexpr, _rexpr),
                                        Expression(LN, _lexpr)
                                    )
                               );
                    }
                    return Expression(
                                PLUS,
                                Expression(
                                    TIMES,
                                    le,
                                    Expression(
                                        TIMES,
                                        _rexpr,
                                        Expression( POWER, _lexpr, Expression(MINUS, _rexpr, 1.0) )
                                    )
                                ),
                                Expression(
                                    TIMES,
                                    re,
                                    Expression(
                                        TIMES,
                                        Expression(POWER, _lexpr, _rexpr),
                                        Expression(LN, _lexpr)
                                    )
                                )
                           );

        default  : return Expression(0.0);
    }
}
//---------------------------------------------------------------------------
void
BinaryExpr::prt(std::ostream& s) const
{
    s << "(" << _lexpr << " " << char(_op) << " " << _rexpr << ")";
}
//---------------------------------------------------------------------------
bool
BinaryExpr::eq(Real v) const
{
    return false;
}
//---------------------------------------------------------------------------
void
BinaryExpr::off(std::string const& n, long const a)
{
    _lexpr.off(n,a);
    _rexpr.off(n,a);
}
//---------------------------------------------------------------------------

///

//---------------------------------------------------------------------------
std::map<ExprNodeType, std::string> TernaryExpr::_opMap;
//---------------------------------------------------------------------------
Real
TernaryExpr::eval(Param& x) const
{
    Real tmp;

    //switch(_opMap[_op])
    switch(_op)
    {
        case PLUS      : tmp = _lexpr.eval(x) + _mexpr.eval(x) + _rexpr.eval(x);
                         return tmp;

        case HILLplus  : tmp = std::pow( (_lexpr.eval(x)/_mexpr.eval(x)), _rexpr.eval(x) );
                         return tmp / ( 1.0 + tmp );

        case HILLminus : tmp = std::pow( (_lexpr.eval(x)/_mexpr.eval(x)), _rexpr.eval(x) );
                         return 1.0 / ( 1.0 + tmp );

        default        : return std::numeric_limits<Real>::signaling_NaN();
    }
}
//---------------------------------------------------------------------------
Real
TernaryExpr::eval(double* y[]) const
{
    Real tmp;

    //switch(_opMap[_op])
    switch(_op)
    {
        case PLUS      : tmp = _lexpr.eval(y) + _mexpr.eval(y) + _rexpr.eval(y);
                         return tmp;

        case HILLplus  : tmp = std::pow( (_lexpr.eval(y)/_mexpr.eval(y)), _rexpr.eval(y) );
                         return tmp / ( 1.0 + tmp );

        case HILLminus : tmp = std::pow( (_lexpr.eval(y)/_mexpr.eval(y)), _rexpr.eval(y) );
                         return 1.0 / ( 1.0 + tmp );

        default        : return std::numeric_limits<Real>::signaling_NaN();
    }
}
//---------------------------------------------------------------------------
Expression
TernaryExpr::df(std::string const& n) const
{
    //switch (_opMap[_op])
    switch (_op)
    {
        case PLUS      : return Expression(PLUS, _lexpr.df(n), _mexpr.df(n), _rexpr.df(n) );

        case HILLplus  : {
                             Expression Hplus  = Expression(HILLplus, _lexpr, _mexpr, _rexpr);
                             Expression tmpp   = Expression(TIMES, Expression(MINUS, 1.0, Hplus), Hplus);
                             Expression term1p = Expression(TIMES,                   Expression(DIVIDE, _rexpr,_lexpr) , _lexpr.df(n));
                             Expression term2p = Expression(TIMES, Expression(MINUS, Expression(DIVIDE, _rexpr,_mexpr)), _mexpr.df(n));
                             Expression term3p = Expression(TIMES, Expression(LN,    Expression(DIVIDE, _lexpr,_mexpr)), _rexpr.df(n));
                             Expression sump   = Expression(PLUS, term1p, term2p, term3p);
                             return Expression(TIMES, tmpp, sump);
                         }

        case HILLminus : {
                             Expression Hp   = Expression(HILLplus,  _lexpr, _mexpr, _rexpr);
                             Expression Hm   = Expression(HILLminus, _lexpr, _mexpr, _rexpr);
                             Expression tmpm = Expression(TIMES, Expression(MINUS, Hm), Hp);
                             Expression term1m = Expression(TIMES,                   Expression(DIVIDE, _rexpr,_lexpr) , _lexpr.df(n));
                             Expression term2m = Expression(TIMES, Expression(MINUS, Expression(DIVIDE, _rexpr,_mexpr)), _mexpr.df(n));
                             Expression term3m = Expression(TIMES, Expression(LN,    Expression(DIVIDE, _lexpr,_mexpr)), _rexpr.df(n));
                             Expression summ = Expression(PLUS, term1m, term2m, term3m);
                             return Expression(TIMES, tmpm, summ);
                         }

        default        : return Expression(0);
    }
}
//---------------------------------------------------------------------------
void
TernaryExpr::prt(std::ostream& s) const
{
    s << _opMap[_op] << "(";
    s << _lexpr << ", " << _mexpr << ", " << _rexpr <<")";
}
//---------------------------------------------------------------------------
bool
TernaryExpr::eq(Real v) const
{
    return false;
}
//---------------------------------------------------------------------------
void
TernaryExpr::off(std::string const& n, long const a)
{
    _lexpr.off(n,a);
    _mexpr.off(n,a);
    _rexpr.off(n,a);
}
//---------------------------------------------------------------------------

///

//---------------------------------------------------------------------------
Real
UserExpr::eval(Param& x) const
{
    return _expr.eval(x);
}
//---------------------------------------------------------------------------
Real
UserExpr::eval(double* y[]) const
{
    return _expr.eval(y);
}
//---------------------------------------------------------------------------
Expression
UserExpr::df(std::string const& n) const
{
    return _expr.df(n);
}
//---------------------------------------------------------------------------
void
UserExpr::prt(std::ostream& s) const
{
    s << "(" << _name << " = " << _expr << ")";
}
//---------------------------------------------------------------------------
bool
UserExpr::eq(Real v) const
{
    return false;
}
//---------------------------------------------------------------------------
void
UserExpr::off(std::string const& n, long const a)
{
    _expr.off(n,a);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Real
PARKIN::sign(Real x)
{
    if ( x > 0.0 ) return 1.0;
    if ( x < 0.0 ) return -1.0;
    return 0.0;
}
//---------------------------------------------------------------------------
