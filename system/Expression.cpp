// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-11-29 td
// last changed:
//
#include "ExprNode.h"
// #include "Expression.h"

using namespace PARKIN;

//----------------------------------------------------------------------------
Expression::Expression() :
    _node( new Literal(0.0) )
{
}
//----------------------------------------------------------------------------
Expression::Expression(Real const& val) :
    _node( new Literal(val) )
{
}
Expression::Expression(std::string const& name) :
    _node( new Variable(name) )
{
}
Expression::Expression(char const* name) :
    _node( new Variable(std::string(name)) )
{
}
//----------------------------------------------------------------------------
Expression::Expression(ExprNodeType const& op,
                       Expression const& e) :
    _node( new UnaryExpr(op, e) )
{
}
//----------------------------------------------------------------------------
Expression::Expression(ExprNodeType const& op,
                       Expression const& e1,
                       Expression const& e2) :
    _node( new BinaryExpr(op, e1, e2) )
{
}
//----------------------------------------------------------------------------
Expression::Expression(ExprNodeType const& op,
                       Expression const& e1,
                       Expression const& e2,
                       Expression const& e3) :
    _node( new TernaryExpr(op, e1, e2, e3) )
{
}
//----------------------------------------------------------------------------
Expression::~Expression()
{
    --(_node->_count);
    if ( _node->_count <= 0) delete _node;
}
//----------------------------------------------------------------------------
Expression::Expression(Expression const& e) :
    _node(e._node)
{
    ++(_node->_count);
}
void
Expression::operator= (Expression const& e)
{
    ++(e._node->_count);
    --(_node->_count);
    if (_node->_count == 0) delete _node;
    _node = e._node;
}
//----------------------------------------------------------------------------
Real
Expression::eval(Param& x) const
{
    return _node -> eval(x);
}
Real
Expression::eval(double* y[]) const
{
    return _node -> eval(y);
}
Expression
Expression::df(std::string const& n) const
{
    return _node -> df(n);
}
void
Expression::prt(std::ostream& os) const
{
    _node -> prt(os);
}
bool
Expression::eq(Real v) const
{
    return _node -> eq(v);
}
void
Expression::off(std::string const& n, long const a)
{
    _node -> off(n, a);
}
//----------------------------------------------------------------------------

///
///

//----------------------------------------------------------------------------
std::ostream& PARKIN::operator<< (std::ostream& os, Expression const& expr)
{
    expr.prt(os);
    return os;
}
//----------------------------------------------------------------------------
