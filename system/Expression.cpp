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
//----------------------------------------------------------------------------
Expression::Expression(ExprNodeType const& op,
                       Real const& val) :
    _node( new UnaryExpr(op, Expression(val)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       std::string const& name) :
    _node( new UnaryExpr(op, Expression(name)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Expression const& e) :
    _node( new UnaryExpr(op, e) )
{
}
//----------------------------------------------------------------------------
Expression::Expression(ExprNodeType const& op,
                       Real const& v1,
                       Real const& v2) :
    _node( new BinaryExpr(op, Expression(v1), Expression(v2)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Real const& v,
                       std::string const& n) :
    _node( new BinaryExpr(op, Expression(v), Expression(n)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       std::string const& n,
                       Real const& v) :
    _node( new BinaryExpr(op, Expression(n), Expression(v)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Real const& v,
                       Expression const& e) :
    _node( new BinaryExpr(op, Expression(v), e) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Expression const& e,
                       Real const& v) :
    _node( new BinaryExpr(op, e, Expression(v)) )
{
}
//----------------------------------------------------------------------------
Expression::Expression(ExprNodeType const& op,
                       std::string const& n1,
                       std::string const& n2) :
    _node( new BinaryExpr(op, Expression(n1), Expression(n2)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       std::string const& n,
                       Expression const& e) :
    _node( new BinaryExpr(op, Expression(n), e) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Expression const& e,
                       std::string const& n) :
    _node( new BinaryExpr(op, e, Expression(n)) )
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
                       Real const& v1,
                       Real const& v2,
                       Real const& v3) :
    _node( new TernaryExpr(op, Expression(v1), Expression(v2), Expression(v3)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       std::string const& n,
                       Real const& v1,
                       Real const& v2) :
    _node( new TernaryExpr(op, Expression(n), Expression(v1), Expression(v2)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Real const& v1,
                       std::string const& n,
                       Real const& v2) :
    _node( new TernaryExpr(op, Expression(v1), Expression(n), Expression(v2)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Real const& v1,
                       Real const& v2,
                       std::string const& n) :
    _node( new TernaryExpr(op, Expression(v1), Expression(v2), Expression(n)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       std::string const& n1,
                       std::string const& n2,
                       Real const& v) :
    _node( new TernaryExpr(op, Expression(n1), Expression(n2), Expression(v)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       std::string const& n1,
                       Real const& v,
                       std::string const& n2) :
    _node( new TernaryExpr(op, Expression(n1), Expression(v), Expression(n2)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Real const& v,
                       std::string const& n1,
                       std::string const& n2) :
    _node( new TernaryExpr(op, Expression(v), Expression(n1), Expression(n2)) )
{
}
//----------------------------------------------------------------------------
Expression::Expression(ExprNodeType const& op,
                       Expression const& e,
                       Real const& v1,
                       Real const& v2) :
    _node( new TernaryExpr(op, e, Expression(v1), Expression(v2)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Real const& v1,
                       Expression const& e,
                       Real const& v2) :
    _node( new TernaryExpr(op, Expression(v1), e, Expression(v2)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Real const& v1,
                       Real const& v2,
                       Expression const& e) :
    _node( new TernaryExpr(op, Expression(v1), Expression(v2), e) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Expression const& e1,
                       Expression const& e2,
                       Real const& v) :
    _node( new TernaryExpr(op, e1, e2, Expression(v)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Expression const& e1,
                       Real const& v,
                       Expression const& e2) :
    _node( new TernaryExpr(op, e1, Expression(v), e2) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Real const& v,
                       Expression const& e1,
                       Expression const& e2) :
    _node( new TernaryExpr(op, Expression(v), e1, e2) )
{
}
//----------------------------------------------------------------------------
Expression::Expression(ExprNodeType const& op,
                       std::string const& n1,
                       std::string const& n2,
                       std::string const& n3) :
    _node( new TernaryExpr(op, Expression(n1), Expression(n2), Expression(n3)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Expression const& e,
                       std::string const& n1,
                       std::string const& n2) :
    _node( new TernaryExpr(op, e, Expression(n1), Expression(n2)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       std::string const& n1,
                       Expression const& e,
                       std::string const& n2) :
    _node( new TernaryExpr(op, Expression(n1), e, Expression(n2)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       std::string const& n1,
                       std::string const& n2,
                       Expression const& e) :
    _node( new TernaryExpr(op, Expression(n1), Expression(n2), e) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Expression const& e1,
                       Expression const& e2,
                       std::string const& n) :
    _node( new TernaryExpr(op, e1, e2, Expression(n)) )
{
}
Expression::Expression(ExprNodeType const& op,
                       Expression const& e1,
                       std::string const& n,
                       Expression const& e2) :
    _node( new TernaryExpr(op, e1, Expression(n), e2) )
{
}
Expression::Expression(ExprNodeType const& op,
                       std::string const& n,
                       Expression const& e1,
                       Expression const& e2) :
    _node( new TernaryExpr(op, Expression(n), e1, e2) )
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
