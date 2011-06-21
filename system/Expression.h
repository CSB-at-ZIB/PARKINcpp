// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-11-23 td
// last changed:
//
#ifndef __EXPRESSION_H
#define __EXPRESSION_H

#include <iostream>
#include <string>
#include <map>


namespace PARKIN
{
    ///
    enum ExprNodeType
    {
          PLUS  = int('+'), MINUS  = int('-'),
          TIMES = int('*'), DIVIDE = int('/'), POWER = int('^')
        , ABS  = 256
        , ARCCOS,   ARCOSH
        // , ARCCOT,   ARCOTH
        // , ARCCSC,   ARCSCH
        // , ARCSEC,   ARSECH
        , ARCSIN,   ARSINH
        , ARCTAN,   ARTANH
        , CEIL
        , COS,      COSH
        // , COT,      COTH
        // , CSC,      CSCH
        , DELTA
        , EXP
        , FLOOR
        , LN,       LOG
        // , SEC,      SECH
        , SIGN
        , SIN,      SINH
        , TAN,      TANH
        , HILLplus
        , HILLminus
    };
    ///

    class ExprNode;         // needed here for declaration of private "_node" pointer as root of Expression

    class Expression
    {
        public:
            typedef std::map <std::string, Real>  Param;            // for declaration of method eval(Param&) below
            typedef Param::const_iterator         ParamIterConst;

            ///

            Expression();

            Expression(Real const& val);
            Expression(std::string const& name);

            //
            Expression(ExprNodeType const&, Real const&);
            Expression(ExprNodeType const&, std::string const&);
            Expression(ExprNodeType const&, Expression const&);

            //
            Expression(ExprNodeType const&, Real const&, Real const&);
            Expression(ExprNodeType const&, Real const&, std::string const&);
            Expression(ExprNodeType const&, std::string const&, Real const&);
            Expression(ExprNodeType const&, Real const&, Expression const&);
            Expression(ExprNodeType const&, Expression const&, Real const&);

            Expression(ExprNodeType const&, std::string const&, std::string const&);
            Expression(ExprNodeType const&, std::string const&, Expression const&);
            Expression(ExprNodeType const&, Expression const&, std::string const&);

            Expression(ExprNodeType const&, Expression const&, Expression const&);

            //
            Expression(ExprNodeType const&, Real const&, Real const&, Real const&);
            Expression(ExprNodeType const&, std::string const&, Real const&, Real const&);
            Expression(ExprNodeType const&, Real const&, std::string const&, Real const&);
            Expression(ExprNodeType const&, Real const&, Real const&, std::string const&);
            Expression(ExprNodeType const&, std::string const&, std::string const&, Real const&);
            Expression(ExprNodeType const&, std::string const&, Real const&, std::string const&);
            Expression(ExprNodeType const&, Real const&, std::string const&, std::string const&);

            Expression(ExprNodeType const&, Expression const&, Real const&, Real const&);
            Expression(ExprNodeType const&, Real const&, Expression const&, Real const&);
            Expression(ExprNodeType const&, Real const&, Real const&, Expression const&);
            Expression(ExprNodeType const&, Expression const&, Expression const&, Real const&);
            Expression(ExprNodeType const&, Expression const&, Real const&, Expression const&);
            Expression(ExprNodeType const&, Real const&, Expression const&, Expression const&);

            Expression(ExprNodeType const&, std::string const&, std::string const&, std::string const&);
            Expression(ExprNodeType const&, Expression const&, std::string const&, std::string const&);
            Expression(ExprNodeType const&, std::string const&, Expression const&, std::string const&);
            Expression(ExprNodeType const&, std::string const&, std::string const&, Expression const&);
            Expression(ExprNodeType const&, Expression const&, Expression const&, std::string const&);
            Expression(ExprNodeType const&, Expression const&, std::string const&, Expression const&);
            Expression(ExprNodeType const&, std::string const&, Expression const&, Expression const&);

            Expression(ExprNodeType const&, Expression const&, Expression const&, Expression const&);

            ///

            ~Expression();

            Expression(Expression const& e);
            void operator= (Expression const& e);

            Real        eval(Param&)             const;
            Real        eval(double*[])          const;
            Expression  df  (std::string const&) const;
            void        prt (std::ostream&)      const;
            bool        eq  (Real)               const;
            void        off (std::string const&, long const);

        private:
            ExprNode* _node;
    };

    //

    std::ostream& operator<< (std::ostream& os, Expression const& expr);

}
#endif // __EXPRESSION_H
