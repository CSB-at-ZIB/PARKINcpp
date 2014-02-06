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
#include <vector>
#include <map>

#include "common/Types.h"

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

            typedef std::vector<std::string>      Names;
            typedef std::vector<Real>               Coeff;
            typedef std::vector<long>              Multix;

            ///

            Expression();

            Expression(Real const& val);
            Expression(std::string const& name);
            Expression(char const* name);

            ///

            // Here, we make use of implicit type conversion: calls such as
            //
            //      Expression(PLUS, "a", 3.8)
            //      Expression(TIMES, Expression(LOG, "x"), "pi")
            //
            // are still possible and valid because of the tacitly perfomed type conversion.

            //
            Expression(ExprNodeType const&, Expression const&);

            //
            Expression(ExprNodeType const&, Expression const&, Expression const&);

            //
            Expression(ExprNodeType const&, Expression const&, Expression const&, Expression const&);

            ///

            Expression(Names const&, Coeff const&, std::vector<Multix> const&);

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
            void        off (Names const&, std::vector<long> const&);

        private:
            ExprNode* _node;
    };

    //

    std::ostream& operator<< (std::ostream& os, Expression const& expr);

}
#endif // __EXPRESSION_H
