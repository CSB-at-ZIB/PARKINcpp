// Copyright (C) 2010 - 2013
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2013-09-25 td
// last changed:
//
#ifndef __PARSER_H
#define __PARSER_H

#include <vector>
#include <string>


namespace PARKIN
{
    // Forward declarations
    class Expression;
    // class Node;

    // Token class
    //--------------------------------------------------------------------------
    class Token
    {
        public:

            // Token type
            enum TokenType
            {
                TypeUnknown = -1,
                TypeEndOfLine = 0,

                // Basic types
                TypeOpenParenthesis, // '('
                TypeCloseParenthesis, // ')'
                // TypeEqual, // '='
                TypePlus, // '+'
                TypeHyphen, // '-'
                TypeAsterisk, // '*'
                TypeForwardSlash, // '/'
                TypeHat, // '^'
                // TypeAmpersand, // '&'
                // TypeComma, // ','
                // TypeSemicolon, // ';'
                TypeIdentifier, // name
                TypeFunction, // math function
                TypeValue, // 0.2, .6, 2.1
            };

            // Constructors
            Token(TokenType type);
            Token(std::string const&  ident, // May throw if can not copy string
                                   int  associativity = 0);
            Token(double value);

            // Information
            TokenType               getType() const;
            std::string const&    getIdentifier() const;
            double                 getValue() const;

            unsigned               getPrecedence() const;
            int                     getAssociativity() const;

            bool operator< (Token const& other) const;

        private:

            TokenType        _type;      // Type of token
            std::string     _ident;     // If type is Identifier
            double          _value;     // If type is Value

            unsigned        _precedence;
            int              _associativity;
    };



    // Parser class
    //--------------------------------------------------------------------------
    class Parser
    {
        public:
            typedef std::vector<Token>::size_type size_type;

            Parser();
            ~Parser();

            Expression parse(std::string const& exstr);

            //

            Token const& operator[] (size_type pos) const;

            //

        private:

            Expression buildExpression(std::vector<Token>& stack);

            void buildTokens(std::string const& exstr);
            //

            std::vector<Token>  _tokens; // Store token and not pointer to token
    };
}

#endif // __PARSER_H

