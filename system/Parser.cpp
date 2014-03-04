// Copyright (C) 2010 - 2013
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2013-09-25 td
// last changed:
//
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "Expression.h"
#include "Parser.h"

using namespace PARKIN;

// Private functions to solve locale problems
// ------------------------------------------
namespace
{
    bool parkin_isalpha(char c)
    {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
    }

    bool parkin_isdigit(char c)
    {
        return (c >= '0' && c <= '9');
    }

    bool parkin_isalnum(char c)
    {
        return parkin_isalpha(c) || parkin_isdigit(c);
    }

    bool parkin_isspace(char c)
    {
        return (c == ' ') || (c == '\t') || (c == '\r') || (c == '\n');
    }
};



// Token
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Constructor
Token::Token(TokenType type) :
        _type(type),
        _ident("./."),
        _value(0.0),
        _precedence(0),
        _associativity(0)
{
    switch( type )
    {
        case TypeComma:
            _ident = ",";
            _precedence = 8;
            break;
        case TypeOpenParenthesis:
            _ident = "(";
            _precedence = 8;
            break;
        case TypeCloseParenthesis:
            _ident = ")";
            _precedence = 8;
            break;
        case TypePlus:
            _ident = "PLUS";
            _precedence = 2;
            _associativity = -1; // left
            break;
        case TypeHyphen:
            _ident = "MINUS";
            _precedence = 2;
            _associativity = -1; // left
            break;
        case TypeAsterisk:
            _ident = "TIMES";
            _precedence = 3;
            _associativity = -1; // left
            break;
        case TypeForwardSlash:
            _ident = "DIVIDE";
            _precedence = 3;
            _associativity = -1; // left
            break;
        case TypeHat:
            _ident = "POWER";
            _precedence = 4;
            _associativity = 1; // right
            break;
        default:
            break;
    }
}
//------------------------------------------------------------------------------
// Construct identifier token
Token::Token(std::string const& ident, int associativity) :
        _type(Token::TypeIdentifier),
        _ident(ident),
        _value(0.0),
        _precedence(0),
        _associativity(associativity)
{
    if (associativity != 0)
    {
        _type = TypeFunction;
        _precedence = 9;
    }
}
//------------------------------------------------------------------------------
// Construct value token
Token::Token(double value) :
        _type(Token::TypeValue),
        _value(value),
        _precedence(0),
        _associativity(0)
{
    _ident.resize(64); // should be enough for a double value string
    std::sprintf(&_ident[0],"%E", value);
}
//------------------------------------------------------------------------------
// Get type
Token::TokenType
Token::getType() const
{
    return _type;
}
//------------------------------------------------------------------------------
// Get identifier
std::string const&
Token::getIdentifier() const
{
    return _ident;
}
//------------------------------------------------------------------------------
// Get value
double
Token::getValue() const
{
    return _value;
}
//------------------------------------------------------------------------------
// Get precedence
unsigned
Token::getPrecedence() const
{
    return _precedence;
}
//------------------------------------------------------------------------------
// Get associativity
int
Token::getAssociativity() const
{
    return _associativity;
}
//------------------------------------------------------------------------------
bool
Token::operator< (Token const& other) const
{
    if ( _associativity < 0 ) // left
    {
        return _precedence <= other._precedence;
    }
    else
    {
        return _precedence < other._precedence;
    }

    return true;
}
//------------------------------------------------------------------------------


// Parser
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Constructor
Parser::Parser() :
    _tokens()
{
    _tokens.clear();
    // if(expr == 0)
    //     throw(NullPointerException("Parser::Parser"));
}
//------------------------------------------------------------------------------
// Destructor
Parser::~Parser()
{
}
//------------------------------------------------------------------------------
// Parse an expression string
Expression
Parser::parse(std::string const& exstr)
{
    std::vector<Token>  rpnStack;
    std::vector<Token>  opStack;
    std::string         rpnStr = "";
    // _tokens.clear();

    buildTokens(exstr);

    // Make sure it is not still empty
    if(_tokens.size() == 0)
    {
        // throw(EmptyExpressionException());
        std::cerr << std::endl;
        std::cerr << "Error: Parser::Parse: Empty expression string.";
        std::cerr << std::endl;
        return exstr;
    }

    // Parse the range
    //return parseRegion(0, _tokens.size() - 1);

///

    // std::cout << std::endl;
    // std::cout << "#tokens = " << _tokens.size();
    // std::cout << std::endl;
    // std::cout << std::endl;
    //
    // for (size_type pos = 0; pos < _tokens.size(); ++pos)
    // {
    //     std::cout << " " << _tokens[pos].getIdentifier();
    // }
    // std::cout << std::endl;
    // std::cout << std::endl;

///

    for (size_type pos = 0; pos < _tokens.size(); ++pos)
    {
        Token::TokenType    op2Type;
        Token               term = _tokens[pos];

        switch( term.getType() )
        {
            case Token::TypePlus:
            case Token::TypeHyphen:
            case Token::TypeAsterisk:
            case Token::TypeForwardSlash:
            case Token::TypeHat:
                op2Type = (!opStack.empty()) ?
                                opStack.back().getType() : Token::TypeUnknown;

                if ( (op2Type == Token::TypePlus ||
                      op2Type == Token::TypeHyphen ||
                      op2Type == Token::TypeAsterisk ||
                      op2Type == Token::TypeForwardSlash ||
                      op2Type == Token::TypeHat) && (term < opStack.back()) )
                {
                    rpnStack.push_back( opStack.back() );
                    rpnStr += " " + opStack.back().getIdentifier();

                    opStack.pop_back();
                }
                opStack.push_back(term);
                break;

            case Token::TypeFunction:
                opStack.push_back(term);
                break;

            case Token::TypeComma:
                while ( opStack.back().getType() != Token::TypeOpenParenthesis )
                {
                    rpnStack.push_back( opStack.back() );
                    rpnStr += " " + opStack.back().getIdentifier();

                    opStack.pop_back();
                }
                break;

            case Token::TypeOpenParenthesis:
                opStack.push_back(term);
                break;

            case Token::TypeCloseParenthesis:
                while ( opStack.back().getType() != Token::TypeOpenParenthesis )
                {
                    rpnStack.push_back( opStack.back() );
                    rpnStr += " " + opStack.back().getIdentifier();

                    opStack.pop_back();
                }
                opStack.pop_back();
                if ( !opStack.empty() &&
                     opStack.back().getType() == Token::TypeFunction )
                {
                    rpnStack.push_back( opStack.back() );
                    rpnStr += " " + opStack.back().getIdentifier();

                    opStack.pop_back();
                }
                break;

            default:
                rpnStack.push_back( term );
                rpnStr += " " + term.getIdentifier();
                break;
        }
    }

    while ( !opStack.empty() )
    {
        rpnStack.push_back( opStack.back() );
        rpnStr += " " + opStack.back().getIdentifier();

        opStack.pop_back();
    }

///
    // std::cout << std::endl;
    // std::cout << "Parser::parse: rpnStr = '" << rpnStr << "'" << std::endl;
///

    return buildExpression(rpnStack);
}


// Get a token
Token const&
Parser::operator[] (Parser::size_type pos) const
{
    return _tokens[pos];
}


// Build expression recursively
Expression
Parser::buildExpression(std::vector<Token>& rpn)
{
    Expression arg1, arg2;

    if ( !rpn.empty() )
    {
        switch ( rpn.back().getType() )
        {
            case Token::TypePlus:
                rpn.pop_back();
                arg2 = buildExpression( rpn );
                rpn.pop_back();
                arg1 = buildExpression( rpn );
                return Expression(PARKIN::PLUS, arg1, arg2);
            case Token::TypeHyphen:
                rpn.pop_back();
                arg2 = buildExpression( rpn );
                rpn.pop_back();
                arg1 = buildExpression( rpn );
                return Expression(PARKIN::MINUS, arg1, arg2);

            case Token::TypeAsterisk:
                rpn.pop_back();
                arg2 = buildExpression( rpn );
                rpn.pop_back();
                arg1 = buildExpression( rpn );
                return Expression(PARKIN::TIMES, arg1, arg2);
            case Token::TypeForwardSlash:
                rpn.pop_back();
                arg2 = buildExpression( rpn );
                rpn.pop_back();
                arg1 = buildExpression( rpn );
                return Expression(PARKIN::DIVIDE, arg1, arg2);

            case Token::TypeHat:
                rpn.pop_back();
                arg2 = buildExpression( rpn );
                rpn.pop_back();
                arg1 = buildExpression( rpn );
                return Expression(PARKIN::POWER, arg1, arg2);

            case Token::TypeFunction:
            {
                std::string fun = rpn.back().getIdentifier();

                std::transform(fun.begin(),fun.end(),
                                 fun.begin(), ::tolower);

                if ( fun == "abs" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::ABS, arg1);
                }
                else if ( fun == "exp" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::EXP, arg1);
                }
                else if ( fun == "cos" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::COS, arg1);
                }
                else if ( fun == "sin" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::SIN, arg1);
                }
                else if ( fun == "tan" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::TAN, arg1);
                }
                else if ( fun == "cosh" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::COSH, arg1);
                }
                else if ( fun == "sinh" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::SINH, arg1);
                }
                else if ( fun == "tanh" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::TANH, arg1);
                }
                else if ( fun == "acos" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::ARCCOS, arg1);
                }
                else if ( fun == "asin" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::ARCSIN, arg1);
                }
                else if ( fun == "atan" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::ARCTAN, arg1);
                }
                else if ( fun == "acosh" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::ARCOSH, arg1);
                }
                else if ( fun == "asinh" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::ARSINH, arg1);
                }
                else if ( fun == "atanh" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::ARTANH, arg1);
                }
                else if ( fun == "ceil" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::CEIL, arg1);
                }
                else if ( fun == "floor" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::FLOOR, arg1);
                }
                else if ( fun == "log" || fun == "ln" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::LN, arg1);
                }
                else if ( fun == "log10" || fun == "lg" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::LOG, arg1);
                }
                else if ( fun == "sign" )
                {
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::SIGN, arg1);
                }
                else if ( fun == "pow" )
                {
                    rpn.pop_back();
                    arg2 = buildExpression( rpn );
                    rpn.pop_back();
                    arg1 = buildExpression( rpn );
                    return Expression(PARKIN::POWER, arg1, arg2);
                }
                else
                {
                    std::cerr << std::endl;
                    std::cerr << "Error: Parser::buildExpression: Unknown function.";
                    std::cerr << std::endl;
                    return Expression();
                }
            }

            case Token::TypeIdentifier:
                return Expression( rpn.back().getIdentifier() );

            case Token::TypeValue:
                return Expression( rpn.back().getValue() );

            default:
                std::cerr << std::endl;
                std::cerr << "Error: Parser::buildExpression: Unknown token syntax.";
                std::cerr << std::endl;
                return Expression();
        }
    }

    return Expression();
}


// Build tokens
void
Parser::buildTokens(std::string const& exstr)
{
    _tokens.clear();

    // Test zero-length expression
    if(exstr.length() == 0)
    {
        // throw(EmptyExpressionException());
        std::cerr << std::endl;
        std::cerr << "Error: Parser::BuildTokens: Empty expression string.";
        std::cerr << std::endl;
        return;
    }

    // Search through list
    std::string::size_type pos;
    bool comment = false;

    for(pos = 0; pos < exstr.length(); ++pos)
    {
///
        // std::cout << " c:'" << exstr[pos] << "'";
///
        // Take action based on character
        switch(exstr[pos])
        {
            // Comment
            case '%':
            case '#':
            {
                comment = true;

                break;
            }

            // Newline ends comment
            case '\r':
            case '\n':
            {
                comment = false;

                break;
            }

            // Open parenthesis
            case '(':
            {
                if(comment == false)
                {
                    _tokens.push_back(Token(Token::TypeOpenParenthesis));
                }
                break;
            }

            // Close parenthesis
            case ')':
            {
                if(comment == false)
                {
                    _tokens.push_back(Token(Token::TypeCloseParenthesis));
                }
                break;
            }

            // // Equal
            // case '=':
            // {
            //     if(comment == false)
            //     {
            //         _tokens.push_back(Token(Token::TypeEqual, 0));
            //     }
            //     break;
            // }

            // Plus
            case '+':
            {
                if(comment == false)
                {
                    _tokens.push_back(Token(Token::TypePlus));
                }
                break;
            }

            // Hyphen
            case '-':
            {
                if(comment == false)
                {
                    _tokens.push_back(Token(Token::TypeHyphen));
                }
                break;
            }

            // Asterisk
            case '*':
            {
                if(comment == false)
                {
                    if ( pos < exstr.length() && exstr[pos+1] == '*' )
                    {
                        _tokens.push_back(Token(Token::TypeHat));
                        ++pos;
                    }
                    else
                    {
                        _tokens.push_back(Token(Token::TypeAsterisk));
                    }
                }
                break;
            }

            // Forward slash
            case '/':
            {
                if(comment == false)
                {
                    _tokens.push_back(Token(Token::TypeForwardSlash));
                }
                break;
            }

            // Hat (exponent)
            case '^':
            {
                if(comment == false)
                {
                    _tokens.push_back(Token(Token::TypeHat));
                }
                break;
            }

            // // Ampersand
            // case '&':
            // {
            //     if(comment == false)
            //     {
            //         _tokens.push_back(Token(Token::TypeAmpersand, 0));
            //     }
            //     break;
            // }

            // Comma
            case ',':
            {
                if(comment == false)
                {
                    _tokens.push_back(Token(Token::TypeComma));
                }
                break;
            }

            // // Semicolon
            // case ';':
            // {
            //     if(comment == false)
            //     {
            //         _tokens.push_back(Token(Token::TypeSemicolon, 0));
            //     }
            //     break;
            // }

            // None of the above, but it may be an identifier or value
            default:
            {
                if(comment == false)
                {
                    // First, test for value
                    if(exstr[pos] == '.' || parkin_isdigit(exstr[pos]))
                    {
                        // We are a value
                        std::string::size_type start = pos;

                        // Digits before period
                        while(parkin_isdigit(exstr[pos]))
                            pos++;

                        // Period
                        if(exstr[pos] == '.')
                            pos++;

                        // Digits after period
                        while(parkin_isdigit(exstr[pos]))
                            pos++;

                        // Scientific notation
                        if (exstr[pos] == 'e' || exstr[pos] == 'E')
                        {
                            pos++;

                            if (exstr[pos] == '+' || exstr[pos] == '-')
                                pos++;

                            while(parkin_isdigit(exstr[pos]))
                                pos++;
                        }


                        // Create token
                        std::string ident = exstr.substr(start, pos - start);
                        char *buf;
                        _tokens.push_back(Token( std::strtod(ident.c_str(),&buf) ));

                        // Move pos back so pos++ will set it right
                        --pos;
                    }
                    else if(exstr[pos] == '_' || parkin_isalpha(exstr[pos]))
                    {
                        // We are an identifier
                        std::string::size_type start = pos;
                        // bool foundname = true; // Found name part

                        // Search for name, then period, etc
                        // An identifier can be multiple parts.  Each part
                        // is formed as an identifier and seperated by a period,
                        // An identifier can not end in a period
                        //
                        // color1.red : 1 identifier token
                        // color1. : color1 is identifier, . begins new token
                        // color1.1red : Not value (part 2 is not right)
                        // while(foundname)
                        {
                            // Part before period
                            while(exstr[pos] == '_' || parkin_isalnum(exstr[pos]))
                                pos++;

                            // // Is there a period
                            // if(exstr[pos] == '.')
                            // {
                            //     pos++;
                            //
                            //     // There is a period, look for the name again
                            //     if(exstr[pos] == '_' || parkin_isalpha(exstr[pos]))
                            //     {
                            //         foundname = true;
                            //     }
                            //     else
                            //     {
                            //         // No name after period
                            //         foundname = false;
                            //
                            //         // Remove period from identifier
                            //         pos--;
                            //     }
                            // }
                            // else
                            // {
                            //     // No period after name, so no new name
                            //     foundname = false;
                            // }
                        }

                        if (exstr[pos] == '(')
                        {
                            _tokens.push_back(Token(exstr.substr(start, pos - start), -1));
                        }
                        else
                        {
                            // Create token
                            _tokens.push_back(Token(exstr.substr(start, pos - start)));
                        }

                        // Move pos back so pos++ will set it right
                        --pos;
                    }
                    else if (parkin_isspace(exstr[pos]))
                    {
                        // Do nothing, just ignore white space, but it still
                        // seperates tokens
                    }
                    else
                    {
                        /*
                        // Unknown token
                        UnknownTokenException e;
                        e.SetStart(pos);
                        e.SetEnd(pos);

                        throw(e);
                        */
                        std::cerr << std::endl;
                        std::cerr << "Error: Parser::BuildTokens: Unknown token." << std::endl;
                        return;
                    }
                }
                break;
            }
        }
    }
    // // Final token
    // _tokens.push_back(Token(Token::TypeEndOfLine));
}


