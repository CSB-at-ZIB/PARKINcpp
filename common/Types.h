// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 09.03.2010 td
// last changed:
//

#ifndef __TYPES_H
#define __TYPES_H

#include <map>
#include <string>
#include <vector>
#include <iostream>

///

namespace PARKIN
{
    typedef double Real;

    static std::map<int,std::string> errStrMap;

    class ErrType
    {
        public:
            ErrType(int ierr = 0, std::string s = "") :
                _ierr(ierr), _msg(s) { }
            ~ErrType() { }

            int getIerr() const
            { return _ierr; }
            void setIerr(int ierr, std::string s = "")
            { _ierr = ierr; _msg = s; }

            std::string getMsg() const
            { return _msg; }
            void setMsg(std::string s)
            { _msg = s; }

            std::ostream& operator<< (std::ostream& os)
            {
                return os << "message " << _ierr << ": " << _msg;
            }

        private:
            int             _ierr;
            std::string     _msg;
    };
}

#endif // __TYPES_H
