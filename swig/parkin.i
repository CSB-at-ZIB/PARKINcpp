// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-11-10 td
// last changed:
//
%module parkin


%include linalg.i
%include nonlin.i
%include odelib.i
%include system.i



%include <typemaps.i>

%{
#include <iostream>
%}


// Suppress warnings about unknown classes std::streambuf and std::ostream
%warnfilter(SWIGWARN_TYPE_UNDEFINED_CLASS) CPyOutbuf;
%warnfilter(SWIGWARN_TYPE_UNDEFINED_CLASS) CPyOstream;
// code from: http://www.nabble.com/Using-std%3A%3Aistream-from-Python-ts7920114.html#a7923253
%inline %{
class CPyOutbuf : public std::streambuf
{
    public:
        CPyOutbuf(PyObject* obj) : m_PyObj(obj) { Py_INCREF(m_PyObj); }
        ~CPyOutbuf() { Py_DECREF(m_PyObj); }
    protected:
        int_type overflow(int_type c)
        {
            // Call to PyGILState_Ensure ensures there is Python
            // thread state created/assigned.
            PyGILState_STATE gstate = PyGILState_Ensure();
            PyObject_CallMethod(m_PyObj, (char*) "write", (char*) "c", c);
            PyGILState_Release(gstate);
            return c;
        }
        std::streamsize xsputn(const char* s, std::streamsize count)
        {
            // Call to PyGILState_Ensure ensures there is Python
            // thread state created/assigned.
            PyGILState_STATE gstate = PyGILState_Ensure();
            PyObject_CallMethod(m_PyObj, (char*) "write", (char*) "s#", s, int(count));
            PyGILState_Release(gstate);
            return count;
        }

        PyObject* m_PyObj;
};

class CPyOstream : public std::ostream
{
    public:
        CPyOstream(PyObject* obj) : std::ostream(&m_Buf), m_Buf(obj) {}
    private:
        CPyOutbuf m_Buf;
};
%}

%pythoncode %{
    def MakeCppOStream(file):
        os = CPyOstream(file)
        os.thisown = False
        return os
        # ls = LogStream_ostream(os)
        # ls.thisown = False
        # return ls

%}

