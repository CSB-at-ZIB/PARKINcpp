// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-11-10 td
// last changed:
//
%module nonlin


%{
#include <nonlin/GaussNewton.h>
#include <nonlin/YeOldeParkinCore.h>
%}


%include linalg.i


// %apply int& INOUT { int& ifail };
// %ignore PARKIN::GaussNewton::setWk;
%ignore PARKIN::Info;

%include "nonlin/UserFunc.h"
%include "nonlin/GaussNewton.h"
%include "nonlin/YeOldeParkinCore.h"
