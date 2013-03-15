// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-11-10 td
// last changed:
//
%module odelib

%{
// #include <common/Types.h>
#include <odelib/FirstOrderODESystem.h>
#include <odelib/ODETrajectory.h>
#include <odelib/ODESolver.h>
// #include <odelib/DOP853.h>
%}

// %include <std_map.i>
// %include <std_vector.i>
// %include <std_string.i>

// %template(Grid)       std::vector< double >;
// %template(Trajectory) std::map< unsigned, PARKIN::ODESolver::Grid >;

%include "common/Types.h"
%include "odelib/FirstOrderODESystem.h"
%include "odelib/ODETrajectory.h"
%include "odelib/ODESolver.h"
// %include "odelib/DOP853.h"

