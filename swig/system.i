// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-11-10 td
// last changed:
//
%module system

%{
#include <common/Types.h>
#include <system/Expression.h>
#include <system/BioSystem.h>
#include <system/BioPAR.h>
#include <system/BioProcessor.h>

#include <sstream>
%}


%include <std_string.i>
%include <std_vector.i>
%include <std_map.i>
%include <std_pair.i>

%template(Param)                std::map< std::string, double >;
%template(StringList)           std::vector< std::string >;                             // for Species, Parameter, Auxiliary
%template(ValueList)            std::vector< double >;                                  // for ODESolver::Grid
%template(VectorList)           std::vector< PARKIN::Vector >;                          // for ODESolver::Grid
%template(ValuePair)            std::pair< double, double >;                            // for MeasurementPoint, see next line
%template(MeasurementPoint)     std::map< std::string, std::pair<double,double> >;      // mapping from species name to pair (meas.value, std.dev)
%template(Trajectory)           std::map< unsigned, std::vector<double> >;              // for ODESolver (forward declaration)
%template(TrajectoryMap)        std::map< std::string, std::vector<double> >;           // for BioProcessor (forward declaration)

%include linalg.i
%include nonlin.i
%include odelib.i

%ignore PARKIN::Expression::operator=;
%ignore PARKIN::BioRHS::operator=;
%ignore PARKIN::operator<<;


%include "common/Types.h"
%include "system/Expression.h"
%include "system/BioRHS.h"
%include "system/BioSystem.h"
%include "system/BioPAR.h"
%include "system/BioProcessor.h"


%extend PARKIN::Expression
{

    std::string __str__()
    {
        std::ostringstream s;
        (*self).prt(s);
        return s.str();
    }

};

%ignore std::vector< PARKIN::Expression >::vector(size_type);
%ignore std::vector< PARKIN::Expression >::resize(size_type);

%template(MeasurementList)      std::vector< PARKIN::MeasurementPoint >;
// %template(ExpressionList)       std::vector< PARKIN::Expression >;
%template(ExpressionMap)        std::map< std::string, PARKIN::Expression >;


