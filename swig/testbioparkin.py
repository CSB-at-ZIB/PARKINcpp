"""This script illustrates the use of the python parkin package"""

__author__ = "Thomas Dierkes (dierkes at zib dot de)"
__date__ = "2010-12-20"
__copyright__ = "Copyright (C) 2010-2011  ZIB - Zuse Institue Berlin, German"

#

import sys		# for sys.exit()
from parkin import *

#

tstart = 2.0
tend   = 8.0

species = StringList()
parameter = StringList()
elist = ExpressionMap()
biosys = BioSystem(tstart, tend)


# set names / identifies of species
species.push_back("x1")
species.push_back("x2")
# biosys.setSpecies(species)

# set names / identifies of parameters
parameter.push_back("a")
parameter.push_back("b")
biosys.setParameters(parameter)

# set expressions for ODE system
elist["x1"] =  Expression(TIMES, parameter[0], species[1])
elist["x2"] =  Expression(TIMES, Expression(MINUS, parameter[1]), species[0])
biosys.setODESystem(elist)

# set the initial value(s) for ODE system
biosys.setInitialValue( species[0],  5.0 )
biosys.setInitialValue( species[1], -3.0 )

# set the time points of taken measurements
meastp = Vector( ValueList([tstart + j*(tend-tstart)/9.0 for j in range(1,10)]) )
biosys.setMeasurementTimePoints( meastp )
#print meastp.t()

# compute trajectories of model (i.e. model values) for a specific set of variables / parameters
var = Param()
var[parameter[0]] = 2.0/3.0
var[parameter[1]] = 3.0
#print var.items()

biosys.computeModel(var,"adaptive")
#sys.exit()

# extract and print out the computed solution vectors
tp = biosys.getOdeTrajectoryTimePoints()
x1 = biosys.getOdeTrajectory(0)
x2 = biosys.getOdeTrajectory(1)
print "t = "
print tp.t()
print species[0]+" = "
print x1.t()
print species[1]+" = "
print x2.t()

#
synscal = Vector()

syndata = biosys.getMeasurements()
synscal.zeros( syndata.nr() )

print
print "  Simulated  Experimental  Data  :   nr() = ", syndata.nr()
print syndata
print

# prepare for solution of inverse ODE problem
rtol    = 1.0e-5
p       = Vector(2)
pscal   = Vector()
#p.zeros(2)
p[0] = 1.42
p[1] = 0.8
pscal.zeros(2)
synscal.zeros( syndata.nr() )

par1 = StringList()
par1.push_back("a")
par1.push_back("b")

prob = BioPAR(biosys)
prob.setParameter(par1)

gn = GaussNewton()
wk = gn.getWk()
iopt = IOpt()

iopt.mode      = 0
iopt.jacgen    = 3
iopt.qrank1    = False
iopt.nonlin    = 4
iopt.norowscal = False
iopt.mprmon    = 2
iopt.mprerr    = 1

gn.setProblem(prob)
gn.initialise(syndata.nr(), p,pscal, syndata,synscal, rtol, iopt,wk)

gn.run()
#gn.analyse()
gn.printCounter()
#sys.exit()

print
print "Inv.Prob. solution = "
print gn.getSolution().t()

