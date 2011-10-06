"""This script illustrates the use of the python parkin package"""

__author__ = "Thomas Dierkes (dierkes at zib dot de)"
__date__ = "2010-12-20"
__copyright__ = "Copyright (C) 2010-2011  ZIB - Zuse Institue Berlin, German"

#

import sys		# for sys.exit()
from parkin import *

#

tstart = 1.0
tend   = 14.0

species = StringList()
parameter = StringList()
emap = ExpressionMap()
biosys = BioSystem(tstart, tend)


# set names / identifies of species
species.push_back("s1")
species.push_back("s2")
# biosys.setSpecies(species)

# set names / identifies of parameters
parameter.push_back("compartment")
parameter.push_back("k1")
biosys.setParameters(parameter)

# set expressions for ODE system
#elist.push_back( Expression(TIMES, parameter[0], species[1]) )
#elist.push_back( Expression(TIMES, Expression(MINUS, parameter[1]), species[0]) )
aux = ExpressionMap()
aux["react1"] = Expression(TIMES, "compartment", Expression(TIMES, "k1", "s1"))
emap["s1"] = Expression(MINUS, aux["react1"])
emap["s2"] = aux["react1"]

biosys.setODESystem(emap)

# set the initial value(s) for ODE system
biosys.setInitialValue( species[0], 1.5e-4 )
biosys.setInitialValue( species[1], 0.0 )


# set the time points of taken measurements
timepoints = [tstart + j*(tend-tstart)/10.0 for j in range(1, 11)]
print "Timepoints: %s" % timepoints
meastp = Vector( ValueList(timepoints) )
biosys.setMeasurementTimePoints( meastp )
print meastp.t()


## set breakpoints of events
#breakpoints = [tstart + j*(tend-tstart)/2.0 for j in range(0,3)]
#print "Breakpoints: %s" % breakpoints
#breaktp = Vector( ValueList(breakpoints) )
#biosys.setBreakpoints( breaktp )
#print breaktp.t()


biosys.setSolverRTol(1.0e-6)
biosys.setSolverATol(1.0e-8)

# compute trajectories of model (i.e. model values) for a specific set of variables / parameters
var = Param()
var["compartment"] = 1.0
var["k1"] = 1.5
#print var.items()

biosys.computeModel(var,"init")

# extract and print out the computed solution vectors
tp = biosys.getOdeTrajectoryTimePoints()
s1 = biosys.getOdeTrajectory(0)
s2 = biosys.getOdeTrajectory(1)
s1sim = biosys.getSimTrajectoryPoints(0)
s2sim = biosys.getSimTrajectoryPoints(1)
print "t = "
print tp.t()
print species[0]+" = "
print s1.t()
print species[1]+" = "
print s2.t()
print "Simulated %s = %s" % (species[0], s1sim.t())
print "Simulated %s = %s" % (species[1], s2sim.t())


#
# prepare for solution of inverse ODE problem
#

mlist     = biosys.getMeasurementList()
invbiosys = BioSystem( tstart, tend )

# invbiosys.setSpecies(species)
invbiosys.setParameters(parameter)

invbiosys.setODESystem(emap)
invbiosys.setInitialValue( species[0], biosys.getInitialValue(species[0]) )
invbiosys.setInitialValue( species[1], biosys.getInitialValue(species[1]) )

invbiosys.setParamValue( parameter[0], var[parameter[0]] )
invbiosys.setParamValue( parameter[1], var[parameter[1]] )

invbiosys.setMeasurementList( meastp, mlist )

# setting solution accuracies
invbiosys.setSolverRTol(1.0e-7)
invbiosys.setSolverATol(1.0e-8)
xtol = 1.0e-4


par1 = Param()
# par1["compartment"] = 2.0
par1["k1"] = 0.8

par1Thres = Param()
# par1Thres["compartment"] = EPMACH
par1Thres["k1"] = EPMACH

iopt = IOpt()

iopt.mode      = 0      # 0: normal run, 1: single step
iopt.jacgen    = 3      # 1: user supplied Jacobian, 2: num.diff., 3: num.diff.(with feedback)
iopt.qrank1    = False  # allow Broyden rank-1 updates if True
iopt.nonlin    = 3      # 1: linear, 2: mildly nonlin., 3: highly nonlin., 4: extremely nonlin.
iopt.norowscal = False  # allow for automatic row scaling of Jacobian if False
iopt.lpos      = False
iopt.itmax     = 45
iopt.mprmon    = 2
iopt.mprerr    = 1


proc = BioProcessor( invbiosys, "nlscon" )
proc.setIOpt( iopt )
proc.setCurrentParamValues( par1 )
proc.setCurrentParamThres( par1Thres )


speThres = proc.getCurrentSpeciesThres()

for species in speThres.keys():
    speThres[species] = EPMACH

proc.setCurrentSpeciesThres( speThres );


## ... and identify!
#proc.identifyParameters(xtol)
#
##sys.exit()
#
#final = proc.getIdentificationResults()
#print
#print "Inv.Prob. solution:"
#print "-------------------"
#for parm in final.keys():
#    var[parm] = final[parm]
#    print "%10s = %e" % ( parm, final[parm] )
#print
#
#
#print
#print "Final Solution (for all system parameters)"
#print "------------------------------------------"
#print "%25s:  %f" % ( parameter[0], var[parameter[0]] )
#print "%25s:  %f" % ( parameter[1], var[parameter[1]] )


## from here start sensitivity computation
#
proc.setCurrentParamValues( par1 )
proc.setCurrentParamThres( par1Thres )


speThres = proc.getCurrentSpeciesThres()

for species in speThres.keys():
    speThres[species] = EPMACH

proc.setCurrentSpeciesThres( speThres )

timepoints = [ 5.56, 10.5 ]
tsttp = Vector( ValueList(timepoints) )
#tsttp.zeros(1)
#tsttp[0] = 10.5

print
print
print "*** prepDetailedSens() ***"
print "rc = %d" % (proc.prepareDetailedSensitivities( tsttp ))
print "***"

sensMat = proc.getSensitivityMatrices()

for j in range(0, len(sensMat)):
    print "Sensitivity for timepoint #%d (t = %f)" % (j+1, tsttp[j])
    print "mat (%d x %d) =" % (sensMat[j].nr(), sensMat[j].nc())
    print sensMat[j]

