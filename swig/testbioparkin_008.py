"""This script illustrates the use of the python parkin package"""
from math import sqrt
from parkin import Param, StringList, BioSystem, ExpressionMap, Expression, TIMES, POWER, PLUS, MINUS, Vector, MeasurementList, IOpt, GaussNewtonWk, GaussNewton, BioPAR, Matrix, ValueList, BioProcessor, EPMACH


__author__ = "Moritz Wade (wade at zib dot de)"
__date__ = "2011-10-07"
__copyright__ = "Copyright (C) 2010-2011  ZIB - Zuse Institute Berlin, Germany"

#

import sys		# for sys.exit()
#from parkin import *

jacgen = 3
args = sys.argv
if len(args) > 1:
    jacgen = int(args[1])
#

tstart = 0.0
tend   = 100.0


parValues, parThres = Param(), Param()

speciesIDs = StringList()
paramIDs = StringList()
biosys = BioSystem(tstart, tend)

aux = ExpressionMap()
emap = ExpressionMap()

# set names / identifies of species
speciesIDs.push_back("C")
speciesIDs.push_back("X")
speciesIDs.push_back("M")
speciesIDs.push_back("Y")
speciesIDs.push_back("Z")
# biosys.setSpecies(species)

# set names / identifies of parameters

paramIDs.push_back("glo_V1p");     parValues["glo_V1p"] = 0.75
paramIDs.push_back("r04_K1");      parValues["r04_K1"]  = 0.1

paramIDs.push_back("r05_V2");      parValues["r05_V2"]  = 0.25
paramIDs.push_back("r05_K2");      parValues["r05_K2"]  = 0.1

# param.push_back("V3");      var["V3"]  = 0.0
paramIDs.push_back("glo_V3p");     parValues["glo_V3p"] = 0.3
paramIDs.push_back("r06_K3");      parValues["r06_K3"]  = 0.2

paramIDs.push_back("r07_V4");      parValues["r07_V4"]  = 0.1
paramIDs.push_back("r07_K4");      parValues["r07_K4"]  = 0.1

paramIDs.push_back("r02_K5");      parValues["r02_K5"]  = 0.02
paramIDs.push_back("glo_K6");      parValues["glo_K6"]  = 0.3

paramIDs.push_back("r01_vi");      parValues["r01_vi"]  = 0.1
paramIDs.push_back("r02_k1");      parValues["r02_k1"]  = 0.5
paramIDs.push_back("r03_kd");      parValues["r03_kd"]  = 0.02
paramIDs.push_back("r11_kd");      parValues["r11_kd"]  = 0.02
paramIDs.push_back("r08_a1");      parValues["r08_a1"]  = 0.05
paramIDs.push_back("r09_a2");      parValues["r09_a2"]  = 0.05
paramIDs.push_back("r10_d1");      parValues["r10_d1"]  = 0.05
paramIDs.push_back("r13_d1");      parValues["r13_d1"]  = 0.05
paramIDs.push_back("r12_vs");      parValues["r12_vs"]  = 0.2
paramIDs.push_back("r10_alpha");   parValues["r10_alpha"] = 0.1
paramIDs.push_back("r11_alpha");   parValues["r11_alpha"] = 0.1


aux["glo_V1"] = Expression(
                TIMES,
                "C",
                Expression(
                    TIMES,
                    "glo_V1p",
                    Expression(
                        POWER,
                        Expression(PLUS, "C", "glo_K6"),
                        -1.0
                    )
                )
            )

aux["glo_V3"] = Expression(TIMES, "glo_V3p", "M")

#//
#
#// Reaction / rule:
#//
#//      r1  :  vi
#//      r2  :  C * k1 * X * pow( C + K5, -1 )
#//      r3  :  C * kd
#//      r4  :  (1 + -M) * V1 * pow( K1 + -M + 1, -1 )
#//      r5  :  M * V2 * pow( K2 + M, -1 )
#//      r6  :  V3 * (1 + -X) * pow( K3 + -X + 1, -1 )
#//      r7  :  V4 * X * pow( K4 + X, -1 )
#//      r8  :  a1 * C * Y
#//      r9  :  a2 * Z
#//      r10 :  alpha * d1 * Z
#//      r11 :  alpha * kd * Z
#//      r12 :  vs
#//      r13 :  d1 * Y

aux["r1"] = Expression("r01_vi")

aux["r2"] = Expression(
                TIMES,
                "C",
                Expression(
                    TIMES,
                    "r02_k1",
                    Expression(
                        TIMES,
                        "X",
                        Expression(
                            POWER,
                            Expression(PLUS, "C", "r02_K5"),
                            -1.0
                        )
                    )
                )
            )

aux["r3"] = Expression(TIMES, "C", "r03_kd")

aux["r4"] = Expression(
                TIMES,
                Expression(PLUS, 1.0, Expression(MINUS, "M")),
                Expression(
                    TIMES,
                    aux["glo_V1"],
                    Expression(
                        POWER,
                        Expression(PLUS, "r04_K1", Expression(PLUS, Expression(MINUS, "M"), 1.0)),
                        -1.0
                    )
                )
            )

aux["r5"] = Expression(
                TIMES,
                "M",
                Expression(
                    TIMES,
                    "r05_V2",
                    Expression(
                        POWER,
                        Expression(PLUS, "r05_K2", "M"),
                        -1.0
                    )
                )
            )

aux["r6"] = Expression(
                TIMES,
                aux["glo_V3"],
                Expression(
                    TIMES,
                    Expression(PLUS, 1.0, Expression(MINUS, "X")),
                    Expression(
                        POWER,
                        Expression(PLUS, "r06_K3", Expression(PLUS, Expression(MINUS, "X"), 1.0)),
                        -1.0
                    )
                )
            )

aux["r7"] = Expression(
                TIMES,
                "r07_V4",
                Expression(
                    TIMES,
                    "X",
                    Expression(
                        POWER,
                        Expression(PLUS, "r07_K4", "X"),
                        -1.0
                    )
                )
            )

aux["r8"] = Expression(TIMES, "r08_a1", Expression(TIMES, "C", "Y"))

aux["r9"] = Expression(TIMES, "r09_a2", "Z")

aux["r10"] = Expression(TIMES, "r10_alpha", Expression(TIMES, "r10_d1", "Z"))

aux["r11"] = Expression(TIMES, "r11_alpha", Expression(TIMES, "r11_kd", "Z"))

aux["r12"] = Expression("r12_vs")

aux["r13"] = Expression(TIMES, "r13_d1", "Y")

#//
#
#// ODEs:
#//
#//      C' = ( 1*r1) + (-1*r2) + (-1*r3) + (-1*r8) + ( 1*r9) + ( 1*r10)
#//         = (vi) + (- C * k1 * X * pow(K5 + C,-1)) + (- C * kd) +
#//           (- a1 * C * Y) + (a2 * Z) + (alpha * d1 * Z)
#//
#//      X' = ( 1*r6) + (-1*r7)
#//         = (V3p * M * (1 - X) * pow(K3 - X + 1, -1)) + (- V4 * X * pow(K4 + X, -1))
#//
#//      M' = ( 1*r4) + (-1*r5)
#//         = (C * V1p * (1 - M) * pow(K6 + C, -1) * pow(K1 - M + 1, -1)) +
#//           (- M * V2 * pow(K2 + X, -1))
#//
#//      Y' = (-1*r8) + ( 1*r9) + ( 1*r11) + ( 1*r12) + (-1*r13)
#//         = (- a1 * C * Y) + (a2 * Z) + (alpha *kd * Z) + (vs) + (- d1 * Y)
#//
#//      Z' = ( 1*r8) + (-1*r9) + (-1*r10) + (-1*r11)
#//           (a1 * C * Y) + (- a2 * Z) + (- alpha * d1 * Z) + (- alpha * kd * Z)

emap["C"] = Expression(
                        PLUS,
                        Expression(TIMES, 1.0, aux["r1"]),
                        Expression(
                            PLUS,
                            Expression(MINUS, aux["r2"]),
                            Expression(
                                PLUS,
                                Expression(MINUS, aux["r3"]),
                                Expression(
                                    PLUS,
                                    Expression(MINUS, aux["r8"]),
                                    Expression(
                                        PLUS,
                                        Expression(TIMES, 1.0, aux["r9"]),
                                        Expression(TIMES, 1.0, aux["r10"])
                                    )
                                )
                            )
                        )
                    )

emap["X"] = Expression(
                        PLUS,
                        Expression(TIMES, 1.0, aux["r6"]),
                        Expression(MINUS, aux["r7"])
                    )

emap["M"] = Expression(
                        PLUS,
                        Expression(TIMES, 1.0, aux["r4"]),
                        Expression(MINUS, aux["r5"])
                    )

emap["Y"] = Expression(
                        PLUS,
                        Expression(MINUS, aux["r8"]),
                        Expression(
                            PLUS,
                            Expression(TIMES, 1.0, aux["r9"]),
                            Expression(
                                PLUS,
                                Expression(TIMES, 1.0, aux["r11"]),
                                Expression(
                                    PLUS,
                                    Expression(TIMES, 1.0, aux["r12"]),
                                    Expression(MINUS, aux["r13"])
                                )
                            )
                        )
                    )

emap["Z"] = Expression(
                        PLUS,
                        Expression(TIMES, 1.0, aux["r8"]),
                        Expression(
                            PLUS,
                            Expression(MINUS, aux["r9"]),
                            Expression(
                                PLUS,
                                Expression(MINUS, aux["r10"]),
                                Expression(MINUS, aux["r11"])
                            )
                        )
                    )



biosys.setODESystem(emap)
#// biosys.setSpecies(species);  // has now a complete different effect:
#                                // resets to an identity 'emap' mapping
biosys.setParameters(paramIDs)


biosys.setInitialValue("C", 0.0)
biosys.setInitialValue("X", 0.0)
biosys.setInitialValue("M", 0.0)
biosys.setInitialValue("Y", 1.0)
biosys.setInitialValue("Z", 1.0)




# set the time points of taken measurements
timepoints = [tstart + j*(tend-tstart)/10.0 for j in range(1, 11)]
timepoints = [10.0, 20.0]
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


biosys.setSolverRTol(1.0e-7)
biosys.setSolverATol(1.0e-9)

# compute trajectories of model (i.e. model values) for a specific set of variables / parameters

biosys.computeModel(parValues) # ,"adaptive")

# extract and print out the computed solution vectors
tp = biosys.getOdeTrajectoryTimePoints()
print "t =\n%s" % tp.t()

for i in xrange(len(speciesIDs)):
#for i in xrange(2):
    speciesTrajectory = biosys.getOdeTrajectory(i)
    print("%s =\n%s" % (speciesIDs[i],speciesTrajectory.t()))


## from here start sensitivity computation
#
proc = BioProcessor( biosys, "nlscon" )
iopt = proc.getIOpt()

iopt.mode      = 0       # 0: normal run, 1: single step
iopt.jacgen    = jacgen  # 1: user supplied Jacobian, 2: num.diff., 3: num.diff.(with feedback)
iopt.qrank1    = False   # allow Broyden rank-1 updates if True
iopt.nonlin    = 3       # 1: linear, 2: mildly nonlin., 3: highly nonlin., 4: extremely nonlin.
iopt.norowscal = False   # allow for automatic row scaling of Jacobian if False
iopt.transf    = 0
iopt.itmax     = 45
iopt.mprmon    = 2
iopt.mprerr    = 1



# if this is done here, we get a -46
for paramID in paramIDs:
    parThres[paramID] = EPMACH

#proc = BioProcessor( biosys, "nlscon" )
proc.setIOpt( iopt )
proc.setCurrentParamValues( parValues )


# if this is done here, everything works
#for paramID in paramIDs:
#    parThres[paramID] = EPMACH

proc.setCurrentParamThres( parThres )




speThres = proc.getCurrentSpeciesThres()

for speciesIDs in speThres.keys():
    speThres[speciesIDs] = EPMACH

proc.setCurrentSpeciesThres( speThres )



timepoints = [ 30.005 ]
tsttp = Vector( ValueList(timepoints) )

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



