###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        To read in and solve an instance of the multiparametric
#                   Linear Complementarity Problem (mpLCP) using the methodology
#                   presented in "Advancing Parametric Optimization," Springer
#                   2021.
#
################################################################################

# Getting Started
from cypari2 import Pari
import sys
import re
import logging
import time
import pyomo.environ as pyo
from pyomo.opt import SolverFactory, TerminationCondition
import multiprocessing
from read_flags import *
from read_problem import *
from create_pyomo_models import *
from matrix_manipulation import *
from lemke import *
from crisscross import *
from inv_region import *
from counter import Counter


# Initialize pari
pari = Pari()
pari.default("nbthreads", multiprocessing.cpu_count())

#turn off pyomo warnings
logging.getLogger('pyomo.core').setLevel(logging.ERROR)

# Declare Variables, etc
numVar          = 0
numParam        = 0
gxInitialized   = False
mIsNumeric      = True
useLemke        = False
useCrissCross   = False
nlpsAsFeasProbs = False
checkFwithEH    = True
checkDimWithF   = True
linearSolver    = "glpk"
nonlinearSolver = ""
gMatrix         = 0
xVar            = 0
paramSpace      = []
startingPoint   = []
epsilon         = 0.0000001
parallelPivots  = False
parallel        = True
feasible        = False
originalGmatrix = []
originalBasis   = list(range(numVar))
Q0 = []
midpoint = []
ellipsoidTol    = 0.01
ellipsoidIter   = 50


# Set parameters using command line flags
if len(sys.argv) > 2:
    linearSolver, nonlinearSolver, parallelPivots, useCrissCross, nlpsAsFeasProbs, checkFwithEH, checkDimWithF = ReadFlags(sys, 
                                                                                                                            logging, 
                                                                                                                            linearSolver, 
                                                                                                                            nonlinearSolver, 
                                                                                                                            parallelPivots, 
                                                                                                                            useCrissCross, 
                                                                                                                            nlpsAsFeasProbs,
                                                                                                                            checkFwithEH,
                                                                                                                            checkDimWithF)

# Read in the problem instance
numVar, numParam, gMatrix, xVar, paramSpace, mIsNumeric = ReadFile(pari, sys, re, numVar, numParam, gMatrix, xVar, paramSpace, gxInitialized, mIsNumeric)
originalGmatrix = [row[:] for row in gMatrix] #deep copy

if mIsNumeric:
    logging.warning("Warning: The data entered consists of an M matrix containing no parameters. While the method implemented here is applicable for this problem, a more efficient procedure exists. See Adelgren and Wiecek's 'A two phase algorithm for the multiparametric linear complementarity problem' (2016). This method may implemented here in a future release, but is not as of now. Continuing ... ")


#Find a point in the relative interior of the parameter space
constraint_expressions = []
constraint_expressions_rel_int = []
for i in range(len(paramSpace)):
    constraint_expressions.append(pariToPyomo(re, str(paramSpace[i][0].Str())) + " <= " + pariToPyomo(re, str(paramSpace[i][1].Str())))
    constraint_expressions_rel_int.append(pariToPyomo(re, str(paramSpace[i][0].Str())) + " + model.x[" + str(numParam + 1) + "] <= " + pariToPyomo(re, str(paramSpace[i][1].Str())))
model, opt = CreatePyomoModel(pyo, numParam + 1, "model.x[" + str(numParam + 1) + "]", pyo.maximize, constraint_expressions_rel_int, linearSolver)
results = opt.solve(model)

if results.solver.termination_condition != "optimal":
    sys.exit("Unable to find point in relative interior of parameter space. Exiting!")
else:
    for v in model.component_data_objects(pyo.Var, active=True):
        startingPoint.append(pyo.value(v))

printMatrix(gMatrix)

#If no nonlinear solver is specified, NLP's will not be solved to optimility. Instead, we will attempt to find a feasible points with strictly positive optimal value using modified versions of the ellipsoid method. We now attempt to find a starting ellipse containing the parameter space.
if nonlinearSolver == "":
    Q0 = [ [0 for j in range(numParam + 2)] for i in range(numParam + 2)]
    lowVals = []
    highVals = []
    for i in range(numParam):
        model, opt = CreatePyomoModel(pyo, numParam + 1, "model.x[" + str(i + 1) + "]", pyo.maximize, constraint_expressions, linearSolver)
        results = opt.solve(model)
        if results.solver.termination_condition != "optimal":
            sys.exit("Unable to find maximum value of parameter " + str(i+1) + " within the given parameter space. Exiting!")
        else:
            highVals.append(pyo.value(model.o))
        model, opt = CreatePyomoModel(pyo, numParam + 1, "model.x[" + str(i + 1) + "]", pyo.minimize, constraint_expressions, linearSolver)
        results = opt.solve(model)
        if results.solver.termination_condition != "optimal":
            sys.exit("Unable to find minimum value of parameter " + str(i+1) + " within the given parameter space. Exiting!")
        else:
            lowVals.append(pyo.value(model.o))
    highVals.append(2*startingPoint[-1])
    highVals.append(1)
    lowVals.append(0)
    lowVals.append(0)
    distanceSq = 0.0
    for i in range(len(Q0)):
        Q0[i][i] = len(Q0)/4*(highVals[i] - lowVals[i])**2
        midpoint.append((highVals[i] + lowVals[i])/2)
    for row in Q0:
        print(row)
    print(midpoint)

# Perform Phase 1
if useCrissCross:
    basis, gMatrix, feasible = CrissCross(pari, logging, numVar, gMatrix, xVar, startingPoint, parallelPivots, epsilon)
else:
    basis = originalBasis

print("Current Basis:")
print(basis)
rgn0 = InvRgn(re, pari, pyo, gMatrix, basis, xVar, startingPoint, epsilon, nonlinearSolver, Q0, midpoint, ellipsoidTol, ellipsoidIter, paramSpace, nlpsAsFeasProbs, checkFwithEH, checkDimWithF)

if feasible:

    dim = rgn0.Dim(re, pari, pyo)
    print("\n\nDim:")
    print(dim)
    fullDim = True
    if dim < len(xVar) - 2:
        fullDim = False
    
    if not fullDim:
        gMatrix = originalGmatrix
        basis = originalBasis
        rgn0 = InvRgn(gMatrix, basis)
        feasible = False
        
manager = multiprocessing.Manager()
discoveredBases = manager.list([basis])
q = multiprocessing.Queue()
q.put(rgn0)
numToProcess = Counter(0)
numToProcess.Increment()
finalPartition = multiprocessing.Queue()
procs = []

if not feasible:   # Traditional Phase 1 is required for determining the starting basis
    sys.exit("Write Code for Phase 1")

# If a feasible starting solution has been found after Phase 1, start Phase 2
if feasible:
    print("\n\nFeasible. Process Region.")
    while numToProcess.Value() > 0:
        if not q.empty():
            rgn = q.get()
            if parallel:
                p = multiprocessing.Process(target=ProcessRegion, args=(rgn, q, discoveredBases, finalPartition, re, pari, pyo, paramSpace, parallelPivots))
                p.start()
                procs.append(p)
            else:
                ProcessRegion(rgn, q, discoveredBases, finalPartition, re, pari, pyo, paramSpace, parallelPivots)
        numToProcess.Decrement()
    
    if parallel:
        for p in procs:
            p.join()
#    while not q.empty
#    processes = []
#    for k in range(len(M)):
#        p = multiprocessing.Process(target=ProcessRow, args=(q, M[k], M[i], k, i, j))
#        processes.append(p)
#        p.start()
#    for p in processes:
#        ret = q.get() 
#        order.append(ret[0])
#        rows.append(ret[1])
#    for p in processes:
#        p.join()
#    reorder = order.copy()
#    for k in order:
#        reorder[order[k]] = k
#    M = [rows[k] for k in reorder]

printMatrix(gMatrix)
print('\n')
printMatrix(paramSpace)
print('\n')
print(mIsNumeric)
#printLatexTableau(basis, gMatrix)

