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
useAprime       = True
linearSolver    = "glpk"
nonlinearSolver = ""
gMatrix         = 0
xVar            = 0
paramSpace      = []
startingPoint   = []
epsilon         = 0.0000001
parallelPivots  = False
parallel        = False #True
feasible        = False
originalGmatrix = []
originalBasis   = list(range(numVar))
Q0 = []
midpoint = []
ellipsoidTol    = 0.01
ellipsoidIter   = 50
outputFilename  = "Solution.txt"


# Set parameters using command line flags
if len(sys.argv) > 2:
    linearSolver, nonlinearSolver, parallelPivots, useCrissCross, nlpsAsFeasProbs, checkFwithEH, checkDimWithF, useAprime = ReadFlags(sys, 
                                                                                                                                        logging, 
                                                                                                                                        linearSolver, 
                                                                                                                                        nonlinearSolver, 
                                                                                                                                        parallelPivots, 
                                                                                                                                        useCrissCross, 
                                                                                                                                        nlpsAsFeasProbs,
                                                                                                                                        checkFwithEH,
                                                                                                                                        checkDimWithF,
                                                                                                                                        useAprime)

# Read in the problem instance
t = time.time()
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
print(feasible)

rgn0 = InvRgn(re, pari, pyo, gMatrix, basis, xVar, startingPoint, epsilon, nonlinearSolver, Q0, midpoint, ellipsoidTol, ellipsoidIter, paramSpace, nlpsAsFeasProbs, checkFwithEH, checkDimWithF, None)

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
    print(numToProcess.Value())
    #while numToProcess.Value() > 0:
    while not q.empty() or q.qsize() > 0:
        print("Size of Q: " + str(q.qsize()))
        if True: #not q.empty():
            rgn = q.get()
            if parallel:
                p = multiprocessing.Process(target=ProcessRegion, args=(rgn, q, discoveredBases, finalPartition, re, pari, pyo, paramSpace, parallelPivots))
                p.start()
                procs.append(p)
            else:
                ProcessRegion(rgn, q, discoveredBases, finalPartition, re, pari, pyo, paramSpace, parallelPivots)
        numToProcess.Decrement()
        print("Size of Q: " + str(q.qsize()))
        print(q.empty())
    
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

totalTime = time.time() - t

printMatrix(gMatrix)
print('\n')
printMatrix(paramSpace)
print('\n')
print(mIsNumeric)
#printLatexTableau(basis, gMatrix)


# Write the solution

outputFile = open(outputFilename, 'w')

print("The problem entered was an instance of mpLCP having the form\n", file = outputFile)
print("\tw - M(x)z = q(x)\n\tw'z = 0\n\tw,z >= 0\n", file = outputFile)

mx = max((len(str(ele.Str())) for row in originalGmatrix for ele in row[numVar:-1]))
print("with M(x) =\n", file = outputFile)
for row in originalGmatrix:
    print("\t[ " + "  ".join(["{:<{mx}}".format(str((-1*ele).Str()),mx=mx) for ele in row[numVar:-1]]) + " ]", file = outputFile)
    
mx = max((len(str(ele.Str())) for row in originalGmatrix for ele in row[2*numVar:]))
print("\nand q(x) =\n", file = outputFile)
for row in originalGmatrix:
    print("\t[ " + "  ".join(["{:<{mx}}".format(str(ele.Str()),mx=mx) for ele in row[2*numVar:]]) + " ]", file = outputFile)
    
mx = max((len(str(ele.Str())) for row in paramSpace for ele in row))
print("\nsubject to the additional restriction that 'x' must satisfy:\n", file = outputFile)
for row in paramSpace:
    print("\t" + " <= ".join(["{:<{mx}}".format(str(ele.Str()),mx=mx) for ele in row]), file = outputFile)
    
print("\n\n\n**************************************************************************************************\n\nThe solution was computed in " + str(round(totalTime, 2)) + " seconds and consists of the following regions.\n\n**************************************************************************************************\n\n", file = outputFile)
    
k = 1
while not finalPartition.empty():
    print("\n\nRegion " + str(k) + ":\n", file = outputFile)
    rgn = finalPartition.get()
    rhs = rgn.RHS()
    basis = rgn.Basis()
    mx = max((len(str(ele.Str())) for row in rhs for ele in row))
    for i in range(len(rhs)):
        var = ""
        if basis[i] < numVar:
            var = "w_" + str(i + 1)
        else:
            var = "z_" + str(i + 1)
        print("\t" + var + " = " + " ".join(["{:<{mx}}".format(str(ele.Str()),mx=mx) for ele in rhs[i]]) + " >= 0 ", file = outputFile)
    k = k + 1
        
print("\n\n\n\nNote: The region descriptions above do not include the 'additional restrictions' listed at the top of this document, although these restrictions do, of course, apply to all regions. Additionally, all omitted variables should be assumed to be zero.", file = outputFile)
        
outputFile.close()

