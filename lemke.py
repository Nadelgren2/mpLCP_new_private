###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        Use Lemke's method for the Linear Complementary Problem 
#                   (LCP) to solve an instance of multiparametric LCP (mpLCP) at
#                   a fixed point in the parameter space. If the fixed point
#                   happens to lie within the relative interior of a full 
#                   dimensional invariancy region, phase 1 of the algorithm
#                   presented in "Advancing Parametric Optimization" can be 
#                   skipped.
#
################################################################################

from matrix_manipulation import *

# Define Functions

# Utilize a simple implementation of Lemke's algorithm as presented in Chapter 
# 11 of the 3rd edition of Bazaraa's "Nonlinear Programming" in attempt to find
# a starting basis.
#
# Input:    flag -- the command line flag
#           message --  the warning message to print
#           warnings -- the warnings environment
def LemkesMethod(pari, logging, numVar, gMatrix, xVar, startingPoint, parallelPivots, useCrissCross):
    pivotRow = 0
    pivotCol = 0
    val = 0
    minVal = 0
    feasible = True
    keepGoing = True
    minValAssigned = False
    basis = list(range(numVar))
    originalGmatrix = [row[:] for row in gMatrix] #deep copy
    
    #add entries to the tableau for a dummy variable
    for i in range(len(gMatrix)):
        gMatrix[i].append(-pari.one())
   
    #Initialization -- Check if initial point is feasible. If not, enter the dummy variable
    for i in range(len(gMatrix)):
        val = pari.substvec(gMatrix[i][2*numVar], xVar[0:-1], startingPoint)
        if val < minVal:
            minVal = val
            pivotRow = i
            feasible = False
    
    if not feasible:
        basis[pivotRow] = 2*numVar + 1
        gMatrix = matrixPivot(gMatrix, pivotRow, -1, parallelPivots)
        for row in gMatrix:
            print(row)
        pivotCol = pivotRow + numVar
        
        #Main Step
        ite = 0
        while keepGoing:
            #Bazaraa Step 1
            ite += 1
            val = 0
            minVal = 0
            minValAssigned = False
            for i in range(len(gMatrix)):
                val = pari.substvec(gMatrix[i][pivotCol], xVar[0:-1], startingPoint)
                print(val)
                if val > 0.0:
                    val = pari.substvec(gMatrix[i][2*numVar], xVar[0:-1], startingPoint)/val
                    if not minValAssigned or val < minVal:
                        minValAssigned = True
                        minVal = val
                        pivotRow = i
            #Bazaraa Step 4
            if not minValAssigned:
                keepGoing = False
                mess = "Continuing to the traditional Phase 1 ... "
                if useCrissCross:
                    mess = "Now attempting to use the Criss Cross Method to find a starting basis ..."
                logging.warning("Lemke's Method exits with ray termination. No starting basis has been found. " + mess)
                gMatrix = originalGmatrix
            #Bazaraa Steps 2 and 3
            else:
                print(minVal)
                print(pivotRow)
                print(basis)
                #Bazaraa Step 3
                if basis[pivotRow] == 2*numVar + 1:
                    gMatrix = matrixPivot(gMatrix, pivotRow, pivotCol, parallelPivots)
                    basis[pivotRow] = pivotCol
                    keepGoing = False
                    feasible = True
                #Bazaraa Step 2
                else:
                    gMatrix = matrixPivot(gMatrix, pivotRow, pivotCol, parallelPivots)
                    for row in gMatrix:
                        print(row)
                    temp = pivotCol
                    if basis[pivotRow] < numVar:
                        pivotCol = basis[pivotRow] + numVar
                    else:
                        pivotCol = basis[pivotRow] - numVar 
                    print(temp)
                    print(pivotCol)
                    basis[pivotRow] = temp
                    print(basis)
                    print("Pivot Col is " + str(pivotCol))
                    keepGoing = True
            if ite > 2:
                keepGoing = False
        
    print("Finish Writing Lemke Code")
    return basis, gMatrix, feasible
