###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        Use the criss cross method for the Linear Complementary 
#                   Problem (LCP) to solve an instance of multiparametric LCP 
#                   (mpLCP) at a fixed point in the parameter space. If the 
#                   fixed point happens to lie within the relative interior of a
#                   full dimensional invariancy region, phase 1 of the algorithm
#                   presented in "Advancing Parametric Optimization" can be 
#                   skipped.
#
################################################################################

from matrix_manipulation import *

# Define Functions

# Print a warning indicating that the criss cross method has received an 'Exit'
# status
#
# Input:    logging --  the logging environment
#           startingPoint   --  a point in the relative interior of the 
#                               parameter space and at which we are attempting
#                               to identify our starting basis
def ExitWarning(logging, startingPoint):
    logging.warnings("The Criss Cross method received an 'exit' status and has therefore failed to 'process' this instance of (mp)LCP with the parameters fixed at " +
                        str(startingPoint[0:-1]) + ". Note that, as the criss cross method is guaranteed to 'process' instances of LCP in which M is a sufficient matrix, it is likely that your instance of mpLCP does not satisfy Assumption 1.1 of 'Advancing Parametric Optimization', namely, that M must be sufficient at every parameter vector in the parameter space. Nevertheless, we continue and attempt to find a starting basis using the standard Phase 1 procedure ... ")

# Utilize a simple implementation of the criss cross algorithm as presented in
#
#   den Hertog, D., Roos, C., & Terlaky, T. (1993). The linear complimentarity
#   problem, sufficient matrices, and the criss-cross method. Linear Algebra
#   and Its Applications, 187, 1-14. 
#
# in attempt to find a starting basis.
#
# Input:    pari    --  the pari environment
#           logging --  the logging environment
#           numVar  --  the number of variables present in the current instance
#           gMatrix --  the tableau representation of the current instance
#           xVar    --  the array containing the pari variables used to 
#                       represent the instance's parameters
#           startingPoint   --  a point in the relative interior of the 
#                               parameter space and at which we are attempting
#                               to identify our starting basis
#           parallelPivots  --  a boolean indicating whether or not principal
#                               pivots carried out on gMatrix should be computed
#                               in parallel
#           epsilon         --  a small value used to avoid numerical issues
#
# Output:   basis   --  a list indicating the basic variables at the starting
#                       solution
#           gMatrix --  an updated tableau representing the solution to the 
#                       mpLCP at the current basis
#           feasible    --  a boolean indicating whether or not the criss cross
#                           method discovered a feasible solution to the (mp)LCP
def CrissCross(pari, logging, numVar, gMatrix, xVar, startingPoint, parallelPivots, epsilon):
    pivotRow = -1
    pivotRow2 = 0
    pivotCol = 0
    pivotCol2 = 0
    val = 0
    val2 = 0
    feasible = True
    keepGoing = True
    pivotFound = False
    basis = list(range(numVar))
    originalGmatrix = [row[:] for row in gMatrix] #deep copy
    
    #Initialization -- Check if initial point is feasible. If not, enter the dummy variable
    while keepGoing:
#        print("current basis: " + str(basis))
        pivotRow = -1
        for i in range(len(gMatrix)):
            val = pari.substvec(gMatrix[i][2*numVar], xVar[0:-1], startingPoint)
#            print("RHS value " + str(i) + ": " + str(val))
            if val < 0.0:
                pivotRow = i
#                print(str(basis[pivotRow]) + " will exit the basis.")
                if basis[pivotRow] < numVar:
                    pivotCol = basis[pivotRow] + numVar
                else:
                    pivotCol = basis[pivotRow] - numVar
                break

        if pivotRow >= 0:
            #Diagonal Pivot Check
            val = pari.substvec(gMatrix[pivotRow][pivotCol], xVar[0:-1], startingPoint)
            if val < -epsilon:
                basis[pivotRow] = pivotCol
                gMatrix = matrixPivot(gMatrix, pivotRow, pivotCol, parallelPivots)
#                print("A diagonal pivot will be performed")
            elif val > epsilon:
                ExitWarning(logging, startingPoint)
                keepGoing = False
            else:
                #Exchange Pivot Check
                pivotFound = False
                for i in range(len(gMatrix)):
                    val = pari.substvec(gMatrix[i][pivotCol], xVar[0:-1], startingPoint)
                    val2 = pari.substvec(gMatrix[pivotRow][i + numVar], xVar[0:-1], startingPoint)
                    if val > 0.0 or val2 < 0.0:
                        if val*val2  >= 0:
                            ExitWarning(logging, startingPoint)
                            keepGoing = False
                        else:
                            pivotRow2 = i
#                            print(str(basis[pivotRow2]) + " will also exit the basis. An exchange pivot will be performed.")
                            if basis[pivotRow2] < numVar:
                                pivotCol2 = basis[pivotRow2] + numVar
                            else:
                                pivotCol2 = basis[pivotRow2] - numVar
                            pivotFound = True
                            basis[pivotRow] = pivotCol
                            basis[pivotRow2] = pivotCol2
                            gMatrix = matrixPivot(gMatrix, pivotRow, pivotCol2, parallelPivots)
                            gMatrix = matrixPivot(gMatrix, pivotRow2, pivotCol, parallelPivots)
                            temp = gMatrix[pivotRow]
                            gMatrix[pivotRow] = gMatrix[pivotRow2]
                            gMatrix[pivotRow2] = temp
                            break
                if not pivotFound:
                    #The instance is not feasible at the given starting point
                    feasible = False
                    keepGoing = False
                    gMatrix = originalGmatrix
        else:
            #A feasible solution and starting basis have been found
            keepGoing = False
#    print("-------------------------------------------------")
    return basis, gMatrix, feasible
