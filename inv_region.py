###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        Define the class to be associated with an invariancy region.
#
################################################################################

from matrix_manipulation import *
from create_pyomo_models import *
from ellipsoid import *
#from read_problem import printMatrix
import time
import copy

# Define Wrapper Functions

# Process an Invariancy Region.
#
# Input:    rgn --  an instance of the InvRgn class to process
#           q   --  a queue containing regions that still need processing
#           discoveredBases --  a list of bases that have already been processed
#                               or already discovered and shown to not need
#                               processing
#           finalPartition  --  a queue containing the regions which form the 
#                               final partition of the parameter space
#
def ProcessRegion(rgn, q, discoveredBases, finalPartition, re, pari, pyo, paramSpace, parallelPivots):
    print("Processing region associated with basis: " + str(rgn.Basis()))
    rgn.Process(q, discoveredBases, finalPartition, re, pari, pyo, paramSpace, parallelPivots)



# Define the Invariancy Region Class

class InvRgn:
    def __init__(   self, 
                    re, 
                    pari, 
                    pyo, 
                    gMatrix, 
                    basis, 
                    xVar, 
                    startingPoint, 
                    epsilon, 
                    nonlinearSolver, 
                    Q0, 
                    midpoint, 
                    ellipsoidTol, 
                    ellipsoidIter, 
                    paramSpace, 
                    nlpsAsFeasProbs, 
                    checkFwithEH,
                    checkDimWithF):
        self.solver     = nonlinearSolver
        self.grads      = [ [] for _ in range(len(basis) + len(paramSpace)) ]
        self.defIneq    = [ pari.zero() for _ in range(len(basis) + len(paramSpace)) ]
        self.eps        = epsilon
        self.tableau    = gMatrix
        self.basis      = basis
        self.xVar       = xVar
        self.startPnt   = startingPoint
        self.rhs        = []
        self.point      = []
        self.f          = []
        self.ellipseQ   = Q0
        self.ellipsePnt = midpoint
        self.ellipseTol = ellipsoidTol
        self.ellipseIt  = ellipsoidIter
        self.nlpFeas    = nlpsAsFeasProbs
        self.fWithEH    = checkFwithEH
        self.dWithF     = checkDimWithF
        self.dim        = 0
        if self.solver == "":
            self.GetIneqAndGradients(pari, paramSpace, True)
        else:
            self.GetIneqAndGradients(pari, paramSpace, False)
        self.BuildZ(pari)
        self.BuildEandH(re, pari, pyo)

    # Processing the invariancy region is finished. Only store the basis and the
    # defining constraints. No need to continue storing the entire tableau
    def Finalize(self):
        for el in self.tableau:
            rhs.append(el.Str())
        del self.tableau
        
    
    # Getters
    def Tableau(self):
        return self.tableau
        
    def Basis(self):
        return self.basis
    
    def RHS(self):
        return self.rhs
        
    def Dim(self, re, pari, pyo):
        if self.dim == 0:
            pt, val, feas = self.IsFullDim(re, pari, pyo)
            if feas and val < -self.eps:
                self.dim = len(self.xVar) - 2
            else:
                self.dim = len(self.xVar) - 3
        return self.dim
        
    # Other Methods
    
    # Solve NlpD wrapper
    def IsFullDim(self, re, pari, pyo):
        retVal = 0
        retPnt = []
        retFeas = False
        if self.solver != "":
            retPnt, retVal, retFeas = self.SolveNlpD_pyo(re, pari, pyo)
        else:
            retPnt, retVal, retFeas = self.SolveNlpD_ell(pari)
        self.point = retPnt
        return retPnt, retVal, retFeas
    
    # Solve NlpD using pyomo
    def SolveNlpD_pyo(self, re, pari, pyo):
        feasible = False
        constraint_expressions = []
        for k in range(len(self.defIneq)):
            if k < len(self.basis) and not self.z[k]:
                constraint_expressions.append(pariToPyomo(re, str(self.defIneq[k])) + " + model.x[" + str(len(self.xVar) - 1) + "] <= 0")
            elif k >= len(self.basis):
                constraint_expressions.append(pariToPyomo(re, str(self.defIneq[k])) + " <= 0")
                
        print("Constraint Expressions for NLP D:")
        print(constraint_expressions)
        
        if not self.nlpFeas:
            model, opt = CreatePyomoModel(pyo, len(self.xVar), "-model.x[" + str(len(self.xVar) - 1) + "]", pyo.minimize, constraint_expressions, self.solver)
        else:
            model, opt = CreatePyomoModel(pyo, len(self.xVar), "0", pyo.minimize, constraint_expressions, self.solver)
            model.x[len(self.xVar) - 1].setlb(self.eps)
        results = opt.solve(model)
        
        print(results.solver.termination_condition)
        
        if results.solver.termination_condition == "optimal" or results.solver.termination_condition == "locallyOptimal":
            if not self.nlpFeas and pyo.value(model.o) < 0.0:
                feasible = True
            elif self.nlpFeas:
                feasible = True
        
        x = [pyo.value(model.x[k+1]) for k in range(len(model.x))]
        print("Solution:")
        print(x)
        
        if self.nlpFeas:
            return x, -1*pyo.value(model.x[len(self.xVar) - 1]), feasible
        else:
            return x, pyo.value(model.o), feasible
        
        
    # Solve NlpD using the ellipsoid method
    def SolveNlpD_ell(self, pari):
        grads = []
        ineq_constraints = []
        for k in range(len(self.defIneq)):
            if k < len(self.basis) and not self.z[k]:
                grads.append(copy.deepcopy(self.grads[k]))
                ineq_constraints.append(self.defIneq[k] + self.xVar[-2])
            elif k >= len(self.basis):
                grads.append(self.grads[k])
                ineq_constraints.append(self.defIneq[k])
        optPoint, objVal, feasible = Ellipsoid(pari, -1*self.xVar[-2], ineq_constraints, self.xVar, grads, self.ellipseQ, self.ellipsePnt, self.ellipseTol, self.ellipseIt)
        print("Optimal point:")
        print(optPoint)
        return optPoint, objVal, feasible
        
    # Process the region
    def Process(self, q, discoveredBases, finalPartition, re, pari, pyo, paramSpace, parallelPivots):
        self.BuildF(re, pari, pyo)
        for i in range(len(self.basis)):
            if self.f[i]:
                newRegions, newBases = self.GetAdjacentRegionsAcross(i, discoveredBases, re, pari, pyo, paramSpace, parallelPivots)
        print("Finish writing code for region processing")
        
    # Find all full dimensional regions across the given 'facet'
    def GetAdjacentRegionsAcross(self, i, discoveredBases, re, pari, pyo, paramSpace, parallelPivots):
        print("Looking for adjacent regions across: " + str(i))
        # line 1
        newRegions = []
        newBases = []
        numVar = len(self.basis)
        newBasis = copy.deepcopy(self.basis)
        iComp = i
        if self.basis[i] < numVar:
            iComp = self.basis[i] + numVar
        newBasis[i] = iComp
        
        #line 2
        if newBasis not in discoveredBases:
            #line 3
            print(self.tableau[i][iComp])
            print(pari.poldegree(self.tableau[i][iComp]))
            if pari.poldegree(self.tableau[i][iComp]) >= 0: # pari degree is negative inf for 0, 0 for constants, and >0 otherwise
                #generate new region using a diagonal pivot
                newRegion = InvRgn(re, 
                                   pari, 
                                   pyo, 
                                   matrixPivot(self.tableau, i, iComp, parallelPivots), 
                                   newBasis, 
                                   self.xVar, 
                                   self.startPnt, 
                                   self.eps, 
                                   self.solver, 
                                   self.ellipseQ, 
                                   self.ellipsePnt,
                                   self.ellipseTol, 
                                   self.ellipseIt, 
                                   paramSpace, 
                                   self.nlpFeas, 
                                   self.fWithEH, 
                                   self.dWithF)
        #printMatrix(newRegion.Tableau())
        print('\n')
        print(self.basis)
        print(newBasis)
        print("Finish writing code for getting adjacent regions")
        return newRegions, newBases
        
    # Determine the basis elements whose associated RHS's in the tableau form
    # (k-1)-dimensional boundaries of the invariancy region
    def BuildF(self, re, pari, pyo):
        for i in range(len(self.basis)):
            print("Checking if constraint " + str(i) + " is a 'facet'.")
            if not self.z[i] and not self.e[i] and not self.f[i]:
                print("Entering if for " + str(i))
                optPoint, objVal, feasible = self.SolveNlpF(re, pari, pyo, i)
        
    # Build the set Z of all RHS entries that are identically zero.
    def BuildZ(self, pari):
        self.z = [False for i in range(len(self.basis))]
        for i in range(len(self.z)):
            if pari.poldegree(self.tableau[i][-1]) < 0: #in pari, the degree of zero is negative infinity. The degree of other constants is 0.
                self.z[i] = True;
        
    # Build the sets:
    #   (1) E   --  all basis elements whose associated RHS's in the 
    #               tableau do not intersect the invariancy region
    #   (2) H[i]--  for each "i" the basis, other elements in the basis whose 
    #               associated RHS's in the tableau contain the RHS associated 
    #               with "i". Generally, this indicates that the RHS associated
    #               with each "j" in H[i] is a (polynomial) multiple of the RHS
    #               associated with "i".
    def BuildEandH(self, re, pari, pyo):
        t = time.time()
        self.e = [False for i in range(len(self.basis))]
        self.f = [False for i in range(len(self.basis))]
        self.h = [ [False for i in range(len(self.basis))] for i in range(len(self.basis)) ]
        for i in range(len(self.basis)):
            print("Processing Basis Element: " + str(self.basis[i]))
            if not self.z[i] and not self.e[i]:
                for j in range(len(self.basis)):
                    if j != i and not self.z[j] and not self.e[j] and not self.h[i][j]:
                        print("\tComparing to Basis Element: " + str(self.basis[j]))
                        optPoint, objVal, feasible = self.SolveNlpH(re, pari, pyo, i, j)
                        if feasible and objVal > -self.eps and objVal < self.eps:
                            self.h[i][j] = True
                            for k in range(len(self.h[j])):
                                if self.h[j][k]:
                                    self.h[i][k] = True
                        elif not feasible or (feasible and objVal > self.eps):
                            self.e[i] = True
                            break
                        elif self.fWithEH and feasible and objVal < -self.eps:
                            strict = True
                            print("Constraint Violations:")
                            for l in range(len(self.defIneq)):
                                if l != i:
                                    if pari.substvec(self.defIneq[l], self.xVar, optPoint) > self.eps:
                                        strict = False
                                        break
                            if strict:
                                self.f[i] = True
            if self.fWithEH and self.dWithF and self.dim < len(self.xVar) - 2 and sum(self.h[i]) == 0 and self.f[i]:
                self.dim = len(self.xVar) - 2
        print("Gradients After Build E and H:")
        print(self.grads)
        print("E")
        print(self.e)
        print("H")
        for row in self.h:
            print(row)
        print("F")
        print(self.f)
        print("dim")
        print(self.dim)
        print("Time:")
        print(time.time() - t)
        
            
    # Solve NlpH wrapper
    def SolveNlpH(self, re, pari, pyo, i, j):
        retVal = 0
        retPnt = []
        retFeas = False
        if self.solver != "":
            retPnt, retVal, retFeas = self.SolveNlpH_pyo(re, pari, pyo, i, j)
        else:
            retPnt, retVal, retFeas = self.SolveNlpH_ell(pari, i, j)
        return retPnt, retVal, retFeas
        
    # Solve NlpH using Pyomo
    def SolveNlpH_pyo(self, re, pari, pyo, i, j):
        feasible = False
        constraint_expressions = []
        constraint_expressions.append(pariToPyomo(re, str(self.defIneq[i])) + " == 0")
        for k in range(len(self.defIneq)):
            if k != i and (k >= len(self.basis) or not self.z[k]):
                if k == j:
                    constraint_expressions.append(pariToPyomo(re, str(self.defIneq[k])) + " + model.x[" + str(len(self.xVar) - 1) + "] <= 0")
                else:
                    constraint_expressions.append(pariToPyomo(re, str(self.defIneq[k])) + " <= 0")
        if not self.nlpFeas:
            model, opt = CreatePyomoModel(pyo, len(self.xVar), "-model.x[" + str(len(self.xVar) - 1) + "]", pyo.minimize, constraint_expressions, self.solver)
        else:
            model, opt = CreatePyomoModel(pyo, len(self.xVar), "0", pyo.minimize, constraint_expressions, self.solver)
            model.x[len(self.xVar) - 1].setlb(self.eps)
        results = opt.solve(model)
        
        print(results.solver.termination_condition)
        
        if results.solver.termination_condition == "optimal" or results.solver.termination_condition == "locallyOptimal":
            if not self.nlpFeas and pyo.value(model.o) < 0.0:
                feasible = True
            elif self.nlpFeas:
                feasible = True
        
        x = [pyo.value(model.x[k+1]) for k in range(len(model.x))]
        print("Solution:")
        print(x)
        
        if self.nlpFeas:
            return x, -1*pyo.value(model.x[len(self.xVar) - 1]), feasible
        else:
            return x, pyo.value(model.o), feasible
        
        
    # Solve NlpH using the modified ellipsoid method
    def SolveNlpH_ell(self, pari, i, j):
        grads = []
        ineq_constraints = []
        eq_constraint = self.defIneq[i]
        for k in range(len(self.defIneq)):
            if k != i and (k >= len(self.basis) or not self.z[k]):
                if k == j:
                    grads.append(copy.deepcopy(self.grads[k]))
                    ineq_constraints.append(self.defIneq[k] + self.xVar[-2])
                    grads[-1][-2] += 1
                    print("J-th Gradient:")
                    print(grads[-1])
                    print("Row of J-th Gradient:")
                    print(str(len(grads) - 1))
                else:
                    grads.append(self.grads[k])
                    ineq_constraints.append(self.defIneq[k])
        grads.append(self.grads[i])
        optPoint, objVal, feasible = EllipsoidEq(pari, -1*self.xVar[-2], eq_constraint, ineq_constraints, self.xVar, grads, self.ellipseQ, self.ellipsePnt, self.ellipseTol, self.ellipseIt)
        print("Optimal point:")
        print(optPoint)
        return optPoint, objVal, feasible
        
    # Solve NlpF wrapper
    def SolveNlpF(self, re, pari, pyo, i):
        retVal = 0
        retPnt = []
        retFeas = False
        if self.solver != "":
            retPnt, retVal, retFeas = self.SolveNlpF_pyo(re, pari, pyo, i)
        else:
            retPnt, retVal, retFeas = self.SolveNlpF_ell(pari, i)
        return retPnt, retVal, retFeas
        
    # Solve NlpF using Pyomo
    def SolveNlpF_pyo(self, re, pari, pyo, i):
        feasible = False
        constraint_expressions = []
        constraint_expressions.append(pariToPyomo(re, str(self.defIneq[i])) + " == 0")
        for k in range(len(self.defIneq)):
            if k != i and (k >= len(self.basis) or (not self.z[k] and not self.h[i][k])):
                constraint_expressions.append(pariToPyomo(re, str(self.defIneq[k])) + " + model.x[" + str(len(self.xVar) - 1) + "] <= 0")
        if not self.nlpFeas:
            model, opt = CreatePyomoModel(pyo, len(self.xVar), "-model.x[" + str(len(self.xVar) - 1) + "]", pyo.minimize, constraint_expressions, self.solver)
        else:
            model, opt = CreatePyomoModel(pyo, len(self.xVar), "0", pyo.minimize, constraint_expressions, self.solver)
            model.x[len(self.xVar) - 1].setlb(self.eps)
        results = opt.solve(model)
        
        print(results.solver.termination_condition)
        
        if results.solver.termination_condition == "optimal" or results.solver.termination_condition == "locallyOptimal":
            if not self.nlpFeas and pyo.value(model.o) < 0.0:
                feasible = True
            elif self.nlpFeas:
                feasible = True
        
        x = [pyo.value(model.x[k+1]) for k in range(len(model.x))]
        print("Solution:")
        print(x)
        
        if self.nlpFeas:
            return x, -1*pyo.value(model.x[len(self.xVar) - 1]), feasible
        else:
            return x, pyo.value(model.o), feasible
            
    # Solve NlpH using the modified ellipsoid method
    def SolveNlpF_ell(self, pari, i):
        grads = []
        ineq_constraints = []
        eq_constraint = self.defIneq[i]
        for k in range(len(self.defIneq)):
            if k != i and (k >= len(self.basis) or not (not self.z[k] and not self.h[i][k])):
                grads.append(copy.deepcopy(self.grads[k]))
                ineq_constraints.append(self.defIneq[k] + self.xVar[-2])
                grads[-1][-2] += 1
                print("J-th Gradient:")
                print(grads[-1])
                print("Row of J-th Gradient:")
                print(str(len(grads) - 1))
        grads.append(self.grads[i])
        optPoint, objVal, feasible = EllipsoidEq(pari, -1*self.xVar[-2], eq_constraint, ineq_constraints, self.xVar, grads, self.ellipseQ, self.ellipsePnt, self.ellipseTol, self.ellipseIt)
        print("Optimal point:")
        print(optPoint)
        return optPoint, objVal, feasible
    
    # Use pari to compute the defining inequalities of the invariancy region 
    # (stored in less-than-or-equal-to form). Then, compute the gradient of each
    # These are used for the ellipsoid method.
    def GetIneqAndGradients(self, pari, paramSpace, storeGrads):
        for i in range(len(self.basis)):
            val = pari.substvec(pari.denominator(self.tableau[i][-1]), self.xVar[0:-1], self.startPnt)
            if val > 0.0:
                self.defIneq[i] = -1*pari.numerator(self.tableau[i][-1])
            else:
                self.defIneq[i] = pari.numerator(self.tableau[i][-1])
            if storeGrads:
                for v in self.xVar:
                    self.grads[i].append(pari.deriv(self.defIneq[i],v))
                print("Gradient of Ineq " + str(i))
                print(self.grads[i])
        for i in range(len(paramSpace)):
            self.defIneq[i + len(self.basis)] = paramSpace[i][0] - paramSpace[i][1]
            if storeGrads:
                for v in self.xVar:
                    self.grads[i + len(self.basis)].append(pari.deriv(self.defIneq[i + len(self.basis)],v))
                print("Gradient of Ineq " + str(i + len(self.basis)))
                print(self.grads[i + len(self.basis)])
