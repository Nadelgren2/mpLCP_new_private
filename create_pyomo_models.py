###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        Create optimization models using Pyomo. Original intention
#                   is to use this code as part of a larger software package
#                   designed to solve multiparametric Linear Complementarity
#                   Problems (mpLCP's).
#
################################################################################


# Define Functions

# Convert a Pari string to a Pyomo string
#
# Input:    re  --  the re environment (for parsing regular expressions)
#           s   --  the Pari string
#
# Output:   sNew    --  the Pyomo string
def pariToPyomo(re, s):
    sNew = re.sub('(._)(\d+)', r'model.x[\2]', s)
    sNew2 = re.sub('\^', r'**', sNew)
    return sNew2

# Build the pyomo model
#
# Input:    pyo --  the pyomo environment
#           nv  --  the number of variables in the current model
#           obj --  the expression representing the objective function
#           sen --  the sense of the objective function (max or min)
#           con --  the vector of expressions representing the constraints
#           sol --  the choice of optimization solver
#
# Output    model   --  the optimization model
#           opt     --  the pyomo object containing the optimization options
def CreatePyomoModel(pyo, nv, obj, sen, con, sol):
    model = pyo.ConcreteModel()
    model.x = pyo.Var(pyo.RangeSet(nv), initialize = 0.0)
    model.o = pyo.Objective(expr = eval(obj), sense = sen)
    model.c = pyo.ConstraintList()
    for exp in con:
        model.c.add(expr = eval(exp))
    opt = pyo.SolverFactory(sol)

    return model, opt

