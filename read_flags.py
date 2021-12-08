###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        Read in and set parameters passed from the command line in 
#                   order to solve an instance of the multiparametric Linear 
#                   Complementarity Problem (mpLCP).
#
################################################################################


# Define Functions

# Print a warning message if a commandline flag is passed with an unrecognized
# value
#
# Input:    flag -- the command line flag
#           message --  the warning message to print
#           logging --  the logging environment
def PrintInvalidParameterMessage(flag, curValue, validValues, logging):
    logging.warning("Invalid value for flag '" + str(flag) + "', ignoring and leaving at default value of " + str(curValue) + ".For future reference, valid values are " + str(validValues) + ".")

# Parse the input flags
#
# Input:    sys --  the variable containing any information passed at the
#                   command line
#           logging --  the logging environment
#           linearSolver    --  a string indicating the solver to be used by
#                               pyomo when solving LP or MIP problems
#           nonlinearSolver --  a string indicating the solver to be used by
#                               pyomo when solving NLP or MINLIP problems
#           parallelPivots  --  a boolean indicating whether or not principal
#                               pivot operations carried out on the problem's
#                               tableau should be computed in parallel
#           useCrissCross   --  a boolean indicating whether or not the Criss
#                               Cross method for LCP should be used in attempt
#                               to find a starting basis
#           nlpsAsFeasProbs --  a boolean indicating whether or not NLP problems
#                               should be posed as feasibility problems rather
#                               than optimization problems. Specifically, the
#                               NLP's solved herein do not need solved to
#                               optimality in order for the mpLCP method to
#                               function correctly. Only a feasible point with
#                               strictly positive (or negative, depending on
#                               the situation) objective value is needed. So,
#                               posing the problem as a feasibility problem by
#                               setting the objective as a constant and adding
#                               a constraint forcing the appropriate strict sign
#                               of the original objective may result in quicker
#                               run times.
#           checkFwithEH    --  a boolean indicating whether or not to check
#                               if any of the defining constraints of a region
#                               can be shown to form (k-1)-dimensional 
#                               boundaries of the region while building the 
#                               special sets E and H.
#           checkDimWithF   --  a boolean indicating whether or not to check
#                               the dimension of an invariancy region when 
#                               building the set F.
#           useAprime       --  a boolean indicating whether or not to use a
#                               "simplified" version of NLP_A when looking for 
#                               adjacent regions using exchange pivots
#           checkSwithEHF   --  a boolean indicating whether or not to check for
#                               the possibility of skipping NLP_S during phase 1
#                               if an appropriate point is found when building
#                               sets E, H, and F
#
# Outputs:  linearSolver
#           nonlinearSolver
#           parallelPivots
#           useCrissCross
#           nlpsAsFeasProbs
#           checkFwithEH
#           checkDimWithF
#           useAprime
#           checkSwithEHF
def ReadFlags(sys, logging, linearSolver, nonlinearSolver, parallelPivots, useCrissCross, nlpsAsFeasProbs, checkFwithEH, checkDimWithF, useAprime, checkSwithEHF):
    # Read the flags
    
    i = 2
    while i < len(sys.argv):
        if sys.argv[i][0] != '-':
            print("Invalid Command Line Argument (missing '-'). Exiting.\n")
            sys.exit("Got " + str(sys.argv[i]))
        else:
            if sys.argv[i] == "-linearSolver":
                i += 1
                linearSolver = sys.argv[i]
                print("Note: Setting linear solver as '" + str(linearSolver) + "'. It is the user's responsibility to ensure that appropriate software is installed and is 'visible' to pyomo. Available solvers and appropriate strings for calling them from pyomo can be found by typing 'pyomo help --solvers' in a Linux terminal.")
            elif sys.argv[i] == "-nonlinearSolver":
                i += 1
                nonlinearSolver = sys.argv[i]
                print("Note: Setting nonlinear solver as '" + str(nonlinearSolver) + "'. It is the user's responsibility to ensure that appropriate software is installed and is 'visible' to pyomo. Available solvers and appropriate strings for calling them from pyomo can be found by typing 'pyomo help --solvers' in a Linux terminal.")
            elif sys.argv[i] == "-parPivot":
                i += 1
                if sys.argv[i].upper() == "T":
                    parallelPivots = True
                elif sys.argv[i].upper() == "F":
                    parallelPivots = False
                else:
                    PrintInvalidParameterMessage("-parPivot", parallelPivot, "T and F", logging);
            elif sys.argv[i] == "-useCrissCross":
                i += 1
                if sys.argv[i].upper() == "T":
                    useCrissCross = True
                elif sys.argv[i].upper() == "F":
                    useCrissCross = False
                else:
                    PrintInvalidParameterMessage("-useCrissCross", useCrissCross, "T and F", logging);
            elif sys.argv[i] == "-nlpFeas":
                i += 1
                if sys.argv[i].upper() == "T":
                    nlpsAsFeasProbs = True
                elif sys.argv[i].upper() == "F":
                    nlpsAsFeasProbs = False
                else:
                    PrintInvalidParameterMessage("-nlpFeas", nlpsAsFeasProbs, "T and F", logging);
            elif sys.argv[i] == "-checkFwithEH":
                i += 1
                if sys.argv[i].upper() == "T":
                    checkFwithEH = True
                elif sys.argv[i].upper() == "F":
                    checkFwithEH = False
                else:
                    PrintInvalidParameterMessage("-checkFwithEH", checkFwithEH, "T and F", logging);
            elif sys.argv[i] == "-checkDimWithF":
                i += 1
                if sys.argv[i].upper() == "T":
                    checkDimWithF = True
                elif sys.argv[i].upper() == "F":
                    checkDimWithF = False
                else:
                    PrintInvalidParameterMessage("-checkDimWithF", checkDimWithF, "T and F", logging);
            elif sys.argv[i] == "-useAp":
                i += 1
                if sys.argv[i].upper() == "T":
                    useAprime = True
                elif sys.argv[i].upper() == "F":
                    useAprime = False
                else:
                    PrintInvalidParameterMessage("-useAp", checkDimWithF, "T and F", logging);
            elif sys.argv[i] == "-checkSwithEHF":
                i += 1
                if sys.argv[i].upper() == "T":
                    checkSwithEHF = True
                elif sys.argv[i].upper() == "F":
                    checkSwithEHF = False
                else:
                    PrintInvalidParameterMessage("-checkSwithEHF", checkDimWithF, "T and F", logging);
        i += 1
    
    return linearSolver, nonlinearSolver, parallelPivots, useCrissCross, nlpsAsFeasProbs, checkFwithEH, checkDimWithF, useAprime, checkSwithEHF
