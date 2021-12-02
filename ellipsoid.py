###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        Approximately "solve" nonlinear programs using a modified
#                   ellipsoid method.
#
################################################################################

import numpy

# Define Functions

# Run the ellipsoid method. The method we use can be found in Chapter 9.7 of 
# Introduction to Operations Research by Ecker and Kupperschmid. However, we 
# modify this method by exiting early if we find a feasible point with strictly 
# negative objective value, as this is enough to guarantee proper functionality
# of the mpLCP method.
#
# Note that we assume the problem is a minimization subject to less than or 
# equal to constraints
#
# Input:    pari    --  the pari environment
#           obj     --  a pari object representing the objective function
#           in_con  --  a vector of pari objects representing the constraints
#           xVar    --  a vector of pari objects representing the variables
#           Q       --  the matrix defining the starting ellipse
#           x       --  the vector defining the starting point
#           tol     --  the termination tolerance
#           maxIter --  the maximum number of iterations
#
# Output    feas    --  a boolean indicating whether or not a feasible point was
#                       found with strictly negative objective value
#           x       --  a feasible point with strictly negative objective value,
#                       when one is found
def Ellipsoid(pari, obj, in_con, xVar, grads, Q, x, tol, maxIter):
    s = tol*2
    index = 0
    i = 0
    feas = False
    x = numpy.array(x)
    Q = numpy.array(Q)
    obj_grad = []
    n = len(Q)
    last_index = 0
    val = 0
    iter_w_no_improve = [0 for _ in range(len(in_con))]
    violation = [0 for _ in range(len(in_con))]
    
    for v in xVar:
        obj_grad.append(pari.deriv(obj,v))

    
    print("---------- Ellipsoid Method ------------------")
    
    print("Starting point:")
    print(x)
    
    # Find the maximum violation in the constraints and set x[-2] to this value. This attempts to find a feasible starting point.
    for l in range(len(in_con)):
        violation[l] = pari.substvec(in_con[l], xVar, x)
    
#    print("Inequality violations:")
#    print(violation)
    
    maxViol = max(violation)
    if maxViol > x[-2]:
        x[-2] = -1*maxViol
        Q[-2][-2] = x[-2]**2
    
    while s > tol and i <= maxIter and not feas:

        i += 1
        p = obj
        j = 0
        alpha = 1

#        print("******** Optimization Step, Iteration: " + str(i) + " *********")
#        print("Original Q:")
#        for row in Q:
#            print(row)
        for k in range(len(in_con)):
            j = (index + k) % len(in_con)
            val = pari.substvec(in_con[j], xVar, x)
            if val > tol:
#                print("Violated Row: " + str(j) + "\tValue: " + str(val))
                p = in_con[j]
                percent_change = (violation[j] - val)/val
#                print("Percent change from last time: ")
#                print(percent_change)
                if percent_change > tol:
                    iter_w_no_improve[j] = 0
                elif violation[j] > 0:
                    iter_w_no_improve[j] += 1
                violation[j] = val
                break
        index = j
        # exit if there is no improvement. It is unlikely that we can satisfy the equality and the inequalities at the same time
        if i > 1 and iter_w_no_improve[index] > 1:
#            print("Exiting due to no improvement")
            break
        if p == obj:
#            print("p is the objective")
            val = pari.substvec(p, xVar, x)
#            print(val)
            if val < 0.0:
                feas = True
        if not feas:
            if p == obj:
                grad = obj_grad
            else:
                grad = grads[index]
            grad_at_x = []
            for d in grad:
                grad_at_x.append(float(pari.substvec(d, xVar, x)))
#            print("Gradient:")
#            print(grad)
#            print("Gradient at x:")
#            print(grad_at_x)
            grad_at_x = numpy.multiply(1.0/numpy.linalg.norm(grad_at_x), grad_at_x)
#            print("Normalized Gradient at x:")
#            print(grad_at_x)
            
            
            num = numpy.matmul(Q,grad_at_x)
            rt = numpy.dot(grad_at_x, num)
            if rt < 0.0:
                break
            den = -1*(rt)**0.5
#            print("Denominator of d:")
#            print(den)
            direction = numpy.multiply(1.0/den, num)
#            print("Direction:")
#            print(direction)
            
            x = numpy.add(x, numpy.multiply(1/(n + 1), direction))
#            print("Updated x:")
#            print(x)
            
            Q = numpy.multiply(n**2/(n**2 - 1), numpy.subtract(Q, numpy.multiply(2/(n + 1), numpy.outer(direction, direction))))
#            print("Updated Q:")
#            for row in Q:
#                print(row)
            
#            print("Not finished writing equality constrained ellipsoid method!!")
#            feas = True
        last_index = index
        index = (index + 1) % len(in_con)
        
    print("-----------------------------------------------------------------------------")
    return x, pari.substvec(obj, xVar, x), feas
    
    
# Run a variant of the ellipsoid method that can handle equality constraints. 
# A basic description of the method we use can be found in:
#
#   Shah, Sharmila, John E. Mitchell, and Michael Kupferschmid. "An ellipsoid 
#   algorithm for equality-constrained nonlinear programs." Computers & 
#   Operations Research 28.1 (2001): 85-92.
#
# A more detailed description of the method can be found in the first author's
# PhD dissertation. Again, we modify this method by exiting early if we find a 
# feasible point with strictly negative objective value, as this is enough to 
# guarantee proper functionality of the mpLCP method.
#
# Note that we assume the problem is a minimization subject to ONE equality 
# constraint and one or more less than or equal to constraints
#
# Input:    pari    --  the pari environment
#           obj     --  a pari object representing the objective function
#           eq_con  --  a pari object representing the equality constraint
#           in_con  --  a vector of pari objects representing the inequality 
#                       constraints
#           xVar    --  a vector of pari objects representing the variables
#           Q       --  the matrix defining the starting ellipse
#           x       --  the vector defining the starting point
#           tol     --  the termination tolerance
#           maxIter --  the maximum number of iterations
#
# Output    feas    --  a boolean indicating whether or not a feasible point was
#                       found with strictly negative objective value
#           x       --  a feasible point with strictly negative objective value,
#                       when one is found
def EllipsoidEq(pari, obj, eq_con, in_con, xVar, grads, Q, x, tol, maxIter):
    s = tol*2
    index = 0
    i = 0
    feas = False
    x = numpy.array(x)
    Q = numpy.array(Q)
    eq_grad = grads[-1]
    obj_grad = []
    n = len(Q)
    last_index = 0
    val = 0
    iter_w_no_improve = [0 for _ in range(len(in_con))]
    violation = [0 for _ in range(len(in_con))]
    
    for v in xVar:
        obj_grad.append(pari.deriv(obj,v))

    
    print("---------- Ellipsoid Method With A Single Equality Constraint ------------------")
    
    print("Starting point:")
    print(x)
    
    # Perform initial projection of x to the flat of the equality constraint.
    eq_at_x = float(pari.substvec(eq_con, xVar, x))
    while eq_at_x > tol or eq_at_x < -1*tol:
        #Project x to the flat of the constraint.
        eq_grad_at_x = []
        for d in eq_grad:
            eq_grad_at_x.append(float(pari.substvec(d, xVar, x)))
        eq_grad_at_x = numpy.array(eq_grad_at_x)
#            print("A:")
#            print(eq_grad_at_x)
        # eq_grad_at_x is matrix A as in above referenced dissertation
        # similarly, eq_at_x is v
        alpha = -1*eq_at_x/(numpy.dot(eq_grad_at_x, eq_grad_at_x))
        x = x + numpy.multiply(alpha, eq_grad_at_x)
        eq_at_x = float(pari.substvec(eq_con, xVar, x))
        
    
#    print("Point after initial projection:")
#    print(x)
    
    # Find the maximum violation in the constraints and set x[-2] to this value. This attempts to find a feasible starting point.
    for l in range(len(in_con)):
        violation[l] = pari.substvec(in_con[l], xVar, x)
    
#    print("Inequality violations:")
#    print(violation)
    
    maxViol = max(violation)
    if maxViol > x[-2]:
        x[-2] = -1*maxViol
        Q[-2][-2] = x[-2]**2
    
    while s > tol and i <= maxIter and not feas:
#        print("Iteration: " + str(i))
#        last_index = index
        i += 1
        p = obj
        j = 0
        alpha = 1
        eq_at_x = float(pari.substvec(eq_con, xVar, x))
#        print("******** Projecting x *********")
#        print("x before projection:")
#        print(x)
#        print("Equality evaluated at x:")
#        print(eq_at_x)
        if i > 0:
            while eq_at_x > tol or eq_at_x < -1*tol:
                #Project x to the flat of the constraint.
                eq_grad_at_x = []
                for d in eq_grad:
                    eq_grad_at_x.append(float(pari.substvec(d, xVar, x)))
                eq_grad_at_x = numpy.array(eq_grad_at_x)
    #            print("A:")
    #            print(eq_grad_at_x)
                # eq_grad_at_x is matrix A as in above referenced dissertation
                # similarly, eq_at_x is v
                alpha = -1*eq_at_x/(numpy.dot(eq_grad_at_x, eq_grad_at_x))
                x = x + numpy.multiply(alpha, eq_grad_at_x)
                eq_at_x = float(pari.substvec(eq_con, xVar, x))
            
#        print("Equality evaluated at x:")
#        print(eq_at_x)
        
#        print("******** Optimization Step, Iteration: " + str(i) + " *********")
#        print("Original Q:")
#        for row in Q:
#            print(row)
        for k in range(len(in_con)):
            j = (index + k) % len(in_con)
            val = pari.substvec(in_con[j], xVar, x)
            if val > tol:
#                print("Violated Row: " + str(j) + "\tValue: " + str(val))
                p = in_con[j]
                percent_change = (violation[j] - val)/val
#                print("Percent change from last time: ")
#                print(percent_change)
                if percent_change > tol:
                    iter_w_no_improve[j] = 0
                elif violation[j] > 0:
                    iter_w_no_improve[j] += 1
                violation[j] = val
                break
        index = j
        # exit if there is no improvement. It is unlikely that we can satisfy the equality and the inequalities at the same time
        if i > 1 and iter_w_no_improve[index] > 1:
#            print("Exiting due to no improvement")
            break
        if p == obj:
#            print("p is the objective")
            val = pari.substvec(p, xVar, x)
#            print(val)
            if val < 0.0:
                feas = True
        if not feas:
            if p == obj:
                grad = obj_grad
            else:
                grad = grads[index]
            grad_at_x = []
            for d in grad:
                grad_at_x.append(float(pari.substvec(d, xVar, x)))
#            print("Gradient:")
#            print(grad)
#            print("Gradient at x:")
#            print(grad_at_x)
            grad_at_x = numpy.multiply(1.0/numpy.linalg.norm(grad_at_x), grad_at_x)
#            print("Normalized Gradient at x:")
#            print(grad_at_x)
            
            
            qat = numpy.matmul(Q,eq_grad_at_x)
#            print("QAt:")
#            print(qat)
            aq = numpy.matmul(eq_grad_at_x, Q)
#            print("AQ:")
#            print(aq)
            aqat = numpy.matmul(aq, eq_grad_at_x)
#            print("AQAt:")
#            print(aqat)
            num = numpy.dot(numpy.subtract(Q, numpy.outer(numpy.multiply(1.0/aqat, qat), aq)), grad_at_x)
#            print("Numerator of d:")
#            print(num)
            rt = numpy.dot(grad_at_x, num)
            if rt < 0.0:
                break
            den = -1*(rt)**0.5
#            print("Denominator of d:")
#            print(den)
            direction = numpy.multiply(1.0/den, num)
#            print("Direction:")
#            print(direction)
            
            x = numpy.add(x, numpy.multiply(1/(n + 1), direction))
#            print("Updated x:")
#            print(x)
            
            Q = numpy.multiply(n**2/(n**2 - 1), numpy.subtract(Q, numpy.multiply(2/(n + 1), numpy.outer(direction, direction))))
#            print("Updated Q:")
#            for row in Q:
#                print(row)
            
#            print("Not finished writing equality constrained ellipsoid method!!")
#            feas = True
        last_index = index
        index = (index + 1) % len(in_con)
        
    print("-----------------------------------------------------------------------------")
    return x, pari.substvec(obj, xVar, x), feas
