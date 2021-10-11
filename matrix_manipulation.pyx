###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        Perform matrix manipulations -- primarily principal pivots.
#                   Original intention is for use as part of a solver for multi-
#                   parametric Linear Complementarity Problems (mpLCP's).
#
################################################################################


from cython.parallel import prange, threadid

# Define Functions

# Perform a pricipal pivot on the given matrix, i.e., perform elementary row 
# operations so that column j of M is an identity vector with a 1 in the i-th 
# position.
#
# Input:    M   --  the matrix to manipulate
#           i   --  the row index
#           j   --  the column index
#
# Output:   M   --  the updated matrix
cpdef object[:,:] matrixPivot(object[:,:] M, int i, int j):
#    oldM = M
#    t = time.time()
#    temp = []
#    temp.append(M[i][j])
    cdef object temp = M[i][j]
#    for a in M[i]:
#        a = a/temp;
#        print(a);
#    for row in range(len(M)):
#        if row != i:
#            temp2 = M[row][j]
#            for col in range(len(M[row])):
#                M[row][col] -= temp2*M[i][col]
#                print(M[row][col])
##    for row in M:
##        print(row)
#    print(time.time() - t)

    for a in prange(len(M[i]), nogil = True):
        M[i][a] = M[i][a]/temp;
        print(M[i][a]);
        print("thread: " + str(threadid()))
    for row in prange(len(M)):
        if row != i:
            temp2 = M[row][j]
            print("thread: " + str(threadid()))
            for col in prange(len(M[row])):
                M[row][col] -= temp2*M[i][col]
                print(M[row][col])
                print("thread: " + str(threadid()))

#    print(M[i])
#    print(Div2(M[i],temp))
#    random_vector = Parallel(n_jobs=-1, backend='threading')(delayed(Div)(k, l) for k in M[i] for l in [3.])
#    print(random_vector)


#    for row in M:
#        print(row)




