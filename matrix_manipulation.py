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

import multiprocessing

# Functions use for parallelization

# Create a rational Pari object as the division of two Pari objects
#
# Input:    a   --  a tuple containing:
#                   num --  the numerator
#                   den --  the denominator
#
# Output:   num/den
def Div(a):
    num, den = a
    return num/den

# Process a single element of a matrix upon which a principal pivoting operation
# is being performed
#
# Input:    a   --  a tuple containing:
#                   const   --  the original value of the element
#                   mult    --  the entry in the same column as const, but in 
#                               but in the pivot row
#                   num     --  the entry in the same row as const, but in the
#                               pivot column
#                   den     --  the entry in the pivot row, pivot column
#
# Output:   const - mult*num/den
def ProcessElement(a):
    const, mult, num, den = a
    return const - mult*num/den
    
# Process a row of the a matrix upon which a principal pivoting operation
# is being performed. Namely:
#   (1) For the pivot row   --  divide every entry by the value in the pivot
#                               column
#   (2) For other rows      --  from every entry, subtract the entry in the 
#                               pivot row, current column times the entry in the
#                               current row, pivot column divided by the entry 
#                               in the pivot row, pivot column
#
# Input:    q           --  the multiprocessing queue
#           row1        --  the pivot row
#           row2        --  the current row
#           row1index   --  the row index of the pivot row
#           row2index   --  the row index of the current row
#           pivotColumn --  the column index of the pivot column
def ProcessRow(q, row1, row2, row1index, row2index, pivotColumn):
    pool = multiprocessing.Pool(processes = None)
    if row1index == row2index:
        v = pool.map_async(Div, ((i, row1[pivotColumn]) for i in row1))
    else:
        v = pool.map_async(ProcessElement, ((row1[i], row2[i], row1[pivotColumn], row2[pivotColumn]) for i in range(len(row1))))
    q.put([row1index,v.get()])


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
def matrixPivot(M, i, j, parallel):
    temp = M[i][j]
    
    if not parallel:
        for k in range(len(M[i])):
            M[i][k] = M[i][k]/temp;
        for row in range(len(M)):
            if row != i:
                temp2 = M[row][j]
                for col in range(len(M[row])):
                    M[row][col] -= temp2*M[i][col]
    else:
        q = multiprocessing.Queue()
        processes = []
        order = []
        rows = []
        for k in range(len(M)):
            p = multiprocessing.Process(target=ProcessRow, args=(q, M[k], M[i], k, i, j))
            processes.append(p)
            p.start()
        for p in processes:
            ret = q.get() 
            order.append(ret[0])
            rows.append(ret[1])
        for p in processes:
            p.join()
        reorder = order.copy()
        for k in order:
            reorder[order[k]] = k
        M = [rows[k] for k in reorder]
        
    return(M)



