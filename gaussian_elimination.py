# @source: https://github.com/lihiSabag/Numerical-Analysis-2023.git

import numpy as np
from colors import bcolors

def gaussianElimination(mat):
    print("\033[94m" + "Gaussian Elimination" + "\033[0m")
    N = len(mat)
    singular_flag = forward_substitution(mat)
    if singular_flag != -1:

        if mat[singular_flag][N]:
            return "Singular Matrix (Inconsistent System)"
        else:
            return "Singular Matrix (May have infinitely many solutions)"

# if matrix is non-singular:
    forward_substitution_to_diagonal(mat)
    print(np.array(mat))
    # get solution to system using backward substitution
    return backward_substitution(mat)


# The function receives an upper triangular matrix and returns a fully ranked matrix
def forward_substitution(mat):
    N = len(mat)
    for k in range(N):
        pivot_row = k
        v_max = abs(mat[k][k])  # Setting the maximum value to the diagonal element itself
        for i in range(k + 1, N):
            if abs(mat[i][k]) > v_max:
                v_max = abs(mat[i][k])
                pivot_row = i

        if not mat[pivot_row][k]:  # Checking if the diagonal element is zero
            return k  # Matrix is singular

        # Swap the current row with the pivot row
        if pivot_row != k:
            # Swap entire rows, including the augmented column
            SaveRowi = mat[k].copy()
            mat[k] = mat[pivot_row]
            mat[pivot_row] = SaveRowi
        # End Partial Pivoting
        for i in range(k + 1, N):
            m = (mat[i][k] / mat[k][k])
            for j in range(k + 1, N + 1):
                mat[i][j] -= (mat[k][j] * m)
                if abs(mat[i][j]) < 1e-10:  # Small values are treated as zeros
                    mat[i][j] = 0

            mat[i][k] = 0  # Ensure lower triangular elements are zeroed out
    return -1

# function to calculate the values of the unknowns
def forward_substitution_to_diagonal(mat):
    N = len(mat)
    for k in range(N - 1, -1, -1):
        scalar = mat[k][k]
        for j in range(N + 1):
            mat[k][j] /= scalar

        for i in range(k - 1, -1, -1):
            scalar = mat[i][k]
            for j in range(N + 1):
                mat[i][j] -= mat[k][j] * scalar

def backward_substitution(mat):
    N = len(mat)
    x = np.zeros(N)  # An array to store solution
    # Start calculating from last equation up to the first
    for i in range(N - 1, -1, -1):
        x[i] = mat[i][N]
        # Initialize j to i+1 since matrix is upper triangular
        for j in range(i + 1, N):
            x[i] -= mat[i][j] * x[j]
        x[i] = (x[i] / mat[i][i])
    return x
def norm(mat):
    size = len(mat)
    max_row = 0
    for row in range(size):
        sum_row = 0
        for col in range(size):
            sum_row += abs(mat[row][col])
        if sum_row > max_row:
            max_row = sum_row
    return max_row
if __name__ == '__main__':

    CB = [
             [-1,  1,     3,   -3,     1,   3],
             [3,  -3,    -4,    2,     3,   8],
             [2,   1,    -5,   -3,     5,   2],
             [-5, -6,     4,    1,     3,  14],
             [3,  -2,    -2,   -3,     5,   6]]
    """norm_A = norm(CB)
    print("\nThe norm of the matrix CB is ", norm_A)"""

    AB = [
             [1,   2,     3,    4,     5],
             [2,   3,     4,    5,     1],
             [8,   8,     8,    8,     1],
             [24, 15,     22,   1,     8],]
    E = [[2, 1/3, -1.2518],
         [1/4, 2, -2.3957]]

    result = gaussianElimination(E)
    if isinstance(result, str):
        print(result)
    else:
        print(bcolors.OKBLUE,"\nSolution for the system:")
        for x in result:
            print("{:.6f}".format(x))