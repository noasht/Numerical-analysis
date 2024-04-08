# @source: https://github.com/lihiSabag/Numerical-Analysis-2023.git

import numpy as np
from matrix_utility import (lowerTriangularMatrix, upperTriangularMatrix, diagonalMatrix, inverseD, norm, inverse,)

"""
Performs Jacobi iterations to solve the line system of equations, Ax=b, 
starting from an initial guess, ``x0``.

Terminates when the change in x is less than ``tol``, or
if ``N`` [default=200] iterations have been exceeded.

Receives 5 parameters:
    1.  a, the NxN matrix that method is being performed on.
    2.  b, vector of solution. 
    3.  X0,  the desired initial guess.
        if x is None, the initial guess will bw determined as a vector of 0's.
    4.  TOL, tolerance- the desired limitation of tolerance of solution's anomaly.
        if tolerance is None, the default value will set as 1e-16.
    5.  N, the maxim number of possible iterations to receive the most exact solution.
        if N is None, the default value will set as 200.

Returns variables:
    1.  x, the estimated solution
"""

def G_jacobi_iterative(A,n):
    L = lowerTriangularMatrix(A, n)
    U = upperTriangularMatrix(A, n)
    D = diagonalMatrix(A, n)
    InversD = inverseD(D,n)
    InversD = InversD * -1
    LUSum = L + U
    try:
        inversLandU = inverse(LUSum)
    except ValueError:
        return 2
    G = np.dot(InversD,inversLandU)
    normG = norm(G)
    return normG


def jacobi_iterative(A, b , X0 , TOL=1e-10, N=200):
    print("\033[94m" + "Jacobi Iterative"+ "\033[0m")
    n = len(A)
    if not n == A.shape[1]:
        raise ValueError("The matrix is not square")
    normG = G_jacobi_iterative(A, n)
    if normG > 1:
        raise ValueError(f"The convergence criterion is not met since the norm is {normG} which is bigger then 1")
    k = 1
    x = np.zeros(n, dtype=np.double)
    print('preforming jacobi algorithm\n')
    print(
        "Iteration" + "\t\t\t".join([" {:>12}".format(var) for var in ["x{}".format(i) for i in range(1, len(A) + 1)]]))
    print("-----------------------------------------------------------------------------------------------")
    while k <= N:
        for i in range(n):
            sigma = 0
            for j in range(n):
                if j != i:
                    sigma += A[i][j] * X0[j]
            x[i] = (b[i] - sigma) / A[i][i]

        print("{:<15} ".format(k) + "\t\t".join(["{:<15} ".format(val) for val in x]))
        # The condition if norm(x - X0, np.inf) < TOL:
        # Checks whether the maximum distance between elements of the vectors x and X0 is less than the threshold TOL.
        # This is used to assess the similarity between the vectors. If the maximum distance is smaller than TOL,
        # it can be inferred that the vectors are sufficiently close to each other, and thus considered similar
        if np.linalg.norm(x - X0, np.inf) < TOL:
            return tuple(x)
        k += 1
        X0 = x.copy()
    print("Maximum number of iterations exceeded")
    return tuple(x)

if __name__ == "__main__":

    """A = np.array([
             [-1,  1,     3,   -3,     1],
             [3,  -3,    -4,    2,     3],
             [2,   1,    -5,   -3,     5],
             [-5, -6,     4,    1,     3],
             [3,  -2,    -2,   -3,     5]])
    b = np.array([3, 8, 2, 14, 6])"""
    A = np.array([[3, -1, 1], [0, 1, -1], [1, 1, -2]])
    b = np.array([4, -1, -3])
    x = np.zeros_like(b, dtype=np.double)

    try:
        solution = jacobi_iterative(A, b, x)
        solution = tuple(map(lambda x: round(x, 2), solution))
        print("\nThe root found using the " + "\033[94m" + "Jacobi Iterative" + "\033[0m" + " is " + "\033[94m" + f"{solution}" + "\033[0m")
    except ValueError as e:
        print(str(e))
