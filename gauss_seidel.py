# @source: https://github.com/lihiSabag/Numerical-Analysis-2023.git
import numpy as np
from numpy.linalg import norm
from colors import bcolors
from matrix_utility import (lowerTriangularMatrix, upperTriangularMatrix, diagonalMatrix, norm, inverse,
                            Function_for_Checking_Dominant_Diagonal_in_Matrix)

def G_gauss_seidel(A,n):
    L = lowerTriangularMatrix(A, n)
    U = upperTriangularMatrix(A, n)
    D = diagonalMatrix(A, n)
    sumLandD = L + D
    try:
        inversLsumD = inverse(sumLandD)
    except ValueError:
        return 2
    inversLsumD = inversLsumD * -1
    G = np.dot(inversLsumD, U)
    normG = norm(G)
    return normG


def gauss_seidel(A, b, X0, TOL=1e-10, N=200):
    print("\033[94m" +"Gauss Seidel"+ "\033[0m" )
    n = len(A)
    if not n == A.shape[1]:
        raise ValueError("The matrix is not square")
    k = 1
    x = np.zeros(n, dtype=np.double)
    """if not Function_for_Checking_Dominant_Diagonal_in_Matrix(A):
        raise ValueError(f"The convergence criterion is not met")"""
    normG = G_gauss_seidel(A, n)
    if normG > 1:
            raise ValueError(f"The convergence criterion is not met since the norm is {normG} which is bigger then 1")
    print('preforming gauss seidel algorithm')
    print("Iteration" + "\t\t\t".join([" {:>12}".format(var) for var in ["x{}".format(i) for i in range(1, len(A) + 1)]]))
    print("-----------------------------------------------------------------------------------------------")
    while k <= N:
        for i in range(n):
            sigma = 0
            for j in range(n):
                if j != i:
                    sigma += A[i][j] * x[j]
            x[i] = (b[i] - sigma) / A[i][i]
        print("{:<15} ".format(k) + "\t\t".join(["{:<15} ".format(val) for val in x]))
        if np.linalg.norm(x - X0, np.inf) < TOL:
            return tuple(x)
        k += 1
        X0 = x.copy()
    print("Maximum number of iterations exceeded")
    return tuple(x)

if __name__ == '__main__':

    """A = np.array([
        [-1, 1, 3, -3, 1],
        [3, -3, -4, 2, 3],
        [2, 1, -5, -3, 5],
        [-5, -6, 4, 1, 3],
        [3, -2, -2, -3, 5]])
    b = np.array([3, 8, 2, 14, 6])"""
    A = np.array([[3, -1, 1], [0, 1, -1], [1, 1, -2]])
    b = np.array([4, -1, -3])
    X0 = np.zeros_like(b)
    try:
        solution =gauss_seidel(A, b, X0)
        solution = tuple(map(lambda x: round(x, 2), solution))
        print("\nThe root found using the " + "\033[94m" + "Gauss Seidel" + "\033[0m" + " is " + "\033[94m" + f"{solution}" + "\033[0m")
    except ValueError as e:
        print(str(e))
