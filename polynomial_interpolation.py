# @source: https://github.com/lihiSabag/Numerical-Analysis-2023.git
from colors import bcolors
from matrix_utility import *
import numpy as np


def polynomialInterpolation(matrix,table_points, x):

    b = [[point[1]] for point in table_points]
    matrixNew = np.hstack((matrix, b))
    print(bcolors.OKBLUE, "The matrix obtained from the points: ", bcolors.ENDC,'\n', np.array(matrix))
    print(bcolors.OKBLUE, "\nb vector: ", bcolors.ENDC,'\n',np.array(b))
    matrixSol = gaussianElimination(matrixNew)
    if matrixSol is not None:
        print(bcolors.OKBLUE, "\nResult Gauss: ", bcolors.ENDC, '\n', np.array(matrixSol))
        result = sum([matrixSol[i] * (x ** i) for i in range(len(matrixSol))])
        print(bcolors.OKBLUE, "\nThe polynom:", bcolors.ENDC)
        print('P(X) = '+'+'.join([ '('+str(matrixSol[i])+') * x^' + str(i) + ' ' for i in range(len(matrixSol))]))
        print(bcolors.OKGREEN, f"\nThe Result of P(X={x}) is:", bcolors.ENDC)
        print(result)
        return result
    return None

def Prerequisite(table_points):
    matrix = [[point[0] ** i for i in range(len(table_points))] for point in table_points]  # Makes the initial matrix
    if not np.linalg.det(matrix):
        print("Singular Matrix")
        return None
    return matrix

if __name__ == '__main__':

    print(bcolors.OKBLUE, "----------------- Interpolation & Extrapolation Methods -----------------", bcolors.ENDC)
    table_points = [(1, 1), (2, 0), (5, 2)]
    x = 3
    matrix = Prerequisite(table_points)
    if matrix is not None:
        print(bcolors.OKBLUE, "Table Points: ", bcolors.ENDC, table_points)
        print(bcolors.OKBLUE, "Finding an approximation to the point: ", bcolors.ENDC, x,'')
        if polynomialInterpolation(matrix, table_points, x) is None:
            print("     Singular Matrix")
        print(bcolors.OKBLUE, "---------------------------------------------------------------------------", bcolors.ENDC)
