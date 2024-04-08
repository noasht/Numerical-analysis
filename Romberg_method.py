import numpy as np
from colors import bcolors
import math

def romberg_integration(func, a, b, n, error):
    """
    Romberg Integration

    Parameters:
    func (function): The function to be integrated.
    a (float): The lower limit of integration.
    b (float): The upper limit of integration.
    n (int): The number of iterations (higher value leads to better accuracy).

    Returns:
    float: The approximate definite integral of the function over [a, b].
    """
    h = b - a
    R = np.zeros((n, n), dtype=float)

    number = 0
    R[0, 0] = 0.5 * h * (func(a) + func(b))
    print(bcolors.OKGREEN, "-> " + str(R[0, 0]), bcolors.ENDC)
    number += 1

    for i in range(1, n):
        h /= 2
        sum_term = 0

        for k in range(1, 2 ** i, 2):
            sum_term += func(a + k * h)

        R[i, 0] = 0.5 * R[i - 1, 0] + h * sum_term
        number += 1

        for j in range(1, i + 1):
            R[i, j] = R[i, j - 1] + (R[i, j - 1] - R[i - 1, j - 1]) / ((4 ** j) - 1)
            number += 1
        print(bcolors.OKGREEN, "-> " + str(R[i, i]), bcolors.ENDC)

        if n > 0 and abs(R[i, i] - R[i - 1, i - 1]) < error:
            difference = abs(R[i, i] - R[i - 1, i - 1])
            print("difference =" + str(difference))
            return R[i, i], i

    return R[n - 1, n - 1], number


def f(x):
    return (math.sin(x ** 2 + 5 * x + 6)) / (2 * math.e ** (-x))


if __name__ == '__main__':
    a = 2.7
    b = 4.1
    n = 30

    error = 0.00001
    integral, index = romberg_integration(f, a, b, n, error)

    print(f" Division into n={index} sections ")
    print(bcolors.OKBLUE, f"Approximate integral in range [{a},{b}] is {integral}", bcolors.ENDC)