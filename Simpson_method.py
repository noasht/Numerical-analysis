# @source: https://github.com/lihiSabag/Numerical-Analysis-2023.git
import math
import sympy as sp
from colors import bcolors
"""
   Simpson's Rule for Numerical Integration

   Parameters:
   f (function): The function to be integrated.
   a (float): The lower limit of integration.
   b (float): The upper limit of integration.
   n (int): The number of subintervals (must be even).

   Returns:
   float: The approximate definite integral of the function over [a, b].
   """
x = sp.symbols('x')
def simpsons_rule(f, a, b, n, Error):
    if n % 2 != 0:
        raise ValueError("Number of subintervals (n) must be even for Simpson's Rule.")
    integralPrev = Simpsons(f, a, b, n)
    n *= 2
    integralNext = Simpsons(f, a, b, n)
    difference = integralNext - integralPrev
    counter = 0
    while difference > Error and counter < 200:
        integralPrev = integralNext
        n *= 2
        integralNext = Simpsons(f, a, b, n)
        difference = integralNext - integralPrev
        counter += 1
    return integralNext


def Simpsons(f, a, b, n):

    h = (b - a) / n
    integral = f(a) + f(b)  # Initialize with endpoints
    for i in range(1, n):
        x_i = a + i * h
        if i % 2 == 0:
            integral += 2 * f(x_i)
        else:
            integral += 4 * f(x_i)
    integral *= h / 3
    return integral


if __name__ == '__main__':
    print("\033[94m" + "Simpsons Rule" + "\033[0m")
    f = lambda x: (math.sin(x ** 2 + 5 * x + 6)) / (2 * math.e ** (-x))
    n = 50
    a=2.7
    b=4.1
    Error = 0.00001
    print( f" Division into n={n} sections ")
    integral = simpsons_rule(f, b, a, n,Error)
    print(bcolors.OKBLUE, f"Numerical Integration of definite integral in range [{a},{b}] is {integral}", bcolors.ENDC)


""""
import math
import sympy as sp
from colors import bcolors
from sympy.utilities.lambdify import lambdify
x = sp.symbols('x')
def simpsons_rule(f, a, b, n):
    Simpson's Rule for Numerical Integration

    Parameters:
    f (function): The function to be integrated.
    a (float): The lower limit of integration.
    b (float): The upper limit of integration.
    n (int): The number of subintervals (must be even).

    Returns:
    float: The approximate definite integral of the function over [a, b].
    if n % 2 != 0:
        return ValueError("Number of subintervals (n) must be even for Simpson's Rule.")

    h = (b - a) / n

    integral = f(a) + f(b)  # Initialize with endpoints

    for i in range(1, n):
        x_i = a + i * h
        if i % 2 == 0:
            integral += 2 * f(x_i)
        else:
            integral += 4 * f(x_i)

    integral *= h / 3

    return integral

def error_s(f_s, a, b, error,n):
    h = (b - a) / n
    x = sp.symbols('x')
    df1 = sp.diff(f_s, x)
    df2 = sp.diff(df1, x)
    df3 = sp.diff(df2, x)
    df4 = sp.diff(df3, x)
    f_tag4 = lambdify(x, df4)
    er=f_tag4(error)
    y = (1/180)*(b-a)*(h**4)
    finallError= y*er
    print("-----------------------------------------------------")
    print(bcolors.OKBLUE,"The Error: ", finallError,bcolors.ENDC)
    print("-----------------------------------------------------")



if __name__ == '__main__':
    f = lambda x: math.e ** (x ** 2)
    n = 30
    a=0
    b=1
    print( f" Division into n={n} sections ")
    integral = simpsons_rule(f, 0, 1, n)
    print(bcolors.OKBLUE, f"Numerical Integration of definite integral in range [{a},{b}] is {integral}", bcolors.ENDC)
    f_s =  math.e ** (x ** 2)
    error_s(f_s, a, b, 1, n)"""