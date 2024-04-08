# @source: https://github.com/lihiSabag/Numerical-Analysis-2023.git
import sympy as sp
def secant_method(polynomial, intervalA, intervalB, TOL=1e-6, N=50):
    x = sp.symbols('x')
    f = sp.sympify(polynomial)
    func = sp.lambdify(x, f)
    Xr0 = intervalA
    Xr = intervalB
    if func(Xr) - func(Xr0) == 0:
        print("method cannot continue. Please select a different range")
        return
    print("{:<10} {:<15} {:<15} {:<15}".format("Iteration", "Xr0", "Xr", "Ans"))
    for i in range(N):
        if func(Xr) - func(Xr0) == 0:
            print("Method cannot continue. Please select a different range")
            return
        Ans = Xr0 - (func(Xr0) * ((Xr - Xr0) / (func(Xr) - func(Xr0))))
        if abs(Ans - Xr) < TOL:
            if intervalA <= Ans <= intervalB:
                return Ans
            else:
                return None # Procedure completed successfully
        print("{:<10} {:<15.6f} {:<15.6f} {:<15.6f}".format(int(i), float(Xr0), float(Xr),float(Ans)))
        Xr0 = Xr
        Xr = Ans
    if intervalA <= Ans <= intervalB:
        return Ans
    else:
        return None


if __name__ == '__main__':
    print("\033[1m" + "Secant method " + "\033[0m")
    x = sp.symbols('x')
    polynomial = x**2 - 5*x + 2
    intervalA = 0
    intervalB = 5
    roots = secant_method(polynomial, intervalA, intervalB)
    if not roots:
        print("\033[0m" +f"No roots were found within the given range [{intervalA}, {intervalB}]"+"\033[94m")
    else:
        print("\nThe Final Roots found using the " + "\033[94m" + "Secant Method" + "\033[0m" + "\033[94m" + "is x = " +
              "\033[94m" +f" {roots} "+ "\033[0m")