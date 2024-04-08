# @source: https://github.com/lihiSabag/Numerical-Analysis-2023.git
import sympy as sp
def newton_raphson(polynomial, intervalA, intervalB, TOL=1e-6, N=50):
    x = sp.symbols('x')
    f = sp.sympify(polynomial)
    df = sp.diff(f, x)
    df_func = sp.lambdify(x, df)
    Xr0 = float((intervalA + intervalB) / 2)
    if df_func(Xr0) == 0:
        raise ValueError(f"Derivative is zero at {Xr0}, method cannot continue.")
    print("\033[94m" +"Newton Raphson" + "\033[0m")
    print("{:<10} {:<15} {:<15} ".format("Iteration", "Xr0", "Xr"))
    for i in range(N):
        if df_func(Xr0) == 0:
            raise ValueError(f"Derivative is zero at {Xr0}, method cannot continue.")
        Xr = Xr0 - f.evalf(subs={x: Xr0}) / df_func(Xr0)
        if abs(Xr - Xr0) < TOL:
            if intervalA <= Xr <= intervalB:
                return Xr
            else:
                return None
        print("{:<10} {:<15.9f} {:<15.9f} ".format(int(i), float(Xr0), float(Xr)))
        Xr0 = Xr
    if intervalA <= Xr <= intervalB:
        return Xr
    else:
        return None

if __name__ == '__main__':
    x = sp.symbols('x')
    f = x**2 + 3*x + 2.25
    intervalA = 0
    intervalB = 6
    try:
        roots = newton_raphson(f, intervalA, intervalB)
        if roots is None:
            print("\033[1m" + "Newton Raphson " +"\033[0m"+ "\033[0m" +f"No roots were found within the given range [{intervalA}, {intervalB}]"+"\033[94m")
        else:
            roots = round(roots, 5)
            print("\nThe Final Root found using the " + "\033[94m" + "Newton Raphson" + "\033[0m" +"\033[94m" +
                  f"in range [{intervalA}, {intervalB}] x =" "{:<15.9f}".format(roots)+ "\033[0m")
    except ValueError as e:
        print(str(e))
