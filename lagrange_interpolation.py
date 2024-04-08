# @source: https://github.com/lihiSabag/Numerical-Analysis-2023.git
def lagrange_interpolation(x_data, y_data, x):
    n = len(x_data)
    result = 0.0
    for i in range(n):
        term = y_data[i]
        for j in range(n):
            if i != j:
                term *= (x - x_data[j]) / (x_data[i] - x_data[j])
        result += term
    return result


if __name__ == '__main__':
    print("\033[94m"+"Lagrange Interpolation" +"\033[0m")
    x_data = [1, 2, 5]
    y_data = [1, 0, 2]
    x_interpolate = 3  # The x-value where you want to interpolate
    y_interpolate = lagrange_interpolation(x_data, y_data, x_interpolate)
    print("\nInterpolated value at" +"\033[94m"+f" x = {x_interpolate}" +"\033[0m"+ " is " + "\033[94m"+f"y = {round(y_interpolate,4)}")