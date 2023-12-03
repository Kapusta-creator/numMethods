def divided_difference(x, y):
    n = len(x)
    f = [[0] * n for i in range(n)]
    for i in range(n):
        f[i][0] = y[i]
    for j in range(1, n):
        for i in range(n - j):
            f[i][j] = (f[i + 1][j - 1] - f[i][j - 1]) / (x[i + j] - x[i])
    return f[0]


def newton_interpolation_kus(x, y, xi):
    n = len(x)
    f = divided_difference(x, y)
    p = f[0]
    for i in range(1, n):
        product = f[i]
        for j in range(i):
            product *= (xi - x[j])
        p += product
    return p


def delta_f(x, y):
    n = len(x)
    f = [[0] * n for i in range(n)]
    for i in range(n):
        f[i][0] = y[i]
    for j in range(1, n):
        for i in range(n - j):
            f[i][j] = (f[i + 1][j - 1] - f[i][j - 1])
    return f


def newton_interpolation_1(x, y, xi, h):
    q = (xi - x[0]) / h
    f = delta_f(x, y)[0]
    p = f[0]
    fact = 1
    for i in range(1, len(f)):
        product = f[i]
        for j in range(i):
            product *= (q - j)
        p += product / fact
        fact *= (i + 1)
    return p


def newton_interpolation_2(x, y, xi, h):
    q = (xi - x[3]) / h
    f = delta_f(x, y)
    p = f[len(f) - 1][0]
    fact = 1
    for i in range(1, len(f)):
        product = f[len(f) - 1 - i][i]
        for j in range(i):
            product *= (q + j)
        p += product / fact
        fact *= (i + 1)
    return p


X = [0, 1, 2, 3]
Y = [1, 2, 4, 1]
print(newton_interpolation_kus(X, Y, 1.5))
print(newton_interpolation_1(X, Y, 1.5, 1))
print(newton_interpolation_2(X, Y, 1.5, 1))
X = [0, 1, 2, 3, 4]
Y = [1, 2, 4, 1, 0]
print(newton_interpolation_1(X, Y, 1.5, 1))
