# Рассчет значения полинома
def calc(f, x):
    res = 0.0
    for i in range(len(f)):
        res += f[i] * (x ** i)
    return res


# Рассчет производной полинома
def d(f):
    df = f[:len(f) - 1:]

    for i in range(len(f) - 1, 0, -1):
        df[i - 1] = f[i] * i
    return df


# Упрощенный расчет производной
def simple_d(f, x0, x1=None, delta=None):
    if delta is None:
        return (calc(f, x1) - calc(f, x0)) / (x1 - x0)
    else:
        return (calc(f, x0) - calc(f, x0 - delta)) / delta


# Метод Ньютона и Ньютона-Бройдена
def newton(f, df, x0, eps, curr_iter=1, c=1.0, iter_cnt=100000):
    x1 = x0 - c * (calc(f, x0) / calc(df, x0))
    while abs((x1 - x0) / (1.0 - (x1 - x0) / (x0 - x1))) > eps:
        if iter_cnt > curr_iter:
            x0 = x1
            curr_iter += 1
        else:
            return None
        x1 = x0 - c * (calc(f, x0) / calc(df, x0))
    return x1


# Упрощенный метод Ньютона
def simple_newton(f, x0, fx0, eps, curr_iter=1, iter_cnt=100000):
    x1 = x0 - calc(f, x0) / fx0
    while abs((x1 - x0) / (1.0 - (x1 - x0) / (x0 - x1))) > eps:
        if iter_cnt > curr_iter:
            x0 = x1
            curr_iter += 1
        else:
            return None
        x1 = x0 - (calc(f, x0) / fx0)
    return x1


# Метод секущих
def seq(f, x0, delta, eps, curr_iter=1, iter_cnt=10000):
    x1 = x0 - calc(f, x0) / simple_d(f, x0, delta=delta)
    while abs((x1 - x0) / (1.0 - (x1 - x0) / (x0 - x1))) > eps:
        if iter_cnt > curr_iter:
            x0, x_last = x1, x0
            curr_iter += 1
        else:
            return None
        x1 = x0 - calc(f, x0) / simple_d(f, x_last, x1=x0)
    return x1


# Метод хорд
def hords(f, a, b, eps, curr_iter=1, iter_cnt=10000):
    if calc(f, a) * calc(d(d(f)), a) < 0:
        a, b = b, a
    x = a - (calc(f, a)*(b - a)) / (calc(f, b) - calc(f, a))
    while abs((x - a) / (1.0 - (x - a) / (a - x))) > eps:
        if iter_cnt > curr_iter:
            a = x
            curr_iter += 1
            x = a - (calc(f, a) * (b - a)) / (calc(f, b) - calc(f, a))
        else:
            return None
    return x


f = [-10, 10, 0, 4]
print("Полином -10 + 10х + 4х^3, x0 = 2\n")
print("Метод Ньютона-Бройдена: x = ", newton(f, d(f), 2, 0.000001, 1, 0.3))
print("Метод Ньютона: x = ", newton(f, d(f), 2, 0.000001, 1))
print("Упрпощенный метод Ньютона: x = ", simple_newton(f, 2, calc(d(f), 2), 0.000001))
print("Метод секущих: delta = 0.0001, x = ", seq(f, 2, 0.0001, 0.000001))
print("Метод хорд: a = 3, b = -1, x = ", hords(f, 3, -1, 0.000001))
print("\n---------------------\n")
f = [3, 0, 1, 1]
print("Полином 3 + x^2 + х^3, x0 = 4\n")
print("Метод Ньютона-Бройдена: x = ", newton(f, d(f), 4, 0.000001, 1, 0.3))
print("Метод Ньютона: x = ", newton(f, d(f), 4, 0.000001, 1))
print("Упрпощенный метод Ньютона: x = ", simple_newton(f, 4, calc(d(f), 4), 0.000001))
print("Метод секущих: delta=10, x =", seq(f, 4, 10, 0.000001))
print("Метод хорд: a = 2, b = -3, x =", hords(f, 2, -3, 0.000001))




