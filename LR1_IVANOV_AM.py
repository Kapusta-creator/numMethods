# Номер 1 (метод дихотомии)
# Рассчет значения полинома
def calc(f, x):
    res = 0.0
    for i in range(len(f)):
        res += f[i] * (x ** i)
    return res


# Метод дихотомии
def dih(f, x1, x2, eps):
    mid = 0
    while abs(x1 - x2) > eps:
        mid = (x1 + x2) / 2
        if calc(f, x1) * calc(f, mid) < 0:
            x2 = mid
        else:
            x1 = mid
    return mid


# Заполнение массива корней для производных полинома и самого полинома
def roots_for_pr(level, first, second, roots_cnt):
    major = 0
    for i in range(level):
        s = abs(first[level][i])
        if s > major:
            major = s
    major += 1.0

    roots_cnt[level] = 0
    for i in range(roots_cnt[level - 1] + 1):
        if i == 0:
            edgeLeft = -major
        else:
            edgeLeft = second[level - 1][i - 1]
        rb = calc(first[level], edgeLeft)
        if rb == 0:
            second[level][roots_cnt[level]] = edgeLeft
            roots_cnt[level] += 1
            continue
        if rb > 0:
            left = 1
        else:
            left = -1
        if i == roots_cnt[level - 1]:
            edgeRight = major
        else:
            edgeRight = second[level - 1][i]

        rb = calc(first[level], edgeRight)
        if rb == 0:
            second[level][roots_cnt[level]] = edgeRight
            roots_cnt[level] += 1
            continue

        if rb > 0:
            right = 1
        else:
            right = -1

        if left == right:
            continue

        if left < 0:
            edgeNegativ = edgeLeft
            edgePositive = edgeRight
        else:
            edgeNegativ = edgeRight
            edgePositive = edgeLeft

        second[level][roots_cnt[level]] = dih(first[level], edgeNegativ, edgePositive, 0.00000001)
        roots_cnt[level] += 1


# Входная функция алгоритма, нахождение производных полинома
def get_roots(polynom):
    first = [[0 for i in range(j + 1)] for j in range(len(polynom))]
    second = [[0 for i in range(j + 1)] for j in range(len(polynom))]
    n = len(polynom) - 1
    roots_cnt = [0 for i in range(n + 1)]
    for i in range(len(polynom)):
        first[n][i] = polynom[i] / polynom[n]
    for i in range(n - 1, -1, -1):
        for j in range(i, -1, -1):
            first[i][j] = first[i + 1][j + 1]*(j + 1)
    for i in range(n - 1, -1, -1):
        for j in range(len(first[i])):
            first[i][j] = first[i][j] / first[i][len(first[i]) - 1]
    second[1][0] = -first[1][0]
    roots_cnt[1] = 1
    for i in range(2, n + 1):
        roots_for_pr(i, first, second, roots_cnt)
    rootsCount = roots_cnt[n]
    roots = [0 for i in range(rootsCount)]
    for i in range(rootsCount):
        roots[i] = second[n][i]
    return roots


p = [-192, -160, 68, 24, -11, 1]  # Полином x^5 - 11x^4 + 24x^3 + 68x^2 - 160x -192
print("Пример работы алгоритма задания 1:\n")
print("Полином x^5 - 11x^4 + 24x^3 + 68x^2 - 160x -192\n")
print(f"Найденные корни: {get_roots(p)}")


# Номер 2(нахождение корня)
def get_sqrt(num, eps):
    xn = num
    m = dict()
    while m.get(xn, 1):
        m[xn] = 0
        next = (xn + num / xn) / 2
        if abs(xn - next) < eps:
            return next
        xn = next
    return -1


print("-" * 20)
print("\nПример работы алгоритма задания 2:\n")
a = 214
print(f"Число a = {a}\n")
print(f"Найденный корень: {get_sqrt(a, 0.0001)}\n")
