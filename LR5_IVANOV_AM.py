from LR3_IVANOV_AM import gaussMethod
import numpy as np

X = [3, 4, 5, 6]
Y = [1, 0, 4, 2]

W = np.array([[X[i] ** j for j in range(len(X))] for i in range(len(X))], dtype=float)
W = np.c_[W, np.array(Y).T]

a = gaussMethod(W)
x = 4.5
value = 0
for i in range(len(a)):
    value += a[i] * x ** i

print(f"Решение глобальным методом интерполяции: {value}\n")


def eval_p1i(x, i):
    return (x - X[i + 1]) / (X[i] - X[i + 1])


def eval_p1i1(x, i):
    return (x - X[i]) / (X[i + 1] - X[i])


def eval_p2i(x, i):
    return (x - X[i + 1]) * (x - X[i + 2]) / ((X[i] - X[i + 1]) * (X[i] - X[i + 2]))


def eval_p2i1(x, i):
    return (x - X[i]) * (x - X[i + 2]) / ((X[i + 1] - X[i]) * (X[i + 1] - X[i + 2]))


def eval_p2i2(x, i):
    return (x - X[i]) * (x - X[i + 1]) / ((X[i + 2] - X[i]) * (X[i + 2] - X[i + 1]))


x = 4.5
i = 1
P1i = eval_p1i(x, i)
P1j = eval_p1i1(x, i)
ans = P1i * Y[i] + P1j * Y[i + 1]
print(f"Решение кусочным-линейным методом {ans}\n")


x = 4.5
L21 = eval_p2i(x, i - 1) * Y[i - 1] + eval_p2i1(x, i - 1) * Y[i] + eval_p2i2(x, i - 1) * Y[i + 1]
L22 = eval_p2i(x, i) * Y[i] + eval_p2i1(x, i) * Y[i + 1] + eval_p2i2(x, i) * Y[i + 2]
ans = (L21 + L22) / 2.0
print(f"Решение кусочным-параболическим методом {ans}\n")
