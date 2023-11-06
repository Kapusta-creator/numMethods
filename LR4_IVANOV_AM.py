import numpy as np
from math import atan

A = np.array(
    [
        [5, 1, 2],
        [1, 4, 1],
        [2, 1, 3]
    ], dtype=float
)

X = np.array(
    [1, 1, 1], dtype=float
)


def get_norm(x):
    return np.sqrt(np.sum(x ** 2))


def iter(a, x, eps=0.001, max_iter=1000):
    lambda_ = 0
    for i in range(max_iter):
        x1 = np.dot(a, x)
        norm1 = get_norm(x1)
        lambda1 = norm1 / get_norm(x)

        x = x1 / norm1
        if abs(lambda1 - lambda_) < eps:
            break
        lambda_ = lambda1

    return lambda_, x


if __name__ == "__main__":
    print("Исходная матрица\n")
    print(A)


    print("Решение методом итераций\n")
    value, vector = iter(A, X)
    print(f"Собственное число {value}\n", f"Собственный вектор {vector}\n", sep="")


def get_max_upper_diag(A):
    max = abs(A[0][1])
    row, col = 0, 1
    for i in range(len(A)):
        for j in range(i + 1, len(A)):
            if abs(A[i][j]) > max:
                max = abs(A[i][j])
                row, col = i, j
    return max, row, col


def rotate(A, eps=0.001, max_iter=1000):
    vec = np.eye(len(A))
    for i in range(max_iter):
        max_a, row, col = get_max_upper_diag(A)

        if max_a < eps:
            break

        if A[row][row] == A[col][col]:
            phi = np.pi / 4
        else:
            phi = 0.5 * atan(2 * A[row][col] / (A[row][row] - A[col][col]))

        H = np.eye(len(A))
        H[row][row] = np.cos(phi)
        H[row][col] = -np.sin(phi)
        H[col][row] = np.sin(phi)
        H[col][col] = np.cos(phi)

        A = H.T @ A @ H
        vec = vec @ H
    lambdas = np.diag(A)

    return lambdas, vec


A = np.array([[5, 1, 2],
              [1, 4, 1],
              [2, 1, 3]])


if __name__ == "__main__":
    print("Решение методом вращения\n")
    lambdas, vec = rotate(A)
    print(f"Собственные числа {lambdas}\n")
    print("Собственные вектора:\n")
    print(*vec.T, sep="\n")

