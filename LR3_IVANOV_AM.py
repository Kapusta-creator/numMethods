import copy
import math

import numpy as np

augMatrix = np.array([[1.0, -1.0, 1.0, -4.0, -2.0],
                      [2.0, 1.0, -5.0, 1.0, 2.0],
                      [8.0, -1.0, -1.0, 2.0, 11.0],
                      [1.0, 6.0, -2.0, -2.0, -7.0]])

augMatrix3 = np.array([[1.0, 4.0, 3.0, 10.0],
                      [2.0, 1.0, -1.0, -1.0],
                      [3.0, -1.0, 1.0, 11.0]])

augMatrix4 = np.array([[0.0, 4.0, 3.0, 10.0],
                      [2.0, 0.0, -1.0, -1.0],
                      [3.0, -1.0, 0.0, 11.0]])

augMatrix6 = np.array([[2, 2, 10, 14], [10, 1, 1, 12], [2, 10, 1, 13]], dtype=float)


def gaussMethod(matrix):
    for i, row in enumerate(matrix):
        a = row[i]
        row /= a
        for lower_row in matrix[i + 1:]:
            lower_row -= lower_row[i] * row

    for i in range(len(matrix) - 1, 0, -1):
        row = matrix[i]
        for upper_row in matrix[:i]:
            upper_row[-1] -= upper_row[i]*row[-1]
            upper_row[i] = 0

    return matrix[:, -1]


def gaussMethodGeneralElem(matrix):
    for i in range(len(matrix)):
        max_row = abs(matrix[i:, i]).argmax() + i
        if max_row != i:
            matrix[[i, max_row]] = matrix[[max_row, i]]
        a = matrix[i, i]
        matrix[i] /= a
        for lower_row in matrix[i + 1:]:
            lower_row -= lower_row[i] * matrix[i]
    for i in range(len(matrix) - 1, 0, -1):
        row = matrix[i]
        for upper_row in matrix[:i]:
            upper_row[-1] -= upper_row[i]*row[-1]
            upper_row[i] = 0

    return matrix[:, -1]


def gaussMethodRect(matrix):
    for k in range(len(matrix)):
        matCopy = copy.deepcopy(matrix)
        a = matrix[k][k]
        for i in range(k, len(matrix)):
            for j in range(len(matrix[0])):
                if i > k and j <= k:
                    matrix[i][j] = 0
                elif i == k:
                    matrix[i][j] /= a
                else:
                    matrix[i][j] = matrix[i][j] - matCopy[k][j] * matCopy[i][k] / a
    for i in range(len(matrix) - 1, 0, -1):
        row = matrix[i]
        for upper_row in matrix[:i]:
            upper_row[-1] -= upper_row[i]*row[-1]
            upper_row[i] = 0

    return matrix[:, -1]


def makeAlphaBeta(matrix):
    alpha = copy.deepcopy(matrix)
    for i in range(len(alpha)):
        for j in range(len(alpha[i])):
            if i == j:
                continue
            if i != j and j < len(alpha[i]) - 1:
                alpha[i][j] = -alpha[i][j] / alpha[i][i]
            else:
                alpha[i][j] = alpha[i][j] / alpha[i][i]
    for i in range(len(alpha)):
        alpha[i][i] = 0
    return alpha[:, 0:len(alpha[0]) - 1], np.column_stack([alpha[:, -1]])


def normaMatrix(matrix):
    summ = 0
    for vec in matrix:
        for i in range(len(vec)):
            summ += vec[i] * vec[i]
    return math.sqrt(summ)


def normaVec(vec):
    summ = 0
    for i in range(len(vec)):
        summ += vec[i] * vec[i]
    return math.sqrt(summ)


def simpleIterMethod(matrix, eps, iterCnt):
    for i in range(len(matrix)):
        max_row = abs(matrix[i:, i]).argmax() + i
        if max_row != i:
            matrix[[i, max_row]] = matrix[[max_row, i]]
    alpha, beta = makeAlphaBeta(matrix)
    x0 = copy.deepcopy(beta)
    x1 = np.dot(alpha, x0) + x0
    currIter = 0
    while normaMatrix(alpha) / (1 - normaMatrix(alpha)) * normaVec(x1 - x0) > eps and currIter < iterCnt:
        x0 = copy.deepcopy(x1)
        x1 = np.dot(alpha, x1) + beta
        currIter += 1
    return x1


def zeydelMethod(matrix, eps, iterCnt):
    for i in range(len(matrix)):
        max_row = abs(matrix[i:, i]).argmax() + i
        if max_row != i:
            matrix[[i, max_row]] = matrix[[max_row, i]]
    alpha, beta = makeAlphaBeta(matrix)
    x0 = copy.deepcopy(beta)
    x1 = copy.deepcopy(x0)
    x1[0] = np.dot(alpha[0, :], x0) + x0[0]
    currIter = 1
    while normaVec(x1 - x0) > eps and currIter < iterCnt:
        x0 = copy.deepcopy(x1)
        x1[currIter % len(alpha)] = np.dot(alpha[currIter % len(alpha), :], x1) + beta[currIter % len(alpha)]
        currIter += 1
    return x1


def solve_LU(a, b):
    lu_matrix = np.matrix(np.zeros([a.shape[0], a.shape[1]]))
    n = a.shape[0]

    for k in range(n):
        for j in range(k, n):
            lu_matrix[k, j] = a[k, j] - lu_matrix[k, :k] * lu_matrix[:k, j]
        for i in range(k + 1, n):
            lu_matrix[i, k] = (a[i, k] - lu_matrix[i, : k] * lu_matrix[: k, k]) / lu_matrix[k, k]

    y = np.matrix(np.zeros([lu_matrix.shape[0], 1]))
    for i in range(y.shape[0]):
        y[i, 0] = b[i, 0] - lu_matrix[i, :i] * y[:i]

    x = np.matrix(np.zeros([lu_matrix.shape[0], 1]))
    for i in range(1, x.shape[0] + 1):
        x[-i, 0] = (y[-i] - lu_matrix[-i, -i:] * x[-i:, 0] )/ lu_matrix[-i, -i]

    return x


print(f"\nМатрица 1 \n{augMatrix}\n")
print("Метод единичного деления, матрица 1\n", gaussMethod(copy.deepcopy(augMatrix)))
print(f"\nМатрица 4\n{augMatrix4}\n")
print("Выбор ведущего элемента, матрица 4\n", gaussMethodGeneralElem(copy.deepcopy(augMatrix4)))
print(f"\nМатрица 3\n{augMatrix3}\n")
print("Метод исключения, матрица 3\n", gaussMethodRect(copy.deepcopy(augMatrix3)))
a3 = np.array([[1.0, 4.0, 3.0],
               [2.0, 1.0, -1.0],
               [3.0, -1.0, 1.0]])
b3 = np.array([[10.0], [-1.0], [11.0]])
print("\n")
print("LU- разложение, матрица 3\n", *solve_LU(copy.deepcopy(a3), copy.deepcopy(b3)))
print(f"\nМатрица 6\n{augMatrix6}\n")
print("Метод простых итераций, матрица 6\n", simpleIterMethod(copy.deepcopy(augMatrix6), 0.000000001, 100000000))
print("\n")
print("Метод Зейделя, матрица 6\n", zeydelMethod(copy.deepcopy(augMatrix6), 0.000000001, 100000000))
print("\n")




