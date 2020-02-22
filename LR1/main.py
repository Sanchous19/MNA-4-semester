import numpy as np
import numpy.linalg as la


def gaussianMethod(n, matrixA, matrixB):
    extendedMatrix = np.column_stack((matrixA, matrixB))
    if la.matrix_rank(extendedMatrix) != la.matrix_rank(matrixA):
        print("Система не имеет решений")
        return
    elif la.matrix_rank(matrixA) != n:
        print("Система имеет бесконечно много решений")
        return
    triangularMatrixA = np.copy(matrixA)
    triangularMatrixB = np.copy(matrixB)
    for i in range(n - 1):
        if triangularMatrixA[i][i] == 0:
            for j in range(i + 1, n):
                if triangularMatrixA[j][i] != 0:
                    temp = triangularMatrixA[i]
                    triangularMatrixA[i] = triangularMatrixA[j]
                    triangularMatrixA[j] = temp
                    break
        for j in range(i + 1, n):
            coefficient = triangularMatrixA[j][i] / triangularMatrixA[i][i]
            triangularMatrixA[j] -= coefficient * triangularMatrixA[i]
            triangularMatrixB[j] -= coefficient * triangularMatrixB[i]
    matrixX = np.zeros(n)
    for i in range(n):
        matrixX[n - 1 - i] = triangularMatrixB[n - 1 - i]
        for j in range(i):
            matrixX[n - 1 - i] -= matrixX[n - 1 - j] * triangularMatrixA[n - 1 - i][n - 1 - j]
        matrixX[n - 1 - i] /= triangularMatrixA[n - 1 - i][n - 1 - i]
    print("Решение:")
    for i in range(n):
        print("x{0}={1}".format(i + 1, matrixX[i]), end=" ")
    print("\n")
    print("Обратная матрица:")
    inverseMatrixA = la.inv(matrixA)
    print(inverseMatrixA)
    normA = la.norm(matrixA, np.inf)
    normB = la.norm(matrixB, np.inf)
    normInverseA = la.norm(inverseMatrixA, np.inf)
    absoluteValue = normInverseA * 0.001
    relativeValue = normA * normInverseA * 0.001 / normB
    print("Оценка абсолютной погрешности решения: {0}".format(absoluteValue))
    print("Оценка относительной погрешности решения: {0}".format(relativeValue))


if __name__ == '__main__':
#     n = 4
#     matrixA = np.array([[3.738, 0.195, 0.275, 0.136],
#                         [0.519, 5.002, 0.405, 0.283],
#                         [0.306, 0.381, 4.812, 0.418],
#                         [0.272, 0.142, 0.314, 3.935]])
#     matrixB = np.array([0.815, 0.191, 0.423, 0.352])
    n = int(input("Введите размерность системы уравнений "))
    print("Система уравнений будет иметь вид:")
    for i in range(n):
        for j in range(n):
            print("a{0}{1}*x{2}".format(i + 1, j + 1, j + 1), end="")
            if j != n - 1:
                print("+", end="")
            else:
                print("=", end="")
        print("b{0}".format(i + 1))
    matrixA = np.zeros((n, n))
    matrixB = np.zeros(n)
    for i in range(n):
        for j in range(n):
            matrixA[i][j] = input("a{0}{1}=".format(i + 1, j + 1))
        matrixB[i] = input("b{0}=".format(i + 1))
    print("Система уравнений будет иметь вид:")
    for i in range(n):
        for j in range(n):
            print("{0}*x{1}".format(matrixA[i][j], j + 1), end="")
            if j != n - 1:
                print("+", end="")
            else:
                print("=", end="")
        print(matrixB[i])
    gaussianMethod(n, matrixA, matrixB)

