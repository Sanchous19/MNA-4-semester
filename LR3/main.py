import numpy as np
import numpy.linalg as la
import math


def squareRootMethod(n, matrixA, matrixB):
    extendedMatrix = np.column_stack((matrixA, matrixB))
    if la.matrix_rank(extendedMatrix) != la.matrix_rank(matrixA):
        print("Система не имеет решений")
        return
    elif la.matrix_rank(matrixA) != n:
        print("Система имеет бесконечно много решений")
        return
    matrixB = np.matmul(np.transpose(matrixA), matrixB)
    matrixA = np.matmul(np.transpose(matrixA), matrixA)
    print("После симметризации система уравнений будет иметь вид:")
    printSystemOfEquations(n, matrixA, matrixB)
    matrixU = np.triu(matrixA)
    for i in range(n):
        for j in range(i):
            matrixU[i][i] -= matrixU[j][i] ** 2
        matrixU[i][i] = math.sqrt(matrixU[i][i])
        for j in range(i + 1, n):
            for k in range(i):
                matrixU[i][j] -= matrixU[k][i] * matrixU[k][j]
            matrixU[i][j] /= matrixU[i][i]
    matrixTransU = np.transpose(matrixU)
    matrixY = np.copy(matrixB)
    for i in range(n):
        for j in range(i):
            matrixY[i] -= matrixY[j] * matrixTransU[i][j]
        matrixY[i] /= matrixTransU[i][i]
    matrixX = np.copy(matrixY)
    for i in range(n):
        for j in range(i):
            matrixX[n - 1 - i] -= matrixX[n - 1 - j] * matrixU[n - 1 - i][n - 1 - j]
        matrixX[n - 1 - i] /= matrixU[n - 1 - i][n - 1 - i]
    print("Решение:")
    for i in range(n):
        print("x{0}={1}".format(i + 1, matrixX[i]), end=" ")
    print("\n")
    determinant = 1
    for i in range(n):
        determinant *= matrixU[i][i] ** 2
    print("Определитель матрицы A равен: |A|={0}".format(determinant), end="\n\n")
    inverseMatrixA = calculateInverseMatrix(matrixU)
    print("Обратная матрица:")
    print(inverseMatrixA)


def calculateInverseMatrix(matrixU):
    matrixTransU = np.transpose(matrixU)
    inverseMatrix = np.empty((n, 0))
    for k in range(n):
        matrixE = np.zeros(n)
        matrixE[k] = 1
        matrixY = np.copy(matrixE)
        for i in range(n):
            for j in range(i):
                matrixY[i] -= matrixY[j] * matrixTransU[i][j]
            matrixY[i] /= matrixTransU[i][i]
        matrixX = np.copy(matrixY)
        for i in range(n):
            for j in range(i):
                matrixX[n - 1 - i] -= matrixX[n - 1 - j] * matrixU[n - 1 - i][n - 1 - j]
            matrixX[n - 1 - i] /= matrixU[n - 1 - i][n - 1 - i]
        inverseMatrix = np.column_stack((inverseMatrix, matrixX))
    return inverseMatrix


def printSystemOfEquations(n, matrixA, matrixB):
    for i in range(n):
        for j in range(n):
            print("{0}*x{1}".format(matrixA[i][j], j + 1), end="")
            if j != n - 1:
                print("+", end="")
            else:
                print("=", end="")
        print(matrixB[i])
    print()


if __name__ == '__main__':
    n = 4
    matrixA = np.array([[3.738, 0.195, 0.275, 0.136],
                        [0.519, 5.002, 0.405, 0.283],
                        [0.306, 0.381, 4.812, 0.418],
                        [0.272, 0.142, 0.314, 3.935]])
    matrixB = np.array([0.815, 0.191, 0.423, 0.352])
    # n = int(input("Введите размерность системы уравнений "))
    # print("Система уравнений будет иметь вид:")
    # for i in range(n):
    #     for j in range(n):
    #         print("a{0}{1}*x{2}".format(i + 1, j + 1, j + 1), end="")
    #         if j != n - 1:
    #             print("+", end="")
    #         else:
    #             print("=", end="")
    #     print("b{0}".format(i + 1))
    # matrixA = np.zeros((n, n))
    # matrixB = np.zeros(n)
    # for i in range(n):
    #     for j in range(n):
    #         matrixA[i][j] = input("a{0}{1}=".format(i + 1, j + 1))
    #     matrixB[i] = input("b{0}=".format(i + 1))
    print("Система уравнений будет иметь вид:")
    printSystemOfEquations(n, matrixA, matrixB)
    squareRootMethod(n, matrixA, matrixB)


