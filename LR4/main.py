import numpy as np
import math


def squareRootMethod(n, matrixA, matrixB):
    if not isSymmetricalMatrix(n, matrixA):
        matrixA = np.matmul(np.transpose(matrixA), matrixA)
    eigenvectors = np.eye(n)
    while not isTrueMatrix(n, matrixA):
        maxElem = 0
        maxI = -1
        maxJ = -1
        for i in range(n):
            for j in range(i + 1, n):
                if abs(matrixA[i][j]) > maxElem:
                    maxElem = abs(matrixA[i][j])
                    maxI = i
                    maxJ = j
        angel = (1 / 2) * math.atan(2 * matrixA[maxI][maxJ] / (matrixA[maxI][maxI] - matrixA[maxJ][maxJ]))
        matrixU = np.eye(n)
        matrixU[maxJ][maxJ] = matrixU[maxI][maxI] = math.cos(angel)
        matrixU[maxI][maxJ] = -math.sin(angel)
        matrixU[maxJ][maxI] = math.sin(angel)
        print(matrixU)
        print()
        matrixA = np.matmul(np.matmul(np.transpose(matrixU), matrixA), matrixU)
        print(matrixA)
        eigenvectors = np.matmul(eigenvectors, matrixU)
    print("Собственные значения матрицы A:")
    print(np.diag(matrixA), end="\n\n")
    print("Собственные векторы матрицы A:")
    print(eigenvectors)


def isTrueMatrix(n, matrixU):
    for i in range(n):
        for j in range(n):
            if i != j and abs(matrixU[i][j]) >= 0.01:
                return False
    return True


def printSystemOfEquations(n, matrixA, matrixB):
    for i in range(n):
        for j in range(n):
            print("{0:+}*x{1}".format(matrixA[i][j], j + 1), end="")
            if j == n - 1:
                print("=", end="")
        print(matrixB[i])
    print()


def isSymmetricalMatrix(n, matrixA):
    for i in range(n):
        for j in range(i + 1, n):
            if matrixA[i][j] != matrixA[j][i]:
                return False
    return True


def main():
    n = 3
    matrixA = np.array([[2, -1, -1],
                        [-1, 3, -1],
                        [-1, -1, 4]], dtype=float)
    matrixB = np.array([0.815, 0.191, 0.423])
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
    # extendedMatrix = np.column_stack((matrixA, matrixB))
    # if la.matrix_rank(extendedMatrix) != la.matrix_rank(matrixA):
    #     print("Система не имеет решений")
    #     return
    # elif la.matrix_rank(matrixA) != n:
    #     print("Система имеет бесконечно много решений")
    #     return
    print("Система уравнений будет иметь вид:")
    printSystemOfEquations(n, matrixA, matrixB)
    squareRootMethod(n, matrixA, matrixB)


if __name__ == '__main__':
    main()


