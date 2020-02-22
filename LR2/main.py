import numpy as np
import numpy.linalg as la
import math


def simpleIterationMethod(n, matrixA, matrixB):
    print("Метод простых итераций")
    for i in range(n):
        coefficient = matrixA[i][i]
        matrixA[i] /= coefficient
        matrixB[i] /= coefficient
    matrixA = np.eye(n) - matrixA
    print("Вид системы для решения методом простых итераций")
    for i in range(n):
        print("x{0}(k+1)=".format(i + 1), end="")
        for j in range(n):
            print("{0}*x{1}(k)".format(matrixA[i][j], j + 1), end="")
            if j != n - 1 and matrixA[i][j + 1] >= 0:
                print("+", end="")
        if matrixB[i] >= 0:
            print("+", end="")
        print(matrixB[i])
    print()
    if la.norm(matrixA, np.inf) >= 1 and la.norm(matrixA, 1) >= 1 and la.norm(matrixA) >= 1:
        print("Система уравнений не решается через итерационный метод", end="\n\n")
        return
    matrixPrevX = np.copy(matrixB)
    matrixX = np.matmul(matrixA, matrixPrevX) + matrixB
    normA = la.norm(matrixA, np.inf)
    normX = la.norm(matrixX - matrixPrevX, np.inf)
    numOfIteration = math.ceil(math.log((0.01 * (1 - normA)) / normX, normA))
    print("Количество итераций, необходимых для достижения точности в 0.01: {0}".format(numOfIteration), end="\n\n")
    numOfIteration -= 1
    for i in range(numOfIteration):
        matrixPrevX = matrixX
        matrixX = np.matmul(matrixA, matrixPrevX) + matrixB
    print("Решение:")
    for i in range(n):
        print("x{0}={1}".format(i + 1, matrixX[i]), end=" ")
    print(end="\n\n")


def seidelMethod(n, matrixA, matrixB):
    print("Метод Зейделя")
    for i in range(n):
        coefficient = matrixA[i][i]
        matrixA[i] /= coefficient
        matrixB[i] /= coefficient
    matrixA = np.eye(n) - matrixA
    print("Вид системы для решения методом Зейделя")
    for i in range(n):
        print("x{0}(k+1)=".format(i + 1), end="")
        isTheFirst = True
        for j in range(n):
            if j != i:
                if not isTheFirst and matrixA[i][j] >= 0:
                    print("+", end="")
                isTheFirst = False
                print("{0}*x{1}".format(matrixA[i][j], j + 1), end="")
                if j > i:
                    print("(k)", end="")
                else:
                    print("(k+1)", end="")
        if matrixB[i] >= 0:
            print("+", end="")
        print(matrixB[i])
    print()
    if la.norm(matrixA, np.inf) >= 1 and la.norm(matrixA, 1) >= 1:
        print("Система уравнений не решается через метод Зейделя", end="\n\n")
        return
    matrixPrevX = np.copy(matrixB)
    matrixX = np.copy(matrixB)
    for i in range(n):
        for j in range(i):
            matrixX[i] += matrixX[j] * matrixA[i][j]
        for j in range(i + 1, n):
            matrixX[i] += matrixPrevX[j] * matrixA[i][j]
    while abs(max(np.fabs(matrixX - matrixPrevX))) > 0.01:
        matrixPrevX = np.copy(matrixX)
        matrixX = np.copy(matrixB)
        for i in range(n):
            for j in range(i):
                matrixX[i] += matrixX[j] * matrixA[i][j]
            for j in range(i, n):
                matrixX[i] += matrixPrevX[j] * matrixA[i][j]
    print("Решение:")
    for i in range(n):
        print("x{0}={1}".format(i + 1, matrixX[i]), end=" ")


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
    for i in range(n):
        for j in range(n):
            print("{0}*x{1}".format(matrixA[i][j], j + 1), end="")
            if j != n - 1 and matrixA[i][j + 1] >= 0:
                print("+", end="")
            else:
                print("=", end="")
        print(matrixB[i])
    print()
    extendedMatrix = np.column_stack((matrixA, matrixB))
    if la.matrix_rank(extendedMatrix) != la.matrix_rank(matrixA):
        print("Система не имеет решений")
    elif la.matrix_rank(matrixA) != n:
        print("Система имеет бесконечно много решений")
    else:
        simpleIterationMethod(n, matrixA, matrixB)
        seidelMethod(n, matrixA, matrixB)
