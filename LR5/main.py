import math
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


# РЕШЕНИЕ УРАВНЕНИЯ
a = 0.3
b = 0.5
eps1 = 0.0001


def f(x):
    if x != 0:
        return x - math.sin(1 / x)
    else:
        return None


def df(x):
    if x != 0:
        return 1 + math.cos(1 / x) / (x ** 2)
    else:
        return None


def chordMethod():
    prevX = a
    nextX = b
    while abs(nextX - prevX) > eps1:
        print(prevX, nextX)
        temp = prevX - (f(prevX) * (nextX - prevX)) / (f(nextX) - f(prevX))
        prevX = nextX
        nextX = temp
    return nextX


def tangentMethod():
    prevX = a
    nextX = (a + b) / 2
    while abs(nextX - prevX) > eps2:
        prevX = nextX
        nextX = prevX - f(prevX) / df(prevX)
    return nextX
# РЕШЕНИЕ УРАВНЕНИЯ


# РЕШЕНИЕ СИСТЕМЫ УРАВНЕНИЙ
x = 0.8
y = 0.4
eps2 = 0.0001


def system(vectorX):
    syst = np.empty(2)
    syst[0] = math.sin(vectorX[0] + vectorX[1]) - 1.1 * vectorX[0]
    syst[1] = 0.8 * vectorX[0] ** 2 + 2 * vectorX[1] ** 2 - 1
    return syst


def coordianatesOfY(x):
    syst = np.empty(3)
    if abs(1.1 * x) > 1:
        syst[0] = None
    else:
        syst[0] = math.asin(1.1 * x) - x
    if 0.5 - 0.4 * x ** 2 < 0:
        syst[1] = None
        syst[2] = None
    else:
        syst[1] = math.sqrt(0.5 - 0.4 * x ** 2)
        syst[2] = -math.sqrt(0.5 - 0.4 * x ** 2)
    return syst


def pronouncedSystem(vectorX):
    syst = np.empty(2)
    syst[0] = math.sin(vectorX[0] + vectorX[1]) / 1.1
    syst[1] = math.sqrt(0.5 - 0.4 * vectorX[0] ** 2)
    return syst


def jacobian(x):
    w = np.empty((2, 2))
    w[0][0] = math.cos(x[0] + x[1]) - 1
    w[0][1] = math.cos(x[0] + x[1])
    w[1][0] = 2 * x[0]
    w[1][1] = 2 * x[1]
    return w


def iterationMethod():
    prevVectorX = np.array([x, y])
    nextVectorX = pronouncedSystem(prevVectorX)
    differenceX = la.norm(nextVectorX - prevVectorX, np.inf)
    numOfIterations = 1
    while differenceX > eps2:
        prevVectorX = nextVectorX
        nextVectorX = pronouncedSystem(prevVectorX)
        differenceX = la.norm(nextVectorX - prevVectorX, np.inf)
        numOfIterations += 1
    return nextVectorX, numOfIterations


def newtonMethod():
    prevVectorX = np.array([x, y])
    nextVectorX = prevVectorX - np.matmul(la.inv(jacobian(prevVectorX)), system(prevVectorX))
    differenceX = la.norm(nextVectorX - prevVectorX, np.inf)
    numOfIterations = 1
    while differenceX > eps2:
        prevVectorX = nextVectorX
        nextVectorX = prevVectorX - np.matmul(la.inv(jacobian(prevVectorX)), system(prevVectorX))
        differenceX = la.norm(nextVectorX - prevVectorX, np.inf)
        numOfIterations += 1
    return nextVectorX, numOfIterations


def modifiedNewtonMethod():
    prevVectorX = np.array([x, y])
    inverseJacobian = la.inv(jacobian(prevVectorX))
    nextVectorX = prevVectorX - np.matmul(inverseJacobian, system(prevVectorX))
    differenceX = la.norm(nextVectorX - prevVectorX, np.inf)
    numOfIterations = 1
    while differenceX > eps2:
        prevVectorX = nextVectorX
        nextVectorX = prevVectorX - np.matmul(inverseJacobian, system(prevVectorX))
        differenceX = la.norm(nextVectorX - prevVectorX, np.inf)
        numOfIterations += 1
    return nextVectorX, numOfIterations
# РЕШЕНИЕ СИСТЕМЫ УРАВНЕНИЙ


def main():
    print("Заданное уравнение имеет вид:")
    print("x - sin(1 / x) = 0", end="\n\n")
    listX = np.arange(a, b, 0.01)
    listY = np.array([f(x) for x in listX])
    plt.plot(listX, listY)
    plt.axhline(0, color='black')
    plt.grid()
    plt.show()
    print("Решение, полученное методом хорд:")
    print(chordMethod(), end="\n\n")
    print("Решение, полученное методом касательных(метод Ньютона):")
    print(tangentMethod(), end="\n\n")

    print("Заданная система имеет вид:")
    print("sin(x + y) - x = 1")
    print("x^2 + y^2 = 1", end="\n\n")
    listX = np.arange(0, 1, 0.01)
    listY = np.empty((3, 0))
    for x in listX:
        listY = np.column_stack((listY, coordianatesOfY(x)))
    plt.plot(listX, listY[0])
    plt.plot(listX, listY[1], 'g')
    plt.plot(listX, listY[2], 'g')
    plt.grid()
    plt.show()
    print("Метод простых итераций:")
    print("x(k + 1) = Ф(x(k))", end="\n\n")
    print("Решение, полученное методом простых итераций:")
    ans, numOfIter1 = iterationMethod()
    print("x = {0}, y = {1}".format(ans[0], ans[1]), end="\n\n")
    print("Метод Ньютона:")
    print("x(k + 1) = x(k) - [J(x(k))]^(-1) * F(x(k))", end="\n\n")
    print("Решение, полученное методом Ньютона:")
    ans, numOfIter2 = newtonMethod()
    print("x = {0}, y = {1}".format(ans[0], ans[1]), end="\n\n")
    print("Модифицированный метод Ньютона:")
    print("x(k + 1) = x(k) - [J(x(0))]^(-1) * F(x(k))", end="\n\n")
    print("Решение, полученное модифицированным методом Ньютона:")
    ans, numOfIter3 = modifiedNewtonMethod()
    print("x = {0}, y = {1}".format(ans[0], ans[1]), end="\n\n")
    print("Количество итераций в методе простых итераций - {0}, в методе Ньютона - {1}, "
          "в модифицированном методе Ньютона - {2}".format(numOfIter1, numOfIter2, numOfIter3))


if __name__ == '__main__':
    main()
