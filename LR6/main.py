import numpy as np
import matplotlib.pyplot as plt


n = 7
listX = [0.231, 0.848, 1.322, 2.224, 2.892, 3.333, 3.789]
listY = [-2.748, -3.225, -3.898, -5.908, -6.506, -7.236, -8.111]
finalDifferences = np.empty((n, n))
dividedDifferences = np.empty((n, n))
piecewiseLinearCoefficients = np.empty((2, n))
piecewiseQuadraticCoefficients = np.empty((3, int((n + 1) / 2)))
cubicSplineCoefficients = np.empty((4, n))


def lagrangePolynomial(x):
    result = 0
    for i in range(n):
        mul = 1
        for j in range(n):
            if i != j:
                mul *= (x - listX[j]) / (listX[i] - listX[j])
        result += listY[i] * mul
    return result


def printLagrangePolynomial():
    print("y=", end="")
    for i in range(n):
        mul = listY[i]
        for j in range(n):
            if i != j:
                mul /= (listX[i] - listX[j])
        print("{0:+}*".format(mul), end="")
        for j in range(n):
            if i != j:
                print("(x{0:+})".format(-listX[j]), end="")
    print(end="\n\n")
    xCoordinates = np.arange(-2, 5, 0.01)
    yCoordinates = np.array([lagrangePolynomial(x) for x in xCoordinates])
    plt.plot(xCoordinates, yCoordinates)
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.grid()
    plt.show()


def calculateFinalDifferences():
    global finalDifferences
    finalDifferences[0] = listY
    for i in range(1, n):
        for j in range(n - i):
            finalDifferences[i][j] = finalDifferences[i - 1][j + 1] - finalDifferences[i - 1][j]
    print("Таблица конечных разностей:")
    for i in range(0, n - 1):
        print("x({0})={1:+.3f}".format(i, listX[i]), end="   ")
        print("y({0})={1:+.3f}".format(i, listY[i]), end="   ")
        for j in range(1, n - i):
            print("/\\({0})y({1})={2:+.3f}".format(j, i, finalDifferences[j][i]), end="   ")
        print()
    print("x({0})={1:+.3f}".format(n - 1, listX[n - 1]), end="   ")
    print("y({0})={1:+.3f}".format(n - 1, listY[n - 1]), end="\n\n")


def calculateDividedDifferences():
    global dividedDifferences
    dividedDifferences[0] = listY
    for i in range(1, n):
        for j in range(n - i):
            dividedDifferences[i][j] = (dividedDifferences[i - 1][j + 1] - dividedDifferences[i - 1][j]) / (listX[j + i] - listX[j])
    print("Таблица разделенных разностей:")
    for i in range(0, n - 1):
        print("x({0})={1:+.3f}".format(i, listX[i]), end="   ")
        print("y({0})={1:+.3f}".format(i, listY[i]), end="   ")
        for j in range(1, n - i):
            print("f({0})y({1})={2:+.3f}".format(j, i, dividedDifferences[j][i]), end="   ")
        print()
    print("x({0})={1:+.3f}".format(n - 1, listX[n - 1]), end="   ")
    print("y({0})={1:+.3f}".format(n - 1, listY[n - 1]), end="\n\n")


def newtonPolynomial(x):
    result = listY[0]
    for i in range(1, n):
        mul = 1
        for j in range(i):
            mul *= (x - listX[j])
        result += dividedDifferences[i][0] * mul
    return result


def printNewtonPolynomial():
    print("y={0}".format(listY[0]), end="")
    for i in range(1, n):
        print("{0:+}*".format(dividedDifferences[i][0]), end="")
        for j in range(i):
            print("(x{0:+})".format(-listX[j]), end="")
    print(end="\n\n")
    xCoordinates = np.arange(-2, 5, 0.01)
    yCoordinates = np.array([newtonPolynomial(x) for x in xCoordinates])
    plt.plot(xCoordinates, yCoordinates)
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.grid()
    plt.show()


def calculatePiecewiseLinearCoefficients():
    for i in range(1, n):
        piecewiseLinearCoefficients[0][i] = (listY[i] - listY[i - 1]) / (listX[i] - listX[i - 1])
        piecewiseLinearCoefficients[1][i] = (listY[i - 1] * listX[i] - listY[i] * listX[i - 1]) / (listX[i] - listX[i - 1])


def piecewiseLinearApproximation(x):
    if x < listX[0]:
        return None
    if x > listX[n - 1]:
        return None
    for i in range(1, n):
        if x <= listX[i]:
            return piecewiseLinearCoefficients[0][i] * x + piecewiseLinearCoefficients[1][i]


def printPiecewiseLinearApproximation():
    calculatePiecewiseLinearCoefficients()
    print("Кусочно-линейная аппроксимация")
    for i in range(1, n):
        print("y ={0:.3f}*x{1:+.3f},   {2} <= x <= {3}".format(piecewiseLinearCoefficients[0][i],
                                                             piecewiseLinearCoefficients[1][i], listX[i - 1], listX[i]))
    print()
    xCoordinates = np.arange(-2, 5, 0.01)
    yCoordinates = np.array([piecewiseLinearApproximation(x) for x in xCoordinates])
    plt.plot(xCoordinates, yCoordinates)
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.grid()
    plt.show()


def calculatePiecewiseQuadraticCoefficients():
    for i in range(1, int((n + 1) / 2)):
        piecewiseQuadraticCoefficients[0][i] = (listY[2 * i] * (listX[2 * i - 1] - listX[2 * i - 2]) +
            listY[2 * i - 1] * (listX[2 * i - 2] - listX[2 * i]) + listY[2 * i - 2] * (listX[2 * i] - listX[2 * i - 1])) / \
            ((listX[2 * i] - listX[2 * i - 1]) * (listX[2 * i] - listX[2 * i - 2]) * (listX[2 * i - 1] - listX[2 * i - 2]))
        piecewiseQuadraticCoefficients[1][i] = (listY[2 * i] - listY[2 * i - 1] - piecewiseQuadraticCoefficients[0][i] *
            (listX[2 * i] ** 2 - listX[2 * i - 1] ** 2)) / (listX[2 * i] - listX[2 * i - 1])
        piecewiseQuadraticCoefficients[2][i] = listY[2 * i] - piecewiseQuadraticCoefficients[0][i] * \
            (listX[2 * i] ** 2) - piecewiseQuadraticCoefficients[1][i] * listX[2 * i]


def piecewiseQuadraticApproximation(x):
    if x < listX[0]:
        return None
    if x > listX[n - 1 if n % 2 == 1 else n - 2]:
        return None
    for i in range(1, int((n + 1) / 2)):
        if x <= listX[2 * i]:
            return piecewiseQuadraticCoefficients[0][i] * (x ** 2) + piecewiseQuadraticCoefficients[1][i] * x + \
                   piecewiseQuadraticCoefficients[2][i]


def printPiecewiseQuadraticApproximation():
    calculatePiecewiseQuadraticCoefficients()
    print("Кусочно-квадратичная аппроксимация")
    for i in range(1, int((n + 1) / 2)):
        print("y = {0:.3f}*x^2{1:+.3f}*x{2:+.3f},   {3} <= x <= {4}".format(piecewiseQuadraticCoefficients[0][i],
               piecewiseQuadraticCoefficients[1][i], piecewiseQuadraticCoefficients[2][i], listX[2 * i - 2], listX[2 * i]))
    print()
    xCoordinates = np.arange(-2, 5, 0.01)
    yCoordinates = np.array([piecewiseQuadraticApproximation(x) for x in xCoordinates])
    plt.plot(xCoordinates, yCoordinates)
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.grid()
    plt.show()


def calculateCubicSplineCoefficients():
    global cubicSplineCoefficients
    h = np.empty(n)
    l = np.empty(n)
    for i in range(1, n):
        h[i] = listX[i] - listX[i - 1]
        l[i] = (listY[i] - listY[i - 1]) / h[i]
    sigma = np.empty(n)
    sigma[1] = (-1.0 / 2) * h[2] / (h[1] + h[2])
    lambd = np.empty(n)
    lambd[1] = (3.0 / 2) * (l[2] - l[1]) / (h[1] + h[2])
    for i in range(3, n):
        sigma[i - 1] = -h[i] / (2 * h[i - 1] + 2 * h[i] + h[i - 1] * sigma[i - 2])
        lambd[i - 1] = (2 * l[i] - 3 * l[i - 1] - h[i - 1] * lambd[i - 2]) / (2 * h[i - 1] + 2 * h[i] + h[i - 1] * sigma[i - 2])
    cubicSplineCoefficients[2][n - 1] = 0
    cubicSplineCoefficients[2][0] = 0
    for i in range(n - 1, 1, -1):
        cubicSplineCoefficients[2][i - 1] = sigma[i - 1] * cubicSplineCoefficients[2][i] + lambd[i - 1]
    for i in range(1, n):
        cubicSplineCoefficients[0][i] = listY[i]
        cubicSplineCoefficients[1][i] = l[i] + (2.0 / 3) * cubicSplineCoefficients[2][i] * h[i] + (1.0 / 3) * h[i] * \
            cubicSplineCoefficients[2][i - 1]
        cubicSplineCoefficients[3][i] = (cubicSplineCoefficients[2][i] - cubicSplineCoefficients[2][i - 1]) / (3 * h[i])


def cubicSpline(x):
    if x < listX[0]:
        return None
    if x > listX[n - 1]:
        return None
    for i in range(1, n):
        if x <= listX[i]:
            return cubicSplineCoefficients[0][i] + cubicSplineCoefficients[1][i] * (x - listX[i]) + \
                        cubicSplineCoefficients[2][i] * ((x - listX[i]) ** 2) + cubicSplineCoefficients[3][i] * \
                        ((x - listX[i]) ** 3)


def printCubicSpline():
    calculateCubicSplineCoefficients()
    print("Кубическая аппроксимация")
    for i in range(1, n):
        print("y = {0:.3f}{1:+.3f}(x{6:+.3f}){2:+.3f}(x{6:+.3f})^2{3:+.3f}(x{6:+.3f})^3,   {4} <= x <= {5}".format(
            cubicSplineCoefficients[0][i], cubicSplineCoefficients[1][i], cubicSplineCoefficients[2][i],
            cubicSplineCoefficients[3][i], listX[i - 1], listX[i], -listX[i - 1]))
    print()
    xCoordinates = np.arange(-2, 5, 0.01)
    yCoordinates = np.array([cubicSpline(x) for x in xCoordinates])
    plt.plot(xCoordinates, yCoordinates)
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.grid()
    plt.show()


def printAllGraphics():
    xCoordinates = np.arange(0, 4, 0.01)
    yCoordinates1 = np.array([lagrangePolynomial(x) for x in xCoordinates])
    yCoordinates2 = np.array([newtonPolynomial(x) for x in xCoordinates])
    yCoordinates3 = np.array([piecewiseLinearApproximation(x) for x in xCoordinates])
    yCoordinates4 = np.array([piecewiseQuadraticApproximation(x) for x in xCoordinates])
    yCoordinates5 = np.array([cubicSpline(x) for x in xCoordinates])
    plt.plot(xCoordinates, yCoordinates1)
    plt.plot(xCoordinates, yCoordinates2)
    plt.plot(xCoordinates, yCoordinates3)
    plt.plot(xCoordinates, yCoordinates4)
    plt.plot(xCoordinates, yCoordinates5)
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.grid()
    plt.show()


def main():
    print("Промежуточные значения:")
    for i in range(n):
        print("f({0})={1}".format(listX[i], listY[i]), end="   ")
    print(end="\n\n")
    print("Интерполяционный многочлен Лагранжа:")
    printLagrangePolynomial()
    print("L(x1+x2)={0}".format(lagrangePolynomial(listX[1] + listX[2])), end="\n\n")
    calculateFinalDifferences()
    calculateDividedDifferences()
    print("Полином Ньютона:")
    printNewtonPolynomial()
    print("N(x1+x2)={0}".format(newtonPolynomial(listX[1] + listX[2])), end="\n\n")
    printPiecewiseLinearApproximation()
    printPiecewiseQuadraticApproximation()
    printCubicSpline()
    printAllGraphics()


if __name__ == '__main__':
    main()
