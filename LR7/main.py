import math
import numpy as np
import matplotlib.pyplot as plt


# Численное интегрирование
a = 1
b = 2.2
eps = 0.001


def f(x):
    if x == 0:
        return None
    else:
        return 1 / (x ** 3 + x)


def integral(x):
    return math.log(abs(x)) - 1 / 2 * math.log(x ** 2 + 1)


def theSecondDerivative(x):
    if x == 0:
        return None
    else:
        return (12 * (x ** 4) + 6 * (x ** 2) + 1) / ((x ** 3 + x) ** 3)


def calculateStepForIntegration():
    maxSecondDerivative = 0.
    x = a
    while x <= b:
        maxSecondDerivative = max(maxSecondDerivative, abs(theSecondDerivative(x)))
        x += 0.01
    n = int(math.sqrt((maxSecondDerivative * (b - a) ** 3) / (12 * eps)) / 4) * 4
    return (b - a) / n


def trapezoidIntegral(h):
    integral = (f(a) + f(b)) / 2
    n = round(abs(b - a) / h)
    x = a
    for i in range(n - 1):
        x += h
        integral += f(x)
    integral *= h
    return integral


def simpsonIntegral(h):
    integral = f(a) + f(b)
    n = int(abs(b - a) / (2 * h))
    x = a + h
    for i in range(n):
        integral += 4 * f(x)
        x += 2 * h
    x = a + 2 * h
    for i in range(n - 1):
        integral += 2 * f(x)
        x += 2 * h
    integral *= h / 3
    return integral


def newtonLeibnizIntegral():
    return integral(b) - integral(a)
# Численное интегрирование


# Задача Коши
a1 = 1
b1 = 2.6
y0 = 1
eps1 = 0.0001
n = 0
listX = np.empty(0)
exactSolutionSequence = np.empty(0)
eulerSequence = np.empty(0)
fourthOrderRungeKuttaSequence = np.empty(0)
adamsSequence = np.empty(0)


def derivative(x, y):
    if x == 0:
        return None
    else:
        return (y ** 2 * math.log(x) - y) / x


def solveFunction(x):
    return 1 / (math.log(x) + 1)


def calculateStepForKoshiProblem():
    global h, n
    n = math.ceil((b1 - a1) / (eps1 ** 0.25))
    n = n if n % 2 == 0 else n + 1
    h = (b - a) / n
    y2 = helperRungeKuttaMethod(h)
    y_2 = helperRungeKuttaMethod(2 * h)
    while (1 / 15) * abs(y2 - y_2) < eps1:
        h *= 2
        y2 = helperRungeKuttaMethod(h)
        y_2 = helperRungeKuttaMethod(2 * h)
    while (1 / 15) * abs(y2 - y_2) >= eps1:
        h /= 2
        y2 = helperRungeKuttaMethod(h)
        y_2 = helperRungeKuttaMethod(2 * h)
    n = int((b1 - a1) / h) + 1
    return h


def helperRungeKuttaMethod(h):
    x, y = a, y0
    for i in range(2):
        k1 = derivative(x, y)
        k2 = derivative(x + h / 2, y + (h / 2) * k1)
        k3 = derivative(x + h / 2, y + (h / 2) * k2)
        k4 = derivative(x + h, y + h * k3)
        y = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        x += h
    return y


def defineListX(h):
    global listX
    listX = np.empty(n)
    x = a1
    for i in range(n):
        listX[i] = x
        x += h


def exactSolutionMethod(h):
    global exactSolutionSequence
    exactSolutionSequence = np.empty(n)
    for i in range(n):
        exactSolutionSequence[i] = solveFunction(listX[i])
    plt.grid()
    plt.plot(listX, exactSolutionSequence)
    plt.show()


def eulerMethod(h):
    global eulerSequence
    eulerSequence = np.empty(n)
    eulerSequence[0] = y0
    for i in range(1, n):
        eulerSequence[i] = eulerSequence[i - 1] + h * derivative(listX[i - 1], eulerSequence[i - 1])
    printTable(eulerSequence, "Метод Эйлера")
    plt.grid()
    plt.plot(listX, eulerSequence)
    plt.show()


def fourthOrderRungeKuttaMethod(h):
    global fourthOrderRungeKuttaSequence
    fourthOrderRungeKuttaSequence = np.empty(n)
    fourthOrderRungeKuttaSequence[0] = y0
    for i in range(1, n):
        x = listX[i - 1]
        y = fourthOrderRungeKuttaSequence[i - 1]
        k1 = derivative(x, y)
        k2 = derivative(x + h / 2, y + (h / 2) * k1)
        k3 = derivative(x + h / 2, y + (h / 2) * k2)
        k4 = derivative(x + h, y + h * k3)
        fourthOrderRungeKuttaSequence[i] = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    plt.grid()
    plt.plot(listX, fourthOrderRungeKuttaSequence)
    plt.show()


def adamsMethod(h):
    global adamsSequence
    adamsSequence = np.empty(n)
    adamsSequence[0] = y0
    adamsSequence[1] = eulerSequence[1]
    for i in range(2, n):
        adamsSequence[i] = adamsSequence[i - 1] + (h / 2) * (3 * derivative(listX[i - 1], adamsSequence[i - 1]) -
            derivative(listX[i - 2], adamsSequence[i - 2]))
    plt.grid()
    plt.plot(listX, adamsSequence)
    plt.show()


def printTable(sequence, nameOfTable):
    maxDiff = 0
    x = 0
    print(nameOfTable)
    print("x           y           y'          /\\=|y-y'|")
    for i in range(n):
        diff = abs(exactSolutionSequence[i] - sequence[i])
        print("{0:.4f}      {1:.4f}      {2:.4f}      {3:.4f}".format(listX[i], exactSolutionSequence[i], sequence[i], diff))
        if maxDiff < diff:
            maxDiff = diff
            x = listX[i]
    print("Максимум модулей отклонений равен {0} при x = {1}".format(maxDiff, x))
    print()


def printAllGraphics():
    plt.grid()
    plt.plot(listX, eulerSequence)
    plt.plot(listX, fourthOrderRungeKuttaSequence)
    plt.plot(listX, adamsSequence)
    plt.show()
# Задача Коши


def main():
    h = calculateStepForIntegration()
    print("Шаг интегрирования h для вычисления интеграла по формуле трапеций с точностью {0}: {1:.3f}".format(eps, h),
          end="\n\n")
    print("Интеграл по формуле трапеций с шагом h: {0}, погрешность равна {1}".format(trapezoidIntegral(h),
          abs(newtonLeibnizIntegral() - trapezoidIntegral(h))), end="\n\n")
    print("Интеграл по формуле трапеций с шагом 2h: {0}, погрешность равна {1}".format(trapezoidIntegral(2 * h),
          abs(newtonLeibnizIntegral() - trapezoidIntegral(2 * h))), end="\n\n")
    print("Интеграл по формуле Симпсона с шагом h: {0}, погрешность равна {1}".format(simpsonIntegral(h),
          abs(newtonLeibnizIntegral() - simpsonIntegral(h))), end="\n\n")
    print("Интеграл по формуле Симпсона с шагом 2h: {0}, погрешность равна {1}".format(simpsonIntegral(2 * h),
          abs(newtonLeibnizIntegral() - simpsonIntegral(2 * h))), end="\n\n")
    print("Интеграл по формуле Ньютона–Лейбница: {0}".format(newtonLeibnizIntegral()), end="\n\n")
    h1 = calculateStepForKoshiProblem()
    print("Шаг интегрирования h для решения задачи Коши методом Рунге-Кутта(4) с точностью {0}: {1}".format(eps1,
          h1), end="\n\n")
    defineListX(h1)
    exactSolutionMethod(h1)
    eulerMethod(h1)
    fourthOrderRungeKuttaMethod(h1)
    adamsMethod(h1)
    printAllGraphics()


if __name__ == '__main__':
    main()
