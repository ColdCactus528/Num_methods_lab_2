import matplotlib.pyplot as plt
import numpy as np
import math

def FunctionLab2(x):
    return math.log(1+1/(1+x*x))

def DerivativeRightDifference(min, maximum, step):
    massY = []
    for i in np.arange(min, maximum-step, step):
        massY.append((FunctionLab2(i+step) - FunctionLab2(i)) / (i+step - i))
    massY.append(massY[-1])
    return massY

def DerivateCentralDifference(min, maximum, step):
    massY = []
    massY.append((-3 * FunctionLab2(min) + 4 * FunctionLab2(min+step) - FunctionLab2(min+2*step)) / (2*step))
    for i in np.arange(min+step, maximum-step, step):
        massY.append((FunctionLab2(i+step) - FunctionLab2(i-step)) / (2 * step))
    massY.append((3 * FunctionLab2(maximum) - 4 * FunctionLab2(maximum - step) + FunctionLab2(maximum - 2 * step)) / (2 * step))

    return massY

def Derivate2CentralDifference2Acc(min, maximum, step):
    massY = []
    massY.append((2*FunctionLab2(min) - 5 * FunctionLab2(min+step) + 4 * FunctionLab2(min+2*step) - FunctionLab2(min+3*step)) / (step*step))
    for i in np.arange(min+step, maximum-step, step):
        massY.append((FunctionLab2(i-step) - 2 * FunctionLab2(i) + FunctionLab2(i+step)) / (step*step))
    massY.append((2 * FunctionLab2(maximum) - 5 * FunctionLab2(maximum - step) + 4 * FunctionLab2(maximum - 2 * step) - FunctionLab2(maximum - 3 * step)) / (step * step))
    return massY

def Derivate2CentralDifference4Acc(min, maximum, step):
    massY = []
    for i in np.arange(min, maximum, step):
        massY.append((-FunctionLab2(i-step*2) + 16 * FunctionLab2(i-step) - 30 * FunctionLab2(i) - FunctionLab2(i+step*2) + 16 * FunctionLab2(i+step)) / (12 * step * step))

    return massY

def TrueDerivate1(x):
    return -2*x / (x*x*x*x + 3 * x*x + 2)

def TrueDerivate2(x):
    return (6*x*x*x*x + 6*x*x - 4) / (x*x*x*x + 3 * x*x + 2) / (x*x*x*x + 3 * x*x + 2)

def Err(trueMassive, massive):
    Mass = []
    for i in range(len(massive)):
        Mass.append(abs(massive[i] - trueMassive[i]))

    return max(Mass)

def LogError(min, maximum, Func, flag):
    step = 0.1
    massErr = []
    massStep = []
    while(step > 1e-4):
        massTD = []
        massX = np.arange(min, maximum, step)
        if flag == 1:
            for i in massX:
                massTD.append(TrueDerivate1(i))
        else:
            for i in massX:
                massTD.append(TrueDerivate2(i))


        massCD = Func(min, maximum, step)
        massStep.append(math.log(step))
        massErr.append(math.log(Err(massTD, massCD)))
        step = step / 10

    return massStep, massErr


step = 0.1
min = -5
maximum = 5
massX = np.linspace(min, maximum, int((maximum - min) / step))
massXDCD = np.linspace(min, maximum, int((maximum - min) / step))
massXD2CD4Acc = np.linspace(min, maximum, int((maximum - min) / step))
massXDRD = np.linspace(min, maximum, int((maximum - min) / step))
massFunc = [FunctionLab2(i) for i in massX]
massDRD = DerivativeRightDifference(min, maximum, step)
massDCD = DerivateCentralDifference(min, maximum, step)
massD2CD2Acc = Derivate2CentralDifference2Acc(min, maximum, step)
massD2CD4Acc = Derivate2CentralDifference4Acc(min, maximum, step)

massStep, massErr = LogError(min, maximum, DerivativeRightDifference, 1)

print(len(massX))
print(len(massXDCD))
print(len(massDCD))
print(len(massXDRD))
print(len(massDRD))
plt.title("19 задание")  # заголовок
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.plot(massX, massFunc, label='функция')
plt.plot(massXDRD, massDRD, label='производная по правой разности')
plt.plot(massXDCD, massDCD, label='производная по центральной разности')
plt.plot(massXDCD, massD2CD2Acc, label='2 производная по центральной разности, точность 2')
plt.plot(massXD2CD4Acc, massD2CD4Acc, label='2 производная по центральной разности, точность 4')
# plt.plot(massStep, massErr, label='логарифм для функции ошибки для первой производной по правой разности')
plt.legend()
plt.show()
