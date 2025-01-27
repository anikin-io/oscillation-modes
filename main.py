import numpy as np
from math import sqrt
import matplotlib.pyplot as plt

N = 10000 # число точек в реализации
dt = 0.01 # шаг интегрирования

#Константы
m = 5 #масса
r = 1 #длина
k = 2 #кф сопротивления
w = pow((k/m), 0.5)
sigma = r/(2*m) #кф затухания

#Инициализация массивов значений N-ой размерности
x = np.zeros(N)
y = np.zeros(N)
sol = []

#Метод Эйлера
def findSol(baseX,baseY,h):
    #начальные условия
    x[0] = baseX
    y[0] = baseY

    for i in range (1,N):
        x[i] = x[i-1] + h * y[i-1]
        y[i] = y[i-1] + h * (-2*sigma*y[i-1] - x[i-1]*pow(w,2))
    return [x,y]

#Вводим стартовые значения
x0,y0 = float(input('Введите стартовое значение по оси X\n')),float(input('Введите стартовое значение по оси Y\n'))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
T=np.arange(0,N*dt,dt)

#Выводим фазовый портрет
sol =findSol(x0,y0,dt)
x_sol,y_sol = sol[0],sol[1]
ax.plot(x_sol,y_sol)
plt.title("Фазовый портрет: Пружинный маятник")

#x(t)
fig2 = plt.figure()
b = fig2.add_subplot()
b.plot(T,x_sol)
plt.xlabel("t")
plt.ylabel("x(t)")

#Исследование расхождений реализаций по среднеквадратическому отклонению
D = np.zeros(20)
dt_mas = np.zeros(20)
j=0
for k in range(0,20):
    x = np.zeros(N)
    y = np.zeros(N)
    sum = 0
    x0 = x0 + 1
    sol_temp =findSol(x0,y0,dt)
    x_temp,y_temp = sol_temp[0],sol_temp[1]
    for i in range(0,N):
        sum = sum + (x_sol[i]-x_temp[i]) * (x_sol[i]-x_temp[i])
    D[k] = sqrt((sum/N))
    dt_mas[k] = x0
fig5 = plt.figure()
e = fig5.add_subplot()
e.plot(dt_mas,D)
plt.xlabel("Начальное условие x0")
plt.ylabel("Среднеквадратическое отклонение")
plt.show()

#Построение графика АКФ
def akf(x_akf, N):
    r = np.zeros(N)
    k_array = np.arange(1 - N / 2, N / 2 + 1, 1, 'int')
    max = 0

    for j in range(N):
        k = k_array[j]
        if k > 0:
            for i in range(1, N - k + 1):
                r[j] += x_akf[i] * x_akf[i + k - 1]
        elif k == 0:
            r[j] = 1
        elif k < 0:
            for i in range(1, N + k):
                r[j] += x_akf[i] * x_akf[(i - k) - 1]

        if r[j] > max:
            max = r[j]

    for j in range(N):
        if k_array[j] != 0:
            r[j] /= max

    plt.plot(k_array, r)
    plt.xlabel("k")
    plt.ylabel("r")
    plt.title("Коррелограмма")
    plt.show()

akf(x_sol,N)
