# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
 
class splineBase:
    def __init__(sz, a, b, c, d, x):
        sz.a = a
        sz.b = b
        sz.c = c
        sz.d = d
        sz.x = x
 
def createSpline(x, y, n):
    splines = [splineBase(0, 0, 0, 0, 0) for _ in range(0, n)]
    for i in range(0, n):
        splines[i].x = x[i]
        splines[i].a = y[i]
    
    splines[0].c = splines[n - 1].c = 0.0
    alpha = [0.0 for _ in range(0, n - 1)]
    beta  = [0.0 for _ in range(0, n - 1)]
 
    for i in range(1, n - 1):
        hi  = x[i] - x[i - 1]
        hi1 = x[i + 1] - x[i]
        A = hi
        C = 2.0 * (hi + hi1)
        B = hi1
        F = 6.0 * ((y[i + 1] - y[i]) / hi1 - (y[i] - y[i - 1]) / hi)
        z = (A * alpha[i - 1] + C)
        alpha[i] = -B / z
        beta[i] = (F - A * beta[i - 1]) / z
  
    for i in range(n - 2, 0, -1):
        splines[i].c = alpha[i] * splines[i + 1].c + beta[i]
    
    for i in range(n - 1, 0, -1):
        hi = x[i] - x[i - 1]
        splines[i].d = (splines[i].c - splines[i - 1].c) / hi
        splines[i].b = hi * (2.0 * splines[i].c + splines[i - 1].c) / 6.0 + (y[i] - y[i - 1]) / hi
    return splines
 
def getInterpollKoef(splines, x):
    if not splines:
        return None
    
    n = len(splines)
    s = splineBase(0, 0, 0, 0, 0)
    
    if x <= splines[0].x: 
        s = splines[0]
    elif x >= splines[n - 1].x: 
        s = splines[n - 1]
    else: 
        i = 0
        j = n - 1
        while i + 1 < j:
            k = i + (j - i) // 2
            if x <= splines[k].x:
                j = k
            else:
                i = k
        s = splines[j]  
    dx = x - s.x
    return s.a + (s.b + (s.c / 2.0 + s.d * dx / 6.0) * dx) * dx;
    
x = []
y= []
print("Введите количество точек n")
n=int(input())

print("Нажмите 1 если хотите заполнить точки случайными значениями")
print("Нажмите 0 если хотите ввести значения вручную")
print("Нажмите 2 если хотите использовать шаблон")
choose = int(input())
if choose == 0:
    for i in range (n):
        print("x", i + 1,"=")
        xi=int(input())
        x.append(xi) 
        print("y", i + 1, "=")
        yi=int(input())
        y.append(yi)
if choose == 1:
    import numpy as np
    import random
    x = np.random.randint(0, 200, n)
    x.sort()
    y = np.random.randint(0, 200, n)
if choose == 2:
    x = [25, 67, 83, 102, 131, 167, 181, 195]
    y = [58, 47, 35, 24, 33, 78, 50, 112]
    print("x  = ", x)
    print("y = ", y)
mini = min(x)
maxi = max(x)

spline = createSpline(x, y, len(x))


w = list()
z = list()
for i in range(mini, maxi+1):
    w.append(i)
    z.append(getInterpollKoef(spline, i))

plt.ylabel('Ось Y')
plt.xlabel('Ось X')
plt.plot(w,z)
plt.scatter(x,y,color = "r")
plt.show()

with open("spl.txt", "w") as my_file:
    for k in z:
        my_file.write('%s\n' % k)
    my_file.close()