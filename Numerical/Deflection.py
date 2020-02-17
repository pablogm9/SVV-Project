# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 13:14:03 2020

@author: Marni
"""

from math import *
import numpy as np
import pandas

x1 = 0.125
x2 = 0.498
x3 = 1.494
La = 1.611
x = 0
dx = 0.01
E = 1000
Ixx = 10
R1y = 2
R2y = 3
R3y = 4
Py = 5
Mx = 10
xa = 0.245

v = []
w = []

q_int = 1

def Macauley(y):
    if y<0:
        y = 0
    return y

def v_prime(x):
    v = -1/(E*Ixx)*(q_int + R1y/2*(Macauley(x-x1)**2) + R2y/2*(Macauley(x-x2)**2) - Py/2*Macauley(x-(x2+xa/2))**2 - Ay/2*Macauley(x-(x2-xa/2))**2 + R3y/2*(Macauley(x-x3)**2) + C1)
    return v

def w_prime(x):
    w = 1/(E*Iyy*(R1x/2*(Macauley(x-x1)**2) + R2x/2*(Macauley(x-x2)**2) - Px/2*Macauley(x-(x2+xa/2))**2 - Ax/2*Macauley(x-(x2-xa/2))**2 + R3x/2*(Macauley(x-x3)**2 + C1)
    return w

def v(x):
    v = -1/(E*Ixx)*(q_int + R1y/6*(Macauley(x-x1)**3) + R2y/6*(Macauley(x-x2)**3) - Py/6*Macauley(x-(x2+xa/2))**3 - Ay/6*Macauley(x-(x2-xa/2))**3 + R3y/6*(Macauley(x-x3)**3) + C1*x + C2)
    return v

def w(x):
    w = -1/(E*Iyy*(R1x/6*(Macauley(x-x1)**3) + R2x/6*(Macauley(x-x2)**3) - Px/6*Macauley(x-(x2+xa/2))**3 - Ax/6*Macauley(x-(x2-xa/2))**3 + R3x/6*(Macauley(x-x3)**3) + C3*x + C4)
    return w
while x<= La:
    v.append(v(x))
    w.append(w(x))
    x+= dx

print(v)

