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
dx = 0.01
E = 1000
Ixx = 10
R1y = 2
R2y = 3
R3y = 4
Py = 5
Mx = 10
xa = 0.245
Pz = 5
R1z = 2
R2z = 3
R3z = 4
Az = 2
C1 = 2
Iyy = 10
x = 0
Izz = 10
Ay = 2
C2 = 1
C3 = 1
C4 = 1

v_tot = []
w_tot = []

q_int = 1

def Macauley(y):
    if y<0:
        y = 0
    return y

def v_prime(x):
    v_prim = -1/(E*Izz)*(q_int + R1y/2*(Macauley(x-x1)**2) + R2y/2*(Macauley(x-x2)**2) - Py/2*Macauley(x-(x2+xa/2))**2 - Ay/2*Macauley(x-(x2-xa/2))**2 + R3y/2*(Macauley(x-x3)**2)) + C1
    return v_prim

def w_prime(x):
    w_prim = -1/(E*Iyy)*(R1z/2*(Macauley(x-x1)**2) + R2z/2*(Macauley(x-x2)**2) - Pz/2*Macauley(x-(x2+xa/2))**2 - Az/2*Macauley(x-(x2-xa/2))**2 + R3z/2*(Macauley(x-x3)**2)) + C1
    return w_prim


def v(x):
    v = -1/(E*Izz)*(q_int + R1y/6*(Macauley(x-x1)**3) + R2y/6*(Macauley(x-x2)**3) - Py/6*Macauley(x-(x2+xa/2))**3 - Ay/6*Macauley(x-(x2-xa/2))**3 + R3y/6*(Macauley(x-x3)**3)) + C1*x + C2
    return v

def w(x):
    w = -1/(E*Iyy)*(R1z/6*(Macauley(x-x1)**3) + R2z/6*(Macauley(x-x2)**3) - Pz/6*Macauley(x-(x2+xa/2))**3 - Az/6*Macauley(x-(x2-xa/2))**3 + R3z/6*(Macauley(x-x3)**3)) + C3*x + C4
    return w

while x<= La:
    v_tot.append(v(x))
    w_tot.append(w(x))
    x+= dx

print(sum(v_tot))

## Verification



