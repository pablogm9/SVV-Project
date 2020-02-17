# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:56:32 2020

@author: edgar
"""

from math import *
import numpy as np

"INPUTS"
La = 1.611 #m
x1 = 0.125
x2 = 0.498
x3 = 1.494 #m
xa = 0.245    #m
theta=30   #degrees
Py = 49200*sin(theta)
Pz = 49200*cos(theta)   #N
zshear = 1   #z coordinate of shear centre
ha = 0.161 #m
G = 1
J = 1

"UNKNOWNS - to be put in a matrix equation"
R1y=1
R2y=1
R3y=1
R1z=1
R2z=1
R3z=1
Ay=1
Az=1
C5=1


"Tools to implement Macauley functions"
def Macaulay(y):
    if y<0:
        mac=0
    else:
        mac=y
    return mac
"Moment around z axis as a function of x"

D=1   
Td=1
"Td = integral from 0 to chord | q(x,z)(z- zshear) dz"
"D=  double integral of ||q dx dx"
S=1
"S=integral from 0 to x |q(x)dx"
DTd=1
"DTd = double integral 0tox ||Td dxdx"

def moments(x):
    Mz = -D + R1y*Macaulay(x-x1) + R2y*Macaulay(x-x2) - Py*Macaulay(x-(x2+xa/2)) - Ay*Macaulay(x-(x2-xa/2)) + R3y*Macaulay(x-x3)
    My = -R1z*Macaulay(x-x1) + Pz*Macaulay(x-(x2+xa/2)) - R2z*Macaulay(x-x2) - R3z*(x-x3)+ Az*Macaulay(x-(x2-xa/2))
    T = -Td - (0-zshear)*R1y*Macaulay(x-x1)**0 - (0-zshear)*R2y*Macaulay(x-x2)**0 - (0-zshear)*R3y*Macaulay(x-x3)**0 + (0-zshear)*Py*Macaulay(x-(x2+xa/2))**0 - (ha/2)*Pz*Macaulay(x-(x2+xa/2))**0 - (ha/2)*Az*Macaulay(x-(x2-xa/2))**0 + (0-zshear)*Ay*Macaulay(x-(x2-xa/2))**0
    Sy = -S + R1y*Macaulay(x-x1)**0 + R2y*Macaulay(x-x2)**0 - Py*Macaulay(x-(x2+xa/2))**0 - Ay*Macaulay(x-(x2-xa/2))**0 + R3y*Macaulay(x-x3)**0
    Sz = -R1z*Macaulay(x-x1)**0 - R2z*Macaulay(x-x2)**0 + Pz*Macaulay(x-(x2+xa/2))**0 + Az*Macaulay(x-(x2-xa/2))**0 - R3z*Macaulay(x-x3)**0
    twist = 1/(G*J)*(-DTd - (0-zshear)*R1y*Macaulay(x-x1)**1 - (0-zshear)*R2y*Macaulay(x-x2)**1 - (0-zshear)*R3y*Macaulay(x-x3)**1 + (0-zshear)*Py*Macaulay(x-(x2+xa/2))**1 - (ha/2)*Pz*Macaulay(x-(x2+xa/2))**1 - (ha/2)*Az*Macaulay(x-(x2-xa/2))**1 + (0-zshear)*Ay*Macaulay(x-(x2-xa/2))**1 ) + C5
    return(Mz,My,T)
    
    
