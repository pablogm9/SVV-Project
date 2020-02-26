# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:11:48 2020

@author: Marni
"""
import numpy as np
from math import *

r = 0.0805
ca = 0.505
t_skin = 0.0011
t_spar = 0.0024
Vy = 1
Izz = 4.56e-5
s1 = sqrt(r**2 + (ca-r)**2)

def f_circ(s):
    r = 0.0805
    return t_skin*Vy/Izz*(r**2*sin(s))

def integral_circ(a, b):
    ds = (b - a)/(1e5+1)
    int_1 = []
    s  = a + ds
    while s < b:
        int_1.append(f_circ(s)*ds)
        s += ds
    return sum(int_1)

def second_intergral_circ(a,b):
    ds = (b - a)/(1e4+1)
    r = 0.0805
    int_1 = [0]
    int_2 = []
    s  = a + ds
    n = 0
    while s < b:
        int_1.append(f_circ(s)*ds*r)
        n+=1
        s += ds
    return sum(int_1)

def f_triang(s):
    s1 = sqrt(r**2 + (ca-r)**2)
    return t_skin*Vy/Izz*(r-(r/s1)*s)

def f_triang2(s):
    return t_skin*Vy/Izz*(r/s1)*s

def integral_triang(a, b):
    ds = (b - a)/(1e3+1)
    int_1 = []
    s  = a + ds
    while s < b:
        int_1.append(f_triang(s)*ds)
        s += ds
    return sum(int_1)

def second_integral_triang(a, b):
    ds = (b - a)/(1e3+1)
    int_1 = []
    s = a + ds
    int_2 = [f_triang(s)*ds*ds]
    s  = a + 2*ds
    while s < b:
        int_2.append(int_2[-1]+ (f_triang(s)*ds**2))
        s += ds
    return sum(int_2)

def integral_triang2(a, b):
    ds = (b - a)/(1e3+1)
    int_1 = []
    s  = a + ds
    while s < b:
        int_1.append(f_triang2(s)*ds)
        s += ds
    return sum(int_1)

def second_integral_triang2(a, b):
    ds = (b - a)/(1e3+1)
    s  = a + ds
    int_2 = [f_triang2(s)*ds**2]
    s  = a + 2*ds
    while s < b:
        int_2.append(int_2[-1]+ (f_triang2(s)*ds**2))
        s += ds
    return sum(int_2)

def f_spar(s):
    return t_spar*Vy/Izz*s

def f_spar2(s):
    return t_spar*Vy/Izz*(-r+s)

def integral_spar(a, b):
    ds = (b - a)/(1e3+1)
    int_1 = []
    s  = a + ds
    while s < b:
        int_1.append(f_spar(s)*ds)
        s += ds
    return sum(int_1)

def second_integral_spar(a, b):
    ds = (b - a)/(1e3+1)
    s  = a + ds
    int_2 = [f_spar(s)*ds**2]
    s  = a + 2*ds
    while s < b:
        int_2.append(int_2[-1]+ (f_spar(s)*ds**2))
        s += ds
    return sum(int_2)

def integral_spar2(a, b):
    ds = (b - a)/(1e3+1)
    int_1 = []
    s  = a + ds
    while s < b:
        int_1.append(f_spar2(s)*ds)
        s += ds
    return sum(int_1)

def second_integral_spar2(a, b):
    ds = (b - a)/(1e3+1)
    s  = a + ds
    int_2 = [f_spar2(s)*ds**2]
    s  = a + 2*ds
    while s < b:
        int_2.append(int_2[-1]+ (f_spar2(s)*ds**2))
        s += ds
    return sum(int_2)