# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:11:48 2020

@author: Marni
"""
r = 0.0805
ca = 0.505
t_skin = 0.0011
t_spar = 0.0024
Vy = 1
Vz = 1
Izz = 4.753851442684436e-06
Iyy = 4.5895702148629556e-05
s1 = sqrt(r**2 + (ca-r)**2)

#--------------------------------------------------------
# Integrator Circular Part

def f_circ(s):
    r = 0.0805
    return -t_skin*Vy/Izz*(r**2*sin(s))

def integral_circ(a, b):
    ds = (b - a)/(1e3+1)
    int_1 = []
    s  = a + ds
    while s < b:
        int_1.append(f_circ(s)*ds)
        s += ds
    return sum(int_1)

def second_integral_circ(a,b):
    ds = (b - a)/(1e3+1)
    r = 0.0805
    s  = a + ds 
    int_2 = [f_circ(s)*ds*ds*r]
    s  = a + ds
    while s < b:
        int_2.append(int_2[-1] + f_circ(s)*ds*ds*r)
        s += ds
    return sum(int_2)


#--------------------------------------------------------
# Integrator Triangular Part top
    
def f_triang(s):
    s1 = sqrt(r**2 + (ca-r)**2)
    return -t_skin*Vy/Izz*(r-(r/s1)*s)

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

#--------------------------------------------------------
# Integrator Triangular Part bottom

def f_triang2(s):
    s1 = sqrt(r**2 + (ca-r)**2)
    return -t_skin*Vy/Izz*-(r/s1)*s

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

#----------------------------------------------------------
# Spar Integrator top
    
def f_spar(s):
    return -t_spar*Vy/Izz*s

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

#----------------------------------------------------------
# Spar Integrator bottom
    
def f_spar2(s):
    return -t_spar*Vy/Izz*(-r+s)

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
    int_2 = [(f_spar2(s)*ds)*ds]
    s  = a + 2*ds
    while s < b:
        int_2.append(int_2[-1]+ (f_spar2(s)*ds**2))
        s += ds
    return sum(int_2)


#############################################################
#
# Integrals for z direction
#
############################################################


#--------------------------------------------------------
# Integrator Circular Part Z

def f_circ_z(s):
    r = 0.0805
    return -t_skin*Vz/Iyy*(r**2*cos(s))

def integral_circ_z(a, b):
    ds = (b - a)/(1e3+1)
    int_1 = []
    s  = a + ds
    while s < b:
        int_1.append(f_circ_z(s)*ds)
        s += ds
    return sum(int_1)

def second_integral_circ_z(a,b):
    ds = (b - a)/(1e3+1)
    r = 0.0805
    s  = a + ds 
    int_2 = [f_circ_z(s)*ds*ds*r]
    s  = a + ds
    while s < b:
        int_2.append(int_2[-1] + f_circ_z(s)*ds*ds*r)
        s += ds
    return sum(int_2)


#--------------------------------------------------------
# Integrator Triangular Part top z
    
def f_triang_z(s):
    return -t_skin*Vz/Iyy*-(ca-r)/s1*s

def integral_triang_z(a, b):
    ds = (b - a)/(1e3+1)
    int_1 = []
    s  = a + ds
    while s < b:
        int_1.append(f_triang_z(s)*ds)
        s += ds
    return sum(int_1)

def second_integral_triang_z(a, b):
    ds = (b - a)/(1e3+1)
    int_1 = []
    s = a + ds
    int_2 = [f_triang_z(s)*ds*ds]
    s  = a + 2*ds
    while s < b:
        int_2.append(int_2[-1]+ (f_triang_z(s)*ds**2))
        s += ds
    return sum(int_2)

#--------------------------------------------------------
# Integrator Triangular Part bottom z

def f_triang2_z(s):
    s1 = sqrt(r**2 + (ca-r)**2)
    return -t_skin*Vz/Iyy*(-(ca-r) + (ca-r)/s1*s)

def integral_triang2_z(a, b):
    ds = (b - a)/(1e3+1)
    int_1 = []
    s  = a + ds
    while s < b:
        int_1.append(f_triang2_z(s)*ds)
        s += ds
    return sum(int_1)

def second_integral_triang2_z(a, b):
    ds = (b - a)/(1e3+1)
    s  = a + ds
    int_2 = [f_triang2_z(s)*ds**2]
    s  = a + 2*ds
    while s < b:
        int_2.append(int_2[-1]+ (f_triang2_z(s)*ds**2))
        s += ds
    return sum(int_2)

#----------------------------------------------------------
# Spar Integrator top z 
    
def f_spar_z(s):
    return -t_spar*Vz/Iyy

def integral_spar_z(a, b):
    ds = (b - a)/(1e3+1)
    int_1 = []
    s  = a + ds
    while s < b:
        int_1.append(f_spar_z(s)*ds)
        s += ds
    return sum(int_1)

def second_integral_spar_z(a, b):
    ds = (b - a)/(1e3+1)
    s  = a + ds
    int_2 = [f_spar_z(s)*ds**2]
    s  = a + 2*ds
    while s < b:
        int_2.append(int_2[-1]+ (f_spar_z(s)*ds**2))
        s += ds
    return sum(int_2)

#----------------------------------------------------------
# Spar Integrator bottom z
    
def f_spar2_z(s):
    return -t_spar*Vz/Iyy

def integral_spar2_z(a, b):
    ds = (b - a)/(1e3+1)
    int_1 = []
    s  = a + ds
    while s < b:
        int_1.append(f_spar2_z(s)*ds)
        s += ds
    return sum(int_1)

def second_integral_spar2_z(a, b):
    ds = (b - a)/(1e3+1)
    s  = a + ds
    int_2 = [f_spar2_z(s)*ds**2]
    s  = a + 2*ds
    while s < b:
        int_2.append(int_2[-1]+ (f_spar2_z(s)*ds**2))
        s += ds
    return sum(int_2)
