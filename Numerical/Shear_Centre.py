from math import *
import matplotlib as plt
import numpy as np

t_stiff = 0.0012
w_stiff = 0.017
h_stiff = 0.013
B = t_stiff*(w_stiff + h_stiff)
t_skin = 0.0011
t_spar = 0.0024 
Vy = 1
ha = 0.161
ca = 0.505
Izz = 4.582456727475651e-06
dz = 1e-5
r = ha/2
s1 = sqrt(r**2 + (ca-r)**2)
s2 = 2*r
ds = 1e-5
dtheta = 1e-5
ds_ = 1e-2
A = pi*r**2/2
x_s = (ca-r)/s1

A_circ = r**2*pi/2
A_triang = s1*r

def qb1(s):
    qb1 = -t_skin/Izz*(r/s1*s**2/2)
    return qb1

def qb2(theta):
    qb2 = -t_skin/Izz*(r**2*cos(theta)) + qb1(s1)
    return qb2

def qb3(s):
    qb3 = -t_skin/Izz*(r*s-s**2/2) + qb1(s1)
    return qb3

def qb4(s):
    return -t_skin/Izz*(-r*s+r/s1*s**2/2) + qb3(2*r) + qb2(-pi/2)

twist_circ = 1/A_circ*( (-1/Izz*r**3*(sin(pi/2) - sin(-pi/2) + qb1(s1)*pi*r)) - (-r*1/Izz*(2*r)**2/6 + r/s1*(2*r)**3/6 + qb1(s1)*2*r))
qs0_circ = twist_circ*(t_skin/(pi*r) + t_spar/(2*r))*2*A_circ   


twist_triang = 1/2/A_triang*((1/Izz*(r/s1*s1**3/6 + (-r*(2*r)**2/2 + r/s1*(2*r)**3/6) + (-r*s1**2/2 + r/s1*s1**3/6))) + qb1(s1)*2*r + qb2(-pi/2)*s1 + qb3(2*r)*s1)

qs0_triang = twist_triang*(t_skin/(2*s1) + t_spar/(2*r))*2*A_triang


def qs1(s):
    return qb1(s) + qs0_triang

def qs2(s):
    return qb2(s) + qs0_circ

def qs3(s):
    return qb3(s) + qs0_triang + qs0_circ

def qs4(s):
    return qb4(s) + qs0_triang


s = 0
n = 0
ksi = 0

while s<=s1:
    ksi += qs1(s)*ds*(ca-r-x_s*n*ds)*ds_
    n += 1
    s += ds

s = 0
n = 0

while s<=s1:
    ksi += qs4(s)*ds*x_s*n*ds*ds_
    n += 1
    s += ds

theta = pi/2
n = 0

while theta>=-pi/2:
    ksi += qs2(theta)*dtheta*ds_*cos(theta)
    n += 1
    theta += -dtheta
    
