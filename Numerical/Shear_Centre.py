# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:56:29 2020

@author: Marni
"""

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
A = pi*r**2/2
x_s = (ca-r)/s1


def qb1(s):
    qb1 = -t_skin/Izz*(r/(ca-r)*s**2/2)
    return qb1

def qb2(theta):
    qb2 = -t_skin/Izz*(r**2*cos(theta))
    return qb2

def qb3(s):
    qb3 = -t_skin/Izz*(r*s-s**2/2)
    return qb3

def qb4(s):
    return -t_skin/Izz*(-r*s+s**2/2)

s = 0

qb1_dist = []
qb2_dist = []
qb3_dist = []
qb4_dist = []

s_tot = []
while s<=s1:
    q = 0
    q += qb1(s)
    if s >= s1/6:
        q += s1/6*r/(ca-r)*B
    if s >= 2*s1/6:
        q += 2*s1/6*r/(ca-r)*B
    if s >= 3*s1/6:
        q += 3*s1/6*r/(ca-r)*B
    if s >= 4*s1/6:
        q += 4*s1/6*r/(ca-r)*B
    if s >= 5*s1/6:
        q += 5*s1/6*r/(ca-r)*B
        
    qb1_dist.append(q)
    s_tot.append(s)
    s+= ds

#plt.pyplot.plot(s_tot, qb1_dist)

s = 0

while s<=s2:
    qb3_dist.append(qb3(s)+qb1(s1))
    s+= ds

theta = pi/2

theta_s = []

while theta >= -pi/2:
    qb2_dist.append(qb2(theta)+qb1(s1))
    theta_s.append(theta)
    theta += -dtheta    

#plt.pyplot.plot(theta_s,qb2_dist)
s = 0

s_tot = []
while s<=s1:
    qb4_dist.append(qb4(s) + qb2(-pi/2) + qb3(s2))
    s_tot.append(s)
    s+=ds

plt.pyplot.plot(s_tot, qb4_dist)

#####
qs01 =  sum(qb2_dist)/len(qb2_dist)

qs02 = sum(qb1_dist)/len(qb2_dist)*s1 + sum(qb4_dist)/len(qb4_dist)*s1 +  (sum(qb4_dist)/len(qb4_dist) - qs01)*2*r

qs1 = []
qs4 = []
qs2 = []

n = 0
ksi = 0

for i in qb1_dist:
    n +=1
    qs1.append(i+qs02)
    ksi += ((ca-r)/s1 * x_s*n)*(i+qs02)/len(qb1_dist)
    
n = 0

for i in qb4_dist:
    n+=1
    qs1.append(i+qs02)
    ksi += ((ca-r)/s1 * x_s*n)*(i+qs02)/len(qb1_dist)
    
n = 0

for i in qb2_dist:
    n+=1
    qs2.append(i+qs01)
    ksi += -(cos(n*dtheta))*(i+qs01)/len(qb2_dist)
    
    

    

    
    
    
    
    
    
    
    