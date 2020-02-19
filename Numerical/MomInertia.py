#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:32:31 2020

@author: paula
"""
from numpy import cos,sin,pi,sqrt, arctan, arange

ha = 16.1
ca = 50.5
t_st = 0.12
w_st = 1.7
h_st = 1.3
A_st = t_st * (w_st + h_st)
t_sk = 0.11
t_sp = 0.24

#Stiffener Spacing in the Triangle
r = ha/2
a = sqrt(r**2+(ca-r)**2)
spacingtriangle = a/5

# Coordinates for Stiffeners
zlist = []
ylist = []

# Coordinates for Triangle Stiffeners (Positive)
alpha = arctan(r/(ca-r))

for i in range(1,5):
    ztriangle = spacingtriangle * i * cos(alpha)
    ytriangle = spacingtriangle * i * sin(alpha)
    zlist.append(ztriangle)
    ylist.append(ytriangle)

# Coordinates for Circle Stiffeners
for i in range(1,4):
    zcircle = spacingtriangle * 5 * cos(alpha) + r * sin(i*pi/4)
    ycircle = r * cos(i*pi/4)
    zlist.append(zcircle)
    ylist.append(ycircle)
    
# Coordinates for Triangle Stiffeners (Negative)
for i in range(1,5):
    ztriangle = spacingtriangle * cos(alpha) * (5 - i)
    ytriangle = - spacingtriangle * (5 - i) * sin(alpha)
    zlist.append(ztriangle)
    ylist.append(ytriangle)

# Moment of Inertia Contribution from Stiffeners
Izz_stlist = []
for i in zlist:
    Izz_st = A_st * i**2
    Izz_stlist.append(Izz_st)

Iyy_stlist = []
for j in ylist:
    Iyy_st = A_st * j**2
    Iyy_stlist.append(Iyy_st)

Izz_stot = sum(Izz_stlist)
Iyy_stot = sum(Iyy_stlist)

# Moment of Inertia Contribution from Thin-Walled Cross-Section
Izz_1 = (1/12) * t_sp * ha**3
Iyy_1 = (1/12) * ha * t_sp**3
J_1 = (1/12) * t_sp * ha * (t_sp**2 + ha**2)
Izz_23 = (1/12) * t_sk * a**3 * (sin(alpha))**2
Iyy_23 = (1/12) * t_sk * a**3 * (cos(alpha))**2

def fz(theta):
    fz = t_sk * r * (r - r * cos(theta))**2
    return fz

def fy(theta):
    fy = t_sk * r * (r * sin(theta))**2
    return fy

flist_z = []
flist_y = []
h = pi/180
angles = arange(0,pi+h,h)
for i in range(len(angles)-1):
    area_z = h * fz((angles[i] + angles[i+1])/2)
    flist_z.append(area_z)
    area_y = h * fy((angles[i] + angles[i+1])/2)
    flist_y.append(area_y)

Izz_4 = sum(flist_z)
Iyy_4 = sum(flist_y)

Izz_cs = Izz_1 + 2 * Izz_23 + Izz_4
Iyy_cs = Iyy_1 + 2 * Iyy_23 + Iyy_4

# Total Moment of Inertia
Izz = (Izz_stot + Izz_cs)*10**(-8) # in m^4
Iyy = (Iyy_stot + Iyy_cs)*10**(-8) # in m^4