#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:53:49 2020

@author: paula
"""
from numpy import cos,sin,pi,sqrt, arctan, arange

ha = 0.161
ca = 0.505
t_st = 0.0012
w_st = 0.017
h_st = 0.013
A_st = t_st * (w_st + h_st)
t_sk = 0.0011
t_sp = 0.0024
nstiff = 11

#Stiffener Spacing in the Triangle
r = ha/2
a = sqrt(r**2+(ca-r)**2)
P = 2*a + pi*r
space_st = P/nstiff

# Coordinates for Stiffeners
zlist = [r]
ylist = [0.0]

# Coordinates for Triangle Stiffeners (Positive)
alpha = arctan(r/(ca-r))
for i in range(4):
    ztriangle = - ca + r + space_st * (1/2 + i) * cos(alpha)
    ytriangle = space_st * (1/2 + i) * sin(alpha)
    zlist.append(ztriangle)
    ylist.append(ytriangle)

# Coordinates for Circle Stiffeners
i = -1
beta = space_st/r
print(beta)
while beta < pi/2:
    zcircle = r * cos(i*beta)
    ycircle = r * sin(i*beta)
    
    beta = beta + space_st/r
    
    zlist.append(zcircle)
    ylist.append(ycircle)
    zlist.append(zcircle)
    ylist.append(-ycircle)
    
    i = i + 1
    print(i)

# Coordinates for Triangle Stiffeners (Negative)
for i in range(4):
    ztriangle = - ca + r + space_st * cos(alpha) * (1/2 + i)
    ytriangle = - space_st * (1/2 + i) * sin(alpha)
    zlist.append(ztriangle)
    ylist.append(ytriangle)

# Centroid of the Cross-Section
A1 = t_sp * ha
A23 = t_sk * a
A4 = pi * r * t_sk
Atot = A1 + 2 * A23 + A4 + nstiff * A_st
zlistcentr = []
for i in zlist:
    zcoordst = A_st * (i)
    zlistcentr.append(zcoordst)
zcoord = (2 * A23 * (-a/2*cos(alpha)) + A4 * ((2*(r-t_sk))/pi)+ sum(zlistcentr))/Atot

# Moment of Inertia Contribution from Stiffeners
Izz_stlist = []
for i in ylist:
    Izz_st = A_st * i**2
    Izz_stlist.append(Izz_st)

Iyy_stlist = []
for j in zlist:
    Iyy_st = A_st * (zcoord - j)**2
    Iyy_stlist.append(Iyy_st)

Izz_stot = sum(Izz_stlist)
Iyy_stot = sum(Iyy_stlist)

# Moment of Inertia Contribution from Thin-Walled Cross-Section
Izz_1 = (1/12) * t_sp * ha**3
Iyy_1 = (1/12) * ha * t_sp**3 + A1 * (zcoord)**2
Izz_23 = (1/12) * t_sk * a**3 * (sin(alpha))**2 + A23 * (a/2*sin(alpha))**2
Iyy_23 = (1/12) * t_sk * a**3 * (cos(alpha))**2 + A23 * (-a/2*cos(alpha)-zcoord)**2

def fz(theta):
    fz = t_sk * r**3 * (cos(theta))**2
    return fz

def fy(theta):
    fy = t_sk * r**3 * (sin(theta))**2
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
Iyy_4 = sum(flist_y) - A4 * ((2*(r-t_sk))/pi)**2 + A4 * (zcoord - (2*(r-t_sk))/pi)**2

Izz_cs = Izz_1 + 2 * Izz_23 + Izz_4
Iyy_cs = Iyy_1 + 2 * Iyy_23 + Iyy_4

# Total Moment of Inertia
Izz = Izz_stot + Izz_cs # in m^4
Iyy = Iyy_stot + Iyy_cs # in m^4
print(zlist)
print(ylist)
print(zcoord)
print(Izz,Iyy)