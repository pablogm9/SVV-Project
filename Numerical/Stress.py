# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:01:46 2020

@author: Marni
"""
from math import *
import matplotlib as plt

ha = 0.161
ca = 0.505
z = 0
z_hat = 0.25
Mx = 1
My = 2
Mz = 3
x = 1

Ixx = 1
Iyy = 1
Izz = 1

def y(z):
    if z<ha/2:
        ys = ((ha/2)**2-(ha/2-z)**2)**0.5
    if z>=ha/2:
        ys = ha/2-(z-ha/2)*(ha/2)/(ca-ha/2)
    return ys

y_tot = []
z_tot = []
   
#while z<=ca:
#    y = profile_y(z)
#    z+=0.001
 #   y_tot.append(y)
 #   z_tot.append(z)
    
#plt.pyplot.plot(z_tot,y_tot)    
    
def stress_xx(z):
    s = My*(z-z_hat)/Iyy - Mz*y(z)/Izz
    return s

def stress_yy(z):
    s = -Mx*z/Ixx + Mz*x/Izz
    return s

def stress_zz(z):
    s = -My*(x)/Iyy + Mx*y(z)/Ixx
    return s

def stress_vm(x,y,z):
    s_vm = (((sigma_xx(z) - sigma_yy(z))**2 + (sigma_zz - sigma_xx)**2 + 6*(Txy**2 + Tyz**2 + Txz**2))/2)**0.5
    return s_vm

s_xx = []
s_yy = []
s_zz = []
zs = []
z = 0

while z<ca:
    s_xx.append(stress_xx(z))
    s_yy.append(stress_yy(z))
    s_zz.append(stress_zz(z))
    zs.append(z)
    z+=0.001
    
#plt.pyplot.plot(zs,s_xx)
#plt.pyplot.plot(zs,s_yy)
#plt.pyplot.plot(zs,s_zz)
    
   
