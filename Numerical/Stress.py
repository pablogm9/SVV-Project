# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:01:46 2020

@author: Marni
"""
from math import *
import matplotlib as plt

#Pseudo code to calculate shear stresses in the cross section
#arrays need to match in order to implement in the main

q_circle1,q_spar1,q_top,q_bottom,q_spar2,q_circle2 = section.get_shearflows()

normal_stresses = np.array([])
shear_stresses = np.array([])
von_Mises = np.array([])
for i in range(new_aerodata.shape[1]):
    
    normal_stresses_crosssection = np.array([])
    shear_stresses_crosssection = np.array([])
    von_Mises_crosssection = np.array([])
    
    V_y = V_y[i]
    W_z = W_z[i]
    T_x = T_x[i]
    #solving for shear flows due to torque
    #effect of stiffeners excluded in calculating qs_0 due to torque
    Y=[T_x[i],
       0]
    A=np.matrix([[2*A1, 2*A2],
                 [1/(2*A1*G)*((pi*ha/2 - 3*w_st)/t_sk + 3*w_st/(t_sk+t_st)) + ha/(2*A2*G*t_sp) + ha/(2*A1*G*t_sp), - 1/(2*A2*G)*((2*sqrt((ha/2)**2+(c-ha/2)**2) - 8*w_st)/t_sk + 8*w_st/(t_sk+t_st)) - ha/(2*A2*G*t_sp) - ha/(2*A1*G*t_sp)]])
    shears=np.linalg.solve(A,Y)
    
    q_circle0=shears[0]
    q_triangle0=shears[1]
    
    for j in range(len(crosssection_z)):
        
        z = crosssection_z[j]
        y = crosssection_y[j]
        
        sigma_xx = M_z[i]*y/Izz + M_y[i]*(z-z_centroid)/Iyy
        
        normal_stresses_crosssection = np.append(normal_stresses_crosssection,sigma_xx)
        #condition (j<5000 etc.) to make sure that calculations are made in the correct part of the aileron
        #shear flow array must match with a coordinate (yz plane) arrays -  z = crosssection_z[j] and y = crosssection_y[j]
        if j<5000:
            tau = (q_circle1 + q_circle0)/tsk)
        if j>5000 and j<10000:
            tau = (q_spar1 - q_circle0 + q_triangle0)/tsp)
#.......
            #......
        if j>100000:
            tau = (q_circle2 + q_circle0)/tsk
        #shear stress is due to both shear forces Vy and Vz, which are included in q_circle and q_spar etc. arrays
        #..and due to torque, which adds q_circle0 and/or q_triangle0 contributions
        
        shear_stresses_crosssection = np.append(shear_stresses_crosssection,tau)
        von_Mises = np.append(von_Mises,sqrt(sigma_xx**2 + 3*tau**2))








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
    
   
