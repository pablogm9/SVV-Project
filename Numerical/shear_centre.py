# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:48:23 2020

@author: Gytha
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
    ztriangle = - spacingtriangle * i * cos(alpha)
    ytriangle = spacingtriangle * i * sin(alpha)
    zlist.append(ztriangle)
    ylist.append(ytriangle)

# Coordinates for Circle Stiffeners
for i in range(1,4):
    zcircle = r * sin(i*pi/4)
    ycircle = r * cos(i*pi/4)
    zlist.append(zcircle)
    ylist.append(ycircle)
    
# Coordinates for Triangle Stiffeners (Negative)
for i in range(1,5):
    ztriangle = - spacingtriangle * cos(alpha) * (5 - i)
    ytriangle = - spacingtriangle * (5 - i) * sin(alpha)
    zlist.append(ztriangle)
    ylist.append(ytriangle)

# Moment of Inertia Contribution from Stiffeners
Izz_stlist = []
for i in ylist:
    Izz_st = A_st * i**2
    Izz_stlist.append(Izz_st)

Iyy_stlist = []
for j in zlist:
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

#print(Izz, Iyy)
print(ylist)



""" SHEAR CENTRE CALCULATION """

import numpy as np
from numpy import math as m

"input (geometry in m)"
h = 0.161
c = 0.505
t_st = 0.0012
w_st = 0.017
h_st = 0.013
A_st = t_st*(w_st + h_st)
t = 0.0011
t_sp = 0.0024
y_st = ylist
A_en_c = np.pi*(h/2)**2/2 #enclosed area circle
A_en_t = h*(c-h/2)/2 #enclosed area triangle
G = 28 #GPa
r = h/2
l = np.sqrt(r**2+(c-h/2)**2) #length of skin
M = -1/Izz #MoI contribution

"q_b_i"
def q_b1(s):
    return M*(t*r**2*(-cos(s)+1))+M*A_st*sum(y_st[4:6])

def q_b2(s):
    return M*(t_sp*s**2/2)

def q_b3(s): 
    return M*(t*(r*s-r/l*s**2/2))+q_b1(np.pi/2)+q_b2(r)+M*A_st*sum(y_st[0:4])

def q_b4(s):
    return M*(t*(-r/l*s**2/2))+q_b3(l)+M*A_st*sum(y_st[7:11])

def q_b5(s):
    return M*(t_sp*s**2/2)+q_b4(l)

def q_b6(s):
    return M*(t*r**2*(-cos(s)+1))+q_b4(l)-q_b5(-r)+M*A_st*y_st[6]

"q_s0"
#integrated over s
q_1 = M*(t*r**2*(-sin(np.pi/2)+np.pi/2+sin(0)))+M*A_st*sum(y_st[4:6])*np.pi/2
q_2 = M*(t_sp*r**3/6)
q_3 = M*(t*(r*l**2/2-r/l*l**3/6))+(q_b1(np.pi/2)+q_b2(r))*l+M*A_st*sum(y_st[0:4])*l
q_4 = M*(t*(r/l*l**3/6))+q_b3(l)*l+M*A_st*sum(y_st[7:11])*l
q_5 = M*(t_sp*r**3/6)+q_b4(l)*r
q_6 = M*(t*r**2*(-sin(0)+sin(-np.pi/2)+np.pi/2))+(q_b4(l)-q_b5(-r))*np.pi/2+M*A_st*y_st[6]*(np.pi/2)

q_s0_c = -(1/t*(q_1+q_6)+1/t_sp*(q_2+q_5))/(1/t*(np.pi*r)+1/t_sp*2*r)
q_s0_t = -(1/t*(q_3+q_4)+1/t_sp*(q_2+q_5))/(1/t*l*2+1/t_sp*2*r)

"moment equilibrium equation"
#integrated over s including distance
qlist = []
q_1s = M*(t*r**2*(-np.pi/2*sin(np.pi/2)-cos(np.pi/2)+(np.pi/2)**2/2+cos(0)))+M*A_st*sum(y_st[4:6])*(np.pi/2)**2/2
qlist.append(q_1s)
q_2s = M*(t_sp*r**4/24)
qlist.append(q_2s)
q_3s = M*(t*(r*l**3/6-r/l*l**4/24))+(q_b1(np.pi/2)+q_b2(r))*l**2/2+M*A_st*sum(y_st[0:4])*l**2/2
qlist.append(q_3s)
q_4s = M*(t*(r/l*l**4/24))+q_b3(l)*l**2/2+M*A_st*sum(y_st[7:11])*l**2/2
qlist.append(q_4s)
q_5s = M*(t_sp*r**4/24)+q_b4(l)*r**2/2
qlist.append(q_5s)
q_6s = M*(t*r**2*(-cos(0)-np.pi/2*sin(-np.pi/2)+cos(np.pi/2)+(np.pi/2)**2/2))+(q_b4(l)-q_b5(-r))*(np.pi/2)**2/2+M*A_st*y_st[6]*(np.pi/2)**2/2
qlist.append(q_6s)

ksi = sum(qlist)+2*A_en_c*q_s0_c+2*A_en_t*q_s0_t
print(ksi)

