#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 10:35:57 2020

@author: pablo
"""


from numpy import cos,sin,pi,sqrt, arctan, arange
from math import *
import matplotlib as plt




def get_MoI(ha,ca,t_st,w_st,h_st,A_st,t_sk,t_sp,nstiff)
    
    
 
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
    print(alpha)
    for i in range(4):
        ztriangle = - ca + r + space_st * (1/2 + i) * cos(alpha)
        ytriangle = space_st * (1/2 + i) * sin(alpha)
        zlist.append(ztriangle)
        ylist.append(ytriangle)
    
    # Coordinates for Circle Stiffeners
    beta = space_st/r
    while beta < pi/2:
        zcircle = r * cos(beta)
        ycircle = r * sin(beta)
        
        beta = beta + space_st/r
        
        zlist.append(zcircle)
        ylist.append(ycircle)
        zlist.append(zcircle)
        ylist.append(-ycircle)
    
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
    Atot = A1 + 2 * A23 + A4 + 11 * A_st
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
    
    
    return Izz,Iyy

def get_J(G,ha,t_sk,w_st,t_sp,t_st,Ca):
    #calculation of J by evaluating shear flows
    #assumption - stiffeners included in analysis when integrating length over thickness
    #stiffeners assumed as horizontal stiffeners, vertical part is disregarded
    A1=(ha/2)**2*pi/2
    A2=sqrt((ha/2)**2 + (ca-ha/2)**2)*ha/2
    randtorque=1
    Y=[randtorque,
    0]
    AA=np.matrix([[2*A1, 2*A2],
                 [1/(2*A1*G)*((pi*ha/2 - 3*w_st)/t_sk + 3*w_st/(t_sk+t_st)) + ha/(2*A2*G*t_sp) + ha/(2*A1*G*t_sp), - 1/(2*A2*G)*((2*sqrt((ha/2)**2+(c-ha/2)**2) - 8*w_st)/t_sk + 8*w_st/(t_sk+t_st)) - ha/(2*A2*G*t_sp) - ha/(2*A1*G*t_sp)]])
    shears=np.linalg.solve(AA,Y)
    q1=shears[0]
    q2=shears[1]
    dodz =  1/(2*A1)*(q1/G*((pi*ha/2 - 3*w_st)/t_sk + 3*w_st/(t_sk+t_st)) + (q1-q2)/G*(ha)/t_sp)
    dodz2=1/(2*A2)*(q2/G*((2*sqrt((ha/2)**2+(c-ha/2)**2) - 8*w_st)/t_sk + 8*w_st/(t_sk+t_st)) + (q2-q1)/G*(ha/t_sp))
    J=randtorque/(G*dodz)
    return J


def get_SC(t_st,w_st,h_st,t_sk,t_sp,V,ha,Ca,Izz,):
    
    r = ha/2
    dz = 1e-5
    Vy = 1
    
    s1 = sqrt(r**2 + (ca-r)**2)
    s2 = 2*r
    ds = 1e-5
    dtheta = 1e-5
    ds_ = 1e-2
    A = pi*r**2/2
    x_s = (ca-r)/s1
    
    A_circ = r**2*pi/2
    A_triang = s1*r
    B = t_st*(w_st + h_st)

    
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
      
    Shear_centre = -ksi + r     
    return Shear_centre
        
    
    
    
    
    
