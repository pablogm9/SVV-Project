#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 10:35:57 2020

@author: pablo
"""


from numpy import cos,sin,pi,sqrt, arctan, arange,matrix,pi,linspace,array,append
from numpy.linalg import solve
from math import *
import matplotlib as plt



def get_MoI(ha,ca,t_st,w_st,h_st,A_st,t_sk,t_sp,nstiff):
    
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
    
    
    return Izz,Iyy,zcoord

def get_J(G,ha,t_sk,w_st,t_sp,t_st,ca):
    #calculation of J by evaluating shear flows
    #assumption - stiffeners included in analysis when integrating length over thickness
    #stiffeners assumed as horizontal stiffeners, vertical part is disregarded
    A1=(ha/2)**2*pi/2
    A2=sqrt((ha/2)**2 + (ca-ha/2)**2)*ha/2
    randtorque=1
    Y=[randtorque,
    0]
    AA=matrix([[2*A1, 2*A2],
                 [1/(2*A1*G)*((pi*ha/2)/t_sk ) + ha/(2*A2*G*t_sp) + ha/(2*A1*G*t_sp), - 1/(2*A2*G)*((2*sqrt((ha/2)**2+(ca-ha/2)**2))/t_sk) - ha/(2*A2*G*t_sp) - ha/(2*A1*G*t_sp)]])
    shears=solve(AA,Y)
    q1=shears[0]
    q2=shears[1]
    dodz =  1/(2*A1)*(q1/G*((pi*ha/2 - 3*w_st)/t_sk + 3*w_st/(t_sk+t_st)) + (q1-q2)/G*(ha)/t_sp)
    dodz2=1/(2*A2)*(q2/G*((2*sqrt((ha/2)**2+(ca-ha/2)**2) - 8*w_st)/t_sk + 8*w_st/(t_sk+t_st)) + (q2-q1)/G*(ha/t_sp))
    J=randtorque/(G*dodz)
    return J


def get_SC(t_skin,t_spar,ha,ca,Izz):
    
    def qb1(theta):     # Shear Flow open section top circular
        return itg.integral_circ(0, theta)

    def qb2(s):     # Spar Top 
        return itg.integral_spar(0, s)

    def qb3(s):     # Triangular part top
        return  itg.integral_triang(0, s) + qb1(pi/2) + qb2(r)

    def qb4(s):     # Triangular part bottom
        return  itg.integral_triang2(0, s) + qb3(s1)

    def qb5(s):     # Spar bottom
        return  itg.integral_spar2(0, s) + qb2(r)

    def qb6(theta):     # Circular part bottom
        return  itg.integral_circ(-pi/2, theta) + qb1(pi/2)

    """

    shears = []
    s_tot = []
    s = 0
    while s<=s1:
        shears.append(qb4(s))
        s_tot.append(s)
        s += ds
    plt.pyplot.plot(s_tot, shears)
    """
    qs0_circ = -((itg.second_intergral_circ(0, pi/2)/t_skin - itg.second_integral_triang(0, r)/t_spar + itg.second_intergral_circ(-pi/2, 0)/t_skin) - itg.integral_spar2(0, r)/t_spar - qb2(r)*r/t_spar +qb1(pi/2)*pi*r/2/t_skin)*(pi*r/t_skin + 2*r/t_spar)**-1
    qs0_triang = -(itg.second_integral_triang(0, s1) + itg.second_integral_triang2(0, s1) + itg.second_integral_spar(0, r) + itg.second_integral_spar2(0, r) + (qb1(pi/2) + qb2(r))*s1/t_skin + qb3(s1)*s1/t_skin + qb2(r)*r/t_spar)*(2*s1/t_skin + 2*r/t_spar)**-1



    def qs1(s):
        return qb1(s) + qs0_circ

    def qs2(s):
        return qb2(s) - qs0_circ + qs0_triang

    def qs3(s):
        return qb3(s) + qs0_triang

    def qs4(s):
        return qb4(s) + qs0_triang

    def qs5(s):
        return qb5(s) - qs0_circ + qs0_triang

    def qs6(s):
        return qb6(s) + qs0_circ


    shears1 = []
    thetas1 = []
    theta = 0
    while theta<=(pi/2):
        shears1.append(qs1(theta))
        thetas1.append(theta)
        theta += dtheta

    shears6 = []
    thetas6 = []
    theta = -pi/2

    while theta<=0:
        shears6.append(qs6(theta))
        thetas6.append(theta + pi/2)
        theta += dtheta

    shears2 = []
    thetas2 = []
    s = 0
    while s<=r:
        shears2.append(qs2(s))
        thetas2.append(s)
        s += dr

    shears5 = []
    thetas5 = []
    s = 0
    while s<=r:
        shears5.append(qs5(s))
        thetas5.append(s)
        s += dr

    shears3 = []
    thetas3 = []
    s = 0
    while s<=s1:
        shears3.append(qs3(s))
        thetas3.append(s)
        s += ds

    shears4 = []
    thetas4 = []
    s = 0
    while s<=s1:
        shears4.append(qs4(s))
        thetas4.append(s)
        s += ds
        '''
    plt.pyplot.plot(thetas1, shears1)
    plt.pyplot.plot(thetas6, shears6)
    plt.pyplot.plot(thetas2, shears2)
    plt.pyplot.plot(thetas5, shears5)
    plt.pyplot.plot(thetas3, shears3)
    plt.pyplot.plot(thetas4, shears4)
    '''
    ksi = 0
    ksi += r*(ca-r)/s1*sum(shears3)/len(shears3)
    ksi += r*(ca-r)/s1*sum(shears4)/len(shears4)
    ksi += r*sum(shears1)/len(shears1)
    ksi += r*sum(shears6)/len(shears3)

    return -ksi
    
    
def discretize_crosssection(n):
    
    # Input n is the number of points PER SECTION of the cross section
    # E.g. if n=1000, the semicircular part of the crossection will have 1000 points
    
    ha = 0.161 #[m]
    Ca = 0.505  #[m]
    
    
    
    circle_1_y = array([])
    circle_1_z = array([])
    
    circle_2_y = array([])
    circle_2_z = array([])
    
    
    spar_1_y = array([])
    spar_1_z = array([])
    
    spar_2_y = array([])
    spar_2_z = array([])
    
    topline_y = array([])
    topline_z = array([])
    
    bottomline_y = array([])
    bottomline_z = array([])

    # Generate circle coords
    r = 0.5*ha
    theta_range_1 = linspace(pi,pi/2,n)
    theta_range_2 = linspace(3*pi/2,pi,n)
    
    for i in range(len(theta_range_1)):
        theta = theta_range_1[i]
        z = - r*cos(theta)
        y = r*sin(theta)
        
        circle_1_y = append(circle_1_y,y)
        circle_1_z = append(circle_1_z,z)
        
    for i in range(len(theta_range_2)):
        theta = theta_range_2[i]
        z = - r*cos(theta)
        y = r*sin(theta)
        
        circle_2_y = append(circle_2_y,y)
        circle_2_z = append(circle_2_z,z)

    
    # Generate spar coords  
    spar_range_1 = linspace(0,r,n)
    spar_range_2 = linspace(-r,0,n)
    
    for i in range(len(spar_range_1)):
        y = spar_range_1[i]
        z = 0
        
        spar_1_y = append(spar_1_y,y)
        spar_1_z = append(spar_1_z,z)  
 

    for i in range(len(spar_range_2)):
        y = spar_range_2[i]
        z = 0
        
        spar_2_y = append(spar_2_y,y)
        spar_2_z = append(spar_2_z,z)  
    
    
    # Generate coords for both lines
    
    run = -(Ca-r)
    rise_top = -r
    rise_bottom = r

    m_top = rise_top/run 
    m_bottom = rise_bottom/run
    
    b_top = r 
    b_bottom = -r
    
    line_range = linspace(0,-(Ca-r),n)
    
    for i in range(len(line_range)):
        
        z = line_range[i]
        
        top_y_coord = m_top*z + b_top
        bottom_y_coord = m_bottom*z + b_bottom
        
        topline_y = append(topline_y,top_y_coord)
        topline_z = append(topline_z,z)
        
        bottomline_y = append(bottomline_y,bottom_y_coord)
        bottomline_z = append(bottomline_z,z)
    

    circle_1_coords = array([circle_1_z,circle_1_y])
    circle_2_coords = array([circle_2_z,circle_2_y])
    spar_1_coords = array([spar_1_z,spar_1_y])
    spar_2_coords = array([spar_2_z,spar_2_y])
    topline_coords = array([topline_z,topline_y])
    bottomline_coords = array([bottomline_z,bottomline_y])

    
    return circle_1_coords,circle_2_coords,spar_1_coords,spar_2_coords,topline_coords,bottomline_coords
    
    

def get_shearflows():
    
    Izz = 4.753851442684436e-06  
    t_skin =  0.0011
    t_spar = 0.0024
    ca = 0.505
    ha = 0.161
    r = ha/2
    dz = 1e-5
    Vy = 1
    dr = r/1e4
    s1 = sqrt(r**2 + (ca-r)**2)
    s2 = 2*r
    ds = s1/1e4
    dtheta = pi/2/1e4
    A = pi*r**2/2
    x_s = (ca-r)/s1
    angle = (0.505-r)/s1
    Vz = 1
    Iyy = 4.5895702148629556e-05
    
    A_circ = r**2*pi/2
    A_triang = (ca-r)*r
    
    
    
    def qb1(theta):
        return Vy*t_skin/Izz*(r**2*-cos(theta))
    
    def qb2(s):
        return Vy*t_spar/Izz*s**2/2
    
    def qb3(s):
        return Vy*t_skin/Izz*(r*s - s**2/2*r/s1) - qb1(0) + qb2(r)
    
    def qb4(s):
        return Vy*t_skin/Izz*(-s**2/2*r/s1) + qb3(s1)
    
    def qb5(s):
        return Vy*t_spar/Izz*(-r*s + s**2/2) + qb2(r)
    
    def qb6(theta):
        return Vy*t_skin/Izz*(r**2*-cos(theta)) - qb1(0)
    
    
    
    def qb1_z(theta):
        return Vz*t_skin/Izz*(r**2*sin(theta))
    
    def qb2_z(s):
        return Vz*t_spar/Izz*s
    
    def qb3_z(s):
        return Vz*t_skin/Izz*angle*s**2/2 + qb1(pi/2) + qb2(r)
    
    def qb4_z(s):
        return Vz*t_skin/Izz*(-(ca-r)*s + angle*s**2/2) + qb3(s1)
    
    def qb5_z(s):
        return Vz*t_spar/Izz*s + qb2(r)
    
    def qb6_z(theta):
        return Vz*t_skin/Izz*(r**2*sin(theta)) + qb1(pi/2)
    
    """
    shears = []
    s_tot = []
    s = 0
    while s<=s1:
        shears.append(qb4(s))
        s_tot.append(s)
        s += ds
    plt.pyplot.plot(s_tot, shears)
    """
    
    qs0_circ = -(Vy/Izz*(-r**3*sin(pi/2) - r**3/6 - (-r**3/2 + r**3/6) + (-r**3*sin(-pi/2))) - qb2(r)*r/t_spar - qb1(0)*pi*r/2/t_skin)*(pi*r/t_skin + 2*r/t_spar)**-1
    qs0_triang = -(Vy/Izz*((r**3/2 - r**4/6/s1) + -r**4/6/s1 + r**3/6 + (-r**3/2 + r**3/6)) - qb1(0)*s1/t_skin + qb2(r)*s1/t_spar + qb3(s1)*s1/t_skin + qb2(r)*r/t_spar)*(2*s1/t_skin + 2*r/t_spar)**-1
    
    qs0_circ_z = -(Vz/Izz*(2*r**3*cos(0) - r**2/2 - r**2/2) - qb2(r)*r/t_spar + qb1(pi/2)*pi*r/2/t_skin)*(pi*r/t_skin + 2*r/t_spar)**-1
    qs0_triang_z = -(Vz/Izz*(s1**3/6 + -(ca-r)*s1**2/2 + angle*s1**3/6 + r**2/2) + qb1(pi/2)*s1/t_skin + qb2(r)*s1/t_spar + qb3(s1)*s1/t_skin + qb2(r)*r/t_spar)*(2*s1/t_skin + 2*r/t_spar)**-1
    
    def qs1(s):
        return qb1(s) + qs0_circ
    
    def qs2(s):
        return qb2(s) - qs0_circ + qs0_triang
    
    def qs3(s):
        return qb3(s) + qs0_triang
    
    def qs4(s):
        return qb4(s) + qs0_triang
    
    def qs5(s):
        return qb5(s) - qs0_circ + qs0_triang
    
    def qs6(s):
        return qb6(s) + qs0_circ
    
    
    
    def qs1_z(s):
        return qb1_z(s) + qs0_circ_z
    
    def qs2_z(s):
        return qb2_z(s) - qs0_circ_z + qs0_triang_z
    
    def qs3_z(s):
        return qb3_z(s) + qs0_triang_z
    
    def qs4_z(s):
        return qb4_z(s) + qs0_triang_z
    
    def qs5_z(s):
        return qb5_z(s) - qs0_circ_z + qs0_triang_z
    
    def qs6_z(s):
        return qb6_z(s) + qs0_circ_z
    
    '''
    shears = []
    s_tot = []
    s = -pi/2
    while s<=0:
        shears.append(qs6_z(s))
        s_tot.append(s)
        s += ds
    plt.pyplot.plot(s_tot, shears)
    '''
    
        
    qsc1 = []
    qsc2 = []
    qsp1 = []
    qsp2 = []
    qtop = []
    qbot = []
     
    s = 0
    while s<=pi/2:
        qsc1.append(qs1_z(s) + qs1(s))
        s += dtheta
    
    s = 0
    
    while s<=r-dr:
        qsp1.append(qs2_z(s) + qs2(s))
        s += dr
    
    s = 0
    
    while s<=s1-ds:
        qtop.append(qs3_z(s) + qs3(s))
        s += ds
     
    s = 0
    
    while s<=s1-ds:
        qbot.append(qs4_z(s) + qs4(s))
        s += ds
    
    s = 0
    
    while s<=r-dr:
        qsp2.append(qs5_z(s) + qs5(s))
        s += dr
     
    s = -pi/2
    
    while s<=0:
        qsc2.append(qs6_z(s) + qs6(s))
        s += dtheta
    

    return qsc1,qsp1,qtop,qbot,qsp2,qsc2






   
    
    
