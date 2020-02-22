# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:56:32 2020

@author: edgar
"""

from math import *
import numpy as np
from matplotlib import pyplot as plt
def matrix_solver(c,L,x1,x2,x3,xa,theta,t_st,t_sk,t_sp,w_st,ha,Py,Pz,d1,d2,d3,G,E,Izz,Iyy,zshear, S,D, Td, DTd1,DTd2,DTd3, FI1, FI2, FI3):
    
    
    A1=pi/2*ha**2      #m^2
    A2=(c-ha/2)*ha/2      #m^2

    
    #calculation of J by evaluating shear flows
    #assumption - stiffeners included in analysis when integrating length over thickness
    #stiffeners assumed as horizontal stiffeners, vertical part is disregarded
    randtorque=1
    Y=[randtorque,
       0]
    A=np.matrix([[2*A1, 2*A2],
                 [1/(2*A1*G)*((pi*ha/2 - 3*w_st)/t_sk + 3*w_st/(t_sk+t_st)) + ha/(2*A2*G*t_sp) + ha/(2*A1*G*t_sp), - 1/(2*A2*G)*((2*sqrt((ha/2)**2+(c-ha/2)**2) - 8*w_st)/t_sk + 8*w_st/(t_sk+t_st)) - ha/(2*A2*G*t_sp) - ha/(2*A1*G*t_sp)]])
    shears=np.linalg.solve(A,Y)
    q1=shears[0]
    q2=shears[1]
    dodz =  1/(2*A1)*(q1/G*((pi*ha/2 - 3*w_st)/t_sk + 3*w_st/(t_sk+t_st)) + (q1-q2)/G*(ha)/t_sp)
    dodz2=1/(2*A2)*(q2/G*((2*sqrt((ha/2)**2+(c-ha/2)**2) - 8*w_st)/t_sk + 8*w_st/(t_sk+t_st)) + (q2-q1)/G*(ha/t_sp))
    J=randtorque/(G*dodz)   #EXCLUDE stiffeners (?)
    
 
    "UNKNOWNS - to be put in a matrix equation, this is the unknown vector - a solution (variable 'RES'). It is a vector of size 13x1"
    #R1y=1
    #R2y=1
    #R3y=1
    #R1z=1
    #R2z=1
    #R3z=1
    #Ay=1
    #Az=1
    #C1=1
    #C2=1
    #C3=1
    #C4=1
    #C5=1
    
    
    "Tools to implement Macauley functions"
    def Macaulay(y):
        if y<0:
            mac=0
        else:
            mac=y
        return mac
    def Macauley(y):
        if y<0:
            mac=0
        else:
            mac=y
        return mac
    
    #NECESSARY INPUTS TO SOLVE THE SYSTEM:
        #ZSHEAR, Izz, Iyy, J
        #AND THE FOLLOWING 5 INTEGRALS
    #D=100  
    #Td=100
    "Td = integral from 0 to chord | q(x,z)(z- zshear) dz. Done for each section by calculating the ersultant torque from aerodynamic load (acts through centroid - causes moment aroiund shear centre"
    "D=  double integral of ||q dx dx - performed by double integrator"
    #S=100
    "S=integral from 0 to x |q(x)dx. Integration to get force in Y direction"
    #DTd=100
    "DTd = double integral 0tox ||Td dxdx, which is just torque at each cross-section times distance x" 
    #q_int=100; #quadruple integral over dx
    #FI=100;  #done by quadruple integrator
    #FI=fourth integral
    #Q=100  #not needed probably
    #Q is total aerodynamic loading sumsum(q(i,j))
    

    
    #AA=np.matrix([[2,2],
    #              [1,2]])
   # U=[8,7]
    #X=np.linalg.solve(AA,U)    
    #print(X)
    "unknown vector X=[R1y,R2y,R3y,R1z,R2z,R3z,Ay,Az,C1,C2,C3,C4,C5]"
    '''RHS vector Y with BCs'''
    Y=[D+Py*Macaulay(L-(x2+xa/2)),  #Mz(L)=0
       0,                             #My(L)=0
       Td-(0-zshear)*Py*Macaulay(L-(x2+xa/2))**0+(ha/2)*Pz*Macaulay(L-(x2+xa/2))**0,    #T(L)=0
       S+Py*Macaulay(L-(x2+xa/2))**0,      #Sy(L)=0
       -Pz*Macaulay(L-(x2+xa/2))**0,       #Sz(L)=0
       d1*cos(theta) +FI1/(Izz*E) - Py/(E*Izz)*1/6*Macauley(x1-(x2+xa/2))**3+(zshear-0)*DTd1/(G*J) - (zshear-0)*(0-zshear)/(G*J)*Py/6*Macauley(x1-(x2+xa/2))**3 + (zshear-0)*ha/2*1/(G*J)*Pz/6*Macauley(x1-(x2+xa/2))**3,      #v(x1) + twist(x1)*(zshear-0) = d1*cos(theta)
       FI2/(Izz*E) - Py/(E*Izz)*1/6*Macauley(x2-(x2+xa/2))**3+(zshear-0)*DTd2/(G*J) - (zshear-0)*(0-zshear)/(G*J)*Py/6*Macauley(x2-(x2+xa/2))**3 + (zshear-0)*ha/2*1/(G*J)*Pz/6*Macauley(x2-(x2+xa/2))**3,                     #v(x2) + twist(x2)*(zshear-0) = 0
       d3*cos(theta) +FI3/(Izz*E) - Py/(E*Izz)*1/6*Macauley(x3-(x2+xa/2))**3+(zshear-0)*DTd3/(G*J) - (zshear-0)*(0-zshear)/(G*J)*Py/6*Macauley(x3-(x2+xa/2))**3 + (zshear-0)*ha/2*1/(G*J)*Pz/6*Macauley(x3-(x2+xa/2))**3,      #v(x3) + twist(x3)*(zshear-0) = d3*cos(theta)
       -d1*sin(theta)+Pz/(E*Iyy)*1/6*Macauley(x1-(x2+xa/2))**3+ DTd1/(G*J)*(ha/2) - (0-zshear)/(G*J)*(ha/2)*Py*Macauley(x1-(x2+xa/2))**1 + ha/(2*G*J)*ha/2*Pz*Macauley(x1-(x2+xa/2))**1,            #w(x1) + twist(x1)*(ha/2) = -d1*sin(theta)
       Pz/(E*Iyy)*1/6*Macauley(x2-(x2+xa/2))**3+ DTd2/(G*J)*(ha/2) - (0-zshear)/(G*J)*(ha/2)*Py*Macauley(x2-(x2+xa/2))**1 + ha/(2*G*J)*ha/2*Pz*Macauley(x2-(x2+xa/2))**1,                           #w(x2) + twist(x2)*(ha/2) = 0
       -d3*sin(theta)+Pz/(E*Iyy)*1/6*Macauley(x3-(x2+xa/2))**3+ DTd3/(G*J)*(ha/2) - (0-zshear)/(G*J)*(ha/2)*Py*Macauley(x3-(x2+xa/2))**1 + ha/(2*G*J)*ha/2*Pz*Macauley(x3-(x2+xa/2))**1,            #w(x3) + twist(x3)*(ha/2) = -d3*sin(theta)
       Pz/(E*Iyy)*1/6*Macauley(x2-xa/2-(x2+xa/2))**3 ,       #w(x2-xa/2)=0 at ACTUATOR 1
       0]          #Ay - Aztan(30) = 0
    '''
    This Y vector includes all boundary conditions - RHS
    The equation becomes XX*RES = Y
    '''
    
    
    XX=np.matrix([[Macaulay(L-x1),Macaulay(L-x2),Macaulay(L-x3),0,0,0,-Macaulay(L-(x2-xa/2)),0,0,0,0,0,0],
                  [0,0,0,-Macaulay(L-x1),-Macaulay(L-x2),-Macaulay(L-x3),0,Macaulay(L-(x2-xa/2)),0,0,0,0,0],
                  [- (0-zshear)*Macaulay(L-x1)**0, -(0-zshear)*Macaulay(L-x2)**0, -(0-zshear)*Macaulay(L-x3)**0,0,0,0,(0-zshear)*Macaulay(L-(x2-xa/2))**0, -(ha/2)*Macaulay(L-(x2-xa/2))**0,0,0,0,0,0],
                  [Macaulay(L-x1)**0,Macaulay(L-x2)**0,Macaulay(L-x3)**0,0,0,0,-Macaulay(L-(x2-xa/2))**0,0,0,0,0,0,0],
                  [0,0,0,-Macaulay(L-x1)**0, -Macaulay(L-x2)**0, -Macaulay(L-x3)**0,0,Macaulay(L-(x2-xa/2))**0, 0,0,0,0,0],
                  [-1/(E*Izz)*1/6*(Macauley(x1-x1)**3)-(zshear-0)*1/(G*J)*(0-zshear)*Macaulay(x1-x1)**1, -1/(E*Izz)*1/6*(Macauley(x1-x2)**3)- (zshear-0)*1/(G*J)*(0-zshear)*Macaulay(x1-x2)**1, 0, 0,0,0, 0,0,x1,1,0,0,1],
                  [-1/(E*Izz)*1/6*(Macauley(x2-x1)**3)-(zshear-0)*1/(G*J)*(0-zshear)*Macaulay(x2-x1)**1, -1/(E*Izz)*1/6*(Macauley(x2-x2)**3)- (zshear-0)*1/(G*J)*(0-zshear)*Macaulay(x2-x2)**1, 0, 0,0,0, 1/(E*Izz)*1/6*Macauley(x2-(x2-xa/2))**3 +(zshear-0)/(G*J)*(0-zshear)*Macaulay(x2-(x2-xa/2))**1,-(zshear-0)/(G*J)*(ha/2)*Macaulay(x2-(x2-xa/2))**1,x2,1,0,0,1],
                  [-1/(E*Izz)*1/6*(Macauley(x3-x1)**3)-(zshear-0)*1/(G*J)*(0-zshear)*Macaulay(x3-x1)**1, -1/(E*Izz)*1/6*(Macauley(x3-x2)**3)- (zshear-0)*1/(G*J)*(0-zshear)*Macaulay(x3-x2)**1, -1/(E*Izz)*1/6*(Macauley(x3-x1)**3)-(zshear-0)*1/(G*J)*(0-zshear)*Macaulay(x3-x1)**1, 0,0,0, 1/(E*Izz)*1/6*Macauley(x3-(x2-xa/2))**3 +(zshear-0)/(G*J)*(0-zshear)*Macaulay(x3-(x2-xa/2))**1,-(zshear-0)/(G*J)*(ha/2)*Macaulay(x3-(x2-xa/2))**1,x3,1,0,0,1],
                  [-ha/(2*G*J)*(0-zshear)*Macaulay(x1-x1)**1, -ha/(2*G*J)*(0-zshear)*Macaulay(x1-x2)**1, -ha/(2*G*J)*(0-zshear)*Macaulay(x1-x3)**1,  -1/(E*Iyy)*1/6*(Macauley(x1-x1)**3), -1/(E*Iyy)*1/6*(Macauley(x1-x2)**3), -1/(E*Iyy)*1/6*(Macauley(x1-x3)**3),    ha/(2*G*J)*(0-zshear)*Macaulay(x1-(x2-xa/2))**1, 1/(E*Iyy)*1/6*Macauley(x1-(x2-xa/2))**3-ha/(2*G*J)*(ha/2)*Macaulay(x1-(x2-xa/2))**1,   0,0,x1,1,1  ],
                  [-ha/(2*G*J)*(0-zshear)*Macaulay(x2-x1)**1, -ha/(2*G*J)*(0-zshear)*Macaulay(x2-x2)**1, -ha/(2*G*J)*(0-zshear)*Macaulay(x2-x3)**1,  -1/(E*Iyy)*1/6*(Macauley(x2-x1)**3), -1/(E*Iyy)*1/6*(Macauley(x2-x2)**3), -1/(E*Iyy)*1/6*(Macauley(x2-x3)**3),    ha/(2*G*J)*(0-zshear)*Macaulay(x2-(x2-xa/2))**1, 1/(E*Iyy)*1/6*Macauley(x2-(x2-xa/2))**3-ha/(2*G*J)*(ha/2)*Macaulay(x2-(x2-xa/2))**1,   0,0,x2,1,1  ],
                  [-ha/(2*G*J)*(0-zshear)*Macaulay(x3-x1)**1, -ha/(2*G*J)*(0-zshear)*Macaulay(x3-x2)**1, -ha/(2*G*J)*(0-zshear)*Macaulay(x3-x3)**1,  -1/(E*Iyy)*1/6*(Macauley(x3-x1)**3), -1/(E*Iyy)*1/6*(Macauley(x3-x2)**3), -1/(E*Iyy)*1/6*(Macauley(x3-x3)**3),    ha/(2*G*J)*(0-zshear)*Macaulay(x3-(x2-xa/2))**1, 1/(E*Iyy)*1/6*Macauley(x3-(x2-xa/2))**3-ha/(2*G*J)*(ha/2)*Macaulay(x3-(x2-xa/2))**1,   0,0,x3,1,1  ],
                  [0,0,0,  -1/(E*Iyy)*1/6*(Macauley(x2-xa/2-x1)**3), -1/(E*Iyy)*1/6*(Macauley(x2-xa/2-x2)**3), -1/(E*Iyy)*1/6*(Macauley(x2-xa/2-x3)**3),   0, 1/(E*Iyy)*1/6*Macauley(x2-xa/2-(x2-xa/2))**3,   0,0,x2-xa/2,1,0],
                  [0,0,0, 0,0,0, 1,-tan(theta), 0,0,0,0,0]])
    
    #RES=[R1y,R2y,R3y,R1z,R2z,R3z,Ay,Az,C1,C2,C3,C4,C5]
    RES=np.linalg.solve(XX,Y)
    return RES
    


"Moment around z,y axis as a function of x, as well as torque, shear forces and twist angle as functions of x"


def moments(x,new_nodes_x,torques, RES):
    
    '''
    Torques is the list of torques caused by the aerodynamic loading at each cross section
    '''
    
    index = np.where(new_nodes_x>x)[0][0]
    Td = torques[0:index][-1]
    
    Mz = -D + R1y*Macaulay(x-x1) + R2y*Macaulay(x-x2) - Py*Macaulay(x-(x2+xa/2)) - Ay*Macaulay(x-(x2-xa/2)) + R3y*Macaulay(x-x3)
    My = -R1z*Macaulay(x-x1) + Pz*Macaulay(x-(x2+xa/2)) - R2z*Macaulay(x-x2) - R3z*Macaulay(x-x3)+ Az*Macaulay(x-(x2-xa/2))
    T_total = -Td - (0-zshear)*R1y*Macaulay(x-x1)**0 - (0-zshear)*R2y*Macaulay(x-x2)**0 - (0-zshear)*R3y*Macaulay(x-x3)**0 + (0-zshear)*Py*Macaulay(x-(x2+xa/2))**0 - (ha/2)*Pz*Macaulay(x-(x2+xa/2))**0 - (ha/2)*Az*Macaulay(x-(x2-xa/2))**0 + (0-zshear)*Ay*Macaulay(x-(x2-xa/2))**0
    Sy = -S + R1y*Macaulay(x-x1)**0 + R2y*Macaulay(x-x2)**0 - Py*Macaulay(x-(x2+xa/2))**0 - Ay*Macaulay(x-(x2-xa/2))**0 + R3y*Macaulay(x-x3)**0
    Sz = -R1z*Macaulay(x-x1)**0 - R2z*Macaulay(x-x2)**0 + Pz*Macaulay(x-(x2+xa/2))**0 + Az*Macaulay(x-(x2-xa/2))**0 - R3z*Macaulay(x-x3)**0
    twist = 1/(G*J)*(-DTd - (0-zshear)*R1y*Macaulay(x-x1)**1 - (0-zshear)*R2y*Macaulay(x-x2)**1 - (0-zshear)*R3y*Macaulay(x-x3)**1 + (0-zshear)*Py*Macaulay(x-(x2+xa/2))**1 - (ha/2)*Pz*Macaulay(x-(x2+xa/2))**1 - (ha/2)*Az*Macaulay(x-(x2-xa/2))**1 + (0-zshear)*Ay*Macaulay(x-(x2-xa/2))**1 ) + C5
    return(Mz,My,T_total,Sy,Sz,twist)
    #deflection functions as functions of x in both y and z directions (v in y direction, w in z direction)


def v_prime(x):
    v_prim = -1/(E*Izz)*(q_int + R1y/2*(Macauley(x-x1)**2) + R2y/2*(Macauley(x-x2)**2) - Py/2*Macauley(x-(x2+xa/2))**2 - Ay/2*Macauley(x-(x2-xa/2))**2 + R3y/2*(Macauley(x-x3)**2)) + C1
    return v_prim

def w_prime(x):
    w_prim = -1/(E*Iyy)*(R1z/2*(Macauley(x-x1)**2) + R2z/2*(Macauley(x-x2)**2) - Pz/2*Macauley(x-(x2+xa/2))**2 - Az/2*Macauley(x-(x2-xa/2))**2 + R3z/2*(Macauley(x-x3)**2)) + C1
    return w_prim

def v(x):
    v = -1/(E*Izz)*(q_int + R1y/6*(Macauley(x-x1)**3) + R2y/6*(Macauley(x-x2)**3) - Py/6*Macauley(x-(x2+xa/2))**3 - Ay/6*Macauley(x-(x2-xa/2))**3 + R3y/6*(Macauley(x-x3)**3)) + C1*x + C2
    return v

def w(x):
    w = -1/(E*Iyy)*(R1z/6*(Macauley(x-x1)**3) + R2z/6*(Macauley(x-x2)**3) - Pz/6*Macauley(x-(x2+xa/2))**3 - Az/6*Macauley(x-(x2-xa/2))**3 + R3z/6*(Macauley(x-x3)**3)) + C3*x + C4
    return w   

#Stress calculations, normal stress in cross section y-z as a function of y and z
#do it for all y and z to graph 
def normalstress(x):
    sigma = -My(x)/Iyy*z + Mz(x)/Izz*y
    return sigma
