#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 19:49:50 2020

@author: pablo
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 15:01:21 2020

@author: pablo
"""



from math import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import time

import Loading
import Interpolation as inter
import Read
#import final_integration as integration
import SectionalProperties as section

# Start timer to measure runtime
start_time = time.time()


# ----------------- PARAMETERS -----------------
Ca = 0.505 #[m]
La = 1.611 #[m]
ha = 0.161 #[m]
z_hingeline = ha/2 #[m]
t_st = 0.0012 #[m]
w_st = 0.017 #[m]
h_st = 0.013 #[m]
A_st = t_st * (w_st + h_st) #[m^2]
t_sk = 0.0011 #[m]
t_sp = 0.0024 #[m]
nstiff = 11 #[-]

x_1 = 0.125
x_2 = 0.498
x_3 = 1.494 #m
x_a = 0.245    #m
x_4 = x_2-0.5*x_a
theta=radians(30)   #degrees 

r = 0.5*ha
A_circ = r**2*pi/2
A_triang = (Ca-r)*r


Py = 49200*sin(theta)
Pz = 49200*cos(theta)   #N
d1=0.00389   #m
d2=0
d3=0.01245    #m
G = 28*10**9    #in Pa
E=73.1*10**9     #in Pa

aircraft = 'F100'




# ----------------- SECTIONAL PROPERTIES  -----------------

Izz, Iyy, z_centroid = section.get_MoI(ha,Ca,t_st,w_st,h_st,A_st,t_sk,t_sp,nstiff) # Moments of Inertia
z_sc = section.get_SC(t_sk, t_sp, ha, Ca, Izz) # Z-location of shear center
J = section.get_J(G, ha, t_sk, w_st, t_sp, t_st, Ca)


# ----------------- AERODYNAMIC DATA -----------------

print('\n')
print('\n')
print('Extracting aerodynamic data from .txt file...')
print('\n')
print('\n')

# Get original aerodynamic data (before interpolation)
original_aerodata,x_coord,z_coord = Read.get_data(aircraft, Ca, La)
original_aerodata = original_aerodata*1000 #Transform aerodynamic data from kPa to Pa
original_shape = original_aerodata.shape # Shape of data in the form (81,41)



# ----------------- INTERPOLATION STEP 0 -----------------
# Obtain preliminary values for interpolation

# Desired resolutions, in [mm]
r_span = 2
r_chord = 5

# Get value of n for both chord and span (= number of nodes - 1)
n_chord = original_shape[0] - 1 # Should be 80
n_span = original_shape[1] - 1 # Should be 40


# Get number of points PER SPLINE in both directions from resolutions
n_points_chord = inter.get_number_of_points(z_coord, r_chord)
n_points_span = inter.get_number_of_points(x_coord, r_span)

n_internal_chord = n_points_chord -2
n_internal_span = n_points_span - 2

total_points_chord = (n_chord + 1) + (n_internal_chord)*n_chord
total_points_span = (n_span + 1) + (n_internal_span)*n_span


# Initialize matrix with interpolated nodes
new_nodes_x = np.array([])
new_nodes_z = np.array([])


# ----------------- INTERPOLATION STEP 1 -----------------
# (Go from 81x41 to 500x41)
print('\n')
print('\n')

# Initialize matrix with intermidiate aerodynamic data (1/2 interpolation steps) 
intermidiate_aerodata = np.zeros((total_points_chord,original_shape[1]))


#Loop through spanwise cross-sections
for i in range(n_span+1):
    
    percentage=round(50*i/n_span,2)
    
    if percentage%5==0:
        print('Interpolating aerodynamic data... '+str(percentage)+'%'+'\n')

    
    # Select individual spanwise cross section from original aerodymic data
    cross_sectional_loading = original_aerodata[:,i]
    
    # Interpolate in chordwise direction
    cross_sectional_interpolants = inter.find_interpolants(z_coord,cross_sectional_loading)

    
    # New cross sectional values and nodes after interpolation
    new_cross_sectional_nodes, new_cross_sectional_loads = inter.new_loading(z_coord,cross_sectional_interpolants,n_points_chord)
    
    # Keep record of new z coordinates only once
    if i==0:
         new_nodes_z = new_cross_sectional_nodes
    
    # Append to matrix
    intermidiate_aerodata[:,i] = new_cross_sectional_loads
    
    
    
    
# ----------------- INTERPOLATION STEP 2 -----------------
# (Go from 500x41 to 500x100)
    
# Initialize matrix with final aerodynamic data (2/2 interpolation steps)     
new_aerodata = np.zeros((total_points_chord,total_points_span))


#Loop through chordwise cross-sections
for i in range(intermidiate_aerodata.shape[0]):
    
    percentage = 50 + 50*i/(intermidiate_aerodata.shape[0]-1)
    
    if percentage%5==0:
        print('Interpolating aerodynamic data... '+str(round(percentage,2))+'%'+'\n')

    
    # Select individual chordwise cross section from intermidiate aerodynamic data
    cross_sectional_loading = intermidiate_aerodata[i,:]
    
    # Interpolate in spanwise direction
    cross_sectional_interpolants = inter.find_interpolants(x_coord,cross_sectional_loading)
    
    # New cross sectional values and nodes after interpolation
    new_cross_sectional_nodes, new_cross_sectional_loads = inter.new_loading(x_coord,cross_sectional_interpolants,n_points_span)

    # Keep record of new x coordinates only once
    if i==0:
         new_nodes_x = new_cross_sectional_nodes
    
    # Append to matrix
    new_aerodata[i,:] = new_cross_sectional_loads
  


# ----------------- COORD. SYSTEM ADJUSTMENTS -----------------
# Flip z-axis direction + translate to match chosen coordinate system 
# By doing this, the positive z-axis points towards the LE and starts
# at the hinge line.

new_nodes_z = -1*new_nodes_z + z_hingeline



# ----------------- RESULTANT LOAD CALCULATIONS -----------------
# Loop through spanwise cross sections, calculate resultant force and centroid.

# Z-coordinates of centroids of each spanwise cross section
centroids_spanwise = np.array([])

# Resultant forces PER UNIT SPAN at each spanwise cross section
resultant_forces = np.array([])

# Torque at each spanwise crosssection
torques = np.array([])
#torques2=np.array([])

# Aerodynamic force generated at each point
aero_loads = np.zeros(new_aerodata.shape)

# Spanwise spacing -- steps in x direction 
delta_x = np.array([])

print('\n')
print('\n')
# Loop through cross sections
for i in range(new_aerodata.shape[1]):
    
    percentage = 100*i/(new_aerodata.shape[1]-1)
    
    if percentage%5==0:                    
        print('Computing resultant forces at grid nodes... '+str(round(percentage,2))+'%\n')

    
    q = new_aerodata[:,i] # Pressure distribution

    # ADD -- Convert q to force by multiplying by d_x,d_z (area in between spanwise crossections)
    if i ==0: 
        delta_span = new_nodes_x[2]-new_nodes_x[1]
        delta_x = np.append(delta_x,delta_span)
    else:
        delta_span = new_nodes_x[i] - new_nodes_x[i-1]
        delta_x = np.append(delta_x,delta_span)
      
    torque = 0    
    resultant = 0
        
    for j in range(new_aerodata.shape[0]):
        
        if j == 0: 
            delta_chord = abs(new_nodes_z[j]-z_hingeline)
        else:
            delta_chord = abs(new_nodes_z[j] - new_nodes_z[j-1])
            
        area = delta_span*delta_chord
        
        load = q[j]*area
        aero_loads[j,i] = load
        
        torque_i = load*(new_nodes_z[j]-z_sc)
        torque += torque_i
        
        resultant += load
    
    centroid = np.sum(np.multiply(aero_loads[:,i],new_nodes_z))/resultant
    centroids_spanwise = np.append(centroids_spanwise,centroid)
    #torques2=np.append(torques2, resultant*centroid)    
    resultant_forces = np.append(resultant_forces,resultant)
    torques = np.append(torques,torque)

          
          


'''

# ----------------- PLOT COUNTOUR MAP AND 3D SURFACE OF INTERPOLATED DATA ----------------- 












'''    
    
# ----------------- COMPUTING INTEGRALS REQUIRED TO SOLVE S.O.E. ----------------- 

# S -- First integral of q w.r.t x
total_resultants = np.cumsum(resultant_forces)
total_resultant_force = total_resultants[-1]

# D -- Second integral of q w.r.t x
# [F_0 * x_0, ..., F_n * x_n]
aero_moment_z_contributions = np.multiply(resultant_forces,new_nodes_x)
aero_moment_z = np.cumsum(aero_moment_z_contributions)
#D = integration.integrate(new_nodes_x,resultant_forces,'double')

#Td is the last torque
Td_contributions = torques
Td = np.cumsum(Td_contributions)
    
#Indicies for DTd1, DTd2, DTd3
index_x1 = np.where(new_nodes_x>x_1)[0][0]
index_x2 = np.where(new_nodes_x>x_2)[0][0]
index_x3 = np.where(new_nodes_x>x_3)[0][0]
index_x4 = np.where(new_nodes_x>x_4)[0][0]

x1_nodes = new_nodes_x[0:index_x1]
x2_nodes = new_nodes_x[0:index_x2]
x3_nodes = new_nodes_x[0:index_x3]
x4_nodes = new_nodes_x[0:index_x4]

# DTd -- First integral of T w.r.t x
DTd_contributions = np.multiply(Td_contributions,delta_x)
DTd = np.cumsum(DTd_contributions)

# DTd1, DTd2, DTd3
DTd1 = DTd[index_x1]
DTd2 = DTd[index_x2]
DTd3 = DTd[index_x3]
DTd4 = DTd[index_x4]

# Third integral of q w.r.t x
aero_deflection_slope_contributions =  np.multiply(aero_moment_z,delta_x)
aero_deflection_slope = np.cumsum(aero_deflection_slope_contributions)


# FI1  -- Quadruple integral of q w.r.t x
aero_deflection_contributions = np.multiply(aero_deflection_slope,delta_x)
aero_deflection = np.cumsum(aero_deflection_contributions)


FI1 = aero_deflection[index_x1]
FI2 = aero_deflection[index_x2]
FI3 = aero_deflection[index_x3]


# ----------------- SOLVING S.O.E. FOR REACTION FORCES AND INTEGRATION CONSTANTS -----------------



#RES=[R1y,R2y,R3y,R1z,R2z,R3z,Ay,Az,C1,C2,C3,C4,C5]
RES = Loading.matrix_solver(Ca,La,x_1,x_2,x_3,x_a,theta,t_st,t_sk,t_sp,w_st,ha,Py,Pz,d1,d2,d3,G,E,Izz,Iyy,z_sc,total_resultant_force,aero_moment_z[-1],Td[-1],DTd1,DTd2,DTd3,DTd4,FI1,FI2,FI3)

R1y = RES[0]
R2y = RES[1]
R3y = RES[2]
R1z = RES[3]
R2z= RES[4] 
R3z= RES[5]
Ay = RES[6]
Az= RES[7]
C1 = RES[8]
C2= RES[9]
C3= RES[10]
C4 = RES[11]
C5 = RES[12]




# ----------------- COMPUTE MOMENTS AND DEFLECTIONS AS A 'FUNCTION' OF X -----------------
print('\n')
print('\n')
M_y = np.array([])
M_z = np.array([])
T_x = np.array([]) #Note that this torque includes the torque due to aero loading and that due to reaction forces
S_y = np.array([])
S_z = np.array([])
V_y = np.array([])
V_y_prime = np.array([])
W_z = np.array([])
W_z_prime = np.array([])
twist = np.array([])

# Loop through spanwise cross sections
for i in range(new_aerodata.shape[1]):
    
    percentage = 100*i/(new_aerodata.shape[1]-1)
                        
    if percentage%5==0:
        print('Computing shear forces, moments, displacements, and twist... '+str(round(percentage,2))+'%\n')

        
    #Calculate required integrals
    D_i = aero_moment_z[i]
    Td_i = Td[i]
    DTd_i = DTd[i]
    S_i = total_resultants[i]
    aero_deflection_slope_i = aero_deflection_slope[i]
    FI_i = aero_deflection[i]
         
    # Moments and torques
    M_z_i = Loading.moment_z(new_nodes_x[i], R1y,  R2y, R3y, Ay, Py, D_i,x_1,x_2,x_3,x_a)
    M_y_i = Loading.moment_y(new_nodes_x[i], R1z, R2z, R3z, Az, Pz,x_1,x_2,x_3,x_a)
    T_x_i = Loading.torque_x(new_nodes_x[i], Td_i, z_sc, R1y, R2y, R3y, Ay, Az, Py, Pz,x_1,x_2,x_3,x_a,ha)
    
    # Shear forces
    S_y_i = Loading.shear_y(new_nodes_x[i], R1y, R2y, R3y, Ay, Py, S_i,x_1,x_2,x_3,x_a)
    S_z_i = Loading.shear_z(new_nodes_x[i], R1z, R2z, R3z, Az, Pz, x_1,x_2,x_3,x_a)
    
    # Deflections
    V_y_i = Loading.v(new_nodes_x[i],E,Izz,FI_i,R1y,R2y,R3y,Py,Ay,C1,C2,x_1,x_2,x_3,x_a)
    V_y_prime_i = Loading.v_prime(new_nodes_x[i], E, Izz, aero_deflection_slope_i, R1y, R2y, R3y, Py, Ay, C1,x_1,x_2,x_3,x_a)
    
    W_z_i = Loading.w(new_nodes_x[i],E,Iyy,R1z,R2z,R3z,Pz,Az,C3,C4,x_1,x_2,x_3,x_a)
    W_z_prime_i = Loading.w_prime(new_nodes_x[i], E, Iyy, R1z, R2z, R3z, Pz, Az, C3,x_1,x_2,x_3,x_a)
    
    # Twist
    
    twist_i = Loading.twist(new_nodes_x[i], G, J, DTd_i, R1y, R2y, R3y, Py, Pz, Ay, Az, C5,x_1,x_2,x_3,x_a,ha,z_sc)

    
    # Appending to arrays
    M_z = np.append(M_z,M_z_i)
    M_y = np.append(M_y,M_y_i)
    T_x = np.append(T_x,T_x_i)
    S_y = np.append(S_y,S_y_i)
    S_z = np.append(S_z,S_z_i)
    V_y = np.append(V_y,V_y_i)
    V_y_prime = np.append(V_y_prime,V_y_prime_i)
    W_z = np.append(W_z,W_z_i)
    W_z_prime = np.append(W_z_prime,W_z_prime_i)
    twist = np.append(twist,twist_i)
    
    


'''
# ----------------- GENERATE DEFLECTION PLOTS -----------------
print(np.max(twist))
plt.figure()
plt.plot(new_nodes_x,T_x)
plt.figure()
plt.plot(new_nodes_x,twist)
plt.show()

'''


# ----------------- STRESS CALCULATIONS -----------------

# Get discretization of cross section -- (z,y) coordiantes of points
circle1,circle2,spar1,spar2,top,bottom = section.discretize_crosssection(1000)

circle_1_z = circle1[0]
circle_1_y = circle1[1]

circle_2_z = circle2[0]
circle_2_y = circle2[1]

spar_1_z = spar1[0]
spar_1_y = spar1[1]

spar_2_z = spar2[0]
spar_2_y = spar2[1]

top_z = top[0]
top_y = top[1]

bottom_z = bottom[0]
bottom_y = bottom[1]

crosssection_z = np.concatenate((circle_1_z,spar_1_z,top_z,bottom_z,spar_2_z,circle_2_z))
crosssection_y = np.concatenate((circle_1_y,spar_1_y,top_y,bottom_y,spar_2_y,circle_2_y))


# Get shear flows due to shear forcesat each points
q_circle1,q_spar1,q_top,q_bottom,q_spar2,q_circle2 = section.get_shearflows()

# Array containing one subarray per cross section. Each subrray contains 60000 points
normal_stresses = np.array([])
shear_stresses = np.array([])
vm_stresses = np.array([])
'''
# Loop through spanwise crossections
for i in range(new_aerodata.shape[1]):
    
    percentage = 100*i/(new_aerodata.shape[1]-1)
    
    if percentage%5==0:
        print('Calculating shear flow and Von Mises stress distributions... '+str(round(percentage,2))+'%\n')
    
    
    normal_stresses_crosssection = np.array([])
    shear_stresses_crosssection = np.array([])
    vm_stresses_crosssection = np.array([])
    
    #V_y = V_y[i]
    #W_z = W_z[i]
    #T_x = T_x[i]
    
    
    # ------- SHEAR STRESSES -------
    
    # Compute shear flows due to torque
    # Effect of stiffeners excluded in calculation
    Y=[T_x[i],0]
    A=np.matrix([[2*A_circ, 2*A_triang],[1/(2*A_circ*G)*((pi*ha/2)/t_sk) + ha/(2*A_triang*G*t_sp) + ha/(2*A_circ*G*t_sp), - 1/(2*A_triang*G)*((2*sqrt((ha/2)**2+(Ca-ha/2)**2))/t_sk) - ha/(2*A_triang*G*t_sp) - ha/(2*A_circ*G*t_sp)]])
    shears=np.linalg.solve(A,Y)
    
    q_circle0=shears[0]
    q_triangle0=shears[1]
    
    # ORDER: Circle 1, Spar 1, Top, Bottom, Spar 2, Circle 2
    
    # Circle 1
    for j in range(len(circle_1_z)):
        z = circle_1_z[j]
        y = circle_1_y[j]
        
        tau = (q_circle1[j]+q_circle0)/t_sk
        shear_stresses_crosssection = np.append(shear_stresses_crosssection,tau)

    # Spar 1
    for j in range(len(spar_1_z)):
        z = spar_1_z[j]
        y = spar_1_y[j]
        
        tau = (q_spar1[j]-q_circle0+q_triangle0)/t_sp
        shear_stresses_crosssection = np.append(shear_stresses_crosssection,tau)
        
    # Top
    for j in range(len(top_z)):
        z = top_z[j]
        y = top_y[j]
        
        tau = (q_top[j]+q_triangle0)/t_sk
        shear_stresses_crosssection = np.append(shear_stresses_crosssection,tau)
       
    # Bottom
    for j in range(len(bottom_z)):
        z = bottom_z[j]
        y = bottom_y[j]
        
        tau = (q_bottom[j]+q_triangle0)/t_sk
        shear_stresses_crosssection = np.append(shear_stresses_crosssection,tau)
        
     # Spar 2
    for j in range(len(spar_2_z)):
        z = spar_2_z[j]
        y = spar_2_y[j]
        
        tau = (q_spar2[j]-q_circle0+q_triangle0)/t_sp
        shear_stresses_crosssection = np.append(shear_stresses_crosssection,tau) 
        
    # Circle 2
    for j in range(len(circle_2_z)):
        z = circle_2_z[j]
        y = circle_2_y[j]
        
        tau = (q_circle2[j]+q_circle0)/t_sk
        shear_stresses_crosssection = np.append(shear_stresses_crosssection,tau)    
    
    
    # ------- NORMAL STRESSES -------
    
    for j in range(len(crosssection_z)):
        
        z = crosssection_z[j]
        y = crosssection_y[j]
        
        sigma_xx = M_z[i]*y/Izz + M_y[i]*(z-z_centroid)/Iyy
        
        normal_stresses_crosssection = np.append(normal_stresses_crosssection,sigma_xx)
        
    
    # ------- VON MISES STRESSES -------
    
    for i in range(len(crosssection_z)):
        vm = sqrt((normal_stresses_crosssection[i] ** 2) * (3*(shear_stresses_crosssection[i] ** 2)))
        vm_stresses_crosssection = np.append(vm_stresses_crosssection,vm)
        
        
    shear_stresses = np.append(shear_stresses,shear_stresses_crosssection)      
    normal_stresses = np.append(normal_stresses,normal_stresses_crosssection) 
    vm_stresses = np.append(vm_stresses,vm_stresses_crosssection)       
 


vm_GPa = 1e-9*vm_stresses


# Print runtime
print('\n')
print('\n')
print('Complete.'+'\n\n'+'Runtime: %f seconds\n' % (time.time()-start_time))

    

# ----------------- PLOT CROSS-SECTION -----------------

# Cross section
plt.plot(circle_1_z,circle_1_y,'k',linewidth=6)
plt.plot(circle_2_z,circle_2_y,'k',linewidth=6)
plt.plot(spar_1_z,spar_1_y,'k',linewidth=6)
plt.plot(spar_2_z,spar_2_y,'k',linewidth=6)
plt.plot(top_z,top_y,'k',linewidth=6)
plt.plot(bottom_z,bottom_y,'k',linewidth=6)

# Plot origin
plt.plot(0,0,'k.')

# Plot axes
plt.axhline(0,linestyle='dashed',color='black',linewidth=1) # x = 0
plt.axvline(0,linestyle='dashed',color='black',linewidth=1) # y = 0


# Plot heatmap of von mises stress



# Set ranges for axes and show plot
plt.xlim(0.1,-0.45)
plt.ylim(-0.1,0.1)
plt.show()
'''












