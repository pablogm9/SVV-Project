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
import Integration
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


x1 = 0.125
x2 = 0.498
x3 = 1.494 #m
xa = 0.245    #m
x4 = x2-0.5*xa
theta=30   #degrees 


Py = 49200*sin(theta)
Pz = 49200*cos(theta)   #N
d1=0.00389   #m
d2=0
d3=0.01245    #m
G = 28*10**9    #in Pa
E=73.1*10**9     #in Pa

aircraft = 'F100'



# ----------------- SECTIONAL PROPERTIES  -----------------

Izz, Iyy = section.get_MoI(ha,Ca,t_st,w_st,h_st,A_st,t_sk,t_sp,nstiff) # Moments of Inertia
z_sc = section.get_SC(t_sk, t_sp, ha, Ca, Izz) # Z-location of shear center



# ----------------- AERODYNAMIC DATA -----------------

# Get original aerodynamic data (before interpolation)
original_aerodata,x_coord,z_coord = Read.get_data(aircraft, Ca, La)
original_aerodata = original_aerodata*1000 #Transform aerodynamic data from kPa to Pa
original_shape = original_aerodata.shape # Shape of data in the form (81,41)



# ----------------- INTERPOLATION STEP 0 -----------------
# Obtain preliminary values for interpolation

# Desired resolutions, in [mm]
r_span = 5
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


# Initialize matrix with intermidiate aerodynamic data (1/2 interpolation steps) 
intermidiate_aerodata = np.zeros((total_points_chord,original_shape[1]))


#Loop through spanwise cross-sections
for i in range(n_span+1):
    
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

centroids_spanwise = np.array([])
resultant_forces = np.array([])
aero_loads = np.zeros(new_aerodata.shape)
torques = np.array([])

for i in range(new_aerodata.shape[1]):
    
    q = new_aerodata[:,i] # Pressure distribution

    # ADD -- Convert q to force by multiplying by d_x,d_z (area in between spanwise crossections)

    if i ==0: 
        delta_span = new_nodes_x[i]
    else:
        delta_span = new_nodes_x[i] - new_nodes_x[i-1]
      
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
    #new_nodes_zn = new_nodes_z[::-1]
    #resultant = Integration.Analytical_Int_1D(inter.find_interpolants(new_nodes_z,aero_loads[:,i]), new_nodes_z).integrator()
    resultant_forces = np.append(resultant_forces,resultant)
    torques = np.append(torques,torque)
    
    
    '''
    numerator = Integration.Analytical_Int_1D(interpolants_qz, new_nodes_z).integrator()
    denominator = Integration.Analytical_Int_1D(interpolants_q, new_nodes_z).integrator()

    z_centroid = numerator/denominator
    resultant = denominator
    
    torque = abs(resultant)*(z_centroid-z_sc)
    
    centroids_spanwise = np.append(centroids_spanwise,z_centroid)
    resultant_forces = np.append(resultant_forces,resultant)
    torques = np.append(torques,torque)
    '''
    
    
# ----------------- COMPUTING INTEGRALS REQUIRED TO SOLVE S.O.E. -----------------

# S
total_resultant = np.sum(resultant_forces) 

# D
aero_moment_z = np.sum(np.multiply(resultant_forces,new_nodes_x))
#aero_moment_z = Integration.Analytical_Int_1D(inter.find_interpolants(new_nodes_x,resultant_forces),new_nodes_x).double_integrator()

#Td is the last torque
torque_L = torques[-1] 


#DTd1 to DTd4
index_x1 = np.where(new_nodes_x>x1)[0][0]
index_x2 = np.where(new_nodes_x>x2)[0][0]
index_x3 = np.where(new_nodes_x>x3)[0][0]
index_x4 = np.where(new_nodes_x>x4)[0][0]


#DTd's
aero_twist1 = np.sum(np.multiply(torques[0:index_x1],new_nodes_x[0:index_x1]))
#aero_twist1 = Integration.Analytical_Int_1D(inter.find_interpolants(new_nodes_x[0:index_x1],torques[0:index_x1]),new_nodes_x[0:index_x1]).double_integrator()

#DTd2
aero_twist2 = np.sum(np.multiply(torques[0:index_x2],new_nodes_x[0:index_x2]))
#aero_twist2 = Integration.Analytical_Int_1D(inter.find_interpolants(new_nodes_x[0:index_x2],torques[0:index_x2]),new_nodes_x[0:index_x2]).double_integrator()

#DTd3
aero_twist3 = np.sum(np.multiply(torques[0:index_x3],new_nodes_x[0:index_x3]))
#aero_twist3 = Integration.Analytical_Int_1D(inter.find_interpolants(new_nodes_x[0:index_x3],torques[0:index_x3]),new_nodes_x[0:index_x3]).double_integrator()

#DTd4
aero_twist4 = np.sum(np.multiply(torques[0:index_x4],new_nodes_x[0:index_x4]))
#aero_twist4 = Integration.Analytical_Int_1D(inter.find_interpolants(new_nodes_x[0:index_x4],torques[0:index_x4]),new_nodes_x[0:index_x4]).double_integrator()


# FI's

FI1 = np.sum(np.multiply(resultant_forces[0:index_x1],np.power(new_nodes_x[0:index_x1],3)))
#FI1 = Integration.Analytical_Int_1D(inter.find_interpolants(new_nodes_x[0:index_x1],torques[0:index_x1]),new_nodes_x[0:index_x1]).quad_integrator()

FI2 = np.sum(np.multiply(resultant_forces[0:index_x2],np.power(new_nodes_x[0:index_x2],3)))
#FI2 = Integration.Analytical_Int_1D(inter.find_interpolants(new_nodes_x[0:index_x2],torques[0:index_x2]),new_nodes_x[0:index_x2]).quad_integrator()

FI3 = np.sum(np.multiply(resultant_forces[0:index_x3],np.power(new_nodes_x[0:index_x3],3)))
#FI3 = Integration.Analytical_Int_1D(inter.find_interpolants(new_nodes_x[0:index_x3],torques[0:index_x3]),new_nodes_x[0:index_x3]).quad_integrator()


# ----------------- SOLVING S.O.E. FOR REACTION FORCES AND INTEGRATION CONSTANTS -----------------


#RES=[R1y,R2y,R3y,R1z,R2z,R3z,Ay,Az,C1,C2,C3,C4,C5]
RES = Loading.matrix_solver(Ca,La,x1,x2,x3,xa,theta,t_st,t_sk,t_sp,w_st,ha,Py,Pz,d1,d2,d3,G,E,Izz,Iyy,z_sc,total_resultant,aero_moment_z,torque_L,aero_twist1,aero_twist2,aero_twist3,aero_twist4,FI1,FI2,FI3)






# Print runtime
print('Runtime: %f seconds' % (time.time()-start_time))
