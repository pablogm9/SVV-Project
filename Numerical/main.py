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


from Numerical import Interpolation as inter
from Numerical import Read
from Numerical import Integration
from Numerical import SectionalProperties as section

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

aircraft = 'F100'





# ----------------- SECTIONAL PROPERTIES  -----------------

Izz, Iyy = section.get_MoI(ha,Ca,t_st,w_st,h_st,A_st,t_sk,t_sp,nstiff) # Moments of Inertia
z_sc = section.get_SC(t_st, w_st, h_st, t_sk, t_sp, ha, Ca, Izz) # Z-location of shear center



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

torques = np.array([])

for i in range(new_aerodata.shape[1]):
    
    q = new_aerodata[:,i]
    qz = np.multiply(q,new_nodes_z)
    
    interpolants_q = inter.find_interpolants(new_nodes_z,q)
    interpolants_qz = inter.find_interpolants(new_nodes_z,qz)
    
    numerator = Integration.Analytical_Int_1D(interpolants_qz, new_nodes_z).integrator()
    denominator = Integration.Analytical_Int_1D(interpolants_q, new_nodes_z).integrator()

    z_centroid = numerator/denominator
    resultant = denominator
    
    torque = abs(resultant)*(z_centroid-z_sc)
    
    centroids_spanwise = np.append(centroids_spanwise,z_centroid)
    resultant_forces = np.append(resultant_forces,resultant)
    torques = np.append(torques,torque)


# ----------------- TORQUE CALCULATIONS -----------------
# Loop through spanwise cross sections, calculate resultant force and centroid.









# Print runtime
print('Runtime: %f seconds' % (time.time()-start_time))