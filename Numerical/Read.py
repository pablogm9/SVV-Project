from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd





def get_data(aircraft, Ca, La)

    '''
    Reads the provided aerodynamic data and returns it in the shape of a numpy array.
    
    Required inputs:
        - AIRCRAFT: Must be one of the following: A320, CRJ700, Do228, F100
        - Ca: Chordwise length, in [m]
        - La: Chrodwise length, in [m]
    
    Outputs:
        - AERODATA = Matrix containing the provided aerodynamic data. 
                 Numpy matrix of dimension 81 by 41. 
        - COORD_X = X-cordinates of data. Numpy array of
                    length 41.
        - COORD_Z = Z-cordinates of data. Numpy array of
                    length 81.
                    
    '''  

    aircraft = str(aircraft).lower()
    
    # Reading aerodata -- write as 81x41 matrix
    ref_file = open("aerodynamicload"+aircraft+".dat", "r")
    lines = ref_file.readlines()
    ref_file.close()
    
    aerodata = np.mat(np.zeros((81, 41)))
    
    
    for line in lines:
        idx = lines.index(line)
        line = line.replace("\n","")
        values = line.split(',')
        values = np.matrix(values)
        
        aerodata[idx,:] = values
      
    
    # Calculating x and z coordinates of data points 
    theta_x = []  # List of theta_z_i angles
    theta_z = []  # List of theta_x_i angles
    coord_x = []  # List of x coordinates
    coord_z = []  # List of z coordinates
    
    
    N_x = 41
    N_z = 81
    
    for i in range(1,N_x+2):
        theta_i = (i-1)*np.pi/N_x
        theta_x.append(theta_i)
        
    
    for i in range(1,N_x+1):
        x_i = -0.5*(   0.5*La*(1-np.cos(theta_x[i-1])) + 0.5*La*(1-np.cos(theta_x[i]))  )
        coord_x.append(x_i)
    
    
    for i in range(1,N_z+2):
        theta_i = (i-1)*np.pi/N_z
        theta_z.append(theta_i) 
        
        
    for i in range(1,N_z+1):
        z_i = -0.5*(   0.5*Ca*(1-np.cos(theta_z[i-1])) + 0.5*Ca*(1-np.cos(theta_z[i]))  )
        coord_z.append(z_i)
      
    
    # Make all x-coordinates positive (positive x-axis starting from root towards tip )
    coord_x = np.abs(coord_x)
    coord_z = np.abs(coord_z)
    
    return aerodata,coord_x,coord_z
    
        

  

 

    
    
    
    
    
    
    
    
    
    