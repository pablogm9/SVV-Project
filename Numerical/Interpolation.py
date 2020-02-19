#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:15:16 2020

@author: pablo


"""


from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


'''

UNIVARIATE INTERPOLATION FUNCTIONS.

Direction of interpolation is denoted by x.

Number of data points = n+1
Number of intervals = n
Note that the indexing starts at zero, i.e. the first node is located at x_0


'''




def find_interpolants(nodes,data):
       
    '''
    Finds the cubic splines at each of the n nodes. Each cubic spline has the format:
        
        a_i(x-x_i)^3 + b_i(x-x_i)^2 + c_i(x-x_i) + d_i
    
    
    Required inputs:
        - NODES: X positions of the (n+1) nodes where the function values are known.
                 Numpy array - length of n+1.
        - DATA: Function values at the specified nodes.
                Numpy array - length of n+1.
    
    Outputs:
        - INTERPOLANTS = Matrix containing cubic spline parameters a_i,b_i,c_i,d_i for 
                         each of the n intervals. Numpy matrix of dimension n by 4.  
                         
    '''    
            
    n = len(nodes)-1

    # Function values at nodes
    f = data
    
    
    # Compute spacings h where h_i = x_i+1 - x_1
    h = np.diff(nodes) 
    
    
    # Initializing matrix
    A = np.zeros((n+1,n+1))
    shape = A.shape
    
    
    # Construct interior rows (all but the 0th and the nth row)
    for i in range(1,n):
        A[i,i-1] = h[i-1]/6 
        A[i,i] = (h[i-1]+h[i])/3 
        A[i,i+1] = 1/6*h[i] 
        
    
    # Enforcing natural boundary conditions on matrix
    A[0,0]=1
    A[n,n]=1
    
    
    # Initializing RHS vector
    b = np.zeros(n+1)
    
    # Constructing interior values of RHS vector
    for i in range(1,n):
        b[i] = (f[i+1]-f[i])/h[i] - (f[i]-f[i-1])/h[i-1]
    
    # Enforcing natural boundary conditions on RHS vector
    b[0] = 0
    b[n] = 0
        
    # Solve system of equations for M_i's
    # Ax = b --> find x, vector of M_i's (denoted by M)
    M = np.linalg.tensorsolve(A,b)
    
    # Initialize results matrix, which contains a_i,b_i,c_i,d_i for 
    # i={0,1,....,n-2,n-1} -- total of n interpolants
    interpolants = np.zeros((n,4))
    
    
    for i in range(0,n):
        a_i = (M[i+1]-M[i])/(6*h[i])
        b_i = M[i]/2
        c_i = (f[i+1]-f[i])/h[i] - h[i]*M[i]/3 - h[i]*M[i+1]/6
        d_i = f[i]
        row = np.array([a_i,b_i,c_i,d_i])
        interpolants[i,:] = row
        
    return interpolants
    
    


def new_loading(nodes, interpolants, resolution):
    
    '''
    Computes the aerodynamic loading at the new resolution 
    using the interpolants from find_interpolants().
    
    
    Required inputs:
        - INTERPOLANTS: Result from find_interpolants(). 
                        Numpy array - dimension n x 4
        - NODES: X positions of the (n+1) nodes where the function values are known.
                 Numpy array - length of n+1.
        - RESOLUTION: Minimum resolution of the new discretized grid. 
                      Must be given in mm. 
    
    Outputs:
        - NEW_NODES = Numpy array containing the new nodes. Length = (n+1)+(n)(n_internal)
                      Thus, the number of new nodes depends on the RESOLUTION.
        - NEW_LOADING = Array containing the loading at each of the
                    discretized stations. 
                    Numpy array of length n_x
    '''
    
    # Convert resolution from mm to m
    resolution = resolution/1000
    
    n = len(nodes)-1
    h = np.diff(nodes)
    max_spacing = np.amax(h)
    
    # Required number of points per spline to obtain the given 
    # resolution at the most critical (i.e. largest) discretization
    # step. By using this number of points per spline, 
    # the actual obtained resolution will always be equal or GREATER than
    # the provided resolution. 
    n_spline = np.ceil(max_spacing/resolution) + 1 
    
    # Number of internal nodes, for reference.
    n_internal = n_spline - 2
    
    
    new_nodes = np.array([])
    new_loading = np.array([])
    
    # Loop through all splines
    for i in range(0,n):
        
        # Get interpolant coefficients
        a_i = interpolants[i,0]
        b_i = interpolants[i,1]
        c_i = interpolants[i,2]
        d_i = interpolants[i,3]
        
        
        # Generate new nodes
        x_b = nodes[i]
        x_f = nodes[i+1]
       
        interval_nodes = np.linspace(x_b,x_f,num=int(n_spline))
        
        if i != 0:
            interval_nodes = interval_nodes[1:int(n_spline)]
            
        new_nodes = np.append(new_nodes,interval_nodes)
    
        
        # Compute new loads at all discretized points
        for i in range(len(interval_nodes)):
            load = a_i*(interval_nodes[i]-x_b)**3 + b_i*(interval_nodes[i]-x_b)**2 + c_i*(interval_nodes[i]-x_b) + d_i
            new_loading = np.append(new_loading,load) 
        
        
    return new_nodes,new_loading
            
        
        
    











