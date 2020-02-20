from math import *
import numpy as np
#This should be used along a chord of a given airfoil
#a,b,c,d should be passed as arrays and n (number of points) as an int
#c_length is the chord length
'''
This class is to be used to analytically integrate between points
INPUTS: matrix: a matrix of a,b,c,d values for the cubic splines between each of the pairs of the coords also given.
                OR a matrix of a,b,c,d,e values for the quartic splines between each of the pairs of the coords given.
                the double_integrator and quad_integrator functions integrate the cubic spline 2 times or 4 times
        coords: an array of the coordinates at which the aerodynamic data was taken for the a,b,c,d cubic splines(array)
OUTPUTS: The sum of the integrated portions for the entire length of the span or the chord over which you have integrated(array)
        OR the sum of the integrated portions for the entire integral that you wish 
        The output of the double integrator function and the quad integrator function both output a single float that 
        is the integrated value
        
'''

class Analytical_Int_1D:
    def __init__(self,matrix,coords):
        self.matrix = matrix
        if len(matrix[0]) == 4:

            self.a = np.array([arr[0] for arr in matrix])
            self.b = np.array([arr[1] for arr in matrix])
            self.c = np.array([arr[2] for arr in matrix])
            self.d = np.array([arr[3] for arr in matrix])

        if len(matrix[0]) == 5:
            self.a = [arr[0] for arr in matrix]
            self.b = [arr[1] for arr in matrix]
            self.c = [arr[2] for arr in matrix]
            self.d = [arr[3] for arr in matrix]
            self.e = [arr[4] for arr in matrix]
        self.coords = coords

    #x_0 and x_1 are the two x values on the chord between which you are integrating
    def integrator(self):
        x_1 = (self.coords[1:])
        x_0 = (self.coords[:-1])
        if len(self.matrix[0]) == 4:
            return sum((self.a*(x_1)**4)/4 + (self.b*(x_1)**3)/3 + (self.c*(x_1)**2)/2 + self.d*(x_1) - (self.a*(x_0)**4)/4 + (self.b*(x_0)**3)/3 + (self.c*(x_0)**2)/2 + self.d*(x_0))
        if len(self.matrix[0]) == 5:
            return sum((self.a*(x_1)**5)/5 + (self.b*(x_1)**4)/4 + (self.c*(x_1)**3)/3 + (self.d*(x_1)**2)/2 + self.e(x_1)-(self.a*(x_0)**5)/5 + (self.b*(x_0)**4)/4 + (self.c*(x_0)**3)/3 + (self.d*(x_0)**2)/2 + self.e(x_0))

    def double_integrator(self):
        x_1 = (self.coords[1:])
        x_0 = (self.coords[:-1])
        self.e = (self.a*x_0**3) + (self.b*x_0**2) + (self.c*x_0) + self.d
        return sum((self.a*.05*(x_1)**5) + ((self.b/12)*(x_1)**4) + ((self.c/6)*(x_1)**3) + ((self.d*0.5)*(x_1)**2) - (self.a*.05*(x_0)**5) + ((self.b/12)*(x_0)**4) + ((self.c/6)*(x_0)**3) + ((self.d*0.5)*(x_0)**2))
    def quad_integrator(self):
        x_1 = (self.coords[1:])
        x_0 = (self.coords[:-1])
        x_d = x_1-x_0
        self.e = (self.a*x_0**3) + (self.b*x_0**2) + (self.c*x_0) + self.d
        self.f = (self.a/4)*x_0**4 + (self.b/3)*x_0**3 + (self.c/2)*x_0**2 + self.d*x_0 + self.e
        self.g = (self.a/20)*x_0**5 + (self.b/12)*x_0**4 + (self.c/6)*x_0**3 + (self.d/2)*x_0**2 + self.e*x_0 + self.f
        self.h = (self.a/120)*x_0**6  + (self.b/60)*x_0**5 + (self.c/24)*x_0**4 + (self.d/6)*x_0**3 + (self.e/2)*x_0**2 + self.f*x_0 + self.g

        return sum(((self.a/840)*(x_1)**7 + (self.b/360)*x_1**6 + (self.c/120)*x_1**5 + (self.d/24)*x_1**4) - ((self.a/840)*(x_0)**7 + (self.b/360)*x_0**6 + (self.c/120)*x_0**5 + (self.d/24)*x_0**4))
    #+ (self.d/24)*x_d**4 + (self.e / 6) * x_d ** 3 + (self.f / 2) * x_d ** 2 + self.g * x_d)