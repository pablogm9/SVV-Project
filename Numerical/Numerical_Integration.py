from math import *
import numpy as np
#This should be used along a chord of a given airfoil
#a,b,c,d should be passed as arrays and n (number of points) as an int
#c_length is the chord length
'''
This class is to be used to analytically integrate between points
INPUTS: matrix: a matrix of a,b,c,d values for the cubic splines between each of the pairs of the coords also given
        coords: an array of the coordinates at which the aerodynamic data was taken for the a,b,c,d cubic splines(array)
OUTPUTS: The sum of the integrated portions for the entire length of the span or the chord over which you have integrated(array)
        
'''

class Analytical_Int_1D:
    def __init__(self,matrix,coords):
        self.matrix = matrix
        if len(matrix[0]) == 4:
            self.a = [arr[0] for arr in matrix]
            self.b = [arr[1] for arr in matrix]
            self.c = [arr[2] for arr in matrix]
            self.d = [arr[3] for arr in matrix]
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
            return sum((self.a*(x_1-x_0)**4)/3 + (self.b*(x_1-x_0)**3)/2 + self.c*(x_1-x_0)**2 + self.d*(x_1-x_0))
        if len(self.matrix[0]) == 5:
            return sum((self.a*(x_1-x_0)**5)/4 + (self.b*(x_1-x_0)**4)/3 + (self.c*(x_1-x_0)**3)/2 + self.d*(x_1-x_0)**2 + self.e(x_1-x_0))

