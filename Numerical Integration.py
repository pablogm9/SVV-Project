from math import *
import numpy as np
#This should be used along a chord of a given airfoil
#a,b,c,d should be passed as arrays and n (number of points) as an int
#c_length is the chord length
class Analytical_Int_1D:
    def __init__(self,matrix,length):
        self.matrix = matrix
        self.a = [arr[0] for arr in matrix]
        self.b = [arr[1] for arr in matrix]
        self.c = [arr[2] for arr in matrix]
        self.d = [arr[3] for arr in matrix]
        self.length = length

    #x_0 and x_1 are the two x values on the chord between which you are integrating
    def integrator(self):
        n_points = len(self.matrix)
        int_dist = self.length/n_points
        x_vals = [0]
        while x_vals[-1] < self.length:
            x_vals.append(x_vals[-1]+int_dist)
        x_1 = np.array(x_vals[1:])
        x_0 = np.array(x_vals[:-1])
        print(x_vals)
        return sum((self.a*(x_1-x_0)**4)/3 + (self.b*(x_1-x_0)**3)/2 + self.c*(x_1-x_0)**2 + self.d*(x_1-x_0))

