"""
 Reflection/Transmission coefficent calculator
 author kaspar Snashall
 TCD
 This program will calculate the reflection and transmission coeffcient for a semi circular potential
"""
import numpy as np
import matplotlib.pyplot as plt
from cmath import *
import sys
from numpy import linalg
#begin definitions and glbal variables

def k(E,V):
    """
    This function will make a list if k values for use in the schrodinger equations
    """
    # ignoring all constants such as h-bar and m, may add in later...
    return sqrt(2*(E-V))     
#print k(0.02,0.06)
def semicircle(R,xval):
    """
    makes the height values for a semi circle
    """
    if xval > R:
        print "Your x value is outside the range of -R to R"
        sys.exit(2)
    return sqrt(R**2-xval**2)
   
def klist(x,R,E):
    """
    this function will make a list of (x) k values for a given R,E
    """
    xvals = np.linspace(-R,R,100) # make a list of x values
    delx = 2*float(R)/x    
    #print xvals
    heights = [] # list of heights
    for value in xvals:
        h1 = semicircle(R,value) # find all the heights of the various walls
        heights.append(h1)
    kvals = [] # list of k vaues
    for height in heights:
        k1 = k(E,height)
        kvals.append(k1) #creates all the kvals
    print "k values", '\n', kvals
    return kvals,delx

def tmatrix(kvals,delx):
    """
    makes a list of T matrices
    """
    tmatrices = []
    for k in kvals:
        t1 = cos(k*delx)
        t2 = complex(0,1)*sin(k*delx)/k
        t3 = complex(0,1)*sin(k*delx)
        t4 = cos(k*delx)
        matrix = np.matrix([[t1, t2], [t3, t4]])
        tmatrices.append(matrix)
    return tmatrices
        
def m1(R,E):
    """
    creates the first k matrix
    """
    k1 = k(E,0)
    ma = e*complex(0,k1*(-R))
    mb =  e*complex(0,-k1*(-R))
    mc =  k1*e*complex(0,k1*(-R))
    md = k1*e*complex(0,-k1*(-R))
    
    return np.matrix([[ma,mb],[mc,md]])
    
def m3(R,E):
    """
    creates the secondary m3 matirx inverse
    """
    k3 = k(E,0)
    ma = e*complex(0,k3*(R))
    mb =  e*complex(0,-k3*(R))
    mc =  k3*e*complex(0,k3*(R))
    md = k3*e*complex(0,-k3*(R))
    A = np.matrix([[-md,mc],[mb,-ma]])
    detA = 1/(ma*md-mc*mc)
    return detA*A
#print tmatrix([1],0.1)
#print m1(1,0.5)
#print m3(1,0.5)
def main():
    R = 1
    x = 10 # number of sections
    E = 0.99
    kvals,delx = klist(x,R,E)
    tmatrices = tmatrix(kvals,delx)
    #print tmatrices
    totalT = np.matrix([[1,1],[1,1]])
    for matrix in tmatrices:
        totalT= totalT*matrix
    A = m3(R,E)
    B = m1(R,E)
    X = A*totalT*B
    print X
    A1 = np.matrix([[complex(1,1)],[complex(1,1)]])
    coeff = X*A1
    print coeff    
if __name__ == '__main__':
    main()
        
