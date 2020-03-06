#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:02:16 2020

@author: BrianTCook
"""

import numpy as np
from scipy.interpolate import griddata

def get_6D_fw(points, values):

    points = np.load()
    values = np.loadtxt()
    
    x = np.linspace(-4., 4., 100)
    y = np.linspace(-4., 4., 100)
    z = np.linspace(-4., 4., 100)
    
    vx = np.linspace(-500., 500., 100)
    vy = np.linspace(-500., 500., 100)
    vz = np.linspace(-500., 500., 100)
    
    X, Y, Z, VX, VY, VZ = np.meshgrid(x,y, z, vx, vy, vz)

    #Ti is 6-D interpolation using method=method
    Ti = griddata(points, values, (X, Y, Z, VX, VY, VZ), method='linear')
        
    return Ti

def simpson(xvals, fvals, x_init, x_final, N_simp): #Simpson's integration rule, \int_{a}^{b} f(x) dx with N sample points

    N_simp = len(fvals)

    a, b = 0, len(fvals)-1

    I = f[a]+f[b]

    odds = [ 4*f[a + k] for k in range(1,N_simp,2) ]
    evens = [ 2*f[a + k] for k in range(2,N_simp,2) ]
    I += sum(odds) + sum(evens)

    I *= h/3.
    
    return I

def get_entropy(grid):
    
    '''
    right now have N x N x N x N x N x N array
    need to integrate s.t. we get a 1 x 1 x 1 x 1 x 1 x 1 array
    the value within it will be the information entropy
    '''
    
    return 0
