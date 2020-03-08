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

def simpson(xvals, fvals): #Simpson's integration rule, \int_{a}^{b} f(x) dx with N sample points

    N_simp = len(fvals)
    h = xvals[1]-xvals[0]

    a, b = 0, len(fvals)-1

    I = fvals[a]+fvals[b]

    odds = [ 4*fvals[a + k] for k in range(1,N_simp,2) ]
    evens = [ 2*fvals[a + k] for k in range(2,N_simp,2) ]
    I += sum(odds) + sum(evens)

    I *= h/3.
    
    return I

def get_entropy(points, values):
    
    '''
    right now have N x N x N x N x N x N array
    need to integrate s.t. we get a 1 x 1 x 1 x 1 x 1 x 1 array
    the value within it will be the information entropy
    '''
    
    N = 100
    
    spatial_vals = np.linspace(-2, 2, N)
    velocity_vals = np.linspace(-500, 500, N)
    
    grid = get_6D_fw(points, values)
    
    sixD_arr = np.multiply(grid, np.log(grid))
    
    fiveD_arr = [ simpson(velocity_vals, sixD_arr[i,j,k,l,m,:])
                  for i in range(N) for j in range(N) 
                  for k in range(N) for l in range(N) 
                  for m in range(N) ]
    
    fiveD_arr = fiveD_arr.reshape(N, N, N, N, N)
    
    fourD_arr = [ simpson(velocity_vals, fiveD_arr[i,j,k,l,:])
                  for i in range(N) for j in range(N) 
                  for k in range(N) for l in range(N) ]
    
    fourD_arr = fourD_arr.reshape(N, N, N, N)
    
    threeD_arr = [ simpson(velocity_vals, fourD_arr[i,j,k,:])
                   for i in range(N) for j in range(N) 
                   for k in range(N) ]
    
    threeD_arr = threeD_arr.reshape(N, N, N)
    
    twoD_arr = [ simpson(spatial_vals, threeD_arr[i,j,:])
                 for i in range(N) for j in range(N)  ]
    
    twoD_arr = fiveD_arr.reshape(N, N)
    
    oneD_arr = [ simpson(spatial_vals, twoD_arr[i,:])
              for i in range(N)  ]
    
    return simpson(spatial_vals, oneD_arr)
