#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:02:16 2020

@author: BrianTCook
"""

import numpy as np
import glob
import re
import time
import matplotlib.pyplot as plt

from scipy.interpolate import griddata

def get_6D_fw(points, values, N):
    
    spatial = list(np.linspace(-4., 4., N))  
    velocity = list(np.linspace(-500., 500., N))
    
    X, Y, Z, VX, VY, VZ = np.meshgrid(*(spatial, spatial, spatial, 
                                        velocity, velocity, velocity))

    #Ti is 6-D interpolation using method=method    
    grid_z0 = griddata(points, values, (X, Y, Z, VX, VY, VZ), method='nearest')
    grid_z1 = griddata(points, values, (X, Y, Z, VX, VY, VZ), method='linear')
    
    #does nearest but fills in nan values at boundary with small number
    fill_value = np.median(grid_z0)  # Whatever you like
    grid_z0[np.isnan(grid_z1)] = fill_value
    
    return grid_z0

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
    
    N = 10
    
    spatial_vals = np.linspace(-4., 4., N)
    velocity_vals = np.linspace(-500., 500., N)
    
    grid = get_6D_fw(points, values, N)
    
    sixD_arr = np.multiply(grid, np.log(grid))
    #sixD_arr = np.nan_to_num(sixD_arr) # 0 * np.log(0) is nan
    
    fiveD_arr = [ simpson(velocity_vals, sixD_arr[i,j,k,l,m,:])
                  for i in range(N) for j in range(N) 
                  for k in range(N) for l in range(N) 
                  for m in range(N) ]
    
    fiveD_arr = np.asarray(fiveD_arr).reshape(N, N, N, N, N)
    
    fourD_arr = [ simpson(velocity_vals, fiveD_arr[i,j,k,l,:])
                  for i in range(N) for j in range(N) 
                  for k in range(N) for l in range(N) ]
    
    fourD_arr = np.asarray(fourD_arr).reshape(N, N, N, N)
    
    threeD_arr = [ simpson(velocity_vals, fourD_arr[i,j,k,:])
                   for i in range(N) for j in range(N) 
                   for k in range(N) ]
    
    threeD_arr = np.asarray(threeD_arr).reshape(N, N, N)
    
    twoD_arr = [ simpson(spatial_vals, threeD_arr[i,j,:])
                 for i in range(N) for j in range(N)  ]
    
    twoD_arr = np.asarray(twoD_arr).reshape(N, N)
    
    oneD_arr = [ simpson(spatial_vals, twoD_arr[i,:])
              for i in range(N)  ]
    
    return simpson(spatial_vals, oneD_arr)

if __name__ in '__main__':
    
    point_files = glob.glob('/Users/BrianTCook/Desktop/Thesis/second_project_gcs/data/enbid_files/*.ascii')
    
    t0 = time.time()
    
    logN_max = 2
    
    xvals_all, yvals_all = [ [] for i in range(logN_max+1) ], [ [] for i in range(logN_max+1) ]
    
    for point_file in point_files:
        
        print('')
        
        ahh = point_file.split('_')
        bhh = [ a.split('.') for a in ahh ]
        chh = [j for i in bhh for j in i]
        
        print(chh[5:11])
        
        #first one is frame, second one is 
        frame, Nclusters  = [ int(s) for s in chh if s.isdigit() ]
        
        print(Nclusters, 2**logN_max)
        
        if Nclusters <= 2**logN_max:
        
            points = np.loadtxt(point_file)[:, 1:]
            values = np.asarray([ 10**(np.random.rand()+3.) for i in range(points.shape[0]) ])
            
            S = get_entropy(points, values)
            
            sim_time = frame/400. * 40. #frame/frame * end_time (Myr)
            
            tf = time.time()
            
            logN = int(np.log2(Nclusters)) #number of clusters
            
            #plotting entropy as a function of time, saving by log2N
            xvals_all[logN].append(sim_time)
            yvals_all[logN].append(S)
        
            tf = time.time()
            print('current time: %.04f minutes'%((tf-t0)/60.))
        
    plt.figure() 
    
    for logN in range(logN_max+1):
        
        xvals, yvals = xvals_all[logN], yvals_all[logN]
        plt.scatter(xvals, yvals, label=r'$\log_{2} N_{\mathrm{clusters}}$ = %i'%(logN))
    
    plt.gca().set_yscale('log')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('Simulation Time (Myr)', fontsize=16)
    plt.ylabel(r'$S$ [kpc km/s]$^{3}$', fontsize=16)
    plt.tight_layout()
    plt.savefig('entropy_evolution.pdf')