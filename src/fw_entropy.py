#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:02:16 2020

@author: BrianTCook
"""

import numpy as np
import glob
import time
import matplotlib.pyplot as plt
import math

from scipy.interpolate import griddata

'''
def nonuniform_bins(list_of_values, Nmax_in_bin):
    
    data = np.sort(list_of_values)
    
    edges = [ min(data), max(data) ]
    
    flag = 0
    
    while flag == 0:
    
        edges_updated = edges
        
        #partitions old edges
        for i in range(len(edges)-1):
            
            N_in_bin = 0
            lo, high = edges[i], edges[i+1]
            
            for datum in data:
                
                if datum >= lo and datum < high:
                    
                    N_in_bin += 1
                
            if N_in_bin > Nmax_in_bin:
                
                new_edge = 0.5 * (lo + high)
                
                if new_edge not in edges_updated:
                
                    edges_updated = np.append( new_edge, edges_updated )
                
        edges_updated = np.sort(edges_updated)
            
        #checks if new edges work
        N_check = 0
    
        #checks new edges
        for i in range(len(edges_updated)-1):
            
            N_in_bin = 0
            lo, high = edges_updated[i], edges_updated[i+1]
            
            for datum in data:
                
                if datum >= lo and datum < high:
                    
                    N_in_bin += 1
                    
            if N_in_bin > Nmax_in_bin:
                
                N_check += 1

        edges = edges_updated

        if N_check == 0:
            
            flag = 1
    
    #want the center of each point, not the edges
    centers = [ 0.5*(edges[i]+edges[i+1]) for i in range(len(edges)-1)  ]
    
    return centers
'''

def lattices_maker(points, uniformity):
    
    '''
    takes 6D phase space coordinates
    tessellates phase space s.t. there are < Npoints_max
    '''
    
    cluster_populations = np.loadtxt('/home/brian/Desktop/second_project_gcs/data/Nstars_in_clusters.txt')

    for cluster in cluster_populations
    
    xvals, yvals, zvals = points[:,0], points[:,1], points[:,2]
    vxvals, vyvals, vzvals = points[:,3], points[:,4], points[:,5]
    
    if uniformity == 'uniform':
    
        Ncells = 12
        
        x_spatial = np.linspace(min(xvals), max(xvals), Ncells)
        y_spatial = np.linspace(min(yvals), max(yvals), Ncells)
        z_spatial = np.linspace(min(zvals), max(zvals), Ncells)
        x_velocity = np.linspace(min(vxvals), max(vxvals), Ncells)
        y_velocity = np.linspace(min(vyvals), max(vyvals), Ncells)
        z_velocity = np.linspace(min(vzvals), max(vzvals), Ncells)
        
    if uniformity == 'non-uniform':
    
        Nmax_in_bin = 50
        
        x_spatial = nonuniform_bins(xvals, Nmax_in_bin)
        y_spatial = nonuniform_bins(yvals, Nmax_in_bin)
        z_spatial = nonuniform_bins(zvals, Nmax_in_bin) 
        x_velocity = nonuniform_bins(vxvals, Nmax_in_bin)
        y_velocity = nonuniform_bins(vyvals, Nmax_in_bin)
        z_velocity = nonuniform_bins(vzvals, Nmax_in_bin) 
    
    #need to make non-uniform samples so that each cell has <= Nmax_in_bin sample points in it

    return x_spatial, y_spatial, z_spatial, x_velocity, y_velocity, z_velocity

def get_6D_fw(points, values, uniformity):
    
    print('unique fw values: ', np.unique(values).shape)

    X, Y, Z, VX, VY, VZ = np.meshgrid(*lattice_maker(points, uniformity))

    #Ti is 6-D interpolation using method=method    
    grid_z0 = griddata(points, values, (X, Y, Z, VX, VY, VZ), method='nearest')    

    print('unique elements in interpolated lattice: ', np.unique(grid_z0).shape)
    
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

def simpson_nonuniform(x, f):
    
    """
    FOUND ON WIKIPEDIA
    
    Simpson rule for irregularly spaced data.

        Parameters
        ----------
        x : list or np.array of floats
                Sampling points for the function values
        f : list or np.array of floats
                Function values at the sampling points

        Returns
        -------
        float : approximation for the integral
    """
    
    N = len(x) - 1
    h = np.diff(x)

    result = 0.0
    for i in range(1, N, 2):
        hph = h[i] + h[i - 1]
        result += f[i] * ( h[i]**3 + h[i - 1]**3
                           + 3. * h[i] * h[i - 1] * hph )\
                     / ( 6 * h[i] * h[i - 1] )
        result += f[i - 1] * ( 2. * h[i - 1]**3 - h[i]**3
                              + 3. * h[i] * h[i - 1]**2)\
                     / ( 6 * h[i - 1] * hph)
        result += f[i + 1] * ( 2. * h[i]**3 - h[i - 1]**3
                              + 3. * h[i - 1] * h[i]**2)\
                     / ( 6 * h[i] * hph )

    if (N + 1) % 2 == 0:
        result += f[N] * ( 2 * h[N - 1]**2
                          + 3. * h[N - 2] * h[N - 1])\
                     / ( 6 * ( h[N - 2] + h[N - 1] ) )
        result += f[N - 1] * ( h[N - 1]**2
                           + 3*h[N - 1]* h[N - 2] )\
                     / ( 6 * h[N - 2] )
        result -= f[N - 2] * h[N - 1]**3\
                     / ( 6 * h[N - 2] * ( h[N - 2] + h[N - 1] ) )
    return result
        

def get_entropy(points, values, uniformity):
    
    '''
    right now have N x N x N x N x N x N array
    need to integrate s.t. we get a 1 x 1 x 1 x 1 x 1 x 1 array
    the value within it will be the information entropy
    '''
    
    x_spatial, y_spatial, z_spatial, x_velocity, y_velocity, z_velocity = lattice_maker(points, uniformity)
    grid = get_6D_fw(points, values, uniformity)
    
    Nx, Ny, Nz = len(x_spatial), len(y_spatial), len(z_spatial)
    Nvx, Nvy, Nvz = len(x_velocity), len(y_velocity), len(z_velocity)
    
    sixD_arr = np.multiply(grid, np.log(grid))
    
    if uniformity == 'uniform':
        
        fiveD_arr = [ simpson(z_velocity, sixD_arr[i,j,k,l,m,:])
                      for i in range(Nx) for j in range(Ny) 
                      for k in range(Nz) for l in range(Nvx) 
                      for m in range(Nvy) ]
        
        fiveD_arr = np.asarray(fiveD_arr).reshape(Nx, Ny, Nz, Nvx, Nvy)
         
        fourD_arr = [ simpson(y_velocity, fiveD_arr[i,j,k,l,:])
                      for i in range(Nx) for j in range(Ny) 
                      for k in range(Nz) for l in range(Nvx)  ]
        
        fourD_arr = np.asarray(fourD_arr).reshape(Nx, Ny, Nz, Nvx)
        
        threeD_arr = [ simpson(x_velocity, fourD_arr[i,j,k,:])
                       for i in range(Nx) for j in range(Ny) 
                       for k in range(Nz) ]
        
        threeD_arr = np.asarray(threeD_arr).reshape(Nx, Ny, Nz)
        
        twoD_arr = [ simpson(z_spatial, threeD_arr[i,j,:])
                     for i in range(Nx) for j in range(Ny) ]
        
        twoD_arr = np.asarray(twoD_arr).reshape(Nx, Ny)
        
        oneD_arr = [ simpson(y_spatial, twoD_arr[i,:])
                  for i in range(Nx)  ]
        
        #returns entropy
        return simpson(x_spatial, oneD_arr)
    
    if uniformity == 'non-uniform':
    
        fiveD_arr = [ simpson_nonuniform(z_velocity, sixD_arr[i,j,k,l,m,:])
                      for i in range(Nx) for j in range(Ny) 
                      for k in range(Nz) for l in range(Nvx) 
                      for m in range(Nvy) ]
        
        fiveD_arr = np.asarray(fiveD_arr).reshape(Nx, Ny, Nz, Nvx, Nvy)
         
        fourD_arr = [ simpson_nonuniform(y_velocity, fiveD_arr[i,j,k,l,:])
                      for i in range(Nx) for j in range(Ny) 
                      for k in range(Nz) for l in range(Nvx)  ]
        
        fourD_arr = np.asarray(fourD_arr).reshape(Nx, Ny, Nz, Nvx)
        
        threeD_arr = [ simpson_nonuniform(x_velocity, fourD_arr[i,j,k,:])
                       for i in range(Nx) for j in range(Ny) 
                       for k in range(Nz) ]
        
        threeD_arr = np.asarray(threeD_arr).reshape(Nx, Ny, Nz)
        
        twoD_arr = [ simpson_nonuniform(z_spatial, threeD_arr[i,j,:])
                     for i in range(Nx) for j in range(Ny) ]
        
        twoD_arr = np.asarray(twoD_arr).reshape(Nx, Ny)
        
        oneD_arr = [ simpson_nonuniform(y_spatial, twoD_arr[i,:])
                  for i in range(Nx)  ]
        
        #returns entropy
        return simpson_nonuniform(x_spatial, oneD_arr)

if __name__ in '__main__':
    
    point_files = glob.glob('/Users/BrianTCook/Desktop/Thesis/second_project_gcs/Enbid-2.0/AMUSE_data/*.ascii')
    
    t0 = time.time()
    
    logN_max = 0
    
    xvals_all, yvals_all = [ [] for i in range(logN_max+1) ], [ [] for i in range(logN_max+1) ]
    
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    
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

            points = np.loadtxt(point_file)
            values = np.loadtxt(point_file + '_output.est')
            
            uniformity = 'uniform'
            S = get_entropy(points, values, uniformity)
            #fw_map, x0, x1, y0, y1 = get_map(points, values)
            #print('number of unique elements in fw_map: ', np.unique(fw_map).shape)
            
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
        plt.scatter(xvals, yvals, label=r'$\log_{2} N_{\mathrm{clusters}}$ = %i'%(logN), s=4)
    
    plt.gca().set_yscale('log')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
    plt.xlabel('Simulation Time (Myr)', fontsize=16)
    plt.ylabel(r'$S$ ([kpc km/s]$^{3}$)', fontsize=16)
    plt.gca().tick_params(labelsize='large')
    plt.tight_layout()
    plt.savefig('entropy_evolution.pdf')