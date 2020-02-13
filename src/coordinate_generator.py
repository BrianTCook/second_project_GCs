#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 16:20:46 2020

@author: BrianTCook
"""

import numpy as np
from scipy.interpolate import CubicSpline

def azimuthal_coords():

    '''
    creating lists of azimuthal coordinates, pulling from U(0, 2*pi)
    '''

    
    try:
        
        data_directory = '~/Desktop/second_project_GCs/data/'
        simulation_OC_phivals = np.loadtxt(data_directory+'simulation_OC_phivals.txt')
        simulation_YMC_phivals = np.loadtxt(data_directory+'simulation_YMC_phivals.txt')
        simulation_GC_phivals = np.loadtxt(data_directory+'simulation_GC_phivals.txt')
        
    except:
        
        simulation_OC_phivals = 2*np.pi * np.random.rand(1537)
        simulation_YMC_phivals = 2*np.pi * np.random.rand(1000)
        simulation_GC_phivals = 2*np.pi * np.random.rand(1000)
        
        np.savetxt('simulation_OC_phivals.txt', simulation_OC_phivals)
        np.savetxt('simulation_YMC_phivals.txt', simulation_OC_phivals)
        np.savetxt('simulation_GsC_phivals.txt', simulation_OC_phivals)
    
    return simulation_OC_phivals, simulation_YMC_phivals, simulation_GC_phivals
    
'''
creating lists of rvals and zvals from GMM of MW distributions
'''

def radial_coords():
    
    try:
        
        data_directory = '~/Desktop/second_project_GCs/data/'
        
        MW_OC_rvals = np.loadtxt(data_directory+'MW_OC_rvals.txt')
        MW_YMC_rvals = np.loadtxt(data_directory+'MW_YMC_rvals.txt')
        MW_GC_rvals = np.loadtxt(data_directory+'MW_GC_rvals.txt')
        
        simulation_OC_rvals = np.loadtxt(data_directory+'simulation_OC_rvals.txt')
        simulation_YMC_rvals = np.loadtxt(data_directory+'simulation_YMC_rvals.txt')
        simulation_GC_rvals = np.loadtxt(data_directory+'simulation_GC_rvals.txt')
        
    except:
        
        #will need to update these from a catalog at some point
        MW_OC_rvals = np.abs( np.random.normal(loc=0., scale=4., size=1537) )
        MW_YMC_rvals = np.abs( np.random.normal(loc=0., scale=4., size=12) )
        MW_GC_rvals = np.abs( np.random.normal(loc=0., scale=4., size=157) )
        
        np.savetxt('MW_OC_rvals.txt', MW_OC_rvals)
        np.savetxt('MW_YMC_rvals.txt', MW_OC_rvals)
        np.savetxt('MW_GsC_rvals.txt', MW_OC_rvals)
        
        #generating histograms of r values for each cluster type
        y_OC, x_OC = np.histogram(MW_OC_rvals)
        x_OC = [ 0.5*(x_OC[i]+x_OC[i+1]) for i in range(len(x_OC)-1) ]
        
        y_YMC, x_YMC = np.histogram(MW_OC_rvals)
        x_YMC = [ 0.5*(x_YMC[i]+x_YMC[i+1]) for i in range(len(x_YMC)-1) ]
        
        y_GC, x_GC = np.histogram(MW_OC_rvals)
        x_GC = [ 0.5*(x_GC[i]+x_GC[i+1]) for i in range(len(x_GC)-1) ]
        
        #fit a cubic spline to each
        OC_spline, YMC_spline, GC_spline = CubicSpline(x_OC, y_OC), CubicSpline(x_YMC, y_YMC), CubicSpline(x_GC, y_GC)
        
        #set up domains for each rvalue to be tried out
        min_x_OC, max_x_OC = min(x_OC), max(x_OC)
        min_x_YMC, max_x_YMC = min(x_YMC), max(x_YMC)
        min_x_GC, max_x_GC = min(x_GC), max(x_GC)
        
        #fill them in using random sampling
        simulation_OC_rvals, simulation_YMC_rvals, simulation_GC_rvals = [], [], []
        
        while len(simulation_OC_zvals) < 1537:
            
            x_proposed = (max_x_OC - min_x_OC)*np.random.random() + min_x_OC
            y_proposed = np.random.random()
            
            if y_proposed <= OC_spline(x_proposed):
                
                simulation_OC_rvals.append(x_proposed)
            
        while len(simulation_YMC_rvals) < 1000:
            
            x_proposed = (max_x_YMC - min_x_YMC)*np.random.random() + min_x_YMC
            y_proposed = np.random.random()
            
            if y_proposed <= YMC_spline(x_proposed):
                
                simulation_YMC_rvals.append(x_proposed)
            
        while len(simulation_GC_rvals) < 1000:
            
            x_proposed = (max_x_GC - min_x_GC)*np.random.random() + min_x_GC
            y_proposed = np.random.random()
            
            if y_proposed <= GC_spline(x_proposed):
                
                simulation_GC_zvals.append(x_proposed)
        
        np.savetxt('simulation_OC_rvals.txt', simulation_OC_rvals)
        np.savetxt('simulation_YMC_rvals.txt', simulation_OC_rvals)
        np.savetxt('simulation_GsC_rvals.txt', simulation_OC_rvals)
        
    return simulation_OC_rvals, simulation_YMC_rvals, simulation_GC_rvals
        
def zed_coords():
    
    try:
        
        data_directory = '~/Desktop/second_project_GCs/data/'
        
        MW_OC_zvals = np.loadtxt(data_directory+'MW_OC_zvals.txt')
        MW_YMC_zvals = np.loadtxt(data_directory+'MW_YMC_zvals.txt')
        MW_GC_zvals = np.loadtxt(data_directory+'MW_GC_zvals.txt')
    
        simulation_OC_zvals = np.loadtxt(data_directory+'simulation_OC_zvals.txt')
        simulation_YMC_zvals = np.loadtxt(data_directory+'simulation_YMC_zvals.txt')
        simulation_GC_zvals = np.loadtxt(data_directory+'simulation_GC_zvals.txt')
        
    except:
        
        #will need to update these from a catalog at some point
        MW_OC_zvals = np.random.normal(loc=0., scale=1., size=1537)
        MW_YMC_zvals = np.random.normal(loc=0., scale=1., size=12)
        MW_GC_zvals = np.random.normal(loc=0., scale=1., size=157)
        
        np.savetxt('MW_OC_zvals.txt', MW_OC_zvals)
        np.savetxt('MW_YMC_zvals.txt', MW_OC_zvals)
        np.savetxt('MW_GsC_zvals.txt', MW_OC_zvals)
        
        #generating normalized histograms of r values for each cluster type
        y_OC, x_OC = np.histogram(MW_OC_zvals, density=True)
        x_OC = [ 0.5*(x_OC[i]+x_OC[i+1]) for i in range(len(x_OC)-1) ]
        
        y_YMC, x_YMC = np.histogram(MW_OC_zvals, density=True)
        x_YMC = [ 0.5*(x_YMC[i]+x_YMC[i+1]) for i in range(len(x_YMC)-1) ]
        
        y_OC, x_OC = np.histogram(MW_OC_zvals, density=True)
        x_GC = [ 0.5*(x_GC[i]+x_GC[i+1]) for i in range(len(x_GC)-1) ]
        
        #fit a cubic spline to each
        OC_spline, YMC_spline, GC_spline = CubicSpline(x_OC, y_OC), CubicSpline(x_YMC, y_YMC), CubicSpline(x_GC, y_GC)
        
        #set up domains for z values to be set up
        min_x_OC, max_x_OC = min(x_OC), max(x_OC)
        min_x_YMC, max_x_YMC = min(x_YMC), max(x_YMC)
        min_x_GC, max_x_GC = min(x_GC), max(x_GC)
        
        #fill them in using random sampling
        simulation_OC_zvals, simulation_YMC_zvals, simulation_GC_zvals = [], [], []
        
        while len(simulation_OC_zvals) < 1537:
            
            x_proposed = (max_x_OC - min_x_OC)*np.random.random() + min_x_OC
            y_proposed = np.random.random()
            
            if y_proposed <= OC_spline(x_proposed):
                
                simulation_OC_zvals.append(x_proposed)
            
        while len(simulation_YMC_zvals) < 1000:
            
            x_proposed = (max_x_YMC - min_x_YMC)*np.random.random() + min_x_YMC
            y_proposed = np.random.random()
            
            if y_proposed <= YMC_spline(x_proposed):
                
                simulation_YMC_zvals.append(x_proposed)
            
        while len(simulation_GC_zvals) < 1000:
            
            x_proposed = (max_x_GC - min_x_GC)*np.random.random() + min_x_GC
            y_proposed = np.random.random()
            
            if y_proposed <= GC_spline(x_proposed):
                
                simulation_GC_zvals.append(x_proposed)
        
        np.savetxt('simulation_OC_zvals.txt', simulation_OC_zvals)
        np.savetxt('simulation_YMC_zvals.txt', simulation_OC_zvals)
        np.savetxt('simulation_GsC_zvals.txt', simulation_OC_zvals)
      
    return simulation_OC_zvals, simulation_YMC_zvals, simulation_GC_zvals

if __name__ in '__main__':
    
    OC_phi, YMC_phi, GC_phi = azimuthal_coords()
    OC_r, YMC_r, GC_r = radial_coords()
    OC_z, YMC_z, GC_z = zed_coords()
    
    print('hello world!')