#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 17:02:05 2020

takes globular cluster coordinates given in Harris (1996) and converts from Cartesian to cylindrical coordinates

@author: BrianTCook
"""

import numpy as np

def harris_coord_converter():
    
    '''
    harris_array loads in a .txt file
    '''
    
    data_directory = '/home/brian/Desktop/second_project_gcs/data/'
    harris_array = np.loadtxt(data_directory+'mwgc.txt')
    
    print(harris_array)
    
    '''
    xvals, yvals, zvals = harris_array[:,0], harris_array[:,1], harris_array[:,2]
    
    rvals = [ np.sqrt(x[i]**2 + y[i]**2) for i in range(len(harris_df.index)) ]
    phivals = [ np.arctan(yvals[i]/xvals[i]) for i in range(len(harris_df.index)) ]
    zvals = zvals
    
    np.savetxt('MW_GC_rvals.txt', rvals)
    np.savetxt('MW_GC_phivals.txt', phivals)
    np.savetxt('MW_GC_zvals.txt', zvals)
    '''
    
    return 1

