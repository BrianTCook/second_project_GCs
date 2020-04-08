#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 11:19:05 2020

@author: BrianTCook
"""

import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def dimensions(code_name, sim_times_unitless, Norbiters):
    
    '''
    creates an array of dimensions [D0(t=0), ..., D0(t=f),
                                    .
                                    .
                                    .
                                    DN(t=0), ..., DN(t=tf)] 
    for the N clusters provided as input
    '''
    
    datadir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
    clusters_init = np.loadtxt(datadir + 'enbid_files/for_enbid_%s_SingleCluster_frame_00000_Norbiters_%s.ascii'%(code_name, str(Norbiters)))
    
    cluster_populations = np.loadtxt(datadir + 'Nstars_in_clusters.txt')
    cluster_populations = list(cluster_populations[:Norbiters])
    
    def domain(m, cluster):
        
        '''
        m is a string denoting which column to take from
        cluster is a pandas dataframe with a cluster's 6D phase space info
        '''
        
        df = pd.DataFrame(data=cluster, columns=['x', 'y', 'z', 'vx', 'vy', 'vz'])
        
        #use 2 standard deviations as bound, avoids outliers
        return np.percentile(df[m], 97.8) - np.percentile(df[m], 2.2)
        
    switch = 5. #decides if direction is on/off
    gadget_flag = int(math.floor(len(sim_times_unitless)/10))
    dimension_array = np.zeros((Norbiters, 10+1)) #see gadget_flag
    
    for k, number_of_stars in enumerate(cluster_populations):
        
        starting_index = int(np.sum( cluster_populations[:k] ))
        ending_index = starting_index + int(number_of_stars)

        time_counter = 0
                
        for j, t in enumerate(sim_times_unitless):
            
            if j%gadget_flag == 0:
                
                clusters_t = np.loadtxt(datadir+'enbid_files/for_enbid_%s_SingleCluster_frame_%s_Norbiters_%s.ascii'%(code_name, str(j).rjust(5, '0'), str(Norbiters)))
            
                cluster_init = clusters_init[starting_index:ending_index, :]                
                cluster_t = clusters_t[starting_index:ending_index, :]
                
                #determine initial size of each cluster in phase space
                
                D = 0
                if domain('x', cluster_t) >= switch * domain('x', cluster_init):
                    D += 1
                if domain('y', cluster_t) >= switch * domain('y', cluster_init):
                    D += 1
                if domain('z', cluster_t) >= switch * domain('z', cluster_init):
                    D += 1
                if domain('vx', cluster_t) >= switch * domain('vx', cluster_init):
                    D += 1
                if domain('vy', cluster_t) >= switch * domain('vy', cluster_init):
                    D += 1
                if domain('vz', cluster_t) >= switch * domain('vz', cluster_init):
                    D += 1
                    
                dimension_array[k, time_counter] = D
                time_counter += 1
        
    return dimension_array

if __name__ in '__main__':
    
    code_name = 'tree'
    tend, dt = 40., 0.1
    
    sim_times_unitless = np.arange(0., tend+dt, dt)
    
    logN_max = 7
    Norbiters_list = [ 2**i for i in range(logN_max) ]
    
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    
    for Norbiters in Norbiters_list:
    
        plt.figure()
        dim_array = dimensions(code_name, sim_times_unitless, Norbiters)
        im = plt.imshow(dim_array, origin='lower', aspect='equal',
                        extent=[min(sim_times_unitless), max(sim_times_unitless), 0, Norbiters],
                        norm=Normalize(vmin=0, vmax=6))
        cbar = plt.colorbar(im)
        cbar.set_label(r'$D$', rotation=270, labelpad=12)
        plt.ylabel('Cluster ID')#, fontsize=16)
        plt.xlabel('Simulation Time (Myr)')#, fontsize=16)
        plt.title('Estimated Manifold Dimensions')#, fontsize=24)
        plt.gca().set_yticks(np.arange(0, Norbiters, math.ceil(Norbiters/4)));
        plt.savefig('dimensions_Norbiters=%i.pdf'%Norbiters)