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
from sklearn.decomposition import PCA
from cluster_table import sort_clusters_by_attribute

def dimensions(code_name, finder_name, sim_times_unitless, Norbiters, cluster_populations):
    
    '''
    creates an array of dimensions [D0(t=0), ..., D0(t=f),
                                    .
                                    .
                                    .
                                    DN(t=0), ..., DN(t=tf)] 
    for the N clusters provided as input
    '''
    
    datadir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
    filename_init = 'tree_data/tree_seba_100Myr/enbid_%s_frame_00000_Norbiters_%s.ascii'%(code_name, str(Norbiters))
    clusters_init = np.loadtxt(datadir+filename_init)
    
    gadget_flag = int(math.floor(len(sim_times_unitless)/50))
    dimension_array = np.zeros((Norbiters, 50+1)) #see gadget_flag
    
    if finder_name == 'naive':
        
        def domain(m, cluster):
            
            '''
            m is a string denoting which column to take from
            cluster is a pandas dataframe with a cluster's 6D phase space info
            '''
            
            df = pd.DataFrame(data=cluster, columns=['x', 'y', 'z', 'vx', 'vy', 'vz'])
            
            #use 2 standard deviations as bound, avoids outliers
            return np.percentile(df[m], 97.8) - np.percentile(df[m], 2.2)
            
        switch = 5. #decides if direction is on/off
    
        time_counter = 0
                    
        for j, t in enumerate(sim_times_unitless):
            
            if j%gadget_flag == 0:
                
                filename_t = 'tree_data/tree_seba_100Myr/enbid_%s_frame_%s_Norbiters_%s.ascii'%(code_name, str(j).rjust(5, '0'), str(Norbiters))
                clusters_t = np.loadtxt(datadir+filename_t)
                    
                for k, number_of_stars in enumerate(cluster_populations):
            
                    starting_index = int(np.sum( cluster_populations[:k] ))
                    ending_index = starting_index + int(number_of_stars)
                
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
    
    if finder_name == 'pca':
        
        time_counter = 0
                    
        for j, t in enumerate(sim_times_unitless):
            
            if j%gadget_flag == 0:
                
                filename_t = 'tree_data/tree_seba_100Myr/enbid_%s_frame_%s_Norbiters_%s.ascii'%(code_name, str(j).rjust(5, '0'), str(Norbiters))
                clusters_t = np.loadtxt(datadir+filename_t)
                
                for k, number_of_stars in enumerate(cluster_populations):
            
                    starting_index = int(np.sum( cluster_populations[:k] ))
                    ending_index = starting_index + int(number_of_stars)
                                
                    cluster_t = clusters_t[starting_index:ending_index, :]
                    
                    pca = PCA()
                    pca.fit(cluster_t)
                    
                    variance_ratios = pca.explained_variance_ratio_
                    
                    D = 0
                    
                    for var_rat in variance_ratios:
                        
                        if var_rat >= 0.01:
                            
                            D += 1
                
                    dimension_array[k, time_counter] = D
                    
                time_counter += 1
                    
    return dimension_array

if __name__ in '__main__':
    
    code_name = 'tree'
    tend, dt = 100., 0.1
    
    sim_times_unitless = np.arange(0., tend+dt, dt)
    
    logN_max = 4
    Norbiters_list = [ 2**i for i in range(logN_max+1) ]
    
    datadir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
    
    cluster_populations_raw = np.loadtxt(datadir+'Nstars_in_clusters.txt')
    indices_dict = sort_clusters_by_attribute('|r|')
    cluster_populations_sorted = [ cluster_populations_raw[indices_dict[i]] for i in range(2**logN_max) ]

    rs = np.loadtxt(datadir+'ICs/dehnen_rvals.txt')
    phis = np.loadtxt(datadir+'ICs/dehnen_phivals.txt')
    zs = np.loadtxt(datadir+'ICs/dehnen_zvals.txt')
    
    vrs = np.loadtxt(datadir+'ICs/bovy_vrvals.txt')
    vphis = np.loadtxt(datadir+'ICs/bovy_vphivals.txt')
    vzs = np.loadtxt(datadir+'ICs/bovy_vzvals.txt')
    
    N = 2**logN_max #total number of initialized clusters I have
    
    #convert from galpy/cylindrical to AMUSE/Cartesian units
    #all in kpc
    xs = [ rs[i] * np.cos(phis[i]) for i in range(N) ]
    ys = [ rs[i] * np.sin(phis[i]) for i in range(N) ]
    zs = [ zs[i] for i in range(N) ]
    
    #all in km/s
    vxs = [ vrs[i] * np.cos(phis[i]) - vphis[i] * np.sin(phis[i]) for i in range(N) ] 
    vys = [ vrs[i] * np.sin(phis[i]) + vphis[i] * np.cos(phis[i]) for i in range(N) ]
    vzs = [ vzs[i] for i in range(N) ]
    
    dists = [ np.sqrt(xs[i]**2 + ys[i]**2 + zs[i]**2) for i in range(N) ]
    dists_sorted = np.sort(dists)

    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    
    dim_finder_names = [ 'pca', 'naive' ] #, 'autoencoding' ]
        
    for finder_name in dim_finder_names:    
        for Norbiters in Norbiters_list:
            
            print(finder_name, Norbiters)
            
            cluster_populations = list(cluster_populations_sorted)[:Norbiters]
            
            plt.figure()
            
            dim_array = dimensions(code_name, finder_name, sim_times_unitless, Norbiters, cluster_populations)
            
            dists_relevant = dists_sorted[:Norbiters]
            dists_min, dists_max = min(dists_relevant), max(dists_relevant)
            dmin, dmax = round(dists_min, 2), round(dists_max, 2)
            
            im = plt.imshow(dim_array, origin='lower', aspect='auto',
                            extent=[min(sim_times_unitless), max(sim_times_unitless), 0, Norbiters],
                            norm=Normalize(vmin=0, vmax=6))
            
            plt.gca().set_yticklabels([dmin, (dmax-dmin)/2., dmax])
            
            cbar = plt.colorbar(im)
            cbar.set_label(r'$D$', rotation=270, labelpad=12)
            plt.ylabel(r'$|\mathbf{r}_{0}|$ (kpc)')#, fontsize=16)
            plt.xlabel('Simulation Time (Myr)')#, fontsize=16)
            plt.title('Estimated Manifold Dimensions (%s)'%(finder_name))#, fontsize=24)
            plt.gca().set_yticks(np.arange(0, Norbiters, math.ceil(Norbiters/4)));
            plt.savefig('dimensions_%s_Norbiters=%i.jpg'%(finder_name, Norbiters))
            plt.close()

    plt.figure()

    for finder_name in dim_finder_names:
        
        if finder_name == 'naive':
        
            finder_color = 'r'
            
        if finder_name == 'pca':
        
            finder_color = 'b'
        
        dim_array = dimensions(code_name, finder_name, sim_times_unitless, max(Norbiters_list), cluster_populations_sorted)
        times_plot = np.linspace(0., tend+dt, len(dim_array[0, :]))
        
        lo, mid, hi = 25, 50, 75
        
        lo_line = [ np.percentile(dim_array[:,i], lo) for i in range(len(times_plot)) ]
        mid_line = [ np.percentile(dim_array[:,i], mid) for i in range(len(times_plot)) ]
        hi_line = [ np.percentile(dim_array[:,i], hi) for i in range(len(times_plot)) ]
        
        plt.plot(times_plot, lo_line, c=finder_color, linewidth=1, linestyle='--')#, label='%i th percentile, %s'%(lo, finder_name))
        plt.plot(times_plot, mid_line, c=finder_color, linewidth=2, label='%s, %i/%i/%i th percentiles'%(finder_name, lo, mid, hi))
        plt.plot(times_plot, hi_line, c=finder_color, linewidth=1, linestyle='--')#, label='%i th percentile, %s'%(hi, finder_name))

    plt.xlabel('Simulation Time (Myr)', fontsize=16)
    plt.ylabel(r'$D$', fontsize=24)
    plt.title(r'Dimension Computation Comparison, $N_{\mathrm{clusters}} = %i$'%(max(Norbiters_list)), fontsize=14)
    plt.legend(loc='best', fontsize=10)
    plt.gca().tick_params(labelsize='large')
    plt.tight_layout()
    plt.savefig('dimension_percentiles_Norbiters=%i.pdf'%(max(Norbiters_list)))
    plt.close()
    