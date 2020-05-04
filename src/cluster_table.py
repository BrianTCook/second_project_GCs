#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 15:14:20 2020

@author: BrianTCook

USE ON MAC NOT ON VIRTUAL MACHINE
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

def sort_clusters_by_attribute(attribute):

    '''
    takes pandas df column name as argument (string)
    outputs a dictionary such that cluster 0 has lowest value, cluster 1 second lowest, and so on
    '''
    
    #data_directory = '/home/s1780638/second_project_gcs/data/'
    #data_directory = '/home/brian/Desktop/second_project_gcs/data/'
    data_directory = '/Users/BrianTCook/Desktop/Thesis/second_project_gcs/data/'
    
    rs = np.loadtxt(data_directory+'ICs/dehnen_rvals.txt')
    phis = np.loadtxt(data_directory+'ICs/dehnen_phivals.txt')
    zs = np.loadtxt(data_directory+'ICs/dehnen_zvals.txt')
    
    vrs = np.loadtxt(data_directory+'ICs/bovy_vrvals.txt')
    vphis = np.loadtxt(data_directory+'ICs/bovy_vphivals.txt')
    vzs = np.loadtxt(data_directory+'ICs/bovy_vzvals.txt')
    
    N = 128 #total number of initialized clusters I have
    
    #convert from galpy/cylindrical to AMUSE/Cartesian units
    #all in kpc
    xs = [ rs[i] * np.cos(phis[i]) for i in range(N) ]
    ys = [ rs[i] * np.sin(phis[i]) for i in range(N) ]
    zs = [ zs[i] for i in range(N) ]
    
    #all in km/s
    vxs = [ vrs[i] * np.cos(phis[i]) - vphis[i] * np.sin(phis[i]) for i in range(N) ] 
    vys = [ vrs[i] * np.sin(phis[i]) + vphis[i] * np.cos(phis[i]) for i in range(N) ]
    vzs = [ vzs[i] for i in range(N) ]
    
    dists = [ round(np.sqrt(xs[i]**2 + ys[i]**2 + zs[i]**2), 3) for i in range(N) ]
    speeds = [ round(np.sqrt(vxs[i]**2 + vys[i]**2 + vzs[i]**2), 2) for i in range(N) ]
    
    Nstars = [ len(np.loadtxt(data_directory+'/star_masses/star_masses_index=%i.txt'%i)) for i in range(N) ]
    masses = [ round(np.sum(np.loadtxt(data_directory+'/star_masses/star_masses_index=%i.txt'%i)), 2) for i in range(N) ]
    radii = np.loadtxt(data_directory+'/ICs/cluster_radii_for_sampling.txt') 
    radii = [ round(r, 2) for r in radii ]
    
    df = pd.DataFrame(list(zip(masses, Nstars, dists, speeds, radii)), columns=['M', 'Nstars', '|r|', '|v|', 'rvir'])
    
    df_sorted_by_r = df.sort_values(by=[attribute])
    
    indices_dict = {}
    
    for df, df_sorted in zip(df.index, df_sorted_by_r.index):
        indices_dict.update( {df : df_sorted} )
        
    return indices_dict

if __name__ in '__main__':
    
    sort_clusters_by_attribute('|r|')

'''
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

fig, axs = plt.subplots(2,2)

mbins = np.logspace(np.log10(min(masses)), np.log10(max(masses)), 40)
nbins = np.logspace(np.log10(min(Nstars)), np.log10(max(Nstars)), 40)

axs[0,0].hist(masses, histtype='step', bins=mbins, normed=True, cumulative=True, color='k')
axs[0,1].hist(Nstars, histtype='step', bins=nbins, normed=True, cumulative=True, color='k')
axs[1,0].hist(dists, histtype='step', bins=40, normed=True, cumulative=True, color='k')
axs[1,1].hist(speeds, histtype='step', bins=40, normed=True, cumulative=True, color='k')

axs[0,0].set_xscale('log')
axs[0,1].set_xscale('log')

axs[0,0].set_xlabel(r'$M_{\mathrm{cluster}}$ ($M_{\odot}$)', fontsize=14)
axs[0,1].set_xlabel(r'$N_{\star}$', fontsize=14)
axs[1,0].set_xlabel(r'$|\mathbf{r} \, (t=0)|$ (kpc)', fontsize=14)
axs[1,1].set_xlabel(r'$|\mathbf{v} \, (t=0)|$ (km/s)', fontsize=14)

axs[0,0].tick_params(labelsize='large')
axs[0,1].tick_params(labelsize='large')
axs[1,0].tick_params(labelsize='large')
axs[1,1].tick_params(labelsize='large')

plt.tight_layout()
plt.savefig('cluster_info.pdf')
'''
