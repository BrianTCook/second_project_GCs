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
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
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
    
    N = 64 #total number of initialized clusters I have
    
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
    radii = [ round(radii[i], 2) for i in range(N) ]
    
    #df = pd.DataFrame(list(zip(masses, Nstars, dists, speeds, radii)), columns=['M', 'Nstars', '|r|', '|v|', 'rvir'])
    df = pd.DataFrame(list(zip(dists, xs, ys, zs, vxs, vys, vzs, Nstars)), columns=['|r|', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'Nstars'])
    df_sorted_by_r = df.sort_values(by=[attribute])
    #df_sorted_by_r = df_sorted_by_r.reset_index(drop=True)
    
    indices_dict = {}
    
    for df, df_sorted in zip(df.index, df_sorted_by_r.index):
        indices_dict.update( {df : df_sorted} )
        
    return indices_dict#, masses, dists

if __name__ in '__main__':
    
    indices_dict = sort_clusters_by_attribute('|r|') #

    '''
    plt.rc('font', family='serif')
    plt.rc('text', usetex=True)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    color_list_xyz = []
    
    N = 64
    
    for i in range(N):
        
        color_i = 'C' + str(math.ceil(np.log2(i+1)))
        inew = indices_dict[i]
        print(i, inew)
        
        if color_i not in color_list_xyz:
            
            min_i, max_i = int(str(math.ceil(np.log2(i+1))))-1, int(str(math.ceil(np.log2(i+1))))
            min_i, max_i = 2**(min_i), 2**(max_i)-1
        
            Nstars = df_sorted_by_r['Nstars'].tolist()
        
            xs = df_sorted_by_r['x'].tolist()
            ys = df_sorted_by_r['y'].tolist()
            zs = df_sorted_by_r['z'].tolist()
            vxs = df_sorted_by_r['vx'].tolist()
            vys = df_sorted_by_r['vy'].tolist()
            vzs = df_sorted_by_r['vz'].tolist()
            
            if i == 0 or i == 1:
                
                ax.scatter(xs[inew], ys[inew], zs[inew], s=2*math.ceil(np.log10(Nstars[inew])), 
                           c=color_i, label=r'Index: %i'%(i))
            
            if i == 2:
                
                ax.scatter(xs[inew], ys[inew], zs[inew], s=2*math.ceil(np.log10(Nstars[inew])), 
                           c=color_i, label=r'Indices: %i, %i'%(i, i+1))
            
            if min_i > 2:
                
                ax.scatter(xs[inew], ys[inew], zs[inew], s=2*math.ceil(np.log10(Nstars[inew])), 
                           c=color_i, label=r'Indices: $%i, \dots, %i$'%(min_i, max_i))
            
            
            ax.quiver(xs[inew], ys[inew], zs[inew], 0.0008*vxs[inew], 0.0008*vys[inew], 0.0008*vzs[inew],
                      color=color_i)
            
            color_list_xyz.append(color_i)
            
        else:
            
            ax.scatter(xs[inew], ys[inew], zs[inew], s=2*math.ceil(np.log10(Nstars[inew])), c=color_i)
            ax.quiver(xs[inew], ys[inew], zs[inew], 0.0008*vxs[inew], 0.0008*vys[inew], 0.0008*vzs[inew],
                      color=color_i)
            
        
    ax.set_xlabel(r'$x$ (kpc)', fontsize=16)
    ax.set_ylabel(r'$y$ (kpc)', fontsize=16)
    ax.set_zlabel(r'$z$ (kpc)', fontsize=16)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_zticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    plt.legend(loc='best', fontsize=8)
    plt.tight_layout()
    plt.savefig('xyz_plane_initial.pdf')
    
    plt.figure()
        
    color_list_xy = []
    
    for i in range(N):
    
        color_i = 'C' + str(math.ceil(np.log2(i+1)))
        inew = indices_dict[i]
        
        if color_i not in color_list_xy:
            
            if i == 0:
        
                plt.scatter(xs[inew], ys[inew], s=2*math.ceil(np.log10(Nstars[inew])), 
                            c=color_i, label=r'Index: 0')
                
            if i == 1:
        
                plt.scatter(xs[inew], ys[inew], s=2*math.ceil(np.log10(Nstars[inew])), 
                            c=color_i, label=r'Index: 1')
                
            if i == 2:
                
                plt.scatter(xs[inew], ys[inew], s=2*math.ceil(np.log10(Nstars[inew])), 
                            c=color_i, label=r'Index: 2, 3')
            
            if i >= 4:
        
                plt.scatter(xs[inew], ys[inew], s=2*math.ceil(np.log10(Nstars[inew])), 
                            c=color_i, label=r'Indices: $%i, \dots, %i$'%(i, 2*i - 1))
                
            plt.arrow(xs[inew], ys[inew], 0.0001*vxs[inew], 0.0001*vys[inew], 
                      alpha = 0.6, color=color_i)
            
            color_list_xy.append(color_i)
            
        else:
            
            plt.scatter(xs[i], ys[i], s=2*math.ceil(np.log10(Nstars[i])), c=color_i)
            plt.arrow(xs[i], ys[i], 0.0002*vxs[i], 0.0002*vys[i], color=color_i)
        
    plt.xlim(-0.2, 0.2)
    plt.ylim(-0.2, 0.2)
    plt.xlabel(r'$x$ (kpc)', fontsize=16)
    plt.ylabel(r'$y$ (kpc)', fontsize=16)
    plt.gca().set_xticks([-0.2, -0.1, 0, 0.1, 0.2])
    plt.gca().set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
    plt.gca().set_aspect('equal')
    plt.gca().tick_params(labelsize='large')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    plt.tight_layout()
    plt.savefig('xy_plane_initial.pdf')

    '''

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

plt.figure()

plt.scatter(masses, radii, s=2, c='k')
plt.xlabel(r'$M_{\mathrm{cluster}}$ ($M_{\odot}$)', fontsize=20)
plt.ylabel(r'$r_{\mathrm{vir}, 0}$ (pc)', fontsize=20)
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().tick_params(labelsize='large')
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig('mass_radius_relation.pdf')
plt.close()
'''
