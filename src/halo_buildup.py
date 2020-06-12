#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 18:07:43 2020

@author: BrianTCook
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:17:37 2020

@author: BrianTCook
"""

import numpy as np
import glob
import time
import matplotlib.pyplot as plt
import pandas as pd
from cluster_table import sort_clusters_by_attribute
pd.options.mode.chained_assignment = None  # default='warn'

def buildup(snapshots, Norbiters, initial_masses):
    
    t0 = time.time()
    
    datadir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
    datadir_AMUSE = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/Enbid-2.0/AMUSE_data/'
    cluster_populations = np.loadtxt(datadir + 'Nstars_in_clusters.txt')
    cluster_radii = np.loadtxt(datadir + '/ICs/cluster_radii_for_sampling.txt') 
    
    indices_dict, df_sorted_by_r = sort_clusters_by_attribute('|r|')
            
    cluster_populations_sorted = [ cluster_populations[indices_dict[i]] for i in range(Norbiters) ]        
    cluster_radii_sorted = [ cluster_radii[indices_dict[i]] for i in range(Norbiters) ]

    true_labels = []
    
    for i in range(Norbiters):
        true_labels += [ i for j in range(int(cluster_populations_sorted[i])) ]
        
    active_orbiters = [ i for i in range(Norbiters) ]
    
    inner_half_clusters = np.zeros(len(snapshots))
    inner_half_field = np.zeros(len(snapshots))
    outer_half_clusters = np.zeros(len(snapshots))
    outer_half_field = np.zeros(len(snapshots))
    
    for k, snapshot in enumerate(snapshots):
        
        print(snapshot)
        print('time: %.02f minutes'%((time.time()-t0)/60.))
        
        data_filename = glob.glob(datadir_AMUSE+'*_%s_Norbiters_%i.ascii'%(snapshot, Norbiters))
        data_3D = np.loadtxt(data_filename[0])[:, :3]
    
        df = pd.DataFrame(data_3D, columns=['x', 'y', 'z'])
        df['labels'] = true_labels

        for cluster_label in active_orbiters:
            
            #print('cluster_label: ', cluster_label)
            #print('current time: %.04f minutes'%((time.time()-t0)/60.))
            
            df_cluster = df.loc[df['labels'] == cluster_label]
            df_cluster['separation distance'] = '' #in parsecs
                    
            #sort by distance from COM of cluster
            xc, yc, zc = np.median(df_cluster['x'].tolist()), np.median(df_cluster['y'].tolist()), np.median(df_cluster['z'].tolist())
        
            for i in df_cluster.index:
                
                dist_sq = (df_cluster.at[i, 'x']-xc)**2 + (df_cluster.at[i, 'y']-yc)**2 + (df_cluster.at[i, 'z']-zc)**2 
                dist_in_pc = np.sqrt(dist_sq) * 1000.
                df_cluster.loc[i, 'separation distance'] = dist_in_pc

            #print('median separation distance, initial radius: %.03f pc, %.03f pc'%(np.median(df_cluster['separation distance'].tolist()), cluster_radii_sorted[cluster_label] ))
            
            df_cluster = df_cluster[ df_cluster['separation distance'] <= 2.*cluster_radii_sorted[cluster_label] ] #outside of original radius
            df_cluster = df_cluster.reset_index(drop=True)
            nstars = len(df_cluster.index)
            
            galactocentric_distance = np.sqrt(xc**2. + yc**2. + zc**2.)
            
            if galactocentric_distance <= 0.5: 
                
                inner_half_clusters[k] += nstars
                inner_half_field[k] += cluster_populations_sorted[cluster_label] - nstars
            
            if galactocentric_distance > 0.5: 
                
                outer_half_clusters[k] += nstars
                outer_half_field[k] += cluster_populations_sorted[cluster_label] - nstars
            
            '''                
            add column saying current label for plot showing N_in birth cluster,
            N_field, and N_adopted by other cluster?
            '''     
                
    return inner_half_clusters, inner_half_field, outer_half_clusters, outer_half_field

if __name__ in '__main__':

    snapshots = [ str(j*10).rjust(5, '0') for j in range(51) ]
    
    initial_masses = 0.
    
    Norbiters = 64
    ihc, ihf, ohc, ohf = buildup(snapshots, Norbiters, initial_masses) #endpoints
    dt = 2. #Myr
    
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')

    sim_times = np.arange(0., len(snapshots)*dt, dt)

    print('last ihc: %i'%(ihc[len(snapshots)-1]))
    print('last ihf: %i'%(ihf[len(snapshots)-1]))
    print('last ohc: %i'%(ohc[len(snapshots)-1]))
    print('last ohf: %i'%(ohf[len(snapshots)-1]))

    plt.semilogy(sim_times, ihc, color='C0', linewidth=1, label=r'Cluster Stars, $|\mathbf{r}| \leq 0.5$ kpc')
    plt.semilogy(sim_times, ohc, color='C1', linewidth=1, label=r'Cluster Stars, $|\mathbf{r}| > 0.5$ kpc')
    plt.semilogy(sim_times, ihf, color='C0', linestyle='--', linewidth=1, label=r'Field Stars, $|\mathbf{r}| \leq 0.5$ kpc')
    plt.semilogy(sim_times, ohf, color='C1', linestyle='--', linewidth=1, label=r'Field Stars, $|\mathbf{r}| > 0.5$ kpc')
        
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10, ncol=1)
    plt.xlim(0., 100.)
    plt.xlabel(r'$t_{\mathrm{sim}}$ (Myr)', fontsize=12)
    plt.ylabel(r'$N_{\star}$', fontsize=12)
    plt.gca().tick_params(labelsize='large')
    plt.tight_layout()
    plt.savefig('buildup_Norbiters_%i.pdf'%(Norbiters))    
    
    print('hello world!')