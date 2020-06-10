#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:17:37 2020

@author: BrianTCook
"""

import numpy as np
import time
import glob
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pylab as pl
from sklearn.cluster import KMeans
from cluster_table import sort_clusters_by_attribute
pd.options.mode.chained_assignment = None  # default='warn'

def sklearn_mapper(true_labels, kmeans_labels):
    
    input_output_dict = {}

    pairings = [ (x, y) for x, y in zip(true_labels, kmeans_labels) ]
    pairings_counted = [ [x,pairings.count(x)] for x in set(pairings) ]

    for counted_pairing in pairings_counted:
        
        pairing, count = counted_pairing
        first, second = pairing
        
        if first not in input_output_dict.keys():
            
            input_output_dict.update( {first: second} )
        
        if count > pairings.count((first, input_output_dict[first])):
            
            input_output_dict.pop(first)
            input_output_dict.update( {first: second} )
    
    return input_output_dict

# http://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit
# for arbitrary axes
def ellipsoid_fit(X):
    x = X[:, 0]
    y = X[:, 1]
    z = X[:, 2]
    D = np.array([x * x + y * y - 2 * z * z,
                 x * x + z * z - 2 * y * y,
                 2 * x * y,
                 2 * x * z,
                 2 * y * z,
                 2 * x,
                 2 * y,
                 2 * z,
                 1 - 0 * x])
    d2 = np.array(x * x + y * y + z * z).T # rhs for LLSQ
    u = np.linalg.solve(D.dot(D.T), D.dot(d2))
    a = np.array([u[0] + 1 * u[1] - 1])
    b = np.array([u[0] - 2 * u[1] - 1])
    c = np.array([u[1] - 2 * u[0] - 1])
    v = np.concatenate([a, b, c, u[2:]], axis=0).flatten()
    A = np.array([[v[0], v[3], v[4], v[6]],
                  [v[3], v[1], v[5], v[7]],
                  [v[4], v[5], v[2], v[8]],
                  [v[6], v[7], v[8], v[9]]])

    center = np.linalg.solve(- A[:3, :3], v[6:9])

    translation_matrix = np.eye(4)
    translation_matrix[3, :3] = center.T

    R = translation_matrix.dot(A).dot(translation_matrix.T)

    evals, evecs = np.linalg.eig(R[:3, :3] / -R[3, 3])
    evecs = evecs.T

    radii = np.sqrt(1. / np.abs(evals))
    radii *= np.sign(evals)

    return center, evecs, radii

def get_kmeans_result(snapshots, Norbiters, initial_masses):
    
    #t0 = time.time()
    
    datadir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
    datadir_AMUSE = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/Enbid-2.0/AMUSE_data/'
    cluster_populations = np.loadtxt(datadir + 'Nstars_in_clusters.txt')
    cluster_radii = np.loadtxt(datadir + '/ICs/cluster_radii_for_sampling.txt') 
    
    indices_dict, cluster_masses = sort_clusters_by_attribute('|r|')
            
    cluster_populations_sorted = [ cluster_populations[indices_dict[i]] for i in range(Norbiters) ]        
    cluster_radii_sorted = [ cluster_radii[indices_dict[i]] for i in range(Norbiters) ]
    cluster_masses_sorted = [ cluster_masses[indices_dict[i]] for i in range(Norbiters) ]

    true_labels = []
    
    for i in range(Norbiters):
        true_labels += [ i for j in range(int(cluster_populations_sorted[i])) ]
    
    delta_max = 0.9
    all_deltas = [ [] for i in range(Norbiters) ]
    endpoints = [ '' for i in range(Norbiters) ]
    active_orbiters = [ i for i in range(Norbiters) ]
    
    star_masses = np.loadtxt(datadir+'tree_data/tree_SingleCluster_masses_Norbiters_64_dt_0.2.txt').T
    
    deltaMs = np.zeros(star_masses[0,:].shape)
    
    for i in range(1, len(deltaMs)):
        deltaMs[i] = np.mean(star_masses[:, i]) - np.mean(star_masses[:, i-1])

    deltaM = np.mean(deltaMs)
    centroids = np.zeros((Norbiters, 6))

    for k, snapshot in enumerate(snapshots):
        
        print(snapshot)
        
        if len(active_orbiters) == 0:
            
            return all_deltas, endpoints
        
        data_filename = glob.glob(datadir_AMUSE+'*_%s_Norbiters_%i.ascii'%(snapshot, Norbiters))
        data_6D = np.loadtxt(data_filename[0])
        
        I6, J6 = data_6D.shape
        data_6D_rescaled = [ ( (data_6D[i,j] - np.amin(data_6D[:,j])) / (np.amax(data_6D[:,j]) - np.amin(data_6D[:,j])) ) 
                             for i in range(I6) for j in range(J6) ]
        
        np_data_6D_rescaled = np.asarray(data_6D_rescaled)
        np_data_6D_rescaled = np.reshape(np_data_6D_rescaled, data_6D.shape)
    
        #apply kmeans clustering
        #kmeans.cluster_centers_
        if k == 0:
            
            kmeans_6D = KMeans(n_clusters=len(active_orbiters), init='k-means++')
        
        if len(active_orbiters) > 0 and k != 0:
        
            kmeans_6D = KMeans(n_clusters=len(active_orbiters), init=centroids)
            
        if len(active_orbiters) == 0:
            
            kmeans_6D = KMeans(n_clusters=1, init=centroids)
        
        kmeans_6D.fit(np_data_6D_rescaled)
        centroids = kmeans_6D.cluster_centers_
        
        print(centroids)
        
        y_kmeans_6D = kmeans_6D.predict(np_data_6D_rescaled)
        
        io_dict_6D = sklearn_mapper(true_labels, y_kmeans_6D)
    
        #k-means clustering with same labelling scheme as true_labels
        y_compare_6D = [ io_dict_6D[y] for y in y_kmeans_6D ]
    
        all_data = np.concatenate((data_6D, np_data_6D_rescaled), axis=1)
    
        df = pd.DataFrame(all_data,
                          columns=['x', 'y', 'z', 'vx', 'vy', 'vz', 
                                   'x (rescaled)', 'y (rescaled)', 'z (rescaled)', 
                                   'vx (rescaled)', 'vy (rescaled)', 'vz (rescaled)'])

        df['labels'] = y_compare_6D
        
        star_masses_truncated = star_masses[:len(df.index), k]
        
        df['masses'] = star_masses_truncated

        for cluster_label in active_orbiters:
            
            #print('cluster_label: ', cluster_label)
            #print('current time: %.04f minutes'%((time.time()-t0)/60.))
            
            df_cluster = df.loc[df['labels'] == cluster_label]
            df_cluster['separation distance'] = '' #in parsecs
            m_cluster_init = cluster_masses_sorted[cluster_label]
                    
            #sort by distance from COM of cluster
            xc, yc, zc = np.median(df_cluster['x'].tolist()), np.median(df_cluster['y'].tolist()), np.median(df_cluster['z'].tolist())
        
            for i in df_cluster.index:
                
                dist_sq = (df_cluster.at[i, 'x']-xc)**2 + (df_cluster.at[i, 'y']-yc)**2 + (df_cluster.at[i, 'z']-zc)**2 
                dist_in_pc = np.sqrt(dist_sq) * 1000.
                df_cluster.loc[i, 'separation distance'] = dist_in_pc
            
            #print('median separation distance, initial radius: %.03f pc, %.03f pc'%(np.median(df_cluster['separation distance'].tolist()), cluster_radii_sorted[cluster_label] ))
            
            df_cluster = df_cluster[ df_cluster['separation distance'] <= 2.*cluster_radii_sorted[cluster_label] ] #outside of original radius
            df_cluster = df_cluster.reset_index(drop=True)
            m_cluster_t = np.sum(df_cluster['masses'].tolist())
            nstars = len(df_cluster.index)
            
            '''                
            add column saying current label for plot showing N_in birth cluster,
            N_field, and N_adopted by other cluster?
            '''     

            if len(df_cluster.index) != 0:

                m_cluster_t = np.sum(df_cluster['masses'].tolist())
                eps = nstars * deltaM / m_cluster_t
                
                delta = 1. - m_cluster_t / m_cluster_init * (1 + eps)
                
                galactocentric_dist = np.sqrt( np.median(df_cluster['x'].tolist())**2. + np.median(df_cluster['y'].tolist())**2. + np.median(df_cluster['z'].tolist())**2. )
                endpoints[cluster_label] = 'final galactocentric distance: %.03f kpc'%(galactocentric_dist)
                
            else:
                
                delta = 1.
                
            all_deltas[cluster_label].append(delta)
            
            if delta >= delta_max:
                
                cluster_ind = active_orbiters.index(cluster_label)
                centroids = np.delete(centroids, cluster_ind, 0) #removes centroid from consideration
                endpoints[cluster_label] = 'disruption time: %.00f Myr'%(k*2.)
                active_orbiters.remove(cluster_label) #removes orbiter label from active ones
                
    return all_deltas, endpoints

if __name__ in '__main__':

    snapshots = [ str(j*10).rjust(5, '0') for j in range(51) ]
    
    initial_masses = 0.
    
    Norbiters = 64
    all_deltas, endpoints = get_kmeans_result(snapshots, Norbiters, initial_masses) #endpoints
    dt = 2. #Myr
    
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    
    colors = pl.cm.rainbow(np.linspace(0,1,Norbiters))
    
    for cluster, deltas in enumerate(all_deltas):
        
        print('cluster: %i, fate: %s'%(cluster, endpoints[cluster]))
        
        y = deltas
        x = np.arange(1e-5, dt*len(deltas), dt) #sim_times
        
        plt.plot(x, y, linewidth=0.5, alpha=0.6, color=colors[cluster], label=str(cluster))
        
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=4, ncol=2)
    plt.axhline(y=0.9, linestyle='--', c='k', linewidth=1)
    plt.xlim(0., 100.)
    plt.ylim(0., 1.)
    plt.xlabel(r'$t_{\mathrm{sim}}$ (Myr)', fontsize=12)
    plt.ylabel(r'$\delta \simeq 1 - M_{\mathrm{cluster}}(t)/M_{\mathrm{cluster}}(t=0)$', fontsize=12)
    plt.tight_layout()
    plt.savefig('mass_loss_Norbiters_%i_kmeans.pdf'%(Norbiters))    
    
    print('hello world!')
    
'''
success plotting

    hits_2D = [ 1 if y1 == y2 else 0 for y1, y2 in zip(true_labels, y_compare_2D) ]
    hits_3D = [ 1 if y1 == y2 else 0 for y1, y2 in zip(true_labels, y_compare_3D) ]
    hits_6D = [ 1 if y1 == y2 else 0 for y1, y2 in zip(true_labels, y_compare_6D) ]
    
    if hits_2D.count(1) > hits_2D.count(0):
        
        success_2D = np.sum(hits_2D)/len(hits_2D)
        
    if hits_2D.count(1) <= hits_2D.count(0):
        
        success_2D = 1. - np.sum(hits_2D)/len(hits_2D)
        
    if hits_3D.count(1) > hits_3D.count(0):
        
        success_3D = np.sum(hits_3D)/len(hits_3D)
        
    if hits_3D.count(1) <= hits_3D.count(0):
        
        success_3D = 1. - np.sum(hits_3D)/len(hits_3D)
        
    if hits_6D.count(1) > hits_6D.count(0):
        
        success_6D = np.sum(hits_6D)/len(hits_6D)
        
    if hits_6D.count(1) <= hits_6D.count(0):
        
        success_6D = 1. - np.sum(hits_6D)/len(hits_2D)


    if plot_2D == True:
        
        plt.figure()    
        plt.scatter(data_3D[:, 1], data_3D[:, 2], c=y_kmeans_2D, s=5, cmap='Dark2')
        plt.annotate('2D coordinates', xy=(0.6, 0.2), xycoords='axes fraction', fontsize=12)
        plt.annotate('Success rate: %.05f'%(success_2D), xy=(0.6, 0.1), xycoords='axes fraction', fontsize=12)
        plt.xlabel(r'$y$ (kpc)', fontsize=20)
        plt.ylabel(r'$z$ (kpc)', fontsize=20)
        
        plt.xlim(-2, 2)
        plt.ylim(-1, 1)
        
        plt.title(r'$k$-means, $N_{clusters} = %i$'%(Norbiters), fontsize=16)
        plt.tight_layout()
        plt.savefig('kmeans_Norbiters_%i_2D.jpg'%(Norbiters))
        plt.close()    
    
    if plot_3D == True:
        
        plt.figure()    
        plt.scatter(data_3D[:, 1], data_3D[:, 2], c=y_kmeans_3D, s=5, cmap='Dark2')
        plt.annotate('3D coordinates', xy=(0.6, 0.2), xycoords='axes fraction', fontsize=12)
        plt.annotate('Success rate: %.05f'%(success_3D), xy=(0.6, 0.1), xycoords='axes fraction', fontsize=12)
        plt.xlabel(r'$y$ (kpc)', fontsize=20)
        plt.ylabel(r'$z$ (kpc)', fontsize=20)
        
        plt.xlim(-2, 2)
        plt.ylim(-1, 1)
        
        plt.title(r'$k$-means, $N_{clusters} = %i$'%(Norbiters), fontsize=16)
        plt.tight_layout()
        plt.savefig('kmeans_Norbiters_%i_3D.jpg'%(Norbiters))
        plt.close()
        
    if plot_6D == True:
        
        plt.figure()    
        plt.scatter(data_6D[:, 1], data_6D[:, 2], c=y_kmeans_6D, s=5, cmap='Dark2')
        plt.annotate('6D coordinates', xy=(0.6, 0.2), xycoords='axes fraction', fontsize=12)
        plt.annotate('Success rate: %.05f'%(success_6D), xy=(0.6, 0.1), xycoords='axes fraction', fontsize=12)
        plt.xlabel(r'$y$ (kpc)', fontsize=20)
        plt.ylabel(r'$z$ (kpc)', fontsize=20)
        
        plt.xlim(-2, 2)
        plt.ylim(-1, 1)
        
        plt.title(r'$k$-means, $N_{clusters} = %i$'%(Norbiters), fontsize=16)
        plt.tight_layout()
        plt.savefig('kmeans_Norbiters_%i_6D.jpg'%(Norbiters))
        plt.close()
    
    return success_2D, success_3D, success_6D
   
for logN in range(logN_max+1):
    
    print('N = %i'%(2**logN))
    if plot_2D == True or plot_3D == True or plot_6D == True:
        
        plt.rc('text', usetex = True)
        plt.rc('font', family = 'serif')
    
    s_2, s_3, s_6 = get_kmeans_result(2**logN, plot_3D, plot_6D)
    s2_vals.append(s_2)
    s3_vals.append(s_3)
    s6_vals.append(s_6)
    
plt.figure()
plt.plot(logN_vals, s2_vals, label='2D distance metric')
plt.plot(logN_vals, s3_vals, label='3D distance metric')
plt.plot(logN_vals, s6_vals, label='6D distance metric')
plt.ylim(0, 1)
plt.legend(loc='best')
plt.title(r'$k$-means Cluster Identification', fontsize=16)
plt.xlabel(r'$\log_{2} N_{\mathrm{clusters}}$', fontsize=14)
plt.ylabel('Labelling Accuracy', fontsize=14)
plt.tight_layout()
plt.savefig('accuracies_kmeans.pdf')
    
plt.figure()
plt.plot(logN_vals, s2_vals, label='2D distance metric')
plt.plot(logN_vals, s3_vals, label='3D distance metric')
plt.plot(logN_vals, s6_vals, label='6D distance metric')
plt.ylim(0, 1)
plt.legend(loc='best')
plt.title(r'$k$-means Cluster Identification', fontsize=16)
plt.xlabel(r'$\log_{2} N_{\mathrm{clusters}}$', fontsize=14)
plt.ylabel('Labelling Accuracy', fontsize=14)
plt.tight_layout()
plt.savefig('accuracies_kmeans.pdf')
'''
