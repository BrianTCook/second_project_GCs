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
import matplotlib.cm as cm
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm, Normalize

from cluster_table import sort_clusters_by_attribute


def lattices_maker(points, uniformity, Norbiters):
    
    '''
    takes 6D phase space coordinates
    tessellates phase space s.t. there are < Npoints_max
    '''
    
    #datadir = '/home/s1780638/second_project_GCs/data/'
    datadir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
    cluster_populations = list( np.loadtxt(datadir + 'Nstars_in_clusters.txt') ) 
    
    #need to give clusters sorted by an attribute, in our case increasing |r|
    #new_index = indices_dict[old_index]
    indices_dict = sort_clusters_by_attribute('|r|')
    
    cluster_populations_sorted = [ int(cluster_populations[indices_dict[i]]) for i in range(Norbiters) ]

    grids = []

    for k, number_of_stars in enumerate(cluster_populations_sorted):
        
        starting_index = int(np.sum( cluster_populations_sorted[:k] ))
        ending_index = starting_index + int(number_of_stars)
    
        xvals, yvals, zvals = points[starting_index:ending_index,0], points[starting_index:ending_index,1], points[starting_index:ending_index,2]
        vxvals, vyvals, vzvals = points[starting_index:ending_index,3], points[starting_index:ending_index,4], points[starting_index:ending_index,5]
        
        if uniformity == 'uniform':
        
            Ncells = 8
            
            x_spatial = np.linspace(min(xvals), max(xvals), Ncells)
            y_spatial = np.linspace(min(yvals), max(yvals), Ncells)
            z_spatial = np.linspace(min(zvals), max(zvals), Ncells)
            x_velocity = np.linspace(min(vxvals), max(vxvals), Ncells)
            y_velocity = np.linspace(min(vyvals), max(vyvals), Ncells)
            z_velocity = np.linspace(min(vzvals), max(vzvals), Ncells)
        
        grids.append([x_spatial, y_spatial, z_spatial, x_velocity, y_velocity, z_velocity])

    return grids

def get_6D_fw(points, values, uniformity, Norbiters):
    
    #datadir = '/home/s1780638/second_project_GCs/data/'
    datadir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
    
    cluster_populations = list( np.loadtxt(datadir + 'Nstars_in_clusters.txt') ) 
    
    #need to give clusters sorted by an attribute, in our case increasing |r|
    #new_index = indices_dict[old_index]
    indices_dict = sort_clusters_by_attribute('|r|')
    
    cluster_populations_sorted = [ int(cluster_populations[indices_dict[i]]) for i in range(Norbiters) ]
    
    grids = lattices_maker(points, uniformity, Norbiters)
    
    fw_values_interpolated = []

    for k, number_of_stars in enumerate(cluster_populations_sorted):

        starting_index = int( np.sum( cluster_populations_sorted[:k] ))
        ending_index = starting_index + int(number_of_stars)
            
        X, Y, Z, VX, VY, VZ = np.meshgrid(*grids[k])
    
        #Ti is 6-D interpolation using method=method    
        grid_z0 = griddata(points[starting_index:ending_index, :], values[starting_index:ending_index], 
                           (X, Y, Z, VX, VY, VZ), method='nearest')    
        
        fw_values_interpolated.append(grid_z0)
    
    return fw_values_interpolated

def simpson_6D(grid, fw_values): #Simpson's integration rule, \int_{a}^{b} f(x) dx with N sample points

    x, y, z, vx, vy, vz = grid
    Ncells = len(x)
    
    hx, hy, hz = (max(x)-min(x))/Ncells, (max(y)-min(y))/Ncells, (max(z)-min(z))/Ncells
    hvx, hvy, hvz = (max(vx)-min(vx))/Ncells, (max(vy)-min(vy))/Ncells, (max(vz)-min(vz))/Ncells    

    C = np.ones((Ncells, Ncells, Ncells, Ncells, Ncells, Ncells))

    for i in range(Ncells):
        if i == 0 or i == (Ncells-1):
            C[i,:,:,:,:,:] *= 1.
        if i%2 == 0 and i != 0:
            C[i,:,:,:,:,:] *= 4.
        if i%2 != 0 and i != (Ncells-1):
            C[i,:,:,:,:,:] *= 2.
            
    for j in range(Ncells):
        if j == 0 or j == (Ncells-1):
            C[:,j,:,:,:,:] *= 1.
        if j%2 == 0 and j != 0:
            C[:,j,:,:,:,:] *= 4.
        if j%2 != 0 and j != (Ncells-1):
            C[:,j,:,:,:,:] *= 2.
        
    for k in range(Ncells):
        if k == 0 or k == (Ncells-1):
            C[:,:,k,:,:,:] *= 1.
        if k%2 == 0 and k != 0:
            C[:,:,k,:,:,:] *= 4.
        if k%2 != 0 and k != (Ncells-1):
            C[:,:,k,:,:,:] *= 2.
            
    for l in range(Ncells):
        if l == 0 or l == (Ncells-1):
            C[:,:,:,l,:,:] *= 1.
        if l%2 == 0 and l != 0:
            C[:,:,:,l,:,:] *= 4.
        if l%2 != 0 and l != (Ncells-1):
            C[:,:,:,l,:,:] *= 2.
            
    for m in range(Ncells):
        if m == 0 or m == (Ncells-1):
            C[:,:,:,:,m,:] *= 1.
        if m%2 == 0 and m != 0:
            C[:,:,:,:,m,:] *= 4.
        if m%2 != 0 and m != (Ncells-1):
            C[:,:,:,:,m,:] *= 2.
            
    for n in range(Ncells):
        if n == 0 or n == (Ncells-1):
            C[:,:,:,:,:,n] *= 1.
        if n%2 == 0 and n != 0:
            C[:,:,:,:,:,n] *= 4.
        if n%2 != 0 and n != (Ncells-1):
            C[:,:,:,:,:,n] *= 2.

    F = fw_values
    
    I_normalize = np.sum( [ C[i,j,k,l,m,n]*F[i,j,k,l,m,n] 
                            for i in range(Ncells) for j in range(Ncells) 
                            for k in range(Ncells) for l in range(Ncells)
                            for m in range(Ncells) for n in range(Ncells) ] ) 

    I_normalize *= (hx*hy*hz*hvx*hvy*hvz)/(3**6.)
    
    F_normalized = fw_values / I_normalize
    
    G = np.multiply(F_normalized, np.log(F_normalized))
        
    I = np.sum( [ C[i,j,k,l,m,n]*G[i,j,k,l,m,n] 
                    for i in range(Ncells) for j in range(Ncells) 
                    for k in range(Ncells) for l in range(Ncells)
                    for m in range(Ncells) for n in range(Ncells) ] ) 
    
    I *= (hx*hy*hz*hvx*hvy*hvz)/(3**6.)
    
    return I
        

def get_entropy(points, values, uniformity, Norbiters):
    
    '''
    right now have N x N x N x N x N x N array
    need to integrate s.t. we get a 1 x 1 x 1 x 1 x 1 x 1 array
    the value within it will be the information entropy
    '''
    
    #datadir = '/home/s1780638/second_project_GCs/data/'
    datadir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
    
    cluster_populations = list( np.loadtxt(datadir + 'Nstars_in_clusters.txt') ) 
    
    #need to give clusters sorted by an attribute, in our case increasing |r|
    #new_index = indices_dict[old_index]
    indices_dict = sort_clusters_by_attribute('|r|')

    cluster_populations_sorted = [ int(cluster_populations[indices_dict[i]]) for i in range(Norbiters) ]
    
    fw_values_interpolated = get_6D_fw(points, values, uniformity, Norbiters)
    grids = lattices_maker(points, uniformity, Norbiters)

    #S = 0 #entropy
    entropies = []

    for k, number_of_stars in enumerate(cluster_populations_sorted):
    
        if k <= 7:
        
            integral = simpson_6D(grids[k], fw_values_interpolated[k])  
            entropies.append(-integral)
            #S -= integral
        
    #returns entropy
    #return S
    return entropies

if __name__ in '__main__':
    
    #datadir = '/home/s1780638/second_project_GCs/data/'
    
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    
    #plot whole picture in (x,z) plane
    fig, axs = plt.subplots(1, 3)
    l = 0
    
    sim_times = np.linspace(0., 100., 51)  
    gadget_flags = [ 0, 25, 50 ]
    
    logN_max = 6
    zvals_all = [ i for i in range(2**3) ] #[ i for i in range(2**logN_max) ]
    zvals_filler = [ 0 for i in range(8, 2**6) ]
    
    zvals_all = zvals_all + zvals_filler
    
    for i, t in enumerate(sim_times):
        
        if i in gadget_flags:
            
            frame = str(i*10).rjust(5, '0')

            print('time: %.02f Myr'%(t))
    
            datadir_AMUSE = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/Enbid-2.0/AMUSE_data/'
            point_files = glob.glob(datadir_AMUSE+'*'+frame+'*.ascii')
            
            t0 = time.time()
            
            xvals_all = [ [ -2 for j in range(2**logN_max) ] for i in range(logN_max+1) ]
            yvals_all = [ [ 0 for j in range(2**logN_max) ] for i in range(logN_max+1) ]

            #plt.xlabel('Simulation Time (Myr)', fontsize=12)
            #plt.ylabel(r'$S/N_{\mathrm{clusters}}$ (entropy units per cluster)', fontsize=12)
            
            '''
            if l == 0:
                axs[l].set_ylabel(r'$S$ (entropy units)', fontsize=16)
            
            if l == 1:
                axs[l].set_xlabel(r'$\log_{2}N_{\mathrm{clusters}}$', fontsize=16)
            '''
                
            axs[l].set_xticks([0, 1, 2, 3, 4, 5, 6])
            axs[l].tick_params(labelsize='large')
            
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
                    #values = np.asarray([ float(val)/np.sum(values) for val in values])
                    
                    uniformity = 'uniform'
                    #S = get_entropy(points, values, uniformity, Nclusters)
                    entropies = get_entropy(points, values, uniformity, Nclusters)
                    #fw_map, x0, x1, y0, y1 = get_map(points, values)
                    #print('number of unique elements in fw_map: ', np.unique(fw_map).shape)
                    
                    sim_time = frame/500. * 100. #frame/frame * end_time (Myr)
                    
                    tf = time.time()
                    
                    logN = int(np.log2(Nclusters)) #number of clusters
                    
                    #plotting entropy as a function of time, saving by log2N
                    #xvals_all[logN].append(sim_time)
                    #yvals_all[logN].append(S/Nclusters)
                    
                    print('logN, len(entropies) are', logN, len(entropies))
                    
                    xvals_all[logN][:len(entropies)] = [ logN for i in range(len(entropies)) ]
                    yvals_all[logN][:len(entropies)] = entropies
                
                    #print('normalized entropy: %.03e (kpc km/s)^3'%(S/Nclusters))
                
                    tf = time.time()
                    print('current time: %.04f minutes'%((tf-t0)/60.))
            
            for logN in range(logN_max+1):
                
                #plt.plot(xvals, yvals, label=r'$\log_{2} N_{\mathrm{clusters}}$ = %i'%(logN), linewidth=1)
                
                xvals, yvals = xvals_all[logN], yvals_all[logN]
                
                cmap = cm.get_cmap('nipy_spectral', 8)
                sc = axs[l].scatter(xvals, yvals, c=zvals_all, s=4, alpha=1.0, cmap=cmap, norm=Normalize(vmin=-.5, vmax=7.5))
                
                if logN == logN_max and l == 2:
            
                    cbar = fig.colorbar(sc)
                    cbar.set_label(r'Index', rotation=270, labelpad=12)
            
            #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
            axs[l].set_title(r'$t_{\mathrm{sim}} = %.02f$ Myr'%(t), fontsize=8)#, xy=(0.1, 0.85 ), xycoords='axes fraction', fontsize=6)
            axs[l].set_xlim(0-0.5, logN_max+0.5)
            axs[l].set_ylim(-20, 20)
            
            l += 1
    
    axs[2].set_yticklabels([])    
    axs[0].set_ylabel(r'$S$ (entropy units)', fontsize=14)   
    axs[1].set_xlabel(r'$\log_{2}N_{\mathrm{clusters}}$', fontsize=14)
    
    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()
    
    plt.suptitle('Differential Entropy', fontsize=18)
    #plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('entropy_evolution_interactions_firsteight.pdf')
    #plt.savefig('entropy_evolution_time.pdf')