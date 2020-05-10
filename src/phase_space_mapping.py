import gzip
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from cluster_table import sort_clusters_by_attribute

import glob

def maps():
    
    datadir_AMUSE = '/Users/BrianTCook/Desktop/Thesis/second_project_gcs/Enbid-2.0/AMUSE_data/'
    
    logN_max = 6
    Norbiters_list = [ 2**6 ] #int(2**i) for i in range(logN_max)]
    
    datadir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
    
    cluster_populations_raw = np.loadtxt(datadir+'Nstars_in_clusters.txt')
    indices_dict = sort_clusters_by_attribute('|r|')
    cluster_populations_sorted = [ cluster_populations_raw[indices_dict[i]] for i in range(2**logN_max) ]
    
    plt.rc('font', family='serif')
    
    for Norbiters in Norbiters_list:
    
        cluster_populations = list(cluster_populations_sorted)[:Norbiters]
        
        sim_times = np.linspace(0., 100., 51)    
        
        for i, t in enumerate(sim_times):
            
            if i == 0 or i == 25 or i == 50:
                
                print('time: %.02f Myr'%(t))
                
                phase_space_data = np.loadtxt(datadir_AMUSE+'enbid_tree_frame_%s_Norbiters_%i.ascii'%(str(i*10).rjust(5, '0'), Norbiters))

                for k, number_of_stars in enumerate(cluster_populations):
            
                    if k < 10 or k > 53:
                        
                        print('Index: %i'%(k))
                    
                        starting_index = int(np.sum( cluster_populations[:k] ))
                        ending_index = starting_index + int(number_of_stars)
                
                        x = phase_space_data[starting_index:ending_index, 0]
                        y = phase_space_data[starting_index:ending_index, 1]
                        z = phase_space_data[starting_index:ending_index, 2]
                        vx = phase_space_data[starting_index:ending_index, 3]
                        vy = phase_space_data[starting_index:ending_index, 4]
                        vz = phase_space_data[starting_index:ending_index, 5]
                        
                        dx = (np.percentile(x, 95) - np.percentile(x, 5))/2.
                        dy = (np.percentile(y, 95) - np.percentile(y, 5))/2.
                        dz = (np.percentile(z, 95) - np.percentile(z, 5))/2.
                        dvx = (np.percentile(vx, 95) - np.percentile(vx, 5))/2.
                        dvy = (np.percentile(vy, 95) - np.percentile(vy, 5))/2.
                        dvz = (np.percentile(vz, 95) - np.percentile(vz, 5))/2.

                        fig, axs = plt.subplots(5, 5)
                        
                        #first column            
                        
                        xy = np.vstack([x,y])
                        z00 = gaussian_kde(xy)(xy)
                        idx = z00.argsort()
                        x, y, colors = x[idx], y[idx], z00[idx]
                        axs[0, 0].scatter(x, y, s=0.5, c=colors, edgecolor='')
                        axs[0, 0].set_ylabel(r'$y$ (kpc)', fontsize=8)
                        axs[0, 0].tick_params(labelsize='xx-small')
                        axs[0, 0].set_title(r'$x = %.03f \pm %.03f$ kpc'%(np.median(x), dx), fontsize=3)

                        
                        xz = np.vstack([x,z])
                        z10 = gaussian_kde(xz)(xz)
                        idx = z10.argsort()
                        x, z, colors = x[idx], z[idx], z10[idx]
                        axs[1, 0].scatter(x, z, s=0.5, c=colors, edgecolor='')
                        axs[1, 0].set_ylabel(r'$z$ (kpc)', fontsize=8)
                        axs[1, 0].tick_params(labelsize='xx-small')

                        xvx = np.vstack([x,vx])
                        z20 = gaussian_kde(xvx)(xvx)
                        idx = z20.argsort()
                        x, vx, colors = x[idx], vx[idx], z20[idx]                        
                        axs[2, 0].scatter(x, vx, s=0.5, c=colors, edgecolor='')
                        axs[2, 0].set_ylabel(r'$v_{x}$ (km/s)', fontsize=8)
                        axs[2, 0].tick_params(labelsize='xx-small')
                        
                        xvy = np.vstack([x,vy])
                        z30 = gaussian_kde(xvy)(xvy)
                        idx = z30.argsort()
                        x, vy, colors = x[idx], vy[idx], z30[idx]                       
                        axs[3, 0].scatter(x, vy, s=0.5, c=colors, edgecolor='')
                        axs[3, 0].set_ylabel(r'$v_{y}$ (km/s)', fontsize=8)
                        axs[3, 0].tick_params(labelsize='xx-small')
                        
                        xvy = np.vstack([x,vz])
                        z40 = gaussian_kde(xvy)(xvy)
                        idx = z40.argsort()
                        x, vz, colors = x[idx], vz[idx], z40[idx]  
                        axs[4, 0].scatter(x, vz, s=0.5, c=colors, edgecolor='')
                        axs[4, 0].set_xlabel(r'$x$ (kpc)', fontsize=8)
                        axs[4, 0].set_ylabel(r'$v_{z}$ (km/s)', fontsize=8)
                        axs[4, 0].tick_params(labelsize='xx-small')
                        
                        #second column
                        
                        axs[0, 1].axis('off')
                        
                        yz = np.vstack([y,z])
                        z11 = gaussian_kde(yz)(yz)
                        idx = z11.argsort()
                        y, z, colors = y[idx], z[idx], z11[idx]  
                        axs[1, 1].scatter(y, z, s=0.5, c=colors, edgecolor='')
                        axs[1, 1].tick_params(labelsize='xx-small')
                        axs[1, 1].set_title(r'$y = %.03f \pm %.03f$ kpc'%(np.median(y), dy), fontsize=3)

                        
                        yvx = np.vstack([y,vx])
                        z21 = gaussian_kde(yvx)(yvx)
                        idx = z21.argsort()
                        y, vx, colors = y[idx], vx[idx], z21[idx]  
                        axs[2, 1].scatter(y, vx, s=0.5, c=colors, edgecolor='')
                        axs[2, 1].tick_params(labelsize='xx-small')
                        
                        yvy = np.vstack([y,vy])
                        z31 = gaussian_kde(yvy)(yvy)
                        idx = z31.argsort()
                        y, vy, colors = y[idx], vy[idx], z31[idx] 
                        axs[3, 1].scatter(y, vy, s=0.5, c=colors, edgecolor='')
                        axs[3, 1].tick_params(labelsize='xx-small')
                        
                        yvz = np.vstack([y,vx])
                        z41 = gaussian_kde(yvz)(yvz)
                        idx = z41.argsort()
                        y, vz, colors = y[idx], vz[idx], z41[idx] 
                        axs[4, 1].scatter(y, vz, s=0.5, c=colors, edgecolor='')
                        axs[4, 1].set_xlabel(r'$y$ (kpc)', fontsize=8)
                        axs[4, 1].tick_params(labelsize='xx-small')
                        
                        #third column
                        
                        axs[0, 2].axis('off')
                        axs[1, 2].axis('off')
                        
                        zvx = np.vstack([z,vx])
                        z22 = gaussian_kde(zvx)(zvx)
                        idx = z22.argsort()
                        z, vx, colors = z[idx], vx[idx], z22[idx] 
                        axs[2, 2].scatter(z, vx, s=0.5, c=colors, edgecolor='')
                        axs[2, 2].tick_params(labelsize='xx-small')
                        axs[2, 2].set_title(r'$z = %.03f \pm %.03f$ kpc'%(np.median(z), dz), fontsize=3)
                        
                        zvy = np.vstack([z,vy])
                        z32 = gaussian_kde(zvy)(zvy)
                        idx = z32.argsort()
                        z, vy, colors = z[idx], vy[idx], z32[idx] 
                        axs[3, 2].scatter(z, vy, s=0.5, c=colors, edgecolor='')
                        #axs[3, 2].set_xlabel=r'$z$', ylabel=r'$v_{y}$')
                        axs[3, 2].tick_params(labelsize='xx-small')
                        
                        zvz = np.vstack([z,vz])
                        z42 = gaussian_kde(zvz)(zvz)
                        idx = z42.argsort()
                        z, vz, colors = z[idx], vz[idx], z42[idx] 
                        axs[4, 2].scatter(z, vz, s=0.5, c=colors, edgecolor='')
                        axs[4, 2].set_xlabel(r'$z$ (kpc)', fontsize=8)
                        axs[4, 2].tick_params(labelsize='xx-small')
                        
                        #fourth column
                        
                        axs[0, 3].axis('off')
                        axs[1, 3].axis('off')                        
                        axs[2, 3].axis('off')
                        
                        vxvy = np.vstack([vx,vy])
                        z33 = gaussian_kde(vxvy)(vxvy)
                        idx = z33.argsort()
                        vx, vy, colors = vx[idx], vy[idx], z33[idx] 
                        axs[3, 3].scatter(vx, vy, s=0.5, c=colors, edgecolor='')
                        axs[3, 3].set_xlabel(r'$v_{x}$ (km/s)', fontsize=8)
                        axs[3, 3].set_ylabel(r'$v_{y}$ (km/s)', fontsize=8)
                        axs[3, 3].tick_params(labelsize='xx-small')
                        axs[3, 3].set_title(r'$v_{x} = %.03f \pm %.03f$ km/s'%(np.median(vx), dvx), fontsize=3)
   
                        
                        vxvz = np.vstack([vx,vz])
                        z43 = gaussian_kde(vxvz)(vxvz)
                        idx = z43.argsort()
                        vx, vz, colors = vx[idx], vz[idx], z43[idx]
                        axs[4, 3].scatter(vx, vz, s=0.5, c=colors, edgecolor='')
                        axs[4, 3].set_ylabel(r'$v_{z}$ (km/s)', fontsize=8)
                        axs[4, 3].tick_params(labelsize='xx-small')
                        
                        #fifth column
                        
                        axs[0, 4].axis('off')
                        axs[1, 4].axis('off')
                        axs[2, 4].axis('off')
                        axs[3, 4].axis('off')
                        
                        vyvz = np.vstack([vy,vz])
                        z44 = gaussian_kde(vyvz)(vyvz)
                        idx = z44.argsort()
                        vx, vz, colors = vy[idx], vz[idx], z44[idx]
                        axs[4, 4].scatter(vy, vz, s=0.5, c=colors, edgecolor='')
                        axs[4, 4].set_xlabel(r'$v_{y}$ (km/s)', fontsize=8) 
                        axs[4, 4].tick_params(labelsize='xx-small')
                        axs[4, 4].set_title(r'$v_{y} = %.03f \pm %.03f$ km/s'%(np.median(vy), dvy), fontsize=3)
                            
                        # Hide x labels and tick labels for top plots and y ticks for right plots.
                        for ax in axs.flat:
                            ax.label_outer()

                        #axs[4, 4].set_ylabel(r'$v_{z} = %.03f \pm %.03f$ km/s'%(np.median(vz), dvz), fontsize=3)
                        #axs[4, 4].yaxis.set_label_position("right")
                        
                        #fig.align_ylabels(axs[:, 0])
                        fig.suptitle(r'Cluster %i ($D_{\mathrm{naive}, \mathrm{pca}} = %i, %i$), $t_{\mathrm{simulation}}$ = %.02f Myr'%(1, 2, k, t), fontsize=12)
                        plt.savefig('phase_space_map_frame_%s_Norbiters_%s_Index_%i.pdf'%(str(i).rjust(5, '0'), str(Norbiters), k))
                        plt.close()
                
                
                plt.figure()
                
                x_total = phase_space_data[:,0]
                z_total = phase_space_data[:,2]
                
                xtotztot = np.vstack([x_total,z_total])
                colortot = gaussian_kde(xtotztot)(xtotztot)
                
                idx = colortot.argsort()
                x_tot, z_tot, colors = x_total[idx], z_total[idx], colortot[idx]
                
                plt.scatter(x_tot, z_tot, s=1, c=colors)
                plt.xlim(-2.0, 2.0)
                plt.ylim(-1.0, 1.0)
                plt.xlabel(r'$x$ (kpc)', fontsize=12)
                plt.ylabel(r'$z$ (kpc)', fontsize=12)
                plt.title('All Clusters, $t_{\mathrm{simulation}}$ = %.02f Myr'%(t), fontsize=16)
                plt.savefig('snapshot_frame__%s_Norbiters_%s_all.pdf'%(str(i).rjust(5, '0'), str(Norbiters)))
                plt.close()
        
    return 0

if __name__ in '__main__':
    
    maps()
