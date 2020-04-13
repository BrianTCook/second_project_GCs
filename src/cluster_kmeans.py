#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:17:37 2020

@author: BrianTCook
"""

import numpy as np
import glob
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

def sklearn_mapper(input_list):
    
    input_output_dict = {}
    output_index = 0
    
    for element in input_list:
        
        if element not in input_output_dict.keys():
            
            input_output_dict.update( {element:output_index} )
            output_index += 1
    
    return input_output_dict
        

def get_kmeans_result(Norbiters, plot_6D, plot_3D):
    
    datadir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
    cluster_populations = np.loadtxt(datadir + 'Nstars_in_clusters.txt')
    cluster_populations = list(cluster_populations[:Norbiters])
    
    labels = []
    
    for logN in range(int(np.log2(Norbiters))+1):
        labels += [ logN for i in range(int(cluster_populations[logN])) ]
    
    data_filename = glob.glob(datadir+'enbid_files/*_00400_Norbiters_%i.ascii'%(Norbiters))
    data_6D = np.loadtxt(data_filename[0])
    data_3D = data_6D[:, 0:3]
    data_2D = data_6D[:, 1:3]

    #apply kmeans clustering
    kmeans_2D = KMeans(n_clusters=Norbiters)
    kmeans_2D.fit(data_2D)
    y_kmeans_2D = kmeans_2D.predict(data_2D)
    io_dict_2D = sklearn_mapper(y_kmeans_2D)
    
    kmeans_3D = KMeans(n_clusters=Norbiters)
    kmeans_3D.fit(data_3D)
    y_kmeans_3D = kmeans_3D.predict(data_3D)
    io_dict_3D = sklearn_mapper(y_kmeans_3D)
    
    kmeans_6D = KMeans(n_clusters=Norbiters)
    kmeans_6D.fit(data_6D)
    y_kmeans_6D = kmeans_6D.predict(data_6D)
    io_dict_6D = sklearn_mapper(y_kmeans_6D)

    y_compare_2D = [ io_dict_2D[y] for y in y_kmeans_2D ]
    y_compare_3D = [ io_dict_3D[y] for y in y_kmeans_3D ]
    y_compare_6D = [ io_dict_6D[y] for y in y_kmeans_6D ]
    
    hits_2D = [ 1 if y1 == y2 else 0 for y1, y2 in zip(labels, y_compare_2D) ]
    hits_3D = [ 1 if y1 == y2 else 0 for y1, y2 in zip(labels, y_compare_3D) ]
    hits_6D = [ 1 if y1 == y2 else 0 for y1, y2 in zip(labels, y_compare_6D) ]
    
    success_2D = np.sum(hits_2D)/len(hits_2D)
    success_3D = np.sum(hits_3D)/len(hits_3D)
    success_6D = np.sum(hits_6D)/len(hits_6D)
    
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
        plt.savefig('kmeans_Norbiters_%i_2D.pdf'%(Norbiters))
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
        plt.savefig('kmeans_Norbiters_%i_3D.pdf'%(Norbiters))
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
        plt.savefig('kmeans_Norbiters_%i_6D.pdf'%(Norbiters))
        plt.close()
    
    return success_2D, success_3D, success_6D

if __name__ in '__main__':
    
    logN_max = 7
    plot_2D, plot_3D, plot_6D = True, True, True

    for logN in range(logN_max):
        
        print('N = %i'%(2**logN))
        if plot_2D == True or plot_3D == True or plot_6D == True:
            
            plt.rc('text', usetex = True)
            plt.rc('font', family = 'serif')
        
        s_2, s_3, s_6 = get_kmeans_result(2**logN, plot_3D, plot_6D)
        print(s_2, s_3, s_6)
    
    print('hello world!')