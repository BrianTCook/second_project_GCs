#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 21:57:22 2020

@author: BrianTCook
"""

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def plotting_things(orbiter_names, code_names, tend, dt):
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr), dt.value_in(units.Myr))
    npanels_x = len(orbiter_names)
    
    '''
    3 by 1, a plot for each orbiter (single star, etc.)
    each panel contains 3 curves, one for each gravity solver
    '''
    
    #energies
    
    plt.rc('font', family='serif')
    
    fig, axs = plt.subplots(1, len(orbiter_names))

    for i, orbiter_name in enumerate(orbiter_names): 
        
        #axs[i].set_ylim(0.8, 1.2)
        
        #if orbiter_name == 'SingleCluster':
                
        axs[i].set_xlabel('Simulation Time (Myr)', fontsize=12)
        
        if orbiter_name == 'SingleStar':
        
            axs[i].set_ylabel(r'$\Delta E/E(t=0)$', fontsize=12)
        
        axs[i].set_title(orbiter_name, fontsize=8)
        
        for code_name in code_names:

            try:
                
                energies = np.loadtxt(code_name + '_' + orbiter_name + '_energies.txt')
                
                e0 = energies[0]
                print('e0 is: %.04e joules'%e0)
                scaled_energies = [ e/e0 - 1. for e in energies ]                
                axs[i].plot(sim_times_unitless, scaled_energies, label=code_name)
                
            except:
                
                print('%s, %s could not be found'%(orbiter_name, code_name))
            
        axs[i].legend(loc='upper right')
   
    plt.tight_layout()         
    plt.savefig('testing_nemesis_energy.pdf')
    plt.close()
    
    #median radial coordinates
    
    fig, axs = plt.subplots(1, len(orbiter_names))

    for i, orbiter_name in enumerate(orbiter_names): 
        
        #if orbiter_name == 'SingleCluster':
                
        axs[i].set_xlabel('Simulation Time (Myr)', fontsize=12)
        
        if orbiter_name == 'SingleStar':
        
            axs[i].set_ylabel('Median Radial Coordinate (kpc)', fontsize=12)
        
        axs[i].set_title(orbiter_name, fontsize=8)
        
        for code_name in code_names:
            
            try:
                median_radial_coords = np.loadtxt(code_name + '_' + orbiter_name + '_median_radial_coords.txt')
                axs[i].plot(sim_times_unitless, median_radial_coords, label=code_name)
                
            except:
                print('%s, %s could not be found'%(orbiter_name, code_name))
            
        axs[i].legend(loc='upper right')
            
    plt.tight_layout()
    plt.savefig('testing_nemesis_radialcoords.pdf')
    plt.close()
    
    #median speeds
    
    fig, axs = plt.subplots(1, len(orbiter_names))

    for i, orbiter_name in enumerate(orbiter_names): 
        
        #if orbiter_name == 'SingleCluster':
                
        axs[i].set_xlabel('Simulation Time (Myr)', fontsize=12)
        
        if orbiter_name == 'SingleStar':
        
            axs[i].set_ylabel('Median Speed (km/s)', fontsize=12)
        
        axs[i].set_title(orbiter_name, fontsize=8)
        
        for code_name in code_names:
            
            try:
                median_speeds = np.loadtxt(code_name + '_' + orbiter_name + '_median_speeds.txt')
                axs[i].plot(sim_times_unitless, median_speeds, label=code_name)                    
            
            except:
                print('%s, %s could not be found'%(orbiter_name, code_name))
            
        axs[i].legend(loc='upper right')
       
    plt.tight_layout() 
    plt.savefig('testing_nemesis_speeds.pdf')
    plt.close()
    
    #clock times
    
    fig, axs = plt.subplots(1, len(orbiter_names))

    for i, orbiter_name in enumerate(orbiter_names): 
        
        #if orbiter_name == 'SingleCluster':
                
        axs[i].set_xlabel('Simulation Time (Myr)', fontsize=12)
            
        if orbiter_name == 'SingleStar':
            
            axs[i].set_ylabel('Clock Time (s)', fontsize=12)
        
        axs[i].set_title(orbiter_name, fontsize=8)
        
        for code_name in code_names:
            
            try:
                clock_times = np.loadtxt(code_name + '_' + orbiter_name + '_clock_times.txt')
                axs[i].semilogy(sim_times_unitless, clock_times, label=code_name)
                axs[i].set_ylim(1e-1, 5e4) #1/10th of a second to ~15 hours
                
            except:
                print('%s, %s could not be found'%(orbiter_name, code_name))
            
        axs[i].legend(loc='upper right')
       
    plt.tight_layout() 
    plt.savefig('testing_nemesis_clocktimes.pdf')
    plt.close()
    
    #center of mass
    
    fig, axs = plt.subplots(1, len(orbiter_names))

    for i, orbiter_name in enumerate(orbiter_names): 
                
        axs[i].set_xlabel('x (kpc)', fontsize=12)
        axs[i].set_ylabel('y (kpc)', fontsize=12)
        axs[i].set_title(orbiter_name, fontsize=8)
        
        for code_name in code_names:
            
            COMs = np.load('COM_data_%s_%s.npy'%(code_name, orbiter_name))
            Norbiters = len(COMs[0, 0, :])
            
            for j in range(Norbiters): #Norbiters
                
                xvals, yvals = COMs[:, 0, j], COMs[:, 1, j]
                axs[i].plot(xvals, yvals, linewidth=1) #label='orbiter %i, %s'%(j, code_name))
                    
        #axs[i].legend(loc='upper right', fontsize=8)
       
    plt.tight_layout() 
    plt.savefig('testing_nemesis_COMs.pdf')
    plt.close()
    
    return 0