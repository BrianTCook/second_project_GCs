#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 21:56:36 2020

@author: BrianTCook
"""

from __future__ import division

from amuse.lab import *
from amuse.couple import bridge
from amuse.support import io

from galpy.potential import MWPotential2014, to_amuse
from galpy.util import bovy_conversion

from gravity_code import gravity_code_setup
from fw_entropy import get_entropy

import numpy as np
import pandas as pd
import gzip
import time
import math

def print_diagnostics(time, simulation_bodies, E_dyn, dE_dyn):
    
    print('------------')
    print('time: ', time)
    print('simulation_bodies.center_of_mass(): ', simulation_bodies.center_of_mass().value_in(units.kpc))
    print('E_dyn: ', E_dyn)
    print('dE_dyn: ', dE_dyn)
    print('------------')

def simulation(code_name, orbiter_name, potential, Mgalaxy, Rgalaxy, sepBinary, 
               rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, Norbiters, tend, dt):
    
    converter_parent = nbody_system.nbody_to_si(Mgalaxy, Rgalaxy)
    converter_sub = nbody_system.nbody_to_si(np.median(masses)|units.MSun, 5.|units.parsec) #masses list is in solar mass units
    
    galaxy_code = to_amuse(potential, t=0.0, tgalpy=0.0, reverse=False, ro=None, vo=None)
    
    #first thing is all particles in one superset, gravity is the code,
    #third thing is the list of orbiter bodies s.t. we can compute COMs independently
    #and plot them with different colors
    
    simulation_bodies, gravity, orbiter_bodies_list, cluster_colors, _ = gravity_code_setup(code_name, orbiter_name, Mgalaxy, Rgalaxy, galaxy_code, sepBinary, 
                                                                                                  rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, Norbiters)

    channel_from_gravity_to_framework = gravity.particles.new_channel_to(simulation_bodies)    
    Ntotal = len(simulation_bodies)
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr)+dt.value_in(units.Myr), dt.value_in(units.Myr))
    sim_times = [ t|units.Myr for t in sim_times_unitless ]
    
    np.savetxt('times_in_Myr_%s_%s_Norbiters_%i.txt'%(code_name, orbiter_name, Norbiters), sim_times_unitless)
    
    delta_energies, clock_times = [], []
    
    body_masses = gravity.particles.mass
    total_mass = body_masses.sum()
    
    '''
    if orbiter_name == 'SingleStar':
            cluster_populations = [1 for i in range(Norbiters) ]
    '''
    
    if orbiter_name == 'SingleCluster':
        
            cluster_populations = np.loadtxt('/home/s1780638/second_project_gcs/data/Nstars_in_clusters.txt')
        
            #need to give clusters sorted by an attribute, in our case increasing |r|
            #new_index = indices_dict[old_index]
            indices_dict = sort_clusters_by_attribute('|r|')
            
            cluster_populations_sorted = [ cluster_populations[indices_dict[i]]
                                           for i in range(len(cluster_populations)) ]
    
    #for 3D numpy array storage
    all_data = np.zeros((len(sim_times), Ntotal, 6))
    mass_data = np.zeros((len(sim_times), Ntotal, 1))    
    COM_data = np.zeros((len(sim_times), Norbiters, 2))

    #for saving in write_set_to_file
    filename = 'data_temp.csv'
    attributes = ('mass', 'x', 'y', 'z', 'vx', 'vy', 'vz')
    
    print('len(sim_times) is', len(sim_times))
    gadget_flag = 1 #int(math.floor(len(sim_times)/20))
    
    t0 = time.time()
    all_data_index = 0
    
    for j, t in enumerate(sim_times):

        clock_times.append(time.time()-t0) #will be in seconds
    
        if j == 0:
            E_dyn_init = gravity.kinetic_energy + gravity.potential_energy
            
        E_dyn = gravity.kinetic_energy + gravity.potential_energy
        dE_dyn = (E_dyn/E_dyn_init) - 1.
        
        delta_energies.append(dE_dyn)
         
        if j%gadget_flag == 0:
            
            print_diagnostics(t, simulation_bodies, E_dyn, dE_dyn)
                        
            '''
            io.write_set_to_file(gravity.particles, filename, 'csv',
                                 attribute_types = (units.MSun, units.kpc, units.kpc, units.kpc, units.kms, units.kms, units.kms),
                                 attribute_names = attributes)
            
            data_t = pd.read_csv(filename, names=list(attributes))
            data_t = data_t.drop([0, 1, 2]) #removes labels units, and unit names
            
            masses = data_t['mass'].tolist()
            mass_data[all_data_index, :len(data_t.index), :] = masses #in solar masses
            
            data_t = data_t.drop(columns=['mass']) #goes from 7D --> 6D
            data_t = data_t.astype(float) #strings to floats

            all_data[all_data_index, :len(data_t.index), :] = data_t.values
            np.savetxt('enbid_%s_frame_%s_Norbiters_%s.ascii'%(code_name, str(j).rjust(5, '0'), str(Norbiters)), data_t.values)
            
            all_data_index += 1
            '''
        
        x, y = data_t['x'].tolist(), data_t['y'].tolist()
        
        #stuff to analyze COM of each star cluster
        for k, number_of_stars in enumerate(cluster_populations_sorted):
            
            starting_index = int(np.sum( cluster_populations_sorted[:k] ))
            ending_index = starting_index + int(number_of_stars)
            
            cluster_masses = body_masses[starting_index:ending_index]
            cluster_total_mass = cluster_masses.sum()
        
            x_COM = np.sum( [ body_masses[i]*x[i]/cluster_total_mass for i in range(starting_index, ending_index) ] ) #in kpc
            y_COM = np.sum( [ body_masses[i]*y[i]/cluster_total_mass for i in range(starting_index, ending_index) ] ) #in kpc
        
            COM_data[j, k, 0] = x_COM
            COM_data[j, k, 1] = y_COM

        gravity.evolve_model(t)
        channel_from_gravity_to_framework.copy()

    try:
        gravity.stop()
    except:
        print('gravity cannot be stopped, mwahahaha')
    
    #things that are not easily extracted from write_set_to_file
    
    #f_all = gzip.GzipFile('all_data_%s_%s_Norbiters_%s.npy.gz'%(code_name, orbiter_name, str(Norbiters)), 'w')
    #f_COM = gzip.GzipFile('COM_data_%s_%s_Norbiters_%s.npy.gz'%(code_name, orbiter_name, str(Norbiters)), 'w')
    #np.save(file=f_all, arr=all_data, allow_pickle=True)
    #np.save(file=f_COM, arr=COM_data, allow_pickle=True)
    
    #f_all.close()
    #f_COM.close()
  
    np.savetxt(code_name + '_' + orbiter_name + '_masses_Norbiters_' + str(Norbiters) + '.txt', mass_data)
    np.savetxt(code_name + '_' + orbiter_name + '_colors_Norbiters_' + str(Norbiters) + '.txt', cluster_colors)
    np.savetxt(code_name + '_' + orbiter_name + '_dE_Norbiters_' + str(Norbiters) + '.txt', delta_energies)
    np.savetxt(code_name + '_' + orbiter_name + '_clock_times.txt', clock_times)
    
    return 0
