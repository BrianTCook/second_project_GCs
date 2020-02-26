#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 21:56:36 2020

@author: BrianTCook
"""

from amuse.lab import *
#from amuse.ext.bridge import bridge
from amuse.couple import bridge

from galpy.potential import MWPotential2014, to_amuse
from galpy.util import bovy_conversion

from gravity_code import gravity_code_setup

import numpy as np
import time

def print_diagnostics(time, simulation_bodies, E_dyn, dE_dyn):
    
    print('------------')
    print('time: ', time)
    print('simulation_bodies.center_of_mass(): ', simulation_bodies.center_of_mass())
    print('E_dyn: ', E_dyn)
    print('dE_dyn: ', dE_dyn)
    print('------------')

def simulation(orbiter_name, code_name, potential, Mgalaxy, Rgalaxy, sepBinary, 
               rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, Norbiters, tend, dt):
    
    converter_parent = nbody_system.nbody_to_si(Mgalaxy, Rgalaxy)
    converter_sub = nbody_system.nbody_to_si(np.median(masses)|units.MSun, 5.|units.parsec) #masses list is in solar mass units
    
    galaxy_code = to_amuse(potential, t=0.0, tgalpy=0.0, reverse=False, ro=None, vo=None)
    
    #first thing is all particles in one superset, gravity is the code,
    #third thing is the list of orbiter bodies s.t. we can compute COMs independently
    #and plot them with different colors
    
    simulation_bodies, gravity, orbiter_bodies_list, cluster_colors = gravity_code_setup(orbiter_name, code_name, Mgalaxy, Rgalaxy, 
                                                                                         galaxy_code, sepBinary, rvals, phivals, zvals, 
                                                                                         vrvals, vphivals, vzvals, masses)
    
    channel_from_gravity_to_framework = gravity.particles.new_channel_to(simulation_bodies)
    
    Ntotal = len(gravity.particles)
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr), dt.value_in(units.Myr))
    sim_times = [ t|units.Myr for t in sim_times_unitless ]
    
    np.savetxt('times_in_Myr_%s_%s_Norbiters=%i.txt'%(code_name, orbiter_name, Norbiters), sim_times_unitless)
    
    delta_energies, median_radial_coords, median_speeds, clock_times = [], [], [], []
    
    body_masses = gravity.particles.mass
    total_mass = body_masses.sum()
    
    if orbiter_name == 'SingleStar':
            cluster_populations = [1 for i in range(Norbiters) ]
    if orbiter_name == 'SingleCluster':
            cluster_populations = np.loadtxt('/home/brian/Desktop/second_project_gcs/data/Nstars_in_clusters.txt')
            cluster_populations = cluster_populations[:Norbiters]
    
    #COM data
    #COM_data = np.zeros((len(sim_times), 2, Norbiters))
    
    #cluster_pop_flag = 0
    
    #for saving in write_set_to_file
    filename = "data_%s_%s_Norbiters=%i.csv"%(code_name, orbiter_name, Norbiters)
    
    t0 = time.time()
    
    for j, t in enumerate(sim_times):
        
        clock_times.append(time.time()-t0) #will be in seconds
    
        if j == 0:
            E_dyn_init = gravity.kinetic_energy + gravity.potential_energy
            
        E_dyn = gravity.kinetic_energy + gravity.potential_energy
        dE_dyn = (E_dyn/E_dyn_init) - 1.
        
        delta_energies.append(dE_dyn)
        
        '''
        
        #create an R^3 matrix to house phase space data for all particles
        phase_space_data = np.zeros((len(sim_times), 6, len(simulation_bodies)))
        
        x = [ xx.value_in(units.kpc) for xx in gravity.particles.x ]
        y = [ yy.value_in(units.kpc) for yy in gravity.particles.y ]    
        z = [ zz.value_in(units.kpc) for zz in gravity.particles.z ]
        
        vx = [ vxx.value_in(units.kms) for vxx in gravity.particles.vx ]
        vy = [ vyy.value_in(units.kms) for vyy in gravity.particles.vy ]
        vz = [ vzz.value_in(units.kms) for vzz in gravity.particles.vz ]
        
        for l, star in enumerate(gravity.particles):
            
            phase_space_data[j,0,l] = x[l]
            phase_space_data[j,1,l] = y[l]
            phase_space_data[j,2,l] = z[l]
            phase_space_data[j,3,l] = vx[l]
            phase_space_data[j,4,l] = vy[l]
            phase_space_data[j,5,l] = vz[l]

        ##NEEDS TO BE ADJUSTED FOR WRITE_SET_TO_FILE##
        
        rvals = [ np.sqrt(x[i]**2 + y[i]**2 + z[i]**2) for i in range(Ntotal) ]
        median_radial_coords.append(np.median(rvals))
        
        speeds = [ np.sqrt(vx[i]**2 + vy[i]**2 + vz[i]**2) for i in range(Ntotal) ]
        median_speeds.append(np.median(speeds))

        for k, number_of_stars in enumerate(cluster_populations):
            
        starting_index = int(np.sum( cluster_populations[:k] ))
        ending_index = starting_index + int(number_of_stars)
        
        if cluster_pop_flag == 0:
            print('starting, ending indices are', starting_index, ending_index)
    
        #should not use total mass!
        
        cluster_masses = body_masses[starting_index:ending_index]
        cluster_total_mass = cluster_masses.sum()
    
        x_COM = np.sum( [ body_masses[i]*x[i]/cluster_total_mass for i in range(starting_index, ending_index) ] ) #in kpc
        y_COM = np.sum( [ body_masses[i]*y[i]/cluster_total_mass for i in range(starting_index, ending_index) ] ) #in kpc
    
        COM_data[j, 0, k] = x_COM
        COM_data[j, 1, k] = y_COM
            
        cluster_pop_flag = 1
        energies.append( energy.value_in(units.J) )
        
        np.save('sixD_data_%s_%s.npy'%(code_name, orbiter_name), phase_space_data)
        np.save('COM_data_%s_%s.npy'%(code_name, orbiter_name), COM_data)
            
        np.savetxt(code_name + '_' + orbiter_name + '_median_radial_coords.txt', median_radial_coords)
        np.savetxt(code_name + '_' + orbiter_name + '_median_speeds.txt', median_speeds)
        np.savetxt(code_name + '_' + orbiter_name + '_clock_times.txt', clock_times)
        
        '''
        
        gravity.evolve_model(t)
        channel_from_gravity_to_framework.copy()
        
        write_set_to_file(simulation_bodies, filename, "csv")
        print_diagnostics(t, simulation_bodies, E_dyn, dE_dyn)

    gravity.stop()
    
    np.savetxt(code_name + '_' + orbiter_name + '_colors_Norbiters=' + str(Norbiters) + '.txt.', cluster_colors)
    np.savetxt(code_name + '_' + orbiter_name + '_dE_Norbiters=' + str(Norbiters) + '.txt.', delta_energies)
    
    return 0
