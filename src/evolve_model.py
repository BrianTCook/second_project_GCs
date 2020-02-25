#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 21:56:36 2020

@author: BrianTCook
"""

from gravity_code import gravity_code_setup

def simulation(orbiter_name, code_name, potential, Mgalaxy, Rgalaxy, sepBinary, 
               rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, tend, dt):
    
    converter_parent = nbody_system.nbody_to_si(Mgalaxy, Rgalaxy)
    converter_sub = nbody_system.nbody_to_si(np.median(masses)|units.MSun, 5.|units.parsec) #masses list is in solar mass units
    
    galaxy_code = to_amuse(potential, t=0.0, tgalpy=0.0, reverse=False, ro=None, vo=None)
    
    #first thing is all particles in one superset, gravity is the code,
    #third thing is the list of orbiter bodies s.t. we can compute COMs independently
    #and plot them with different colors
    
    simulation_bodies, gravity, orbiter_bodies_list, cluster_colors = gravity_code_setup(orbiter_name, code_name, Mgalaxy, Rgalaxy, 
                                                                                         galaxy_code, sepBinary, rvals, phivals, zvals, 
                                                                                         vrvals, vphivals, vzvals, masses)
    
    Ntotal = len(gravity.particles)
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr), dt.value_in(units.Myr))
    sim_times = [ t|units.Myr for t in sim_times_unitless ]
    
    energies, median_radial_coords, median_speeds, clock_times = [], [], [], []
    
    t0 = time.time()
    
    #create an R^3 matrix to house phase space data for all particles
    phase_space_data = np.zeros((len(sim_times), 6, len(simulation_bodies)))
    
    body_masses = gravity.particles.mass
    total_mass = body_masses.sum()
    
    #COM data
    COM_data = np.zeros((len(sim_times), 2, Norbiters))
    
    if orbiter_name == 'SingleStar':
            cluster_populations = [1 for i in range(Norbiters) ]
    if orbiter_name == 'SingleCluster':
            cluster_populations = np.loadtxt('/home/brian/Desktop/second_project_gcs/data/Nstars_in_clusters.txt')
            cluster_populations = cluster_populations[:Norbiters]
    
    cluster_pop_flag = 0
    
    for j, t in enumerate(sim_times):
        
        clock_times.append(time.time()-t0) #will be in seconds
        
        energy = gravity.particles.kinetic_energy() + gravity.particles.potential_energy(G=constants.G)
        energies.append( energy.value_in(units.J) )
        
        x = [ xx.value_in(units.kpc) for xx in gravity.particles.x ]
        y = [ yy.value_in(units.kpc) for yy in gravity.particles.y ]    
        z = [ zz.value_in(units.kpc) for zz in gravity.particles.z ]
        
        vx = [ vxx.value_in(units.kms) for vxx in gravity.particles.vx ]
        vy = [ vyy.value_in(units.kms) for vyy in gravity.particles.vy ]
        vz = [ vzz.value_in(units.kms) for vzz in gravity.particles.vz ]
        
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
        
        for l, star in enumerate(gravity.particles):
            
            phase_space_data[j,0,l] = x[l]
            phase_space_data[j,1,l] = y[l]
            phase_space_data[j,2,l] = z[l]
            phase_space_data[j,3,l] = vx[l]
            phase_space_data[j,4,l] = vy[l]
            phase_space_data[j,5,l] = vz[l]
        
        rvals = [ np.sqrt(x[i]**2 + y[i]**2 + z[i]**2) for i in range(Ntotal) ]
        median_radial_coords.append(np.median(rvals))
        
        speeds = [ np.sqrt(vx[i]**2 + vy[i]**2 + vz[i]**2) for i in range(Ntotal) ]
        median_speeds.append(np.median(speeds))
                
        gravity.evolve_model(t)
        
    try:
        gravity.stop()
    except:
        'gravity cannot be stopped!'

    np.save('time_data_%s_%s.npy'%(code_name, orbiter_name), sim_times_unitless)
    np.save('sixD_data_%s_%s.npy'%(code_name, orbiter_name), phase_space_data)
    np.save('COM_data_%s_%s.npy'%(code_name, orbiter_name), COM_data)
        
    np.savetxt(code_name + '_' + orbiter_name + '_colors.txt', cluster_colors)
    np.savetxt(code_name + '_' + orbiter_name + '_energies.txt', energies)
    np.savetxt(code_name + '_' + orbiter_name + '_median_radial_coords.txt', median_radial_coords)
    np.savetxt(code_name + '_' + orbiter_name + '_median_speeds.txt', median_speeds)
    np.savetxt(code_name + '_' + orbiter_name + '_clock_times.txt', clock_times)
    
    return 0
