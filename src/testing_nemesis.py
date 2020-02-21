#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:05:09 2020

@author: BrianTCook
"""

from amuse.lab import *
#from amuse.ext.bridge import bridge
from amuse.couple import bridge

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from galpy.df import quasiisothermaldf
from galpy.potential import MWPotential2014, to_amuse
from galpy.util import bovy_conversion
from galpy.actionAngle import actionAngleStaeckel

from phase_space_mapping import maps
from cluster_maker import star_cluster

import random
import numpy as np
import time
import os

random.seed(73)

from nemesis import Nemesis, HierarchicalParticles
from nemesis_supplement import getxv, parent_worker, sub_worker, py_worker, smaller_nbody_power_of_two, distance_function, radius

#Circumvent a problem with using too many threads on OpenMPI
os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "yes"

def orbiter(orbiter_name, code_name, Mgalaxy, Rgalaxy, sepBinary, 
            rvals, phivals, zvals, masses, index):

    converter_parent = nbody_system.nbody_to_si(Mgalaxy, Rgalaxy)
    converter_sub = nbody_system.nbody_to_si(np.median(masses)|units.MSun, 5.|units.parsec) #masses list is in solar mass units
    
    '''
    takes in R, Z value
    returns VR, Vphi, VZ values
    should get appropriate 6D initial phase space conditions
    '''
    
    Rcoord, phicoord, Zcoord = rvals[index], phivals[index], zvals[index]
    vr_init, vphi_init, vz_init = vrvals[index], vphivals[index], vzvals[index]
    
    #convert from galpy/cylindrical to AMUSE/Cartesian units
    x_init = Rcoord*np.cos(phicoord) | units.kpc
    y_init = Rcoord*np.sin(phicoord) | units.kpc
    z_init = Zcoord | units.kpc
    
    #vphi = R \dot{\phi}? assuming yes for now
    vx_init = (vr_init*np.cos(phicoord) - vphi_init*np.sin(phicoord)) | units.kms
    vy_init = (vr_init*np.sin(phicoord) + vphi_init*np.cos(phicoord)) | units.kms
    vz_init = vz_init | units.kms
    
    if orbiter_name == 'SingleStar':
        
        bodies = Particles(1)
        
        bodies[0].mass = 1|units.MSun
        
        #right place in phase space
        bodies[0].x = x_init
        bodies[0].y = y_init
        bodies[0].z = z_init
        bodies[0].vx = vx_init
        bodies[0].vy = vy_init
        bodies[0].vz = vz_init
        
        #sub_worker in Nemesis, should not matter for SingleStar
        if code_name == 'Nbody':
            
            code = Hermite(converter_sub) #Mercury acts weird for SingleStar
            code.particles.add_particles(bodies)
            code.commit_particles()
            
        if code_name == 'tree':
        
            code = BHTree(converter_sub)
            code.particles.add_particles(bodies)
            code.commit_particles()
            
        if code_name == 'nemesis':
            
            parts = HierarchicalParticles(bodies)
            
            dt = smaller_nbody_power_of_two(0.1 | units.Myr, converter_parent)
            dt_nemesis = dt
            dt_bridge = 0.01 * dt
            dt_param = 0.1
            
            print
            
            nemesis = Nemesis( parent_worker, sub_worker, py_worker)
            nemesis.timestep = dt
            nemesis.distfunc = distance_function
            nemesis.threshold = dt_nemesis
            nemesis.radius = radius
            
            nemesis.commit_parameters()
            nemesis.particles.add_particles(parts)
            nemesis.commit_particles()
            
            code = nemesis
        
        return bodies, code
        
    if orbiter_name == 'SingleCluster':
        
        bodies, code, _ = star_cluster(rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, index, code_name)
        
        return bodies, code

def gravity_code_setup(orbiter_name, code_name, Mgalaxy, Rgalaxy, galaxy_code, sepBinary, 
                       rvals, phivals, zvals, vrvals, vphivals, vzvals, masses):
    
    '''
    will need to ask SPZ if he meant for field, orbiter to be separate in non
    Nemesis gravity solvers?
    '''
    
    Norbiters = len(rvals)
    
    converter_parent = nbody_system.nbody_to_si(Mgalaxy, Rgalaxy)
    converter_sub = nbody_system.nbody_to_si(np.median(masses)|units.MSun, 5.|units.parsec) #masses list is in solar mass units
    
    list_of_orbiters = [ orbiter(orbiter_name, code_name, Mgalaxy, Rgalaxy, sepBinary,
                                     rvals, phivals, zvals, masses, i) for i in range(Norbiters) ]
    
    orbiter_bodies_list = [ list_of_orbiters[i][0] for i in range(Norbiters) ] 
    orbiter_codes_list = [ list_of_orbiters[i][1] for i in range(Norbiters) ]
    
    cluster_colors = []
    
    for i, orbiter_code in enumerate(orbiter_codes_list):   

        stars = orbiter_code.particles.copy()
        cluster_color = np.random.rand(3,)
        
        cluster_colors.append([cluster_color]*len(stars))

        channel = stars.new_channel_to(orbiter_code.particles)
        channel.copy_attributes(['x','y','z','vx','vy','vz'])

    cluster_colors = [ j for i in cluster_colors for j in i ]

    if code_name != 'nemesis':
        
        gravity = bridge.Bridge(use_threading=False)

        for i in range(Norbiters):
            
            cluster_code = orbiter_codes_list[i]
            other_clusters = orbiter_codes_list[:i] + orbiter_codes_list[i+1:]
            other_things = tuple(other_clusters) + (galaxy_code,)

            #bridges each cluster with the bulge, not the other way around though
            gravity.add_system(cluster_code, other_things)    
            
        return gravity.particles, gravity, orbiter_bodies_list
            
    if code_name == 'nemesis':
        
        all_bodies = Particles(0)
        
        for i in range(Norbiters):
            all_bodies.add_particles(orbiter_bodies_list[i])
            
        #don't use orbiter_codes_list
        nemesis_parts = HierarchicalParticles(all_bodies)
        
        '''
        need add_subsystem and assign_subsystem in HierarchicalParticles I think
        '''
        
        dt = smaller_nbody_power_of_two(1. | units.Myr, converter_parent)
        dt_nemesis = dt
        print('dt_nemesis: ', dt.in_(units.Myr))
        dt_bridge = 0.01*dt
        dt_param = 0.1
        
        nemesis = Nemesis( parent_worker, sub_worker, py_worker)
        nemesis.timestep = dt
        nemesis.distfunc = distance_function
        nemesis.threshold = dt_nemesis
        nemesis.radius = radius
        
        nemesis.commit_parameters()
        nemesis.particles.add_particles(nemesis_parts)
        nemesis.commit_particles()

        print('nemesis.particles.compound_particles: ', nemesis.particles.compound_particles)        
        print('nemesis.subcodes are: ', nemesis.subcodes)

        gravity = bridge.Bridge(use_threading=False)
        gravity.add_system(nemesis, (galaxy_code,))
        gravity.timestep = dt_bridge
    
        return all_bodies, gravity, orbiter_bodies_list, cluster_colors

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
                    
    
    orbiter_colors = [ [np.random.random(3,)]*len(orbiter_bodies) for orbiter_bodies in orbiter_bodies_list]
    orbiter_colors = [ j for i in orbiter_colors for j in i ] #concatenate the list above
    
    Ntotal = len(gravity.particles)
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr), dt.value_in(units.Myr))
    sim_times = [ t|units.Myr for t in sim_times_unitless ]
    
    energies, median_radial_coords, median_speeds, clock_times = [], [], [], []
    
    t0 = time.time()
    
    #create an R^3 matrix to house phase space data for all particles
    phase_space_data = np.zeros((len(sim_times), 6, len(simulation_bodies)))
    
    for j, t in enumerate(sim_times):
        
        clock_times.append(time.time()-t0) #will be in seconds
        
        energy = gravity.kinetic_energy.value_in(units.J) + gravity.potential_energy.value_in(units.J)
        energies.append(energy)
        
        x = [ xx.value_in(units.kpc) for xx in gravity.particles.x ]
        y = [ yy.value_in(units.kpc) for yy in gravity.particles.y ]
        z = [ zz.value_in(units.kpc) for zz in gravity.particles.z ]
        
        vx = [ vxx.value_in(units.kms) for vxx in gravity.particles.vx ]
        vy = [ vyy.value_in(units.kms) for vyy in gravity.particles.vy ]
        vz = [ vzz.value_in(units.kms) for vzz in gravity.particles.vz ]
        
        for k, star in enumerate(gravity.particles):
            
            phase_space_data[j,0,k] = x[k]
            phase_space_data[j,1,k] = y[k]
            phase_space_data[j,2,k] = z[k]
            phase_space_data[j,3,k] = vx[k]
            phase_space_data[j,4,k] = vy[k]
            phase_space_data[j,5,k] = vz[k]
        
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
    np.savetxt(code_name + '_' + orbiter_name + '_colors.txt', cluster_colors)
    np.savetxt(code_name + '_' + orbiter_name + '_energies.txt', energies)
    np.savetxt(code_name + '_' + orbiter_name + '_median_radial_coords.txt', median_radial_coords)
    np.savetxt(code_name + '_' + orbiter_name + '_median_speeds.txt', median_speeds)
    np.savetxt(code_name + '_' + orbiter_name + '_clock_times.txt', clock_times)
    
    return 0

def plotting_things(orbiter_names, code_names, tend, dt):
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr), dt.value_in(units.Myr))
    npanels_x = len(orbiter_names)
    
    '''
    3 by 1, a plot for each orbiter (single star, etc.)
    each panel contains 3 curves, one for each gravity solver
    '''
    
    #energies
    
    plt.rc('font', family='serif')
    
    fig, axs = plt.subplots(1, 3)

    for i, orbiter_name in enumerate(orbiter_names): 
        
        #axs[i].set_ylim(0.8, 1.2)
        
        if orbiter_name == 'SingleCluster':
                
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
                print('oh no!')
            
        axs[i].legend(loc='upper right')
   
    plt.tight_layout()         
    plt.savefig('testing_nemesis_energy.pdf')
    plt.close()
    
    #median radial coordinates
    
    fig, axs = plt.subplots(1, 3)

    for i, orbiter_name in enumerate(orbiter_names): 
        
        if orbiter_name == 'SingleCluster':
                
            axs[i].set_xlabel('Simulation Time (Myr)', fontsize=12)
        
        if orbiter_name == 'SingleStar':
        
            axs[i].set_ylabel('Median Radial Coordinate (kpc)', fontsize=12)
        
        axs[i].set_title(orbiter_name, fontsize=8)
        
        for code_name in code_names:
            
            try:
                median_radial_coords = np.loadtxt(code_name + '_' + orbiter_name + '_median_radial_coords.txt')
                axs[i].plot(sim_times_unitless, median_radial_coords, label=code_name)
            except:
                print('oh no!')
            
        axs[i].legend(loc='upper right')
            
    plt.tight_layout()
    plt.savefig('testing_nemesis_radialcoords.pdf')
    plt.close()
    
    #median speeds
    
    fig, axs = plt.subplots(1, 3)

    for i, orbiter_name in enumerate(orbiter_names): 
        
        if orbiter_name == 'SingleCluster':
                
            axs[i].set_xlabel('Simulation Time (Myr)', fontsize=12)
        
        if orbiter_name == 'SingleStar':
        
            axs[i].set_ylabel('Median Speed (km/s)', fontsize=12)
        
        axs[i].set_title(orbiter_name, fontsize=8)
        
        for code_name in code_names:
            
            try:
                median_speeds = np.loadtxt(code_name + '_' + orbiter_name + '_median_speeds.txt')
                axs[i].plot(sim_times_unitless, median_speeds, label=code_name)                    
            except:
                print('oh no!')
            
        axs[i].legend(loc='upper right')
       
    plt.tight_layout() 
    plt.savefig('testing_nemesis_speeds.pdf')
    plt.close()
    
    #clock times
    
    fig, axs = plt.subplots(1, 3)

    for i, orbiter_name in enumerate(orbiter_names): 
        
        if orbiter_name == 'SingleCluster':
                
            axs[i].set_xlabel('Simulation Time (Myr)', fontsize=12)
            
        if orbiter_name == 'SingleStar':
            
            axs[i].set_ylabel('Clock Time (s)', fontsize=12)
        
        axs[i].set_title(orbiter_name, fontsize=8)
        
        for code_name in code_names:
            
            try:
                clock_times = np.loadtxt(code_name + '_' + orbiter_name + '_clock_times.txt')
                axs[i].semilogy(sim_times_unitless, clock_times, label=code_name)
                axs[i].set_ylim(1e-1, 5e3) #1/10th of a second to ~1.5 hours
            except:
                print('oh no!')
            
        axs[i].legend(loc='upper right')
       
    plt.tight_layout() 
    plt.savefig('testing_nemesis_clocktimes.pdf')
    plt.close()
    
    return 0

if __name__ in '__main__':
    
    potential = MWPotential2014 #galpy
    
    sepBinary = 20.|units.parsec #not necessary if not doing binary cluster part
    tend, dt = 100.|units.Myr, 1.|units.Myr
    dt_param = 0.1 #for nemesis
    
    #uses a galpy function to evaluate the enclosed mass
    Mgalaxy, Rgalaxy = float(6.8e10)|units.MSun, 2.6|units.kpc #disk mass for MWPotential2014, Bovy(2015)
    
    Norbiters = 1
    
    rvals = np.loadtxt('/home/brian/Desktop/second_project_gcs/data/dehnen_rvals.txt')
    phivals = np.loadtxt('/home/brian/Desktop/second_project_gcs/data/dehnen_phivals.txt')
    zvals = np.loadtxt('/home/brian/Desktop/second_project_gcs/data/dehnen_zvals.txt')
    
    vrvals = np.loadtxt('/home/brian/Desktop/second_project_gcs/data/bovy_vrvals.txt')
    vphivals = np.loadtxt('/home/brian/Desktop/second_project_gcs/data/bovy_vphivals.txt')
    vzvals = np.loadtxt('/home/brian/Desktop/second_project_gcs/data/bovy_vzvals.txt')
    
    masses = np.loadtxt('/home/brian/Desktop/second_project_gcs/data/cluster_masses_for_sampling.txt')
    
    rvals = rvals[:Norbiters]
    phivals = phivals[:Norbiters]
    zvals = zvals[:Norbiters]  
    masses = masses[:Norbiters]

    orbiter_names = [ 'SingleStar', 'SingleCluster' ]
    code_names = ['tree', 'Nbody' ]#, 'nemesis'
    
    t0 = time.time()
    
    for orbiter_name in orbiter_names:
        for code_name in code_names:
            
            print('\\\\\\\\\\\\\\\\\\\\\\\\')
            print(orbiter_name, code_name)
            print('\\\\\\\\\\\\\\\\\\\\\\\\')
            
            simulation(orbiter_name, code_name, potential, Mgalaxy, Rgalaxy, 
                       sepBinary, rvals, phivals, zvals, vrvals, vphivals, vzvals, 
                       masses, tend, dt)
            
            print('current time: %.03f minutes'%((time.time()-t0)/60.))
            
            if orbiter_name == 'SingleCluster' or orbiter_name == 'SingleStar': #should be the same for all three and captures all gravity info
                
                maps(orbiter_name, code_name)
                continue
            
            print('current time: %.03f minutes'%((time.time()-t0)/60.))
            
    plotting_things(orbiter_names, code_names, tend, dt)
