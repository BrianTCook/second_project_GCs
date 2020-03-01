#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 21:55:03 2020

@author: BrianTCook
"""

from amuse.lab import *
#from amuse.ext.bridge import bridge
from amuse.couple import bridge

from galpy.df import quasiisothermaldf
from galpy.potential import MWPotential2014, to_amuse
from galpy.util import bovy_conversion
from galpy.actionAngle import actionAngleStaeckel

from cluster_maker import orbiter
from nemesis import *
from nemesis_supplement import *

import numpy as np
import time

def gravity_code_setup(code_name, orbiter_name, Mgalaxy, Rgalaxy, galaxy_code, sepBinary, 
                       rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, Norbiters):
    
    '''
    will need to ask SPZ if he meant for field, orbiter to be separate in non
    Nemesis gravity solvers?
    '''
    
    converter_parent = nbody_system.nbody_to_si(Mgalaxy, Rgalaxy)
    converter_sub = nbody_system.nbody_to_si(np.median(masses)|units.MSun, 5.|units.parsec) #masses list is in solar mass units
    
    list_of_orbiters = [ orbiter(code_name, orbiter_name, Mgalaxy, Rgalaxy, sepBinary,
                                     rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, i) for i in range(Norbiters) ]
    
    orbiter_bodies_list = [ list_of_orbiters[i][0] for i in range(Norbiters) ] 
    orbiter_codes_list = [ list_of_orbiters[i][1] for i in range(Norbiters) ]
    
    print(len(list_of_orbiters))
    
    cluster_colors = []
    
    for i, orbiter_code in enumerate(orbiter_codes_list):   

        stars = orbiter_code.particles.copy()
        cluster_color = np.random.rand(3,)
        
        cluster_colors.append([cluster_color]*len(stars))

        channel = stars.new_channel_to(orbiter_code.particles)
        channel.copy_attributes(['mass','x','y','z','vx','vy','vz'])

    cluster_colors = [ j for i in cluster_colors for j in i ] #concatenate the list

    if code_name != 'nemesis':
        
        gravity = bridge.Bridge(use_threading=False)

        for i in range(Norbiters):
            
            cluster_code = orbiter_codes_list[i]
            other_clusters = orbiter_codes_list[:i] + orbiter_codes_list[i+1:]
            other_things = tuple(other_clusters) + (galaxy_code,)

            #bridges each cluster with the bulge, not the other way around though
            gravity.add_system(cluster_code, other_things)    
            
    if code_name == 'nemesis':
        
        all_bodies = Particles(0)
        
        for i in range(Norbiters):
            all_bodies.add_particles(orbiter_bodies_list[i])
        
        '''
        need add_subsystem and assign_subsystem in HierarchicalParticles I think
        '''
                
        dt = smaller_nbody_power_of_two(0.1 | units.Myr, converter_parent)
        dt_nemesis = dt
        dt_bridge = 0.01 * dt
        dt_param = 0.1
        
        nemesis = Nemesis(parent_worker, sub_worker, py_worker, all_bodies)
        nemesis.timestep = dt
        nemesis.distfunc = distance_function
        nemesis.threshold = dt_nemesis
        nemesis.radius = radius
        
        nemesis.commit_parameters()
        nemesis.particles.assign_subsystem(all_bodies, HierarchicalParticles(all_bodies)[0])
        print('nemesis.particles are', nemesis.particles)
        nemesis.commit_particles()
        
        code = nemesis
        
        print('nemesis.particles.compound_particles: ', nemesis.particles.compound_particles)        
        print('nemesis.subcodes are: ', nemesis.subcodes)

        gravity = bridge.Bridge(use_threading=False)
        gravity.add_system(nemesis, (galaxy_code,))
        gravity.timestep = dt_bridge
    
    return gravity.particles, gravity, orbiter_bodies_list, cluster_colors
