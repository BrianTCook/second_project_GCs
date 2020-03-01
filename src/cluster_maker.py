#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:05:36 2020

@author: BrianTCook
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 13:29:47 2020

@author: BrianTCook
"""

from amuse.lab import *
from amuse.ic.kingmodel import new_king_model
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.community.mercury.interface import Mercury

from galpy.df import quasiisothermaldf
from galpy.potential import MWPotential2014
from galpy.util import bovy_conversion
from galpy.actionAngle import actionAngleStaeckel

import numpy as np
np.random.seed(73)

from nemesis import Nemesis, HierarchicalParticles
from nemesis_supplement import getxv, parent_worker, sub_worker, py_worker, smaller_nbody_power_of_two, distance_function, radius

def make_king_model_cluster(Rcoord, Zcoord, phicoord, vr_init, vphi_init, vz_init, 
                            W0, Mcluster, star_masses, Mgalaxy, Rgalaxy, code_name, parameters=[]):

    '''
    sets up a cluster with mass M and radius R
    which nbodycode would you like to use?
    '''
            
    converter = nbody_system.nbody_to_si(star_masses.sum(), 5|units.parsec)
    bodies = new_king_model(len(star_masses), W0, convert_nbody=converter)
    bodies.mass = star_masses

    '''
    takes in R, Z value
    returns VR, Vphi, VZ values
    should get appropriate 6D initial phase space conditions
    '''
    
    #convert from galpy/cylindrical to AMUSE/Cartesian units
    x_init = Rcoord*np.cos(phicoord) | units.kpc
    y_init = Rcoord*np.sin(phicoord) | units.kpc
    z_init = Zcoord | units.kpc
    
    #vphi = R \dot{\phi}? assuming yes for now
    vx_init = (vr_init*np.cos(phicoord) - vphi_init*np.sin(phicoord)) | units.kms
    vy_init = (vr_init*np.sin(phicoord) + vphi_init*np.cos(phicoord)) | units.kms
    vz_init = vz_init | units.kms
    
    pos_vec, vel_vec = (x_init, y_init, z_init), (vx_init, vy_init, vz_init)
    
    #initializing phase space coordinates
    bodies.x += x_init
    bodies.y += y_init
    bodies.z += z_init
    bodies.vx += vx_init
    bodies.vy += vy_init
    bodies.vz += vz_init
    
    #sub_worker in Nemesis
    if code_name == 'Nbody':
        
        code = Hermite(converter)
        code.particles.add_particles(bodies)
        code.commit_particles()
        
    if code_name == 'tree' or code_name == 'nemesis':
    
        code = BHTree(converter)
        code.particles.add_particles(bodies)
        code.commit_particles()
       
    return bodies, code
        
    '''
        
    if code_name == 'nemesis':

        ends up not being used?
        
        #uses a galpy function to evaluate the enclosed mass
        Mgalaxy, Rgalaxy = float(6.8e10)|units.MSun, 2.6|units.kpc #disk mass for MWPotential2014, Bovy(2015)
        converter_parent = nbody_system.nbody_to_si(Mgalaxy, Rgalaxy)
        dt = smaller_nbody_power_of_two(0.1 | units.Myr, converter_parent)
        dt_nemesis = dt
        dt_bridge = 0.01 * dt
        dt_param = 0.1
        
        nemesis = Nemesis(parent_worker(bodies), sub_worker(), py_worker())
        nemesis.timestep = dt
        nemesis.distfunc = distance_function
        nemesis.threshold = dt_nemesis
        nemesis.radius = radius
        
        nemesis.commit_parameters()
        #nemesis.particles.assign_subsystem(bodies, HierarchicalParticles(bodies)[0])
        print('nemesis.particles are', nemesis.particles)
        nemesis.commit_particles()
        
        code = nemesis
        for name,value in parameters:
            setattr(code.parameters, name, value)

        return bodies, _

    '''

def star_cluster(rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, index, Mgalaxy, Rgalaxy, code_name):
    
    '''
    takes 3 random numbers and generates open cluster
    with appropriate ICs in 6D phase space
    '''
    
    #limit to within 1 kpc of the galactic center
    Rcoord, phicoord, Zcoord = rvals[index], phivals[index], zvals[index]
    vr_init, vphi_init, vz_init = vrvals[index], vphivals[index], vzvals[index]
    
    data_directory = '/home/brian/Desktop/second_project_gcs/data/'
    star_masses = np.loadtxt(data_directory+'star_masses/star_masses_index=%i.txt'%index)
    Mcluster = np.sum( star_masses )
    
    star_masses, Mcluster = star_masses|units.MSun, Mcluster|units.MSun
    
    W0 = 1.5
    
    #just for converter, right?
    converter_sub = nbody_system.nbody_to_si(Mcluster, 10|units.parsec)
    
    print('~~~~~~~~~~~~~~~~')
    print('index is: %i'%(index))
    print('r, phi, z: %.04f kpc, %.04f radians, %.04f kpc'%(Rcoord, phicoord, Zcoord))
    print('vr, vphi, vz: %.04f km/s, %.04f km/s, %.04f km/s'%(vr_init, vphi_init, vz_init))
    print('~~~~~~~~~~~~~~~~')
    
    bodies, code = make_king_model_cluster(Rcoord, Zcoord, phicoord, vr_init, vphi_init, vz_init, 
                                           W0, Mcluster, star_masses, Mgalaxy, Rgalaxy, code_name, parameters=[])
    
    return bodies, code, converter_sub

def orbiter(code_name, orbiter_name, Mgalaxy, Rgalaxy, sepBinary, 
            rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, index):

    data_directory = '/home/brian/Desktop/second_project_gcs/data/'
    star_masses = np.loadtxt(data_directory+'/star_masses/star_masses_index=%i.txt'%index)
    
    converter_parent = nbody_system.nbody_to_si(Mgalaxy, Rgalaxy)
    converter_sub = nbody_system.nbody_to_si(np.median(star_masses)|units.MSun, 5.|units.parsec) #masses list is in solar mass units
    
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
           
        '''
            
        does not make sense to have nemesis, SingleStar
        
        if code_name == 'nemesis':
            
            dt = smaller_nbody_power_of_two(0.1 | units.Myr, converter_parent)
            dt_nemesis = dt
            dt_bridge = 0.01 * dt
            dt_param = 0.1
            
            nemesis = Nemesis(parent_worker(bodies), sub_worker(), py_worker())
            nemesis.timestep = dt
            nemesis.distfunc = distance_function
            nemesis.threshold = dt_nemesis
            nemesis.radius = radius
            
            nemesis.commit_parameters()
            #nemesis.particles.assign_subsystem(bodies, HierarchicalParticles(bodies)[0])
            print('nemesis.particles are', nemesis.particles)
            nemesis.commit_particles()
            
            code = nemesis
            
        '''
        
        return bodies, code
        
    if orbiter_name == 'SingleCluster':
        
        bodies, code, _ = star_cluster(rvals, phivals, zvals, vrvals, vphivals, vzvals, 
                                       masses, index, Mgalaxy, Rgalaxy, code_name)
        
        return bodies, code
