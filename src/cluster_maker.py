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

from cluster_table import sort_clusters_by_attribute

import numpy as np
np.random.seed(73)

def make_king_model_cluster(r_init, phi_init, z_init, vr_init, vphi_init, vz_init, 
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
    x_init = r_init*np.cos(phi_init) | units.kpc
    y_init = r_init*np.sin(phi_init) | units.kpc
    z_init = z_init | units.kpc
    
    #vphi = R \dot{\phi}? assuming yes for now
    vx_init = (vr_init*np.cos(phi_init) - vphi_init*np.sin(phi_init)) | units.kms
    vy_init = (vr_init*np.sin(phi_init) + vphi_init*np.cos(phi_init)) | units.kms
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
    
        #doing nemesis solver here is not necessary
        
        code = BHTree(converter)
        code.particles.add_particles(bodies)
        code.commit_particles()
       
    return bodies, code

def star_cluster(rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, index, Mgalaxy, Rgalaxy, code_name):
    
    '''
    takes 3 random numbers and generates open cluster
    with appropriate ICs in 6D phase space
    '''
    
    #limit to within 1 kpc of the galactic center
    r_init, phi_init, z_init = rvals[index], phivals[index], zvals[index]
    vr_init, vphi_init, vz_init = vrvals[index], vphivals[index], vzvals[index]
    
    #data_directory = '/home/s1780638/second_project_gcs/data/'
    data_directory = 'home/brian/Desktop/second_project_gcs/data/'
    
    star_masses = np.loadtxt(data_directory+'star_masses/star_masses_index=%i.txt'%index)
    Mcluster = np.sum( star_masses )
    
    star_masses, Mcluster = star_masses|units.MSun, Mcluster|units.MSun
    
    W0 = 1.5
    
    #just for converter, right?
    converter_sub = nbody_system.nbody_to_si(Mcluster, 10|units.parsec)
    
    bodies, code = make_king_model_cluster(r_init, phi_init, z_init, vr_init, vphi_init, vz_init, 
                                           W0, Mcluster, star_masses, Mgalaxy, Rgalaxy, code_name, parameters=[])
    
    return bodies, code, converter_sub

def orbiter(code_name, orbiter_name, Mgalaxy, Rgalaxy, sepBinary, 
            rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, index):

    #need to give clusters sorted by an attribute, in our case increasing |r|
    #sorting needs to happen once, either here or in star_cluster function
    indices_dict = sort_clusters_by_attribute('|r|')
    index_sorted = indices_dict[index]
    
    #data_directory = '/home/s1780638/second_project_gcs/data/'
    data_directory = 'home/brian/Desktop/second_project_gcs/data/'
    star_masses = np.loadtxt(data_directory+'/star_masses/star_masses_index=%i.txt'%index_sorted)
    
    converter_parent = nbody_system.nbody_to_si(Mgalaxy, Rgalaxy)
    converter_sub = nbody_system.nbody_to_si(np.median(star_masses)|units.MSun, 5.|units.parsec) #masses list is in solar mass units
    
    '''
    takes in R, Z value
    returns VR, Vphi, VZ values
    should get appropriate 6D initial phase space conditions
    '''
    
    r_init, phi_init, z_init = rvals[index_sorted], phivals[index_sorted], zvals[index_sorted]
    vr_init, vphi_init, vz_init = vrvals[index_sorted], vphivals[index_sorted], vzvals[index_sorted]
    
    #convert from galpy/cylindrical to AMUSE/Cartesian units
    x_init = r_init*np.cos(phi_init) | units.kpc
    y_init = r_init*np.sin(phi_init) | units.kpc
    z_init = z_init | units.kpc
    
    #vphi = R \dot{\phi}? assuming yes for now
    vx_init = (vr_init*np.cos(phi_init) - vphi_init*np.sin(phi_init)) | units.kms
    vy_init = (vr_init*np.sin(phi_init) + vphi_init*np.cos(phi_init)) | units.kms
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
            
        if code_name == 'tree' or code_name == 'nemesis':
        
            #definition of gravity_code takes care of nemesis
            
            code = BHTree(converter_sub)
            code.particles.add_particles(bodies)
            code.commit_particles()
        
        return bodies, code
        
    if orbiter_name == 'SingleCluster':
        
        bodies, code, _ = star_cluster(rvals, phivals, zvals, vrvals, vphivals, vzvals, 
                                       masses, index_sorted, Mgalaxy, Rgalaxy, code_name)
        
        return bodies, code
