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

def make_king_model_cluster(Rcoord, Zcoord, phicoord, vr_init, vphi_init, vz_init, 
                            W0, Mcluster, code_name, parameters=[]):

    '''
    sets up a cluster with mass M and radius R
    which nbodycode would you like to use?
    '''

    mZams_flag, Nstars = 0, 5
    Mmin_star, Mmax_star = 0.1, 100.    
    
    while mZams_flag == 0:
        
        if abs(Rcoord - 10.) < 1e-3:
    
            #don't want to go through while loop if cluster is just for sub_worker
            converter = nbody_system.nbody_to_si(Mcluster, Rcluster)
            bodies = new_king_model(Nstars, W0, convert_nbody=converter)
            bodies.mass = [Mcluster]
            mZams_flag = 1
        
        mZams = new_salpeter_mass_distribution(Nstars, Mmin_star|units.MSun, Mmax_star|units.MSun)
        print('Mcluster is', Mcluster)
        print('mZams.sum() is', mZams.sum())
        mass_difference_ratio = (Mcluster - mZams.sum())/Mcluster
        
        if np.abs(mass_difference_ratio) > 0.01:
            Nstars -= 1
        if np.abs(mass_difference_ratio) < -0.01:
            Nstars += 1
        else:
            print('Mclusters, Nstars are', Mcluster, Nstars)
            converter = nbody_system.nbody_to_si(Mcluster, 1.|units.parsec)
            bodies = new_king_model(Nstars, W0, convert_nbody=converter)
            bodies.mass = mZams
            mZams_flag = 1

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
        
    if code_name == 'tree':
    
        code = BHTree(converter)
        code.particles.add_particles(bodies)
        code.commit_particles()
        
    if code_name == 'nemesis':
        
        '''
        ends up not being used?
        '''
        
        parts = HierarchicalParticles(bodies)
        
        dt = smaller_nbody_power_of_two(0.1 | units.Myr, converter_parent)
        dt_nemesis = dt
        dt_bridge = 0.01 * dt
        dt_param = 0.1
        
        nemesis = Nemesis(parent_worker, sub_worker, py_worker)
        nemesis.timestep = dt
        nemesis.distfunc = distance_function
        nemesis.threshold = dt_nemesis
        nemesis.radius = radius
        
        nemesis.commit_parameters()
        nemesis.particles.add_particles(parts)
        nemesis.commit_particles()
        
        code = nemesis
        
    for name,value in parameters:
        setattr(code.parameters, name, value)
    
    return bodies, code

def star_cluster(rvals, phivals, zvals, masses, index, code_name):
    
    '''
    takes 3 random numbers and generates open cluster
    with appropriate ICs in 6D phase space
    '''
    
    #limit to within 100 pc of the galactic center
    Rcoord, phicoord, Zcoord = rvals[index], phivals[index], zvals[index]
    
    Mcluster = masses[index]|units.MSun
    
    #using Staeckel
    aAS = actionAngleStaeckel(pot=MWPotential2014, delta=0.45, c=True)
    qdfS = quasiisothermaldf(1./3., 0.2, 0.1, 1., 1., pot=MWPotential2014, aA=aAS, cutcounter=True)
    vr_init, vphi_init, vz_init = qdfS.sampleV(Rcoord, Zcoord, n=1)[0,:]
    
    #220 km/s at 8 kpc, convert back to km/s
    to_kms = bovy_conversion.velocity_in_kpcGyr(220., 8.) * 0.9785
    vr_init *= to_kms
    vphi_init *= to_kms
    vz_init *= to_kms
    
    W0 = 1.5
    
    #just for converter, right?
    converter_sub = nbody_system.nbody_to_si(Mcluster, 10|units.parsec)
    
    bodies, code = make_king_model_cluster(Rcoord, Zcoord, phicoord, vr_init, vphi_init, vz_init, 
                                           W0, Mcluster, code_name, parameters=[])
    
    return bodies, code, converter_sub
