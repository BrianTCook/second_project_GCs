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

def make_king_model_cluster(Rcoord, Zcoord, phicoord, vr_init, vphi_init, vz_init, 
                            W0, Mcluster, Mstars, code_name, parameters=[]):

    '''
    sets up a cluster with mass M and radius R
    which nbodycode would you like to use?

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
    
    #initializing the masses
    bodies.mass = Mstars|units.Msun
    
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

def star_cluster(rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, star_masses, index, code_name):
    
    '''
    takes 3 random numbers and generates open cluster
    with appropriate ICs in 6D phase space
    '''
    
    #limit to within 1 kpc of the galactic center
    Rcoord, phicoord, Zcoord = rvals[index], phivals[index], zvals[index]
    vr_init, vphi_init, vz_init = vrvals[index], vphivals[index], vzvals[index]
    
    Mcluster = masses[index]|units.MSun
    Mstars = star_masses[index]|units.Msun
    
    W0 = 1.5
    
    #just for converter, right?
    converter_sub = nbody_system.nbody_to_si(Mcluster, 10|units.parsec)
    
    print('~~~~~~~~~~~~~~~~')
    print('index is: %i'%(index))
    print('r, phi, z: %.04f kpc, %.04f radians, %.04f kpc'%(Rcoord, phicoord, Zcoord))
    print('vr, vphi, vz: %.04f km/s, %.04f km/s, %.04f km/s'%(vr_init, vphi_init, vz_init))
    print('~~~~~~~~~~~~~~~~')
    
    bodies, code = make_king_model_cluster(Rcoord, Zcoord, phicoord, vr_init, vphi_init, vz_init, 
                                           W0, Mcluster, Mstars, code_name, parameters=[])
    
    return bodies, code, converter_sub
