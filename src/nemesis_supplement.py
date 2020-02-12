#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 13:51:30 2020

@author: BrianTCook
"""

def getxv(converter, M1, a, e, ma=0):
    
    '''
    Get initial phase space coordinates (position and velocity) for an object around a central body
    
    converter - AMUSE unit converter
    M1        - Mass of the central body in AMUSE mass units
    a         - Semi-major axis of the orbit in AMUSE length units
    e         - Eccentricity of the orbit
    ma        - Mean anomaly of the orbit
    
    Returns: (x, v), the position and velocity vector of the orbit in AMUSE length and AMUSE length / time units
    '''
    
    kepler = Kepler(converter)
    kepler.initialize_code()
    kepler.initialize_from_elements(M1, a, e, mean_anomaly=ma) #Intoducing the attributes of the orbit
    
    x = quantities.as_vector_quantity(kepler.get_separation_vector()) #Grabbing the initial position and velocity
    v = quantities.as_vector_quantity(kepler.get_velocity_vector())
    
    kepler.stop()
    
    return x, v

def parent_worker():
    code = Mercury(converter_parent)
    code.parameters.epsilon_squared=0.| units.kpc**2
    code.parameters.end_time_accuracy_factor=0.
    code.parameters.dt_param=0.1
    return code

def sub_worker(parts):
    code = BHTree(converter_sub)
    code.parameters.inttype_parameter=code.inttypes.SHARED4
    return code

def py_worker():
    code=CalculateFieldForParticles(gravity_constant = constants.G)
    return code

'''
also for nemesis
'''

def smaller_nbody_power_of_two(dt, conv):

    nbdt = conv.to_nbody(dt).value_in(nbody_system.time)
    idt = np.floor(np.log2(nbdt))

    return conv.to_si( 2**idt | nbody_system.time)

def distance_function(ipart, jpart, eta=0.1/2., _G=constants.G):
    
    dx = ipart.x-jpart.x
    dy = ipart.y-jpart.y
    dz = ipart.z-jpart.z

    dr = np.sqrt(dx**2 + dy**2 + dz**2)
    dr3 = dr**1.5
    mu = _G*(ipart.mass + jpart.mass)

    tau = eta/2./2.**0.5*(dr3/mu)**0.5 #need an explanation for this!

    return tau

'''
eta = should be dt_param or dt_param/2. but it's not being defined for whatever reason
'''


def radius(sys, eta=0.1, _G=constants.G):

    #variable shouldn't be named radius
    ra = ((_G*sys.total_mass()*dt**2/eta**2)**(1/3.))
    ra = ra*((len(sys)+1)/2.)**0.75
    return 3.*ra #is this the roche-lobe radius?