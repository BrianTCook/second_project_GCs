#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:05:09 2020

@author: BrianTCook
"""

from amuse.lab import *
#from amuse.ext.bridge import bridge
from amuse.couple import bridge

from galpy.df import quasiisothermaldf
from galpy.potential import MWPotential2014, to_amuse
from galpy.util import bovy_conversion
from galpy.actionAngle import actionAngleStaeckel

#other scripts
from phase_space_mapping import maps
from cluster_maker import star_cluster, orbiter
from gravity_code import gravity_code_setup
from simulation_script import simulation
from create_plots import plotting_things
from convert_numpy_to_ascii import convert_numpy
from nemesis import Nemesis, HierarchicalParticles
from nemesis_supplement import getxv, parent_worker, sub_worker, py_worker, smaller_nbody_power_of_two, distance_function, radius


import random
import numpy as np
import time
import os

random.seed(73)

#Circumvent a problem with using too many threads on OpenMPI
os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "yes"

if __name__ in '__main__':
    
    potential = MWPotential2014 #galpy
    
    sepBinary = 20.|units.parsec #not necessary if not doing binary cluster part
    tend, dt = 4.|units.Myr, 0.1|units.Myr
    dt_param = 0.1 #for nemesis
    
    #uses a galpy function to evaluate the enclosed mass
    Mgalaxy, Rgalaxy = float(6.8e10)|units.MSun, 2.6|units.kpc #disk mass for MWPotential2014, Bovy(2015)
    
    data_directory = '/home/brian/Desktop/second_project_gcs/data/'
    
    rvals = np.loadtxt(data_directory+'ICs/dehnen_rvals.txt')
    phivals = np.loadtxt(data_directory+'ICs/dehnen_phivals.txt')
    zvals = np.loadtxt(data_directory+'ICs/dehnen_zvals.txt')
    
    vrvals = np.loadtxt(data_directory+'ICs/bovy_vrvals.txt')
    vphivals = np.loadtxt(data_directory+'ICs/bovy_vphivals.txt')
    vzvals = np.loadtxt(data_directory+'ICs/bovy_vzvals.txt')
    
    masses = np.loadtxt(data_directory+'ICs/cluster_masses_for_sampling.txt')

    Norbiters_list = [ 1, 2 ] #need to make into a list at some point
    orbiter_names = [ 'SingleStar', 'SingleCluster' ]
    code_names = ['tree', 'Nbody' ]#, 'nemesis'

    t0 = time.time()
    
    for orbiter_name in orbiter_names:
        for code_name in code_names:
            for Norbiters in Norbiters_list:
                
                rvals = rvals[:Norbiters]
                phivals = phivals[:Norbiters]
                zvals = zvals[:Norbiters]  
                masses = masses[:Norbiters]
                        
                print('\\\\\\\\\\\\\\\\\\\\\\\\')
                print(code_name, orbiter_name)
                print('\\\\\\\\\\\\\\\\\\\\\\\\')
                
                simulation(code_name, orbiter_name, potential, Mgalaxy, Rgalaxy, 
                           sepBinary, rvals, phivals, zvals, vrvals, vphivals, vzvals, 
                           masses, Norbiters, tend, dt)
                
                '''
                filename = "data_%s_%s_Norbiters=%i.csv"%(code_name, orbiter_name, Norbiters)
                
                bodies = read_set_from_file(filename, "csv",
                                            attribute_types = (units.MSun, units.kpc, units.kpc, units.kpc, units.kms, units.kms, units.kms),
                                            attribute_names = ('mass', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
                
                print('bodies: %s, %s '%(code_name, orbiter_name))
                print(bodies.x)
                
                #print('current time: %.03f minutes'%((time.time()-t0)/60.))
                '''
              
                maps(code_name, orbiter_name, Norbiters)
                #print('current time: %.03f minutes'%((time.time()-t0)/60.))
        
    plotting_things(code_names, orbiter_names, Norbiters_list, tend, dt)
    convert_numpy(code_names, orbiter_names, Norbiters_list)
