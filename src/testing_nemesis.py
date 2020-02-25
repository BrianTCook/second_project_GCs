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

#other scripts
from phase_space_mapping import maps
from cluster_maker import star_cluster, orbiter
from gravity_code import gravity_code_setup
from evolve_model import simulation
from create_plots import plotting_things
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
    tend, dt = 100.|units.Myr, 1.|units.Myr
    dt_param = 0.1 #for nemesis
    
    #uses a galpy function to evaluate the enclosed mass
    Mgalaxy, Rgalaxy = float(6.8e10)|units.MSun, 2.6|units.kpc #disk mass for MWPotential2014, Bovy(2015)
    
    Norbiters = 4
    
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
            
            #maps(orbiter_name, code_name)
            
            print('current time: %.03f minutes'%((time.time()-t0)/60.))
            
    #plotting_things(orbiter_names, code_names, tend, dt)
