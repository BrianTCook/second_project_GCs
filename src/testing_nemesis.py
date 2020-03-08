#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:05:09 2020

@author: BrianTCook
"""

from amuse.lab import *
from amuse.couple import bridge

from galpy.df import quasiisothermaldf
from galpy.potential import MWPotential2014, to_amuse
from galpy.util import bovy_conversion
from galpy.actionAngle import actionAngleStaeckel

#other scripts
from phase_space_mapping import maps
from simulation_script import simulation
from create_plots import plotting_things
from convert_numpy_to_ascii import convert_numpy

import random
import numpy as np

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import time
import os

random.seed(73)

#Circumvent a problem with using too many threads on OpenMPI
os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "yes"

if __name__ in '__main__':
    
    potential = MWPotential2014 #galpy
    
    sepBinary = 20.|units.parsec #not necessary if not doing binary cluster part
    tend, dt = 0.2|units.Myr, 0.1|units.Myr
    #dt_param = 0.2 #for nemesis
    
    #uses a galpy function to evaluate the enclosed mass
    Mgalaxy, Rgalaxy = float(6.8e10)|units.MSun, 2.6|units.kpc #disk mass for MWPotential2014, Bovy(2015)
    
    data_directory = '/home/brian/Desktop/second_project_gcs/data/'
    
    rvals_all = np.loadtxt(data_directory+'ICs/dehnen_rvals.txt')
    phivals_all = np.loadtxt(data_directory+'ICs/dehnen_phivals.txt')
    zvals_all = np.loadtxt(data_directory+'ICs/dehnen_zvals.txt')
    
    vrvals = np.loadtxt(data_directory+'ICs/bovy_vrvals.txt')
    vphivals = np.loadtxt(data_directory+'ICs/bovy_vphivals.txt')
    vzvals = np.loadtxt(data_directory+'ICs/bovy_vzvals.txt')
    
    masses_all = np.loadtxt(data_directory+'ICs/cluster_masses_for_sampling.txt')

    logN_max = 2
    Norbiters_list = [ 2**i for i in range(logN_max) ] #need to make into a list at some point
    orbiter_names = [ 'SingleCluster' ] #,, 'SingleStar',  'BinaryCluster' 
    code_names = [ 'nemesis', 'tree', 'Nbody' ]

    t0 = time.time()
    
    plt.figure()

    plt.xlabel(r'$\log_{2} N_{\mathrm{clusters}}$', fontsize=20)
    plt.ylabel(r'Clock Time (minutes)', fontsize=20)
    
    for orbiter_name in orbiter_names:
        for code_name in code_names:
            
            yvals = []
            
            for Norbiters in Norbiters_list:
                
                rvals = rvals_all[:Norbiters]
                phivals = phivals_all[:Norbiters]
                zvals = zvals_all[:Norbiters]  
                masses = masses_all[:Norbiters]
                        
                print('\\\\\\\\\\\\\\\\\\\\\\\\')
                print(code_name, orbiter_name)
                print('\\\\\\\\\\\\\\\\\\\\\\\\')
                
                t_init = time.time()
                
                simulation(code_name, orbiter_name, potential, Mgalaxy, Rgalaxy, 
                           sepBinary, rvals, phivals, zvals, vrvals, vphivals, vzvals, 
                           masses, Norbiters, tend, dt)
              
                t_final = time.time()
                
                yvals.append((t_final-t_init)/60.)
                
                #maps(code_name, orbiter_name, Norbiters)
                
            plt.scatter(range(logN_max), yvals, label=code_name)
                
    
    plt.legend(loc='upper left', fontsize=12)
    plt.annotate(r'$t_{\mathrm{end}} = 0.2$ Myr', xy=(0.7, 0.25), xycoords='axes fraction', fontsize=14)
    plt.annotate(r'$\Delta t = 0.1$ Myr', xy=(0.7, 0.15), xycoords='axes fraction', fontsize=14)
    
    plt.gca().set_yscale('log')
    plt.gca().tick_params(labelsize='large')
    
    plt.tight_layout() 
    plt.savefig('clock_vs_Norbiters.pdf')
            
    #plotting_things(code_names, orbiter_names, Norbiters_list, tend, dt)
    #convert_numpy(code_names, orbiter_names, Norbiters_list)
    #entropy_stuff(code_names, orbiter_names, Norbiters_list)
