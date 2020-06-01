#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:05:09 2020

@author: BrianTCook
"""

import os
import sys
import time

sys.path.append(os.getcwd())

from amuse.lab import *
from amuse.couple import bridge

import random
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

#other scripts
from phase_space_mapping import maps
from simulation_script import simulation
from create_plots import plotting_things
from convert_numpy_to_ascii import convert_numpy

from galpy.df import quasiisothermaldf
from galpy.potential import MWPotential2014, to_amuse
from galpy.util import bovy_conversion
from galpy.actionAngle import actionAngleStaeckel

random.seed(73)

#Circumvent a problem with using too many threads on OpenMPI
#os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "yes"

if __name__ in '__main__':
    
    potential = MWPotential2014 #galpy
    
    sepBinary = 20.|units.parsec #not necessary if not doing binary cluster part
    tend, dt = 100.|units.Myr, 0.2|units.Myr
            
    #dt_param = 0.2 #for nemesis
    
    #uses a galpy function to evaluate the enclosed mass
    Mgalaxy, Rgalaxy = float(6.8e10)|units.MSun, 2.6|units.kpc #disk mass for MWPotential2014, Bovy(2015)
    
    data_directory = '/home/brian/Desktop/second_project_gcs/data/' #'/home/s1780638/second_project_gcs/data/'
    
    rvals_all = np.loadtxt(data_directory+'ICs/dehnen_rvals.txt')
    phivals_all = np.loadtxt(data_directory+'ICs/dehnen_phivals.txt')
    zvals_all = np.loadtxt(data_directory+'ICs/dehnen_zvals.txt')
    
    vrvals_all = np.loadtxt(data_directory+'ICs/bovy_vrvals.txt')
    vphivals_all = np.loadtxt(data_directory+'ICs/bovy_vphivals.txt')
    vzvals_all = np.loadtxt(data_directory+'ICs/bovy_vzvals.txt')
    
    masses_all = np.loadtxt(data_directory+'ICs/cluster_masses_for_sampling.txt')
    radii_all = np.loadtxt(data_directory+'ICs/cluster_radii_for_sampling.txt')

    logN_max = 6
    Norbiters_list = [ 64 ] #[ 2**i for i in range(logN_max+1) ]
    orbiter_names = [ 'SingleCluster' ] #,, 'SingleStar',  'BinaryCluster' 
    code_names = [ 'tree' ] #, 'nemesis' ]  , 'Nbody'

    t0 = time.time()
    
    #plt.figure()

    #plt.xlabel(r'$\log_{2} N_{\mathrm{clusters}}$', fontsize=20)
    #plt.ylabel(r'Clock Time (minutes)', fontsize=20)
    
    for orbiter_name in orbiter_names:
        for code_name in code_names:
            
            #Nvals, yvals = [], []
            
            for Norbiters in Norbiters_list:
                
                print('\\\\\\\\\\\\\\\\\\\\\\\\')
                print(code_name, orbiter_name)
                print('\\\\\\\\\\\\\\\\\\\\\\\\')
                
                t_init = time.time()

                simulation(code_name, orbiter_name, potential, Mgalaxy, Rgalaxy, 
                           sepBinary, rvals_all, phivals_all, zvals_all, vrvals_all, vphivals_all, vzvals_all, 
                           masses_all, radii_all, Norbiters, tend, dt)
                
                print('time is: %.03f minutes'%((time.time()-t0)/60.))
                print('time to run last simulation: %.03f minutes'%((time.time()-t_init)/60.))
                
                #t_final = time.time()
                
                #Nvals.append(math.log(Norbiters, 2))
                #yvals.append((t_final-t_init)/60.)
    
                #plt.scatter(Nvals, yvals, label=code_name)
      
    '''
    plt.gca().set_yscale('log')
    plt.legend(loc='upper left', fontsize=12)
    plt.annotate(r'$t_{\mathrm{end}} = 0.3$ Myr', xy=(0.7, 0.25), xycoords='axes fraction', fontsize=14)
    plt.annotate(r'$\Delta t = 0.1$ Myr', xy=(0.7, 0.15), xycoords='axes fraction', fontsize=14)

    plt.gca().tick_params(labelsize='large')
    
    plt.tight_layout() 
    plt.savefig('clock_vs_Norbiters.pdf')  
    '''      
    
    #plotting_things(code_names, orbiter_names, Norbiters_list, tend, dt)
    #convert_numpy(code_names, orbiter_names, Norbiters_list)
    #entropy_stuff(code_names, orbiter_names, Norbiters_list)

sys.exit()
