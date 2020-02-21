#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 14:36:18 2020

@author: BrianTCook
"""

from amuse.lab import *
from amuse.ic.salpeter import new_salpeter_mass_distribution

import numpy as np
import time

t0 = time.time()

masses = np.loadtxt('/home/brian/Desktop/second_project_gcs/data/cluster_masses_for_sampling.txt')
print('number of masses to figure out: %i'%(len(masses)))

star_masses = []

for Mcluster in masses:
    
    Mcluster = Mcluster|units.MSun

    Nstars, Mmin_star, Mmax_star = 100, 0.1, 100.
    mZams_flag = 0
    
    while mZams_flag == 0:
        
        mZams = new_salpeter_mass_distribution(Nstars, Mmin_star|units.MSun, Mmax_star|units.MSun, random=np.random)
        mass_difference_ratio = (Mcluster - mZams.sum())/Mcluster
        
        if mass_difference_ratio > 0.01:
            Nstars += 1
            
        if mass_difference_ratio < -0.01:
            Nstars -= 1
            
        if np.abs(mass_difference_ratio) <= 0.01:            
            mZams_flag = 1
            
    star_masses.append(list(mZams.value_in(units.MSun)))
    
    i = len(star_masses)
    
    if i%10 == 0:
        print('number of masses computed: %i'%(len(star_masses)))
        print('time: %.04f minutes'%((time.time()-t0)/60.))  
        np.savetxt('star_masses.txt', star_masses)
