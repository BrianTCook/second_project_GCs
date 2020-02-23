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


Nclusters = 200
star_masses = [ [] for i in range(Nclusters) ]

j = 0 #counter

while j <= Nclusters:
    
    Mcluster = masses[j]|units.MSun

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
            
    star_masses = mZams.value_in(units.MSun)
    np.savetxt('star_masses_index=%i.txt'%j, star_masses)
        
    if j%10 == 0:
        print('j is', j)
        print('time is %.04 minutes'%((time.time()-t0)/60.))
    
    if j >= Nclusters:
        break
