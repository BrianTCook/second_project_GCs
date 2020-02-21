#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 14:36:18 2020

@author: BrianTCook
"""

from amuse.lab import *
from amuse.ic.kingmodel import new_king_model
from amuse.ic.salpeter import new_salpeter_mass_distribution

masses = np.loadtxt('/home/brian/Desktop/second_project_gcs/data/cluster_masses_for_sampling.txt')

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
            
    star_masses.append(mZams)
    
np.savetxt('/home/brian/Desktop/second_project_gcs_data/star_masses.txt', star_masses)