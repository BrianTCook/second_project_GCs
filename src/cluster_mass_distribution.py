#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 17:40:05 2020

creates a list of 10000 AMUSE masses
[100, 1000000] MSuns
sampled from appropriate power law
saves them as a .txt file

@author: BrianTCook
"""

from amuse.lab import *
import numpy as np

if __name__ in '__main__':
    
    mZams = new_powerlaw_mass_distribution(10000,
                                           float(int(1e2))|units.MSun,
                                           float(int(1e6))|units.MSun,
                                           alpha=-2.)
    
    np.savetxt('cluster_masses_for_sampling.txt', mZams.value_in(units.MSun))
