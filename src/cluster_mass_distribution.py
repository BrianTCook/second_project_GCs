#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 17:40:05 2020

@author: BrianTCook
"""

from amuse.lab import *
from amuse.ic.kingmodel import new_king_model
from amuse.community.mercury.interface import Mercury

if __name__ in '__main__':
    
    mZams = new_powerlaw_mass_distribution(10000,
                                           0.1|units.MSun,
                                           100.|units.Msun,
                                           alpha=-2.)
    
    np.savetxt('cluster_masses_for_sampling.txt', mZams)