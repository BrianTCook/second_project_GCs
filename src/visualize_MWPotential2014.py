#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:52:16 2020

@author: BrianTCook
"""

from galpy.potential import MWPotential2014, plotPotentials
import matplotlib.pyplot as plt

ahh = plotPotentials(MWPotential2014, rmin=0.01, rmax=0.9, nrs=51, 
               zmin=-0.49, zmax=0.49, nzs=51, justcontours=True, ncontours=21,
	       aspect=None)

plt.savefig('MWPotential2014_plot.pdf')
