#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 23:32:22 2020

@author: BrianTCook
"""

import numpy as np
import glob
import matplotlib.pyplot as plt

dt_values = [ 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1. ]
tend = 100.

energy_values = glob.glob('/Users/BrianTCook/Desktop/Thesis/second_project_gcs/data/tree_data/energy_check_data/*.txt')
print(energy_values)

plt.rc('font', family='serif')
plt.rc('text', usetex=True)

plt.figure()

for i, dt in enumerate(dt_values):
    
    energies = np.loadtxt(energy_values[i])  
    sim_times = np.linspace(0., tend+dt, len(energies))
    plt.plot(sim_times, energies, linewidth=1, label=r'$\Delta t$ = %.02f Myr'%(dt))
    
plt.ylim(-0.15, 0.05)
plt.axhline(y=0, linewidth=1, linestyle='--', color='k')
plt.xlabel(r'$t_{\mathrm{sim}}$ (Myr)', fontsize=16)
plt.ylabel(r'$\Delta E / E(t=0)$', fontsize=16)
plt.gca().tick_params(labelsize='large')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
plt.annotate(r'$\log_{2}N_{\mathrm{clusters}}=3$', xy=(0.1,0.1), xycoords='axes fraction', fontsize=14)
plt.title('Energy Conservation, Tree Code', fontsize=16)
plt.tight_layout()
plt.savefig('energy_conservation_plot.pdf')