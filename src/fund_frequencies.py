#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 15:31:44 2020

@author: BrianTCook
"""

import numpy as np
import glob
import matplotlib.pyplot as plt
from superfreq import SuperFreq, find_frequencies

datadir_AMUSE = '/Users/BrianTCook/Desktop/Thesis/second_project_gcs/Enbid-2.0/AMUSE_data/'

t = np.linspace(0., 100., 51) #0 to 100 Myr only saving every 2 Myr
sf = SuperFreq(t)

logN_max = 6
simulation_tags = [ 2**i for i in range(1, logN_max+1) ]

plt.rc('text', usetex = True)
plt.rc('font', family = 'serif')
    
plt.figure()

nstars_clusterzero = 1096
nstars_clusterone = 300

end_index = nstars_clusterzero + nstars_clusterone

for i, tag in enumerate(simulation_tags):
    
    print('Nclusters = %i'%(tag))
    
    ws = glob.glob(datadir_AMUSE + '*' + '_' + str(tag) + '.ascii')

    sample_array = np.loadtxt(ws[0])
    nstars, ndim = sample_array.shape

    print(nstars, ndim)

    pst_arr = np.zeros((nstars_clusterone, len(t), ndim))

    for j, wfile in enumerate(ws):
        
        time_slice = np.loadtxt(wfile)
        pst_arr[:,j,:] = time_slice[nstars_clusterzero:end_index]

    all_freqs = []

    for k in range(nstars_clusterone):
        
        w = pst_arr[k, :, :]
        
        fs = [ (w[:,l] * 1j*(0.001022)*w[:,l+ndim//2]) for l in range(ndim//2) ]
        result = sf.find_fundamental_frequencies(fs)
        
        all_freqs += list(result.fund_freqs)
        
    all_freqs = np.sort(all_freqs)
        
    plt.plot(all_freqs, [ i/len(all_freqs) for i in range(len(all_freqs)) ],
             linewidth=0.5, label=r'$\log_{2}N_{\mathrm{clusters}} = %i$'%(np.log2(tag)))


plt.title('Cluster 1, Fundamental Frequencies', fontsize=16)
plt.xlabel(r'Fundamental Frequency (Myr$^{-1}$)', fontsize=12)
plt.ylabel('Cumulative Distribution Function', fontsize=12)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
plt.tight_layout()
plt.savefig('fund_freqs_clusterone.pdf')