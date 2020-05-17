#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 15:31:44 2020

@author: BrianTCook
"""

import numpy as np
import glob
from superfreq import SuperFreq

#datadir_AMUSE = '/home/brian/Desktop/second_project_gcs/Enbid-2.0/AMUSE_data/'
datadir_AMUSE = '/Users/BrianTCook/Desktop/Thesis/second_project_gcs/Enbid-2.0/AMUSE_data/'

t = np.linspace(0., 100., 51) #0 to 100 Myr only saving every 2 Myr
sf = SuperFreq(t)

logN_max = 6
simulation_tags = [ str(2**i) for i in range(logN_max+1) ]

phase_space_time_arrays = []

for i, tag in enumerate(simulation_tags):

    ws = glob.glob(datadir_AMUSE + '*' + '_' + tag + '.ascii')

    sample_array = np.loadtxt(ws[0])
    nstars, ndim = sample_array.shape

    print(nstars, ndim)

    pst_arr = np.zeros((nstars, ndim, len(t)))

    for j, wfile in enumerate(ws):
        
        time_slice = np.loadtxt(wfile)
        pst_arr[:,:,j] = time_slice
        
    fs = [ (pst_arr[:,k] * 1j*pst_arr[:,k+ndim//2]) for k in range(ndim//2) ]
    freqs, tbl, ix = sf.find_fundamental_frequencies(fs)


