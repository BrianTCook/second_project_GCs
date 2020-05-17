#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 15:31:44 2020

@author: BrianTCook
"""

import numpy as np
import glob
from superfreq import SuperFreq

datadir_AMUSE = '/home/brian/Desktop/second_project_gcs/Enbid-2.0/AMUSE_data/'

t = np.linspace(0., 100., 51) #0 to 100 Myr only saving every 2 Myr
ws = glob.glob(datadir_AMUSE + '*.ascii')

for wfile in ws:
    
    w = np.loadtxt(wfile)
    ntimes, ndim = w.shape

sf = SuperFreq(t)

fs = [(w[:,i] * 1j*w[:,i+ndim//2]) for i in range(ndim//2)]
freqs, tbl, ix = sf.find_fundamental_frequencies(fs)


