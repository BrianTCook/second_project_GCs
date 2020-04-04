#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 13:36:39 2020

@author: BrianTCook
"""

import numpy as np
import glob
    
point_files = glob.glob('/Users/BrianTCook/Desktop/Thesis/second_project_gcs/data/enbid_files/*.ascii')

for point_file in point_files:
    points = np.loadtxt(point_file)
    print(points.shape)
    
    