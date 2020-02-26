#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 14:15:40 2020

@author: BrianTCook
"""

from amuse.lab import *

code_name, orbiter_name, Norbiters = 'tree', 'SingleCluster', 4

data_directory = '/home/brian/Desktop/second_project_gcs/data/'
data_file_name = 'data_%s_%s_Norbiters=%i.csv'%(code_name, orbiter_name, Norbiters)

many_snapshots = read_set_from_file(data_file_name, 'txt')

for snapshot in many_snapshots.history:
    print('Nbodies in this snapshot is: ', len(snapshot))