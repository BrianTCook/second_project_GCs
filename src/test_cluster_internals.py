from amuse.lab import *
import numpy as np
import glob

datadir_AMUSE = '/home/brian/Desktop/second_project_gcs/Enbid-2.0/AMUSE_data/'
datadir = '/home/brian/Desktop/second_project_gcs/data/'
ascii_files = 'enbid_*00000*16.ascii' #sufficient number to check
Norbiters = 16

filestrs = glob.glob(datadir_AMUSE+ascii_files)

r = np.loadtxt(datadir + 'ICs/dehnen_rvals.txt')
phi = np.loadtxt(datadir + 'ICs/dehnen_phivals.txt')
z = np.loadtxt(datadir + 'ICs/dehnen_zvals.txt')
vr = np.loadtxt(datadir + 'ICs/bovy_vrvals.txt')
vphi = np.loadtxt(datadir + 'ICs/bovy_vphivals.txt')
vz = np.loadtxt(datadir + 'ICs/bovy_vzvals.txt')

for filestr in filestrs:
	print(filestr)
