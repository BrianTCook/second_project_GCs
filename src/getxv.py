from amuse.lab import *
from testing_nemesis import getxv

Mcluster = 5e6|units.MSun
Rcluster = 10.|units.parsec

converter = nbody_system.nbody_to_si(Mcluster, Rcluster)
sepBinary = 50.|units.parsec

dBinary, vBinary = getxv(converter, Mcluster, sepBinary, e=0)
print(dBinary)
print(vBinary)

print('hello world!')
