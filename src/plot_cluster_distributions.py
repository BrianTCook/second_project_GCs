import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from coordinate_generator import radial_coords, zed_coords

r_OC, r_YMC, r_GC = radial_coords()
z_OC, z_YMC, z_GC = zed_coords()

plt.figure()
plt.hist(r_OC, density=True, cumulative=True, histtype='step', label='Open Clusters')
plt.hist(r_YMC, density=True, cumulative=True, histtype='step', label='Young Massive Clusters')
plt.hist(r_GC, density=True, cumulative=True, histtype='step', label='Globular Clusters')
plt.xlabel(r'$r$ (kpc)', fontsize=16)
plt.ylabel(r'Cumulative Density Function', fontsize=16)
plt.legend(loc='upper right')
plt.savefig('radial_coord_distribution.pdf')
plt.close()

plt.figure()
plt.hist(z_OC, density=True, cumulative=True, histtype='step', label='Open Clusters')
plt.hist(z_YMC, density=True, cumulative=True, histtype='step', label='Young Massive Clusters')
plt.hist(z_GC, density=True, cumulative=True, histtype='step', label='Globular Clusters')
plt.xlabel(r'$z$ (kpc)', fontsize=16)
plt.ylabel(r'Cumulative Density Function', fontsize=16)
plt.legend(loc='upper right')
plt.savefig('zed_coord_distribution.pdf')
plt.close()
