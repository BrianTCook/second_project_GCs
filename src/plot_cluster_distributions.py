import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from coordinate_generator import radial_coords, zed_coords

r_OC, r_YMC, r_GC = radial_coords()
z_OC, z_YMC, z_GC = zed_coords()

plt.figure()
plt.plot(r_OC, label='Open Clusters')
plt.plot(r_YMC, label='Young Massive Clusters')
plt.plot(r_GC, label='Globular Clusters')
plt.xlabel(r'$r$ (kpc)', fontsize=16)
plt.ylabel(r'Probability Density', fontsize=16)
plt.legend(loc='upper right')
plt.savefig('radial_coord_distribution.pdf')
plt.close()

plt.figure()
plt.plot(z_OC, label='Open Clusters')
plt.plot(z_YMC, label='Young Massive Clusters')
plt.plot(z_GC, label='Globular Clusters')
plt.xlabel(r'$z$ (kpc)', fontsize=16)
plt.ylabel(r'Probability Density', fontsize=16)
plt.legend(loc='upper right')
plt.savefig('zed_coord_distribution.pdf')
plt.close()
