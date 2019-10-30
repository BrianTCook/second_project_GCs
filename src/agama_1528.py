import agama
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

#create three components of a composite galaxy potential
pot_bulge = agama.Potential(type='Sersic', mass = 1, scaleRadius=1, sersicIndex=4, axisRatioZ=0.6)
pot_disc = agama.Potential(type='Disk', mass=4, scaleRadius=3, scaleHeight=0.5)
pot_halo = agama.Potential(type='NFW', mass=25, scaleRadius=10)
pot = agama.Potential(pot_bulge, pot_disc, pot_halo)

#represent the density profile as a collection of partciles
snap = pot_bulge.sample(100000)

#create a potential from this N-body snapshot
pot_nbody = agama.Potential(type='Multipole', particles=snap, symmetry='Axisymmetric')

#choose the grid in radius to plot the profiles
r = np.linspace(0., 25., 250)
xyz = np.column_stack((r, r*0, r*0))

#circular velocity as a function of radius: total...
vcirc_total = np.sqrt(-r*pot.force(xyz)[:,0])

plt.figure()
plt.plot(r, vcirc_total, label='Total')
#...and for each potential component separately

for p in pot:
	plt.plot(r, np.sqrt(-r*p.force(xyz)[:,0]),
label=p.name() )
plt.legend(loc='lower right')
plt.show()
plt.savefig('fig1.png')
plt.close()

plt.figure()
for p in pot:
	plt.plot(r, p.density(xyz), label=p.name())
plt.plot(r, pot_nbody.density(xyz), label='from Nbody')
plt.yscale('log')
plt.legend()
plt.show()
plt.savefig('fig2.png')
plt.close()
