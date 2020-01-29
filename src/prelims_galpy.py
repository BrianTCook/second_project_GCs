from amuse.lab import *
from amuse.ext.bridge import bridge
from amuse.couple.bridge import CalculateFieldForParticles
from amuse.ic.kingmodel import new_king_model

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from galpy.df import quasiisothermaldf
from galpy.potential import MWPotential2014, evaluaterforces, evaluatezforces, vcirc, to_amuse
from galpy.util import bovy_conversion, bovy_plot

from astropy import units as u

Rcoord = 5. * np.random.random() #kiloparsecs
Zcoord = 1. * np.random.random() #kiloparsecs
phicoord = 2 * np.pi * np.random.random() * u.radian

galaxy_code = to_amuse(MWPotential2014, t=0.0, tgalpy=0.0, reverse=False, ro=None, vo=None)

Mcluster, Rcluster = 5e6|units.MSun, 10.|units.parsec
converter = nbody_system.nbody_to_si(Mcluster, Rcluster)

N, W0 = 100, 1.5
bodies = new_king_model(N, W0, convert_nbody=converter)

bodies.x += 10.|units.kpc
bodies.vy += 220.|units.km/units.s

code = BHTree(converter)
code.particles.add_particles(bodies)

gravity = bridge()
gravity.add_system(code, (galaxy_code,))

channel_from_bodies_to_code = bodies.new_channel_to(code.particles)
channel_from_code_to_bodies = code.particles.new_channel_to(bodies)

tend, dt = 100., 1.
times = np.arange(0., tend, dt)

for i, t in enumerate(times):
	t = t|units.Myr
	gravity.evolve_model(t)

	print('t is', t)

channel_from_code_to_bodies.copy()
#gravity.stop()

plt.figure()
plt.scatter(bodies.x.value_in(units.kpc), bodies.y.value_in(units.kpc), s=4, c='k')
plt.xlabel('X (kpc)')
plt.ylabel('Y (kpc)')
plt.savefig('bodies_in_MWPotential2014.png')

