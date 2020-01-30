from amuse.lab import *
from amuse.ext.bridge import bridge
from amuse.couple.bridge import CalculateFieldForParticles
from amuse.ic.kingmodel import new_king_model
from amuse.io import write_set_to_file, read_set_from_file

from galpy.df import quasiisothermaldf
from galpy.potential import MWPotential2014, to_amuse
from galpy.util import bovy_conversion
from galpy.actionAngle import actionAngleStaeckel

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import random
import numpy as np
import time
import os

random.seed(73)

from nemesis import Nemesis, HierarchicalParticles

#Circumvent a problem with using too many threads on OpenMPI
os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "yes"

Rmax, Zmax = 5., 1. #in kiloparsecs
Mgal, Rgal = 1.6e10|units.MSun, 1000.|units.parsec
Nstars, W0, Mcluster, Rcluster = 40, 1.5, 100.|units.MSun, 1.|units.parsec
Rinit = 1000.|units.parsec
parameters = [('epsilon_squared', 0.01|(units.parsec**2))]
t_end, dt = 100.|units.Myr, 1.|units.Myr

def make_king_model_cluster(N, W0, M, R, parameters=[]):

    '''
    sets up a cluster with mass M and radius R
    which nbodycode would you like to use?
    '''
    
    converter = nbody_system.nbody_to_si(M,R)
    bodies = new_king_model(N, W0, convert_nbody=converter)
    code = BHTree(converter) #sub_worker in Nemesis    

    for name,value in parameters:
        setattr(code.parameters, name, value)
    code.particles.add_particles(bodies)
    
    return bodies, code
    
'''
worker functions are for nemesis
'''

def parent_worker():
    converter_parent = nbody_system.nbody_to_si(Mgal, Rgal)
    code = Hermite(converter_parent)
    code.parameters.epsilon_squared=0.| units.kpc**2
    code.parameters.end_time_accuracy_factor=0.
    code.parameters.dt_param=0.1
    #print code.parameters.dt_dia.in_(units.yr)
    return code

def sub_worker(parts):
    converter_sub = nbody_system.nbody_to_si(Mcluster, Rcluster)
    code = BHTree(converter_sub)
    return code

def py_worker():
    code=CalculateFieldForParticles(gravity_constant = constants.G)
    return code

'''
also for nemesis
'''

def smaller_nbody_power_of_two(dt, conv):

    nbdt = conv.to_nbody(dt).value_in(nbody_system.time)
    idt = np.floor(np.log2(nbdt))

    return conv.to_si( 2**idt | nbody_system.time)

dt_param = 0.1

def distance_function(ipart, jpart, eta=dt_param/2., _G=constants.G):
    
	dx = ipart.x-jpart.x
	dy = ipart.y-jpart.y
	dz = ipart.z-jpart.z

	dr = np.sqrt(dx**2 + dy**2 + dz**2)
	dr3 = dr**1.5
	mu = _G*(ipart.mass + jpart.mass)

	tau = eta/2./2.**0.5*(dr3/mu)**0.5 #need an explanation for this!

	return tau

def radius(sys, eta=dt_param, _G=constants.G):

	#variable shouldn't be named radius
	ra = ((_G*sys.total_mass()*dt**2/eta**2)**(1/3.))
	ra = ra*((len(sys)+1)/2.)**0.75
	return 100.*ra

converter = nbody_system.nbody_to_si(Mcluster, Rcluster)
stars_all = new_king_model(Nstars, W0, convert_nbody=converter)

Rcoord = Rmax * np.random.random()
Zcoord = Zmax * np.random.random()
phicoord = 2*np.pi * np.random.random()

#using Staeckel, whatever that means
aAS = actionAngleStaeckel(pot=MWPotential2014, delta=0.45, c=True)
qdfS = quasiisothermaldf(1./3., 0.2, 0.1, 1., 1., pot=MWPotential2014, aA=aAS, cutcounter=True)
vr_init, vphi_init, vz_init = qdfS.sampleV(Rcoord, Zcoord, n=1)[0,:]

#convert from galpy/cylindrical to AMUSE/Cartesian units
x_init = (Rcoord*np.cos(phicoord)) | units.kpc
y_init = Rcoord*np.sin(phicoord) | units.kpc
z_init = Zcoord | units.kpc

vx_init = (vr_init*np.cos(phicoord) - Rcoord*vphi_init*np.sin(phicoord)) | units.kms
vy_init = (vr_init*np.sin(phicoord) + Rcoord*vphi_init*np.cos(phicoord)) | units.kms
vz_init = vz_init | units.kms

print('x_init is', x_init)
print('y_init is', y_init)
print('z_init is', z_init)
print('vx_init is', vx_init)
print('vy_init is', vy_init)
print('vz_init is', vz_init)


for star in stars_all:
    
    star.x += x_init
    star.y += y_init
    star.z += z_init
    star.vx += vx_init
    star.vy += vy_init
    star.vz += vz_init

parts = HierarchicalParticles(stars_all)

converter_parent = nbody_system.nbody_to_si(Mgal, Rgal)

dt = smaller_nbody_power_of_two(0.1 | units.Myr, converter_parent)
dt_nemesis = dt
dt_bridge = 0.01 * dt

nemesis = Nemesis( parent_worker, sub_worker, py_worker)
nemesis.timestep = dt
nemesis.distfunc = distance_function
nemesis.threshold = dt_nemesis
nemesis.radius = radius

nemesis.commit_parameters()
nemesis.particles.add_particles(parts)
nemesis.commit_particles()

channel_to_nemesis = stars_all.new_channel_to(nemesis.particles.all())

#MWPotential2014
galaxy_code = to_amuse(MWPotential2014, t=0., tgalpy=0., reverse=False, ro=None, vo=None)

gravity = bridge()
gravity.add_system(nemesis, (galaxy_code,) )
gravity.timestep = dt_bridge

channel_from_bodies_to_code = stars_all.new_channel_to(gravity.particles)
channel_from_code_to_bodies = gravity.particles.new_channel_to(stars_all)

gravity_solver_str = 'nemesis'

sim_times_unitless = np.arange(0., t_end.value_in(units.Myr), dt.value_in(units.Myr))
sim_times = [ t|units.Myr for t in sim_times_unitless]

t0 = time.time()
clock_times = []
mean_radial_coords = []
mean_speeds = []

#total number of stars in the simulation
Ntotal = len(stars_all)

#create an R^3 matrix to house phase space data for all particles
phase_space_data = np.zeros((len(sim_times), 6, len(gravity.particles)))

for j, t in enumerate(sim_times):

	clock_time = time.time()-t0

	if j%5 == 0:
		print('j=%i, time = %.02f seconds'%(j, clock_time))

	clock_times.append(clock_time)

	#for figures 3 through 6, November 24

	x = [ xx.value_in(units.parsec) for xx in gravity.particles.x ]
	y = [ yy.value_in(units.parsec) for yy in gravity.particles.y ]
	z = [ zz.value_in(units.parsec) for zz in gravity.particles.z ]

	vx = [ vxx.value_in(units.kms) for vxx in gravity.particles.vx ]
	vy = [ vyy.value_in(units.kms) for vyy in gravity.particles.vy ]
	vz = [ vzz.value_in(units.kms) for vzz in gravity.particles.vz ]

	for k, star in enumerate(gravity.particles):

		phase_space_data[j, 0, k] = x[k]
		phase_space_data[j, 1, k] = y[k]
		phase_space_data[j, 2, k] = z[k]
		phase_space_data[j, 3, k] = vx[k]
		phase_space_data[j, 4, k] = vy[k]
		phase_space_data[j, 5, k] = vz[k] 

	xmean, ymean, zmean = np.sum(x)/Ntotal, np.sum(y)/Ntotal, np.sum(z)/Ntotal
	mean_rval = np.sqrt(xmean**2 + ymean**2 + zmean**2)
	mean_radial_coords.append(mean_rval)

	vxmean, vymean, vzmean = np.sum(vx)/Ntotal, np.sum(vy)/Ntotal, np.sum(vz)/Ntotal
	mean_speed = np.sqrt(vxmean**2 + vymean**2 + vzmean**2)
	mean_speeds.append(mean_speed)

	#for .gif of orbits

	plt.figure()
	plt.scatter(x, y, marker='*', s=2, color='k')

	Rinit_in_pc = Rinit.value_in(units.parsec)
	#plt.xlim(-3*Rinit_in_pc, 3*Rinit_in_pc)
	#plt.ylim(-3*Rinit_in_pc, 3*Rinit_in_pc)
	plt.annotate(gravity_solver_str, xy=(0.8, 0.8), xycoords='axes fraction')

	plt.xlabel('$x$ (pc)', fontsize=12)
	plt.ylabel('$y$ (pc)', fontsize=12)
	plt.title('Time: %.02f Myr'%(t.value_in(units.Myr)))
	plt.tight_layout()
	plt.savefig('frame_%s_%s_%i.png'%(str(j).rjust(4, '0'), gravity_solver_str, 1))
	plt.close()        

	gravity.evolve_model(t)

channel_from_code_to_bodies.copy()

'''
print(gravity_solver_str)
print(mean_radial_coords)
gravity_solver_info.append([gravity_solver_str, clock_times,
			    mean_radial_coords, mean_speeds])

np.save('time_data_Nclusters=%i_%s.npy'%(Nclusters, gravity_solver_str), sim_times_unitless)
np.save('sixD_data_Nclusters=%i_%s.npy'%(Nclusters, gravity_solver_str), phase_space_data)

plt.figure()

for gs, clock, radii, speeds in gravity_solver_info:

	#clock time versus simulation time
	plt.plot(sim_times_unitless, clock, label=gs)

plt.xlabel('Simulation time (Myr)')
plt.ylabel('Clock time (seconds)')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('clocktime.png')
plt.close()

plt.figure()

for gs, clock, radii, speeds in gravity_solver_info:

	#<radial coordinate>
	plt.plot(sim_times_unitless, mean_radial_coords, label=gs)

plt.xlabel('Simulation time (Myr)')
plt.ylabel('Average distance from origin (pc)')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('radial_coords.png')
plt.close()

plt.figure()

for gs, clock, radii, speeds in gravity_solver_info:

	#<velocity>
	plt.semilogy(sim_times_unitless, mean_speeds, label=gs)

plt.xlabel('Simulation time (Myr)')
plt.ylabel('Average speed (km/s)')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('speeds.png')
plt.close()
'''
