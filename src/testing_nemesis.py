from amuse.lab import *
from amuse.ext.bridge import bridge
from amuse.couple.bridge import CalculateFieldForParticles
from amuse.ic.kingmodel import new_king_model
from amuse.io import write_set_to_file, read_set_from_file
from amuse.units import quantities

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from galpy.df import quasiisothermaldf
from galpy.potential import MWPotential2014, evaluateRforces, evaluatezforces, to_amuse
from galpy.util import bovy_conversion

import random
import numpy as np
import time
import os

random.seed(73)

from nemesis import Nemesis, HierarchicalParticles

#Circumvent a problem with using too many threads on OpenMPI
os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "yes"
    
def getxv(converter, M1, a, e, ma=0):
    
    '''
    Get initial phase space coordinates (position and velocity) for an object around a central body
    
    converter - AMUSE unit converter
    M1        - Mass of the central body in AMUSE mass units
    a         - Semi-major axis of the orbit in AMUSE length units
    e         - Eccentricity of the orbit
    ma        - Mean anomaly of the orbit
    
    Returns: (x, v), the position and velocity vector of the orbit in AMUSE length and AMUSE length / time units
    '''
    
    kepler = Kepler(converter)
    kepler.initialize_code()
    kepler.initialize_from_elements(M1, a, e, mean_anomaly=ma) #Intoducing the attributes of the orbit
    
    x = quantities.as_vector_quantity(kepler.get_separation_vector()) #Grabbing the initial position and velocity
    v = quantities.as_vector_quantity(kepler.get_velocity_vector())
    
    kepler.stop()
    
    return x, v

def make_king_model_cluster(Nstars, W0, Mcluster, Rcluster, code_name, parameters=[]):

    '''
    sets up a cluster with mass M and radius R
    which nbodycode would you like to use?
    '''
    
    converter = nbody_system.nbody_to_si(Mcluster, Rcluster)
    bodies = new_king_model(Nstars, W0, convert_nbody=converter)
    
    #sub_worker in Nemesis
    if code_name == 'Nbody':
        
        code = Hermite(converter)
        
    if code_name == 'tree':
    
        code = BHTree(converter)
        
    for name,value in parameters:
        setattr(code.parameters, name, value)
    code.particles.add_particles(bodies)
    
    return bodies, code

def orbiter_not_nemesis(orbiter_name, code_name, Rmax, Zmax,
                        Nstars, W0, Mcluster, Rcluster, dBinary):

    converter = nbody_system.nbody_to_si(Mcluster, Rcluster)
    
    '''
    takes in R, Z value
    returns VR, Vphi, VZ values
    should get appropriate 6D initial phase space conditions
    '''
    
    Rcoord = Rmax * np.random.random()
    Zcoord = Zmax * np.random.random()
    phicoord = 2*np.pi * np.random.random()
    
    #need to fix, but want to make sure it's working
    #vels_init = quasiisothermaldf.sampleV(Rcoord, Zcoord, n=1)
    vr_init, vphi_init, vz_init = 50.*np.random.random(), 50.*np.random.random(), 50.*np.random.random() #vs[0,:]
    
    #convert from galpy/cylindrical to AMUSE/Cartesian units
    x_init = (Rcoord*np.cos(phicoord))*1000. | units.parsec
    y_init = Rcoord*np.sin(phicoord)*1000. | units.parsec
    z_init = Zcoord*1000. | units.parsec
    
    vx_init = (vr_init*np.cos(phicoord) - Rcoord*vphi_init*np.sin(phicoord)) | units.kms
    vy_init = (vr_init*np.sin(phicoord) + Rcoord*vphi_init*np.cos(phicoord)) | units.kms
    vz_init = vz_init | units.kms
    
    if orbiter_name == 'SingleStar':
        
        bodies = Particles(1)
        bodies[0].mass = 1|units.MSun
        
        #right place in phase space
        print(bodies[0])
        bodies[0].x = [ x_init|units.kpc ]
        bodies[0].y = [ y_init|units.kpc ]
        bodies[0].z = [ z_init|units.kpc ]
        bodies[0].vx = [ vx_init|units.kms ]
        bodies[0].vy = [ vy_init|units.kms ]
        bodies[0].vz = [ vz_init|units.kms ]
        
        #sub_worker in Nemesis, should not matter for SingleStar
        if code_name == 'Nbody':
            
            code = Hermite(converter)
            
        if code_name == 'tree':
        
            code = BHTree(converter)
        
        return bodies, code
        
    if orbiter_name == 'SingleCluster':
        
        bodies, code = make_king_model_cluster(Nstars, W0, Mcluster, Rcluster, code_name)
        
        '''
        need to initialize initial phase space coordinates with AGAMA or galpy
        '''
    
        for body in bodies:
            #right place in phase space
            body.x += x_init
            body.y += y_init
            body.z += z_init
            body.vx += vx_init
            body.vy += vy_init
            body.vz += vz_init
                
        return bodies, code
        
    if orbiter_name == 'BinaryCluster':
        
        bodies_one, code_one = make_king_model_cluster(Nstars, W0, Mcluster, Rcluster, code_name)
        bodies_two, code_two = make_king_model_cluster(Nstars, W0, Mcluster, Rcluster, code_name)
        
        for body in bodies_one:
            #right place in phase space
            body.x += x_init
            body.y += y_init
            body.z += z_init
            body.vx += vx_init
            body.vy += vy_init
            body.vz += vz_init
    
        for body in bodies_two:
            #right place in phase space
            body.x += x_init
            body.y += y_init
            body.z += z_init
            body.vx += vx_init
            body.vy += vy_init
            body.vz += vz_init
    
        '''
        need to initialize initial phase space coordinates with AGAMA or galpy
        '''
            
        mass_one, mass_two = bodies_one.mass.sum(), bodies_two.mass.sum()
        total_mass = mass_one + mass_two
        
        dBinary, vBinary = getxv(converter, total_mass, dBinary, e=0)
        
        for body in bodies_one:
            body.position += dBinary * mass_one/total_mass
            body.velocity += vBinary * mass_one/total_mass
            
        for body in bodies_two:
            body.position -= dBinary * mass_two/total_mass
            body.velocity -= vBinary * mass_two/total_mass
        
        bodies = Particles(0)
        bodies.add_particles(bodies_one)
        bodies.add_particles(bodies_two)
        
        return bodies, code_one, code_two #need to be different so they're bridged

def orbiter_nemesis(orbiter_name, code_name):
    
    return None, None

def gravity_code_setup(orbiter_name, code_name, galaxy_code, Rmax, Zmax,
                       Nstars, W0, Mcluster, Rcluster, dBinary):
    
    '''
    will need to ask SPZ if he meant for field, orbiter to be separate in non
    Nemesis gravity solvers?
    '''
    
    if code_name != 'Nemesis':
        
        gravity = bridge()
        
        if orbiter_name != 'BinaryCluster':
            
            orbiter_bodies, orbiter_code = orbiter_not_nemesis(orbiter_name, code_name, Rmax, Zmax,
                                                               Nstars, W0, Mcluster, Rcluster, dBinary)
            
        if orbiter_name == 'BinaryCluster':
            
            orbiter_bodies, orbiter_code_one, orbiter_code_two = orbiter_not_nemesis(orbiter_name, code_name, Rmax, Zmax,
                                                                                     Nstars, W0, Mcluster, Rcluster, dBinary)
    
        #gravity.particles.add_particles(orbiter_bodies)
    
        #bridges clusters properly, independent of how many there are because
        #orbiter_code is not defined in binary case
        try:
            gravity.add_system(orbiter_code, (galaxy_code,))
        except:
            gravity.add_system(orbiter_code_one, (galaxy_code,))
            gravity.add_system(orbiter_code_two, (galaxy_code,))
            gravity.add_system(orbiter_code_one, (orbiter_code_two,))
            gravity.add_system(orbiter_code_two, (orbiter_code_one,))

    return orbiter_bodies, gravity

def simulation(orbiter_name, code_name, potential, Rmax, Zmax,  
               Nstars, W0, Mcluster, Rcluster, dBinary, tend, dt):
    
    galaxy_code = to_amuse(potential, t=0.0, tgalpy=0.0, reverse=False, ro=None, vo=None)
    bodies, gravity = gravity_code_setup(orbiter_name, code_name, galaxy_code, Rmax, Zmax, Nstars, W0, Mcluster, Rcluster, dBinary)
    
    channel_from_bodies_to_code = bodies.new_channel_to(gravity.particles)
    channel_from_code_to_bodies = gravity.particles.new_channel_to(bodies)
    
    Ntotal = len(gravity.particles)
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr), dt.value_in(units.Myr))
    sim_times = [ t|units.Myr for t in sim_times_unitless ]
    
    energies, mean_radial_coords, mean_speeds, clock_times = [], [], [], []
    
    t0 = time.time()
    
    for j, t in enumerate(sim_times):

        clock_times.append(time.time()-t0) #will be in seconds
        
        energy = gravity.kinetic_energy.value_in(units.J) + gravity.potential_energy.value_in(units.J)
        energies.append(energy)
        
        x = [ xx.value_in(units.parsec) for xx in gravity.particles.x ]
        y = [ yy.value_in(units.parsec) for yy in gravity.particles.y ]
        z = [ zz.value_in(units.parsec) for zz in gravity.particles.z ]
        
        vx = [ vxx.value_in(units.kms) for vxx in gravity.particles.vx ]
        vy = [ vyy.value_in(units.kms) for vyy in gravity.particles.vy ]
        vz = [ vzz.value_in(units.kms) for vzz in gravity.particles.vz ]
        
        xmean, ymean, zmean = np.sum(x)/Ntotal, np.sum(y)/Ntotal, np.sum(z)/Ntotal
        mean_rval = np.sqrt(xmean**2 + ymean**2 + zmean**2)
        mean_radial_coords.append(mean_rval)
        
        vxmean, vymean, vzmean = np.sum(vx)/Ntotal, np.sum(vy)/Ntotal, np.sum(vz)/Ntotal
        mean_speed = np.sqrt(vxmean**2 + vymean**2 + vzmean**2)
        mean_speeds.append(mean_speed)
                
        gravity.evolve_model(t)

        #filename = code_name + '_' + orbiter_name + '_data.hdf5'
        #write_set_to_file(gravity.particles, filename, "hdf5")
        
    channel_from_code_to_bodies.copy()
    #gravity.stop()
    
    np.savetxt(code_name + '_' + orbiter_name + '_energies.txt', energies)
    np.savetxt(code_name + '_' + orbiter_name + '_mean_radial_coords.txt', mean_radial_coords)
    np.savetxt(code_name + '_' + orbiter_name + '_mean_speeds.txt', mean_speeds)
    np.savetxt(code_name + '_' + orbiter_name + '_clock_times.txt', clock_times)
    
    return 0

def plotting_things(orbiter_names, code_names, tend, dt):
    
    #filename = code_name + '_' + orbiter_name + '_data.hdf5'
    #star_data = read_set_from_file(filename)
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr), dt.value_in(units.Myr))
    
    npanels_x = len(orbiter_names)
    
    '''
    3 by 1, a plot for each orbiter (single star, etc.)
    each panel contains 3 curves, one for each gravity solver
    '''
    
    #energies
    
    fig, axs = plt.subplots(1, 3)

    for i, orbiter_name in enumerate(orbiter_names): 
        
        axs[i].set_title(orbiter_name)
        
        for code_name in code_names:
            
            energies = np.loadtxt(code_name + '_' + orbiter_name + '_energies.txt')
            axs[i].plot(sim_times_unitless, energies, label=code_name)
            
        axs[i].legend(loc='best')
            
    plt.savefig('testing_nemesis_energy.png')
    plt.close()
    
    #mean radial coordinates
    
    fig, axs = plt.subplots(1, 3)

    for i, orbiter_name in enumerate(orbiter_names): 
        
        axs[i].set_title(orbiter_name)
        
        for code_name in code_names:
            
            mean_radial_coords = np.loadtxt(code_name + '_' + orbiter_name + '_mean_radial_coords.txt')
            axs[i].plot(sim_times_unitless, mean_radial_coords, label=code_name)
            
        axs[i].legend(loc='best')
            
    plt.savefig('testing_nemesis_radialcoords.png')
    plt.close()
    
    #mean speeds
    
    fig, axs = plt.subplots(1, 3)

    for i, orbiter_name in enumerate(orbiter_names): 
        
        axs[i].set_title(orbiter_name)
        
        for code_name in code_names:
            
            mean_speeds = np.loadtxt(code_name + '_' + orbiter_name + '_mean_speeds.txt')
            axs[i].plot(sim_times_unitless, mean_speeds, label=code_name)
            
        axs[i].legend(loc='best')
            
    plt.savefig('testing_nemesis_speeds.png')
    plt.close()
    
    #clock times
    
    fig, axs = plt.subplots(1, 3)

    for i, orbiter_name in enumerate(orbiter_names): 
        
        axs[i].set_title(orbiter_name)
        
        for code_name in code_names:
            
            clock_times = np.loadtxt(code_name + '_' + orbiter_name + '_clock_times.txt')
            axs[i].plot(sim_times_unitless, clock_times, label=code_name)
            
        axs[i].legend(loc='best')
            
    plt.savefig('testing_nemesis_clocktimes.png')
    plt.close()
    
    return 0

if __name__ in '__main__':
    
    potential = MWPotential2014
    Rmax, Zmax = 5., 1. #in kpc
    Nstars, W0 = 40, 1.5 #cluster parameters
    Mcluster, Rcluster = 5e6|units.MSun, 10|units.parsec
    dBinary = 10.|units.parsec
    tend, dt = 100.|units.Myr, 1.|units.Myr
    
    orbiter_names = [ 'SingleStar', 'SingleCluster', 'BinaryCluster' ]
    code_names = [ 'tree', 'Nbody' ] # 'Nemesis'
    
    for orbiter_name in orbiter_names:
        for code_name in code_names:
            	simulation(orbiter_name, code_name, potential, Rmax, Zmax, 
            	           Nstars, W0, Mcluster, Rcluster, dBinary, tend, dt)
            
    plotting_things(orbiter_names, code_names, tend, dt)
