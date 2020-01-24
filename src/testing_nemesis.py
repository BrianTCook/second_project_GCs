from amuse.lab import *
from amuse.ext.bridge import bridge
from amuse.couple.bridge import CalculateFieldForParticles
from amuse.ic.kingmodel import new_king_model
from amuse.io import write_set_to_file, read_set_from_file

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

class GalacticCenterGravityCode(object):
    def __init__(self, R, M, alpha):
        self.radius=R
        self.mass=M
        self.alpha=alpha

    def get_gravity_at_point(self, eps, x, y, z):
        r2 = x**2+y**2+z**2
        r = r2**0.5
        m = self.mass*(r/self.radius)**self.alpha
        fr = constants.G*m/r2
        ax=-fr*x/r
        ay=-fr*y/r
        az=-fr*z/r
        return ax,ay,az
    
    def circular_velocity(self, r):
        m=self.mass*(r/self.radius)**self.alpha
        vc=(constants.G*m/r)**0.5
        return vc
    
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
    
    converter = nbody_system.nbody_to_si(M,R)
    bodies = new_king_model(N, W0, convert_nbody=converter)
    
    #sub_worker in Nemesis
    if code_name == 'Nbody':
        
        code = Hermite(converter)
        
    if code_name == 'tree':
    
        code = BHTree(converter)

    for name,value in parameters:
        setattr(code.parameters, name, value)
    code.particles.add_particles(bodies)
    
    return bodies, code

def orbiter_not_nemesis(orbiter_name, code_name, Nstars, W0, Mcluster, Rcluster, dBinary):
    
    converter = nbody_system.nbody_to_si(Mcluster, Rcluster)
    
    if orbiter_name == 'SingleStar':
        
        bodies = Particles(1)
        bodies.mass = [ 1|units.MSun ]
        
        #sub_worker in Nemesis, should not matter for SingleStar
        if code_name == 'Nbody':
            
            code = Hermite(converter)
            
        if code_name == 'tree':
        
            code = BHTree(converter)
        
        return bodies, code
        
    if orbiter_name == 'SingleCluster':
        
        bodies, code = make_king_model_cluster(Nstars, W0, Mcluster, Rcluster, code_name)
        
        '''
        need to initialize initial phase space coordinates with AGAMA
        '''
        
        return bodies, code
        
    if orbiter_name == 'BinaryCluster':
        
        bodies_one, code_one = make_king_model_cluster(Nstars, W0, Mcluster, Rcluster, code_name)
        bodies_one, code_two = make_king_model_cluster(Nstars, W0, Mcluster, Rcluster, code_name)
        
        '''
        need to initialize initial phase space coordinates with AGAMA
        '''
        
        mass_one, mass_two = bodies_one.mass.sum(), bodies_two.mass.sum()
        total_mass = mass_one + mass_two
        
        dBinary, vBinary = getxv(converter, total_mass, dBinary, eccentricity=0)
        
        for star in bodies_one:
            star.position += dBinary * mass_one/total_mass
            star.velocity += dBinary * mass_one/total_mass
            
        for star in bodies_two:
            star.position -= dBinary * mass_two/total_mass
            star.velocity -= dBinary * mass_two/total_mass
        
        bodies = Particles(0)
        bodies.add_particles(bodies_one)
        bodies.add_particles(bodies_two)
        
        return bodies, code_one, code_two #need to be different so they're bridged

def orbiter_nemesis(orbiter_name, code_name):
    
    return None, None

def gravity_code_setup(orbiter_name, code_name, galaxy_code):
    
    '''
    will need to ask SPZ if he meant for field, orbiter to be separate in non
    Nemesis gravity solvers?
    '''
    
    if code_name != 'Nemesis':
        
        gravity = bridge()
        
        if code_name != 'BinaryCluster':
            
            orbiter_bodies, orbiter_code = orbiter_not_nemesis(orbiter_name)
            
        if code_name == 'BinaryCluster':
            
            orbiter_bodies, orbiter_code_one, orbiter_code_two = orbiter_not_nemesis(orbiter_name)
    
        gravity.particles.add_particles(orbiter_bodies)
    
        #bridges clusters properly, independent of how many there are because
        #orbiter_code is not defined in binary case
        try:
            gravity.add_system(orbiter_code, galaxy_code)
        except:
            gravity.add_system(orbiter_code_one, galaxy_code)
            gravity.add_system(orbiter_code_two, galaxy_code)
            gravity.add_system(orbiter_code_one, orbiter_code_two)
            gravity.add_system(orbiter_code_two, orbiter_code_one)

    return gravity

def simulation(orbiter_name, code_name, Mgal, Rgal, alpha, 
               Nstars, W0, Mcluster, Rcluster, dBinary, tend, dt):
    
    galaxy_code = GalacticCenterGravityCode(Rgal, Mgal, alpha)
    gravity = gravity_code_setup(orbiter_name, code_name, galaxy_code)
    
    Ntotal = len(gravity.particles)
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr), dt.value_in(units.Myr))
    sim_times = [ t|units.Myr for t in sim_times_unitless ]
    
    energies, mean_radial_coords, mean_speeds, clock_times = [], [], [], []
    
    for j, t in enumerate(sim_times):

        energy = gravity.kinetic_energy + gravity.potential_energy
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

        filename = code_name + '_' + orbiter_name + '_data.hdf5'
	    write_set_to_file(gravity.particles, filename, "hdf5")
        
    gravity.stop()
    
    np.savetxt(code_name + '_' + orbiter_name + '_energies.txt', energies)
    np.savetxt(code_name + '_' + orbiter_name + '_mean_radial_coords.txt', mean_radial_coords)
    np.savetxt(code_name + '_' + orbiter_name + '_mean_speeds.txt', mean_speeds)
    np.savetxt(code_name + '_' + orbiter_name + '_clock_times.txt', clock_times)
    
    return 0

def plotting_things(orbiter_names, code_names, tend, dt):
    
    filename = code_name + '_' + orbiter_name + '_data.hdf5'
    star_data = read_set_from_file(filename)
    
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr), dt.value_in(units.Myr))
    
    npanels_x = len(orbiter_names)
    
    '''
    for i, orbiter_name in enumerate(orbiter_names):
        plt.subplot(i)
        
        for code_name in code_names:
            
            energies = np.loadtxt(code_name + '_' + orbiter_name + '_energies.txt')
            mean_radial_coords = np.savetxt(code_name + '_' + orbiter_name + '_mean_radial_coords.txt')
            mean_speeds = np.savetxt(code_name + '_' + orbiter_name + '_mean_speeds.txt')
            clock_times = np.savetxt(code_name + '_' + orbiter_name + '_clock_times.txt')
            plt.plot()
    
    '''
    
    return 0

if __name__ in '__main__':
    
    Mgal, Rgal, alpha = 1.6e10|units.MSun, 1000.|units.parsec, 1.2
    Nstars, W0cluster, Mcluster, Rcluster = 40, 1.5, 100.|units.MSun, 1.|units.parsec
    dBinary = 10.|units.parsec
    tend, dt = 100.|units.Myr, 1.|units.Myr
    
    orbiter_names = [ 'SingleStar', 'SingleCluster', 'BinaryCluster' ]
    code_names = [ 'Nbody', 'tree' ] # 'Nemesis'
    
    for orbiter_name in orbiter_names:
        for code_name in code_names:
            simulation(orbiter_name, code_name)
            
    plotting_things(orbiter_names, code_names)