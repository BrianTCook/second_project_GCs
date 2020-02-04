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
from galpy.actionAngle import actionAngleStaeckel

from phase_space_mapping import maps

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
        code.particles.add_particles(bodies)
        
    if code_name == 'tree':
    
        code = BHTree(converter)
        code.particles.add_particles(bodies)
        
    if code_name == 'nemesis':
        
        parts = HierarchicalParticles(bodies)

        converter_parent = nbody_system.nbody_to_si(Mcluster, Rcluster)
        
        dt = smaller_nbody_power_of_two(0.1 | units.Myr, converter_parent)
        dt_nemesis = dt
        dt_bridge = 0.01 * dt
        
        nemesis = Nemesis(parent_worker, sub_worker, py_worker)
        nemesis.timestep = dt
        nemesis.distfunc = distance_function
        nemesis.threshold = dt_nemesis
        nemesis.radius = radius
        
        nemesis.commit_parameters()
        nemesis.particles.add_particles(parts)
        nemesis.commit_particles()
        
        code = nemesis
        
    for name,value in parameters:
        setattr(code.parameters, name, value)
    
    return bodies, code

def parent_worker():
    converter_parent = nbody_system.nbody_to_si(Mcluster, Rcluster)
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

def orbiter(orbiter_name, code_name, Rcoord, Zcoord, phicoord,
            vr_init, vphi_init, vz_init, Nstars, W0, Mcluster, Rcluster, sepBinary):

    converter = nbody_system.nbody_to_si(Mcluster, Rcluster)
    
    '''
    takes in R, Z value
    returns VR, Vphi, VZ values
    should get appropriate 6D initial phase space conditions
    '''
    
    #convert from galpy/cylindrical to AMUSE/Cartesian units
    x_init = Rcoord*np.cos(phicoord) | units.kpc
    y_init = Rcoord*np.sin(phicoord) | units.kpc
    z_init = Zcoord | units.kpc
    
    #vphi = R \dot{\phi}? assuming yes for now
    vx_init = (vr_init*np.cos(phicoord) - vphi_init*np.sin(phicoord)) | units.kms
    vy_init = (vr_init*np.sin(phicoord) + vphi_init*np.cos(phicoord)) | units.kms
    vz_init = vz_init | units.kms
    
    pos_vec, vel_vec = (x_init, y_init, z_init)|units.kpc, (vx_init, vy_init, vz_init)|units.kms
    
    print('-----------------------------')
    print('-----------------------------')
    print('-----------------------------')
    print('orbiter_name: ', orbiter_name)
    print('code_name: ', code_name)
    print('initial spatial coordinates (Cartesian): ', pos_vec)
    print('initial velocity coordinates: ', vel_vec)
    print('-----------------------------')
    print('-----------------------------')
    print('-----------------------------')
    
    if orbiter_name == 'SingleStar':
        
        bodies = Particles(1)
        
        bodies[0].mass = 1|units.MSun
        
        #right place in phase space
        bodies[0].x = x_init
        bodies[0].y = y_init
        bodies[0].z = z_init
        bodies[0].vx = vx_init
        bodies[0].vy = vy_init
        bodies[0].vz = vz_init
        
        #sub_worker in Nemesis, should not matter for SingleStar
        if code_name == 'Nbody':
            
            code = Hermite(converter)
            code.particles.add_particles(bodies)
            
        if code_name == 'tree':
        
            code = BHTree(converter)
            code.particles.add_particles(bodies)
            
        if code_name == 'nemesis':
            
            parts = HierarchicalParticles(bodies)

            converter_parent = nbody_system.nbody_to_si(Mcluster, Rcluster)
            
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
            
            code = nemesis
        
        return bodies, code
        
    if orbiter_name == 'SingleCluster':
        
        bodies, code = make_king_model_cluster(Nstars, W0, Mcluster, Rcluster, code_name)
        
        '''
        need to initialize initial phase space coordinates with AGAMA or galpy
        '''
        
        stars = code.particles.copy()
        
        #initializing phase space coordinates
        stars.x += x_init
        stars.y += y_init
        stars.z += z_init
        stars.vx += vx_init
        stars.vy += vy_init
        stars.vz += vz_init
        
        return stars, code
        
    if orbiter_name == 'BinaryCluster':
        
        bodies_one, code_one = make_king_model_cluster(Nstars, W0, Mcluster, Rcluster, code_name)
        bodies_two, code_two = make_king_model_cluster(Nstars, W0, Mcluster, Rcluster, code_name)
        
        stars_one = code_one.particles.copy()
        stars_two = code_two.particles.copy()
        
        #initializing phase space coordinates
        stars_one.x += x_init
        stars_one.y += y_init
        stars_one.z += z_init
        stars_one.vx += vx_init
        stars_one.vy += vy_init
        stars_one.vz += vz_init
        
        stars_two.x += x_init
        stars_two.y += y_init
        stars_two.z += z_init
        stars_two.vx += vx_init
        stars_two.vy += vy_init
        stars_two.vz += vz_init
            
        #initialize binary system
        mass_one, mass_two = bodies_one.mass.sum(), bodies_two.mass.sum()
        total_mass = mass_one + mass_two
        
        dBinary, vBinary = getxv(converter, total_mass, sepBinary, e=0)
        print('dBinary is', dBinary)
        print('vBinary is', vBinary)
        
        stars_one.position += dBinary * mass_one/total_mass
        stars_one.velocity += vBinary * mass_one/total_mass

        stars_two.position -= dBinary * mass_two/total_mass
        stars_two.velocity -= vBinary * mass_two/total_mass
        
        all_stars = Particles(0)
        all_stars.add_particles(stars_one)
        all_stars.add_particles(stars_two)
        
        return all_stars, code_one, code_two #need to be different so they're bridged

def gravity_code_setup(orbiter_name, code_name, galaxy_code, Rcoord, Zcoord, phicoord,
                       vr_init, vphi_init, vz_init, Nstars, W0, Mcluster, Rcluster, sepBinary):
    
    '''
    will need to ask SPZ if he meant for field, orbiter to be separate in non
    Nemesis gravity solvers?
    '''
    
    if code_name != 'nemesis':
        
        gravity = bridge()
        
        if orbiter_name != 'BinaryCluster':
            
            orbiter_bodies, orbiter_code = orbiter(orbiter_name, code_name, Rcoord, Zcoord, phicoord,
                                                   vr_init, vphi_init, vz_init, Nstars, W0, Mcluster, Rcluster, sepBinary)
            
            gravity.add_system(orbiter_code, (galaxy_code,))

        if orbiter_name == 'BinaryCluster':
            
            orbiter_bodies, orbiter_code_one, orbiter_code_two = orbiter(orbiter_name, code_name, Rcoord, Zcoord, phicoord,
                                                                         vr_init, vphi_init, vz_init, Nstars, W0, Mcluster, Rcluster, sepBinary)
    
            gravity.add_system(orbiter_code_one, (galaxy_code,))
            gravity.add_system(orbiter_code_two, (galaxy_code,))
            gravity.add_system(orbiter_code_one, (orbiter_code_two,))
            gravity.add_system(orbiter_code_two, (orbiter_code_one,))
            
    if code_name == 'nemesis':
        
        gravity = bridge()
        
        #just don't use orbiter_code here, just replace it with nemesis
        if orbiter_name != 'BinaryCluster':
            
            orbiter_bodies, orbiter_code = orbiter(orbiter_name, code_name, Rcoord, Zcoord, phicoord,
                                                   vr_init, vphi_init, vz_init, Nstars, W0, Mcluster, Rcluster, sepBinary)
            
        if orbiter_name == 'BinaryCluster':
            
            orbiter_bodies, orbiter_code_one, orbiter_code_two = orbiter(orbiter_name, code_name, Rcoord, Zcoord, phicoord,
                                                                         vr_init, vphi_init, vz_init, Nstars, W0, Mcluster, Rcluster, sepBinary)
            
        
        parts = HierarchicalParticles(orbiter_bodies)

        converter_parent = nbody_system.nbody_to_si(Mcluster, Rcluster)
        dt = smaller_nbody_power_of_two(0.1 | units.Myr, converter_parent)
        dt_nemesis = dt
        dt_bridge = 0.1 * dt
        
        nemesis = Nemesis( parent_worker, sub_worker, py_worker)
        nemesis.timestep = dt
        nemesis.distfunc = distance_function
        nemesis.threshold = dt_nemesis
        nemesis.radius = radius
        
        nemesis.commit_parameters()
        nemesis.particles.add_particles(parts)
        nemesis.commit_particles()

        gravity.add_system(nemesis, (galaxy_code,))

    return orbiter_bodies, gravity

def simulation(orbiter_name, code_name, potential, Rcoord, Zcoord, phicoord,  
               vr_init, vphi_init, vz_init, Nstars, W0, Mcluster, Rcluster, sepBinary, tend, dt):
    
    galaxy_code = to_amuse(potential, t=0.0, tgalpy=0.0, reverse=False, ro=None, vo=None)
    
    bodies, gravity = gravity_code_setup(orbiter_name, code_name, galaxy_code, Rcoord, Zcoord, phicoord, vr_init, vphi_init, vz_init, Nstars, W0, Mcluster, Rcluster, sepBinary)
    
    stars = gravity.particles.copy()
    
    channel = stars.new_channel_to(gravity.particles)
    channel.copy_attributes(['x','y','z','vx','vy','vz'])
    
    Ntotal = len(bodies)
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr), dt.value_in(units.Myr))
    sim_times = [ t|units.Myr for t in sim_times_unitless ]
    
    energies, mean_radial_coords, mean_speeds, clock_times = [], [], [], []
    
    t0 = time.time()
    
    #create an R^3 matrix to house phase space data for all particles
    phase_space_data = np.zeros((len(sim_times), 6, len(bodies)))
    
    for j, t in enumerate(sim_times):

        clock_times.append(time.time()-t0) #will be in seconds
        
        energy = gravity.kinetic_energy.value_in(units.J) + gravity.potential_energy.value_in(units.J)
        energies.append(energy)
        
        x = [ xx.value_in(units.kpc) for xx in gravity.particles.x ]
        y = [ yy.value_in(units.kpc) for yy in gravity.particles.y ]
        z = [ zz.value_in(units.kpc) for zz in gravity.particles.z ]
        
        vx = [ vxx.value_in(units.kms) for vxx in gravity.particles.vx ]
        vy = [ vyy.value_in(units.kms) for vyy in gravity.particles.vy ]
        vz = [ vzz.value_in(units.kms) for vzz in gravity.particles.vz ]
        
        for k, star in enumerate(gravity.particles):
            
            phase_space_data[j,0,k] = x[k]
            phase_space_data[j,1,k] = y[k]
            phase_space_data[j,2,k] = z[k]
            phase_space_data[j,3,k] = vx[k]
            phase_space_data[j,4,k] = vy[k]
            phase_space_data[j,5,k] = vz[k]
        
        xmean, ymean, zmean = np.sum(x)/Ntotal, np.sum(y)/Ntotal, np.sum(z)/Ntotal
        mean_rval = np.sqrt(xmean**2 + ymean**2 + zmean**2)
        mean_radial_coords.append(mean_rval)
        
        vxmean, vymean, vzmean = np.sum(vx)/Ntotal, np.sum(vy)/Ntotal, np.sum(vz)/Ntotal
        mean_speed = np.sqrt(vxmean**2 + vymean**2 + vzmean**2)
        mean_speeds.append(mean_speed)
                
        gravity.evolve_model(t)

        #filename = code_name + '_' + orbiter_name + '_data.hdf5'
        #write_set_to_file(gravity.particles, filename, "hdf5")
     
    try:
        gravity.stop()
    except:
        'gravity cannot be stopped!'
    
    np.save('time_data_%s_%s.npy'%(code_name, orbiter_name), sim_times_unitless)
    np.save('sixD_data_%s_%s.npy'%(code_name, orbiter_name), phase_space_data)
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
            
            try:
                energies = np.loadtxt(code_name + '_' + orbiter_name + '_energies.txt')
                axs[i].plot(sim_times_unitless, energies, label=code_name)
            except:
                print('oh no!')
            
        axs[i].legend(loc='best')
            
    plt.savefig('testing_nemesis_energy.png')
    plt.close()
    
    #mean radial coordinates
    
    fig, axs = plt.subplots(1, 3)

    for i, orbiter_name in enumerate(orbiter_names): 
        
        axs[i].set_title(orbiter_name)
        
        for code_name in code_names:
            
            try:
                mean_radial_coords = np.loadtxt(code_name + '_' + orbiter_name + '_mean_radial_coords.txt')
                axs[i].plot(sim_times_unitless, mean_radial_coords, label=code_name)
            except:
                print('oh no!')
            
        axs[i].legend(loc='best')
            
    plt.savefig('testing_nemesis_radialcoords.png')
    plt.close()
    
    #mean speeds
    
    fig, axs = plt.subplots(1, 3)

    for i, orbiter_name in enumerate(orbiter_names): 
        
        axs[i].set_title(orbiter_name)
        
        for code_name in code_names:
            
            try:
                mean_speeds = np.loadtxt(code_name + '_' + orbiter_name + '_mean_speeds.txt')
                axs[i].plot(sim_times_unitless, mean_speeds, label=code_name)
            except:
                print('oh no!')
            
        axs[i].legend(loc='best')
            
    plt.savefig('testing_nemesis_speeds.png')
    plt.close()
    
    #clock times
    
    fig, axs = plt.subplots(1, 3)

    for i, orbiter_name in enumerate(orbiter_names): 
        
        axs[i].set_title(orbiter_name)
        
        for code_name in code_names:
            
            try:
                clock_times = np.loadtxt(code_name + '_' + orbiter_name + '_clock_times.txt')
                axs[i].plot(sim_times_unitless, clock_times, label=code_name)
            except:
                print('oh no!')
            
        axs[i].legend(loc='best')
            
    plt.savefig('testing_nemesis_clocktimes.png')
    plt.close()
    
    return 0

if __name__ in '__main__':
    
    potential = MWPotential2014
    Rmin, Rmax = 1., 5. #in kpc
    Zmin, Zmax = 0.5, 2. #in kpc
    
    Rcoord = (Rmax-Rmin)*np.random.random() + Rmin
    Zcoord = (Zmax-Zmin)*np.random.random() + Zmin
    phicoord = 2*np.pi*np.random.random()
    
    #using Staeckel, whatever that means
    aAS = actionAngleStaeckel(pot=MWPotential2014, delta=0.45, c=True)
    qdfS = quasiisothermaldf(1./3., 0.2, 0.1, 1., 1., pot=MWPotential2014, aA=aAS, cutcounter=True)
    vr_init, vphi_init, vz_init = qdfS.sampleV(Rcoord, Zcoord, n=1)[0,:]
    
    #220 km/s at 8 kpc, convert back to km/s
    to_kms = bovy_conversion.velocity_in_kpcGyr(220., 8.) * 0.9785
    vr_init *= to_kms
    vphi_init *= to_kms
    vz_init *= to_kms
    
    Nstars, W0 = 25, 1.5 #cluster parameters
    Mcluster, Rcluster = float(Nstars)|units.MSun, 20.|units.parsec
    sepBinary = 80.|units.parsec
    tend, dt = 20.|units.Myr, 1.|units.Myr
    
    orbiter_names = [ 'SingleStar', 'SingleCluster', 'BinaryCluster' ]
    code_names = [ 'tree', 'Nbody', 'nemesis' ]
    
    for orbiter_name in orbiter_names:
        for code_name in code_names:
            
            print('')
            print(orbiter_name, code_name)
            print('')
            
            if orbiter_name == 'SingleStar' and code_name == 'nemesis':
                continue
            
            simulation(orbiter_name, code_name, potential, Rcoord, Zcoord, phicoord, 
                       vr_init, vphi_init, vz_init, Nstars, W0, Mcluster, Rcluster, sepBinary, tend, dt)
            maps(orbiter_name, code_name)
            
    plotting_things(orbiter_names, code_names, tend, dt)
