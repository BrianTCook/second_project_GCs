from amuse.lab import *
from amuse.ext.bridge import bridge
from amuse.couple.bridge import CalculateFieldForParticles
from amuse.ic.kingmodel import new_king_model
#from amuse.community.phiGRAPE.interface import PhiGRAPE

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import numpy as np
import time
import os

#from globular_clusters_in_galaxy import parent_worker, sub_worker, py_worker
from nemesis import Nemesis, HierarchicalParticles, system_type

#Circumvent a problem with using too many threads on OpenMPI
#os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "yes"

#chapter 7 AMUSE textbook

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

def rand():
    
    return 2.*np.random.rand() - 1.
    
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

def gravity_code_setup(gravity_solver_str, galaxy_code,
		       cluster_codes, cluster_bodies_list):

    if gravity_solver_str == 'Brute':

        gravity = bridge()

        for i, cluster_code in enumerate(cluster_codes):  
                
        	other_clusters = cluster_codes[:i] + cluster_codes[i+1:]
                other_things = tuple(other_clusters) + (galaxy_code,)

                #bridges each cluster with the bulge, not the other way around though
                gravity.add_system(cluster_code, other_things)   

	return gravity     

    if gravity_solver_str == 'Nemesis':

        stars_all = Particles(0)
        
        for cluster_bodies in cluster_bodies_list:
            stars_all.add_particles(cluster_bodies)

        parts = HierarchicalParticles(stars_all)
        
        converter_parent = nbody_system.nbody_to_si(Mgal, Rgal)
        dt = smaller_nbody_power_of_two(0.1 | units.Myr, converter_parent)
        dt_nemesis = dt
        dt_bridge = 0.01 * dt
        
        nemesis = Nemesis( parent_worker, sub_worker, py_worker)
        nemesis.timestep = dt
        nemesis.distfunc = dt_bridge #timestep
        nemesis.threshold = dt_nemesis
        nemesis.radius = 4.|units.parsec #radius
        nemesis.commit_parameters()
        nemesis.particles.add_particles(parts)
        nemesis.commit_particles()

        channel_to_nemesis = stars_all.new_channel_to(nemesis.particles.all())

        #gravity = bridge.Bridge(use_threading=False)
        gravity = bridge()
        gravity.add_system(nemesis, (galaxy_code,) )
        gravity.timestep = dt_bridge
        
    	return gravity

def main(Rgal, Mgal, alpha, gravity_solvers, Nclusters, Nstars, W0, M,
     R, Rinit, parameters, t_end, dt):

    #set up clusters
    cluster_bodies_and_codes = [ make_king_model_cluster(Nstars, W0, M, R, parameters) 
                                 for i in range(Nclusters) ] 
    
    cluster_bodies_list = [ cbc[0] for cbc in cluster_bodies_and_codes ]
    cluster_codes = [ cbc[1] for cbc in cluster_bodies_and_codes ]

    bodies = Particles(0)
    for cluster_bodies in cluster_bodies_list:
	bodies.add_particles( cluster_bodies )

    star_colors = []    

    #just the galactic bulge
    galaxy_code = GalacticCenterGravityCode(Rgal, Mgal, alpha)

    for i, cluster_code in enumerate(cluster_codes):   

        stars = cluster_code.particles.copy()
        cluster_color = np.random.rand(3,)

        for j in range(len(stars)):

            star_colors.append(cluster_color)

        xrand, yrand, zrand = rand(), rand(), rand()
        vxrand, vyrand, vzrand = rand(), rand(), rand()

        stars.x += xrand*Rinit # x in (-R_init, R_init)
        stars.y += yrand*Rinit
        stars.z += zrand*Rinit

        R = Rinit*np.sqrt(xrand**2 + yrand**2 + zrand**2)

        stars.vx = vxrand*galaxy_code.circular_velocity(R)
        stars.vy = vyrand*galaxy_code.circular_velocity(R)
        stars.vz = vzrand*galaxy_code.circular_velocity(R)

        channel = stars.new_channel_to(cluster_code.particles)
        channel.copy_attributes(['x','y','z','vx','vy','vz'])

    gravity_solver_info = []

    for gravity_solver_str in gravity_solvers:

        gravity = gravity_code_setup(gravity_solver_str, galaxy_code,
                                     cluster_codes, cluster_bodies_list)
        
        sim_times_unitless = np.arange(0., t_end.value_in(units.Myr), dt.value_in(units.Myr))
        sim_times = [ t|units.Myr for t in sim_times_unitless]
    
        t0 = time.time()
        clock_times = []
        mean_radial_coords = []
        mean_speeds = []
    
        #total number of stars in the simulation
        Ntotal = Nclusters*Nstars
    
        for i, t in enumerate(sim_times):
    
            clock_time = time.time()-t0
    
            if i%5 == 0:
                print('i=%i, time = %.02f seconds'%(i, clock_time))
    
            clock_times.append(clock_time)
   
            #for figures 3 through 6, November 24

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
    
            
            #for .gif of orbits
            
            plt.figure()
            plt.scatter(x, y, marker='*', s=2, color=star_colors)
    
            Rinit_in_pc = Rinit.value_in(units.parsec)
            plt.xlim(-4*Rinit_in_pc, 4*Rinit_in_pc)
            plt.ylim(-4*Rinit_in_pc, 4*Rinit_in_pc)
            plt.annotate(gravity_solver_str, xy=(0.8, 0.8), xycoords='axes fraction')
    
            plt.xlabel('$x$ (pc)', fontsize=12)
            plt.ylabel('$y$ (pc)', fontsize=12)
            plt.title('Time: %.02f Myr'%(t.value_in(units.Myr)))
            plt.tight_layout()
            plt.savefig('frame_%s_%s.png'%(str(i).rjust(4, '0'), gravity_solver_str))
            plt.close()        

            gravity.evolve_model(t, timestep=dt)

        gravity_solver_info.append([gravity_solver_str, clock_times,
                                    mean_radial_coords, mean_speeds])
    
    cluster_code.stop()
    
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

if __name__ == '__main__':

    Mgal, Rgal, alpha = 1.6e10|units.MSun, 1000.|units.parsec, 1.2
    Nclusters = 4
    Nstars, W0cluster, Mcluster, Rcluster = 20, 1.5, 100.|units.MSun, 1.|units.parsec
    Rinit = 1000.|units.parsec
    parameters = [('epsilon_squared', 0.01|(units.parsec**2))]
    t_end, dt = 50.|units.Myr, 1.|units.Myr

    gravity_solvers = ['Brute' ] #'Nemesis',

    main(Rgal, Mgal, alpha, gravity_solvers, Nclusters, Nstars, W0cluster,
         Mcluster, Rcluster, Rinit, parameters, t_end, dt)
