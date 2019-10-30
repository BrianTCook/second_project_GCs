from amuse.lab import*
from amuse.ext.bridge import bridge
from amuse.ic.kingmodel import new_king_model
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np


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

def make_king_model_cluster(nbodycode, N, W0, M, R, parameters=[]):

	'''
	sets up a cluster with mass M and radius R
	which nbodycode would you like to use?
	'''
	
	converter = nbody_system.nbody_to_si(M,R)
	bodies = new_king_model(N, W0, convert_nbody=converter)
	code = nbodycode(converter)

	for name,value in parameters:
		setattr(code.parameters, name, value)
	code.particles.add_particles(bodies)
	return code
	
def main(Rgal, Mgal, alpha, nbodycode, N, W0, M,
         R, Rinit, parameters, t_end, dt):

	galaxy_code = GalacticCenterGravityCode(Rgal, Mgal, alpha)
	cluster_code = make_king_model_cluster(nbodycode, N, W0,
					       M, R, parameters)

	gravity = bridge()
	gravity.add_system(cluster_code, (galaxy_code,))

	stars = cluster_code.particles.copy()
	stars.x += Rinit
	stars.vy = 0.8*galaxy_code.circular_velocity(Rinit)
	channel = stars.new_channel_to(cluster_code.particles)
	channel.copy_attributes(['x','y','z','vx','vy','vz'])

	times = np.arange(0., t_end.value_in(units.Myr), dt.value_in(units.Myr))
	times = [ t|units.Myr for t in times]

	for i, t in enumerate(times):

		x = gravity.particles.x.value_in(units.parsec)
		y = gravity.particles.y.value_in(units.parsec)

		plt.figure()
		plt.scatter(x, y, marker='.', color='k')

		plt.xlim(np.mean(x)-3*np.std(x), np.mean(x)+3.*np.std(x))
		plt.ylim(np.mean(y)-3*np.std(y), np.mean(y)+3.*np.std(y))

		plt.xlabel('$x$ (pc)', fontsize=12)
		plt.ylabel('$y$ (pc)', fontsize=12)
		plt.title('Time: %.02f Myr'%(t.value_in(units.Myr)))
		plt.savefig('frame_%s.png'%(str(i).rjust(4, '0')))
		plt.close()		

		gravity.evolve_model(t, timestep=dt)

	cluster_code.stop()

if __name__ == '__main__':
	
	Mgal, Rgal, alpha = 1.6e10|units.MSun, 1000.|units.parsec, 1.2
	nbodycode = BHTree
	Ncluster, W0cluster, Mcluster, Rcluster = 100, 1.5, 100.|units.MSun, 1.|units.parsec
	Rinit = 800.|units.parsec
	parameters = [('epsilon_squared', 0.01|(units.parsec**2))]
	t_end, dt = 100.|units.Myr, 1.|units.Myr

	main(Rgal, Mgal, alpha, nbodycode, Ncluster, W0cluster,
	     Mcluster, Rcluster, Rinit, parameters, t_end, dt)
