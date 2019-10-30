from amuse.lab import*
import numpy as np

#chapter 7 AMUSE textbook

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
	

galaxy_code = GalacticCenterGravityCode(Rgal, Mgal, alpha)
cluster_code = make_king_model_cluster

gravity = bridge.Bridge()
gravity.add_system(cluster_code, (galaxy_code,))
gravity.add_system(cluster1, (cluster2, cluster3, bulge, bar, disk))
gravity.timestep = 0.1|units.Myr

if __name__ in '__main__':
	Mgal, Rgal = 1.6e10|units.MSun, 1000.|units.parsec
	alpha = 1.2

