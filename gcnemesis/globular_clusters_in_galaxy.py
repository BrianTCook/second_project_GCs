import numpy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from amuse.lab import *

from amuse.ext.solarsystem import new_solar_system
from amuse.datamodel import Particle,Particles,ParticlesOverlay
from amuse.units import units,constants,nbody_system
from amuse.couple.bridge import CalculateFieldForParticles
from amuse.units.quantities import zero


from amuse.ext.basicgraph import UnionFind

from amuse.community.twobody.twobody import TwoBody
from amuse.community.hermite0.interface import Hermite
from amuse.community.mercury.interface import Mercury
from amuse.community.ph4.interface import ph4
#from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.huayno.interface import Huayno
from amuse.community.bhtree.interface import BHTree
from amuse.community.mi6.interface import MI6

from galpy.potential import MWPotential2014, to_amuse

from nemesis import Nemesis,HierarchicalParticles,system_type
from amuse.couple import bridge
from galactic_model import IntegrateOrbit

import logging
#logging.basicConfig(level=logging.DEBUG)

    
def smaller_nbody_power_of_two(dt, conv):
  nbdt=conv.to_nbody(dt).value_in(nbody_system.time)
  idt=numpy.floor(numpy.log2(nbdt))
  return conv.to_si( 2**idt | nbody_system.time)

markerstyles=["+","p","o","x"]
linestyles=["-",":","-.","--"]
colors=["r","g","b","y","c"]

from amuse.ext.galactics_model import new_galactics_model

def make_galaxy_model(N, M_galaxy, R_galaxy, Mmin, Mmax):
    converter = nbody_system.nbody_to_si(M_galaxy, R_galaxy)
    n_halo = 0
    n_bulge = 0
    n_disk = N
    """
    galaxy = new_galactics_model(n_halo,
                                 converter,
                                 do_scale=True,
                                 bulge_number_of_particles=n_bulge,
                                 disk_number_of_particles=n_disk)
    """
    galaxy = new_plummer_model(n_disk, converter)
    
    from amuse.lab import new_powerlaw_mass_distribution
    masses = new_powerlaw_mass_distribution(len(galaxy), Mmin, Mmax, -2.0)
    galaxy.mass = masses
    galaxy.radius = 10 | units.parsec
    galaxy.King_W0 = 3
    return galaxy

def make_plummer_galaxy_model(N, M_galaxy, R_galaxy):
    converter = nbody_system.nbody_to_si(M_galaxy, R_galaxy)
    model = new_plummer_model(N, converter)
    model.radius = 10 | units.parsec
    model.King_W0 = 3
    model.mass = 20 | units.MSun
    return model

def initialize_globular_clusters(cluster_population, nstars):

    mmean = new_kroupa_mass_distribution(1000).sum()/1000.
    stars = Particles(0)
    for ci in cluster_population:
        nstars = max(1, int(ci.mass/mmean))
        masses = new_kroupa_mass_distribution(nstars)        
        converter = nbody_system.nbody_to_si(masses.sum(), ci.radius)
        bodies = new_king_model(nstars, ci.King_W0, converter)
        bodies.parent = ci
        bodies.mass = masses
        bodies.position += ci.position
        bodies.velocity += ci.velocity
        bodies.name = "star"
        ci.mass = bodies.mass.sum()
        bodies.scale_to_standard(converter)
        stars.add_particles(bodies)
    return stars

def plot_figure(stars, nemesis, nstep, LL, filename):
    import matplotlib.cm as cm
    ss=nemesis.particles.all()
    f=pyplot.figure( figsize=(8,8))  
    ax = f.gca()
    circles=[]
    s = 6*ss.mass/ss.mass.max()
    m = ss.mass/ss.mass.max()
    #pyplot.plot(xx[-1],yy[-1],'b.',markersize=2, mew=1.5)
    if 'xy' in filename:
      for si in range(len(ss)):
        pyplot.plot(ss[si].x.value_in(units.kpc),ss[si].y.value_in(units.kpc),'b.',markersize=s[si], mew=s[si], alpha=0.5)
    else:
      for si in range(len(ss)):
        pyplot.plot(ss[si].x.value_in(units.kpc),ss[si].z.value_in(units.kpc),'b.',markersize=s[si], mew=s[si], alpha=0.5)

    for p in nemesis.particles:
        c='k'
        ls='solid'
        code_colors=dict(ph4='b',Mercury='r',Huayno='g',Hermite='y')
        code_ls=dict(ph4='dotted',Mercury='dashed',Huayno='dashdot',Hermite='solid')
        if nemesis.subcodes.has_key(p):
          c=code_colors[nemesis.subcodes[p].__class__.__name__] 
          ls=code_ls[nemesis.subcodes[p].__class__.__name__]
        if 'xy' in filename:
          x=p.x.value_in(units.kpc)
          y=p.y.value_in(units.kpc)
        else:
          x=p.x.value_in(units.kpc)
          y=p.z.value_in(units.kpc)
        r=p.radius.value_in(units.kpc)
        circles.append( pyplot.Circle((x,y),r,color=c,lw=0.8,ls=ls,fill=False) )
    for c in circles:
      ax.add_artist(c)  
                                
    pyplot.xlim(-1.2*LL,1.2*LL)
    pyplot.ylim(-1.2*LL,1.2*LL)
    pyplot.xlabel("kpc")
    #pyplot.savefig(filename+'%6.6i.png'%nstep,bbox_inches='tight')
    pyplot.show()
    f.clear()
    pyplot.close(f)
  
def globular_clusters(N=10, L=6.| units.kpc, dv=1.0 | units.kms):

  M_galaxy = 1.0e+10 | units.MSun
  R_galaxy = 4.5 | units.kpc
  Mmin = 5 | units.MSun
  Mmax = 20 | units.MSun
  cluster_population = make_galaxy_model(N, M_galaxy, R_galaxy, Mmin, Mmax)
  stars = initialize_globular_clusters(cluster_population, N)
  print("Stars:", len(stars), stars.mass.sum().in_(units.MSun), stars.mass.max().in_(units.MSun), stars.mass.mean().in_(units.MSun))

  stellar = SeBa()
  stellar.particles.add_particles(stars)
  channel_from_stellar = stellar.particles.new_channel_to(stars)
  
  conv=nbody_system.nbody_to_si(M_galaxy, R_galaxy)
  conv_sub=nbody_system.nbody_to_si(1000.| units.MSun, 10.| units.parsec)

  dt=smaller_nbody_power_of_two(0.1 | units.Myr, conv)
  dt_nemesis=dt
  print(dt.in_(units.Myr))
  dt_bridge = 0.01*dt

  #dt_param=0.02
  dt_param=0.1
  LL=L.value_in(units.kpc)

  #converter = nbody_system.nbody_to_si(M_galaxy, R_Galaxy)

  OS= 20 |(units.kms/units.kpc)
  OB= 40 |(units.kms/units.kpc)
  A= 1300 |(units.kms**2/units.kpc)
  M= 1.4e10 |units.MSun
  m=2    

  phi_bar, phi_sp= -0.34906, -0.34906
  galaxy= IntegrateOrbit(t_end= 10|units.Myr,
                         dt_bridge= dt_bridge, 
                         phase_bar= phi_bar, phase_spiral= phi_sp, 
                         omega_spiral= OS, omega_bar= OB, 
                         amplitude= A, m=m, mass_bar= M )
  MWG = to_amuse(MWPotential2014, t=0.0, tgalpy=0.0, reverse=False, ro=None, vo=None) #galaxy.galaxy()
  
  def radius(sys,eta=dt_param,_G=constants.G):
    radius=((_G*sys.total_mass()*dt**2/eta**2)**(1./3.))
    radius = radius*((len(sys)+1)/2.)**0.75
    #print "radius=", radius.in_(units.parsec)
    return 100*radius

  def timestep(ipart,jpart, eta=dt_param/2,_G=constants.G):
    dx=ipart.x-jpart.x  
    dy=ipart.y-jpart.y
    dz=ipart.z-jpart.z
    dr2=dx**2+dy**2+dz**2
    dr=dr2**0.5
    dr3=dr*dr2
    mu=_G*(ipart.mass+jpart.mass)
    tau=eta/2./2.**0.5*(dr3/mu)**0.5
    return tau
    
  numpy.random.seed(7654304)

  parts=HierarchicalParticles(stars)

#  ss=new_solar_system()[[0,5,6,7,8]]
#  parts.assign_subsystem(ss,parts[0])

  def parent_worker():
    code=Hermite(conv)
    code.parameters.epsilon_squared=0.| units.kpc**2
    code.parameters.end_time_accuracy_factor=0.
    #code.parameters.dt_param=0.001
    code.parameters.dt_param=0.1
    print(code.parameters.dt_dia.in_(units.yr))
    return code
  
  
  def sub_worker(parts):
    mode=system_type(parts)
    if mode=="twobody":
      #code=TwoBody(conv_sub)
      code=ph4(conv_sub)
    elif mode=="solarsystem":
      #code=Mercury(conv_sub)
      code=Huayno(conv_sub)
    elif mode=="nbody":
      code=BHTree(conv_sub)
      #code.parameters.inttype_parameter=code.inttypes.SHARED4
    return code
      
  def py_worker():
    code=CalculateFieldForParticles(gravity_constant = constants.G)
    return code
      
  nemesis=Nemesis( parent_worker, sub_worker, py_worker)
  nemesis.timestep=dt
  nemesis.distfunc=timestep
  nemesis.threshold=dt_nemesis
  nemesis.radius=radius
  nemesis.commit_parameters()
  nemesis.particles.add_particles(parts)
  nemesis.commit_particles()

  channel_to_nemesis = stars.new_channel_to(nemesis.particles.all())

  gravity = bridge.Bridge(use_threading=False)
  gravity.add_system(nemesis, (MWG,) )
  gravity.timestep = dt_bridge
  
  tend=100.0 | units.Myr
  t=0|units.yr
  dtdiag=dt*2
  
  time=[0.]

  allparts=nemesis.particles.all()
  E=allparts.potential_energy()+allparts.kinetic_energy()
  E1=nemesis.potential_energy+nemesis.kinetic_energy

  com=allparts.center_of_mass()
  mom=allparts.total_momentum()
  ang=allparts.total_angular_momentum()
    
  E0=E
  A0=(ang[0]**2+ang[1]**2+ang[2]**2)**0.5 
  P0=mom[0].value_in(units.MSun*units.kms)
  totalE=[ 0.]
  totalA=[ 0.]
  totalP=[ 0.]
  totaldM=[ 0. ] | units.MSun
  
  ss=nemesis.particles.all()
  x=(ss.x).value_in(units.kpc)
  xx=[x]
  y=(ss.y).value_in(units.kpc)
  yy=[y]
  z=(ss.z).value_in(units.kpc)
  zz=[z]
  
  nstep=0
  while t< tend-dtdiag/2:
    t+=dtdiag
    gravity.evolve_model(t)  
    #nemesis.evolve_model(t)  
    print(t.in_(units.Myr))
    print(len(nemesis.particles))

    m = stars.mass.sum()
    stellar.evolve_model(t)
    channel_from_stellar.copy()
    dm = m-stars.mass.sum()
    totaldM.append(totaldM[-1]+dm)
    print("dM=", totaldM[-1].in_(units.MSun), dm.in_(units.MSun))
    channel_to_nemesis.copy()
    
    time.append( t.value_in(units.yr) )

    #filename = 'GGC_i%6.6i.h5'%nstep
    #write_set_to_file(stars, filename, "hdf4")
    
    allparts=nemesis.particles.all()
    E=allparts.potential_energy()+allparts.kinetic_energy()
    E1=nemesis.potential_energy+nemesis.kinetic_energy
    
    ang=allparts.total_angular_momentum()
    mom=allparts.total_momentum()
    A=(ang[0]**2+ang[1]**2+ang[2]**2)**0.5
    P=mom[0].value_in(units.MSun*units.kms)
    totalE.append(abs((E0-E)/E0))
    totalA.append(abs((A0-A)/A0))
    totalP.append(abs(P0-P))
    print("dE=", totalE[-1],(E-E1)/E)
#    print allparts.potential_energy(),nemesis.potential_energy
  
    ss=nemesis.particles.all()
    x=(ss.x).value_in(units.kpc)
    y=(ss.y).value_in(units.kpc)
    y=(ss.z).value_in(units.kpc)
    lowm=numpy.where( ss.mass.value_in(units.MSun) < 1)[0]
    highm=numpy.where( ss.mass.value_in(units.MSun) >= 1)[0]
    print("N=", len(ss), len(highm), len(lowm))
    
    xcm=nemesis.particles.x.value_in(units.kpc)
    ycm=nemesis.particles.y.value_in(units.kpc)
    r=(nemesis.particles.radius).value_in(units.kpc)
    
    xx.append(x)
    yy.append(y)
    zz.append(z)
    key=ss.key

    plot_figure(stars, nemesis, nstep, LL, 'xy_i')
    plot_figure(stars, nemesis, nstep, LL, 'xz_i')

    nstep+=1

    
  time=numpy.array(time)
  totalE=numpy.array(totalE)
  totalA=numpy.array(totalA)
  totalP=numpy.array(totalP)
  xx=numpy.array(xx)
  yy=numpy.array(yy)

  f=pyplot.figure( figsize=(8,8))  
  pyplot.semilogy(time,totalE,'r')
  pyplot.semilogy(time,totalA,'g')
  pyplot.semilogy(time,totalP,'b')
  pyplot.savefig('dEdA.png')  
    
  f=pyplot.figure( figsize=(8,8))  
  for i in range(len(xx[0])):
    pyplot.plot(xx[:,i],yy[:,i],colors[i%len(colors)]+linestyles[i%len(linestyles)])
  pyplot.xlim(-LL,LL)
  pyplot.ylim(-LL,LL)
  pyplot.savefig('all-xy.png')  


if __name__=="__main__":
  globular_clusters()

  
