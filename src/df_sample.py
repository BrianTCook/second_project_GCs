from amuse.lab import *
from galpy.df import quasiisothermaldf, dehnendf
from galpy.potential import MWPotential2014, to_amuse
from galpy.util import bovy_conversion
from galpy.actionAngle import actionAngleStaeckel

import numpy as np

aAS = actionAngleStaeckel(pot=MWPotential2014, delta=0.45, c=True)
qdfS = quasiisothermaldf(1./3., 0.2, 0.1, 1., 1., pot=MWPotential2014, aA=aAS, cutcounter=True)
dfc = dehnendf(beta=0.)

'''
sample from Dehnen distribution function to get cylindrical coordinates
where |r| <= 1.0 kpc
'''

nsamples = 4000
o = dfc.sample(n=nsamples, returnOrbit=True, nphi=1, rrange=[0.0, 1.0])

rvals_spherical = [ e.R() for e in o ]
thetavals_spherical = [ np.arccos(1 - 2*np.random.random()) for i in range(nsamples) ]

rvals_cylindrical = [ rvals_spherical[i]*np.sin(thetavals_spherical[i]) for i in range(nsamples) ]
phivals_cylindrical = 2*np.pi * np.random.random(nsamples)
zvals_cylindrical = [ rvals_spherical[i]*np.cos(thetavals_spherical[i]) for i in range(nsamples) ]

    
#using Staeckel
aAS = actionAngleStaeckel(pot=MWPotential2014, delta=0.45, c=True)
qdfS = quasiisothermaldf(1./3., 0.2, 0.1, 1., 1., pot=MWPotential2014, aA=aAS, cutcounter=True)

#220 km/s at 8 kpc, convert back to km/s
to_kms = bovy_conversion.velocity_in_kpcGyr(220., 8.) * 0.9785

for r, z in zip(rvals_cylindrical, zvals_cylindrical):
    
    vr_init, vphi_init, vz_init = qdfS.sampleV(Rcoord, Zcoord, n=1)[0,:]
    
    vrvals_cylindrical.append( vr_init * to_kms )
    vphivals_cylindrical.append( vphi_init * to_kms )
    vzvals_cylindrical.append( vz_init * to_kms )

np.savetxt('dehnen_rvals.txt', rvals_cylindrical)
np.savetxt('dehnen_phivals.txt', phivals_cylindrical)
np.savetxt('dehnen_zvals.txt', zvals_cylindrical)

np.savetxt('bovy_vrvals.txt', vrvals_cylindrical)
np.savetxt('bovy_vphivals.txt', vphivals_cylindrical)
np.savetxt('bovy_vzvals.txt', vzvals_cylindrical)
