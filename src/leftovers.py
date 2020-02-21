#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 12:13:07 2020

leftover stuff, not sure if i will need it again

@author: BrianTCook
"""

'''
if orbiter_name == 'BinaryCluster':
    
    list_of_orbiters = [ orbiter(orbiter_name, code_name, Mgalaxy, Rgalaxy, sepBinary,
                                 rvals, phivals, zvals, masses, i) for i in range(Norbiters) ]
    
    orbiter_bodies, orbiter_code_one, orbiter_code_two = list_of_orbiters[0]
'''

'''
if orbiter_name == 'BinaryCluster':
    
    list_of_orbiters = [ orbiter(orbiter_name, code_name, Mgalaxy, Rgalaxy, sepBinary,
                                 rvals, phivals, zvals, masses, i) for i in range(Norbiters) ]
    
    for i in range(Norbiters):
        
        orbiter_bodies, orbiter_code_one, orbiter_code_two = list_of_orbiters[0]
        
        gravity.add_system(orbiter_code_one, (orbiter_code_two, galaxy_code))
        gravity.add_system(orbiter_code_two, (orbiter_code_one, galaxy_code))
'''
        
if orbiter_name == 'BinaryCluster':
        
        bodies_one, code_one, _ = star_cluster(rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, index, code_name)
        bodies_two, code_two, _ = star_cluster(rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, index, code_name)
            
        #initialize binary system
        mass_one, mass_two = bodies_one.mass.sum(), bodies_two.mass.sum()
        total_mass = mass_one + mass_two
        
        dBinary, vBinary = getxv(converter, total_mass, sepBinary, e=0)
        print('dBinary is', dBinary)
        print('position adjustment for one: ', dBinary * mass_one/total_mass)
        print('position adjustment for two: ', -dBinary * mass_two/total_mass)
        print('----')
        print('vBinary is', vBinary)
        print('velocity adjustment for one: ', vBinary * mass_one/total_mass)
        print('velocity adjustment for two: ', -vBinary * mass_two/total_mass)
        
        bodies_one.position += dBinary * mass_one/total_mass
        bodies_one.velocity += vBinary * mass_one/total_mass

        bodies_two.position -= dBinary * mass_two/total_mass
        bodies_two.velocity -= vBinary * mass_two/total_mass
        
        all_bodies = Particles(0)
        all_bodies.add_particles(bodies_one)
        all_bodies.add_particles(bodies_two)
        
        return all_bodies, code_one, code_two #need to be different so they're bridged
