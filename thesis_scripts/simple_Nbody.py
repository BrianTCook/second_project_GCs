'''
m, x, xdot
are lists of masses, positions, and velocities, respectively
t, t_end, and dt are defined to have sensible time units for the N-body problem at hand

np.add/np.subtract from the NumPy library compute list additions/subtractions element-wise

'''

while t < t_end:

    for i in range(N):
        
        F_i = [0., 0., 0.] #initialize force vector as a list of length 3
        
        for j in range(N):
            
            if i != j:
            
                sep_ij = np.subtract(x[j], x[i]) #compute the separation vector between the ith and jth particle
                F_i += G*m[i]*m_[j] / mag(sep_ij)**3 * sep_ij #summand in eq. (2.1)
                
        xdot[i] += np.add(xdot[i], 1/m[i] * F_i * dt) #element wise addition
        x[i] += xdot[i] * dt #element wise addition
        
    t += dt
    