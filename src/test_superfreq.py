import numpy as np 
import matplotlib.pyplot as plt
from superfreq import SuperFreq

gravity_solver_str = 'Brute'
Nclusters_list = [1, 10]

plt.figure()

for Nclusters in Nclusters_list:

    data_dir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
    filename_time = 'time_data_Nclusters=%i_%s.npy'%(Nclusters, gravity_solver_str)
    filename_phase = 'sixD_data_Nclusters=%i_%s.npy'%(Nclusters, gravity_solver_str)
    
    sim_times = np.load(data_dir+filename_time)
    phase_space_data = np.load(data_dir+filename_phase)
    ntimes, ndim, nparticles = phase_space_data.shape
    print(ntimes, ndim, nparticles)
    
    #fundamental frequency stuff
    freqs_all = []
    for k in range(nparticles):
        t, w = sim_times, phase_space_data[:, :, k]
        sf = SuperFreq(t)
        fs = [ (w[:,i] * 1j * w[:,i+ndim//2]) for i in range(ndim//2) ]
        superfreq_result_object = sf.find_fundamental_frequencies(fs)
        freqs_all.append( superfreq_result_object.fund_freqs )
    
    freqs_all = [j for i in freqs_all for j in i] #concatenate the list
    freqs_all = np.sort(freqs_all)
    
    yvals, xvals = np.histogram(freqs_all)
    xvals_CDF = [ 0.5*(xvals[i]+xvals[i+1]) for i in range(len(xvals)-1) ] 
    
    ysum = sum(yvals)
    yvals = [ y/ysum for y in yvals ] #normalization
    
    yvals_CDF = np.cumsum(yvals)
    
    plt.plot(xvals_CDF, yvals_CDF, label='CDF, Nclusters = %i'%(Nclusters))
    
plt.xlabel('Frequency (action-angle variable units)', fontsize=12)
plt.legend(loc='best')
plt.savefig('fundamental_freqs_CDF_%s.png'%(gravity_solver_str))
