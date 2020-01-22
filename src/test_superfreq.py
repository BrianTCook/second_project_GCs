import numpy as np 
import matplotlib.pyplot as plt
from superfreq import SuperFreq

gravity_solver_str = 'Brute'
Nclusters = 1

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

print(freqs_all)
freqs_all = [j for i in freqs_all for j in i]
freqs_all = np.sort(freqs_all)
print(freqs_all)

X_cdf = np.linspace(min(freqs_all), max(freqs_all), len(freqs_all))
dX_cdf = X_cdf[1] - X_cdf[0]

Y_cdf = freqs_all / ( (freqs_all*dX_cdf).sum() )
cumulated_Y = np.cumsum(Y_cdf * dX_cdf)

plt.figure()
plt.hist(freqs_all, bins=100)
#plt.plot(X_cdf, cumulated_Y, 'k', label='CDF')
plt.xlabel('Frequency (units?)', fontsize=12)
plt.savefig('fundamental_freqs_Nclusters=%i_%s.png'%(Nclusters, gravity_solver_str))
