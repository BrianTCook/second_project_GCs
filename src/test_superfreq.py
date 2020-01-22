import numpy as np 
import matplotlib.pyplot as plt
from superfreq import SuperFreq

gravity_solver_str = 'Brute'
data_dir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
filename_time = 'time_data_Nclusters=1_%s.npy'%(gravity_solver_str)
filename_phase = 'sixD_data_Nclusters=1_%s.npy'%(gravity_solver_str)

phase_space_data = np.load(data_dir+filename)
ntimes, ndim, nparticles = phase_space_data.shape
print(ntimes, ndim, nparticles)

#fundamental frequency stuff
freqs_all = []
for k in range(nparticles):
    t, w = sim_times, phase_space_data[:, :, k]
    sf = SuperFreq(t)
    fs = [ (w[:,i] * 1j * w[:,i+ndim//2]) for i in range(ndim//2) ]
    freqs, tbl, ix = sf.find_fundamental_frequencies(fs)
    freqs_all.append(freqs)

freqs_all = np.sort(freqs_all)

X_cdf = np.linspace(min(freqs_all), max(freqs_all), 100)
dX_cdf = X_cdf[1] - X_cdf[0]

Y_cdf = freqs_all / ( (freqs_all*dX_cdf).sum() )
cumulated_Y = np.cumsum(Y_cdf * dX_cdf)

plt.figure()
plt.plot(X_cdf, cumulated_Y, 'k', label='CDF')
plt.xlabel('Frequency (units?)', fontsize=12)
plt.savefig('fundamental_freqs_%s.png'%(gravity_solver_str))
