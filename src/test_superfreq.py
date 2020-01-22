import numpy as np #known to work
from superfreq import SuperFreq #not so much
print('hello world!')

#fundamental frequency stuff
freqs_all = []
for k in range(len(gravity.particles)):
	t, w = sim_times, phase_space_data[:, :, k]
	ntimes, ndim = w.shape
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
