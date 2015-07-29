import batman
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib import rc
import matplotlib.gridspec as gridspec

print(batman.__file__)

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
matplotlib.rcParams.update({'font.size':14})

params = batman.TransitParams()	        #object to store transit parameters
params.t0 = 0. 				#time of periastron passage (for eccentric orbits), OR
					#mid-transit time (for circular orbits)
params.per = 1.58			#orbital period	
params.rp = 0.1				#planet radius (in units of stellar radii)
params.a = 15.				#semi-major axis (in units of stellar radii)
params.inc = 87.			#orbital inclination (in degrees)	
params.ecc = 0.				#eccentricity	
params.w = 90.				#longitude of periastron (in degrees)
params.limb_dark = "nonlinear"          #limb darkening model
params.u = [0.5, 0.1, 0.1, -0.1]       	#limb darkening coefficients
   
t = np.linspace(-0.025, 0.025, 1000)  	#times at which to calculate light curve	
m = batman.TransitModel(params, t)      #initializes model

flux = m.light_curve(params)		#calculates light curve

m.fac = 1.0e-4
flux_precise = m.light_curve(params)

gs = gridspec.GridSpec(2, 1, height_ratios =[3,1], hspace=0.05)

ax = plt.subplot(gs[0, 0])
ax.plot(t*24., flux, color='k')
ax.set_xticks([])
plt.ylabel("Relative flux")
plt.ylim((0.989, 1.001))
plt.xlim((-0.5, 0.5))

ax2 = plt.subplot(gs[1,0])
ax2.plot(t*24., (flux - flux_precise)*1.0e6, color='k')
plt.xlabel("Time from central transit (hours)")
plt.ylabel("Error (ppm)")
plt.xlim((-0.5, 0.5))
plt.ylim((-1.1,0.1))

plt.savefig("f2.eps", dpi=300)
plt.show()

