import batman
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
matplotlib.rcParams.update({'font.size':14})

params = batman.TransitParams()
params.t0 = 0.                               #time of periastron passage (for eccentric orbits), OR
                                             #mid-transit time (for circular orbits)
params.per = 1.                              #orbital period
params.rp = 0.1                              #planet radius (in stellar radii)
params.a = 12.                               #semi-major axis (in stellar radii)
params.inc = 87.                             #orbital inclination
params.ecc = 0.                              #eccentricity
params.w = 90.                               #longitude of periastron  #FIXME check if this makes sense
params.u = [0., 0.7, 0., -0.3]               #limb darkening coefficients
params.limb_dark = "nonlinear"       	     #limb darkening model

t = np.linspace(-0.02, 0.02, 1000)           #times at which to calculate light curve

m = batman.TransitModel(params, t)

flux = m.LightCurve(params)

params.rp = 0.11                             #updates the planet radius
new_flux = m.LightCurve(params)              #recalculates model light curve

plt.plot(t*24., flux, label = "rp = 0.1")
plt.plot(t*24., new_flux, label = "rp = 0.11")
plt.xlabel("Time from central transit (hours)")
plt.ylabel("Flux")
plt.legend()
plt.savefig("new_rp.png")
plt.show()


ld_options = ["uniform", "linear", "quadratic", "nonlinear"]
ld_coefficients = [[], [0.3], [0.1, 0.3], [0.5, 0.1, 0.1, -0.1]]

plt.figure()

for i in range(4):
	params.limb_dark = ld_options[i]             #specifies the limb darkening profile
       	params.u = ld_coefficients[i]		     #updates limb darkening coefficients
	m = batman.TransitModel(params, t)  	     #initializes model with new limb darkening
	flux = m.LightCurve(params)		     #calculates light curve
	plt.plot(t*24., flux, label = ld_options[i])

plt.xlim((-0.3, 0.3))
plt.legend()
plt.xlabel("Time from central transit (hours)")
plt.ylabel("Relative flux")
plt.savefig("lightcurves.png")
plt.show()

