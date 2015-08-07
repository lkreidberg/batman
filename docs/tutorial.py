# The batman package: fast computation of exoplanet transit light curves
# Copyright (C) 2015 Laura Kreidberg	 
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import batman
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
matplotlib.rcParams.update({'font.size':14})

params = batman.TransitParams()	        #object to store transit parameters
params.t0 = 0. 				#time of inferior conjunction 
params.per = 1.				#orbital period	
params.rp = 0.1				#planet radius (in units of stellar radii)
params.a = 15.				#semi-major axis (in units of stellar radii)
params.inc = 87.			#orbital inclination (in degrees)	
params.ecc = 0.				#eccentricity	
params.w = 90.				#longitude of periastron (in degrees)
params.limb_dark = "nonlinear"          #limb darkening model
params.u = [0.5, 0.1, 0.1, -0.1]       	#limb darkening coefficients
   
"""t = np.linspace(-0.025, 0.025, 1000)  	#times at which to calculate light curve	
m = batman.TransitModel(params, t)      #initializes model

flux = m.light_curve(params)		#calculates light curve

radii = np.linspace(0.09, 0.11, 20)
for r in radii:
	params.rp = r				#updates planet radius
	new_flux = m.light_curve(params)        #recalculates light curve
	plt.plot(t, new_flux)

plt.xlabel("Time from central transit (days)")
plt.ylabel("Relative flux")
plt.ylim((0.987, 1.001))

plt.savefig("change_rp.png")


ld_options = ["uniform", "linear", "quadratic", "nonlinear"]
ld_coefficients = [[], [0.3], [0.1, 0.3], [0.5, 0.1, 0.1, -0.1]]

plt.figure()

for i in range(4):
	params.limb_dark = ld_options[i]             	#specifies the limb darkening profile
	params.u = ld_coefficients[i]	         	#updates limb darkening coefficients
	m = batman.TransitModel(params, t)	        #initializes the model
	flux = m.light_curve(params)		        #calculates light curve
	plt.plot(t, flux, label = ld_options[i])

plt.xlim((-0.025, 0.025))
plt.ylim((0.987, 1.001))
plt.legend()
plt.xlabel("Time (days)")
plt.ylabel("Relative flux")
plt.savefig("lightcurves.png")
#plt.show()

m = batman.TransitModel(params, t, max_err = 0.5)
plt.clf()
m.calc_err(plot = True) 

#m = batman.TransitModel(params, t, nthreads = 4)"""

params.fp = 0.001
params.t_secondary = 0.5
t = np.linspace(0.48, 0.52, 1000)
m = batman.TransitModel(params, t, transittype="secondary")	       
flux = m.light_curve(params)

plt.plot(t, flux)
plt.ylim((0.9995, 1.0015))
plt.xlim((0.48, 0.52))
plt.gca().set_yticklabels([0.9995, 1.000, 1.0005, 1.001, 1.0015])
plt.xlabel("Time (days)") 
plt.ylabel("Relative flux")
plt.savefig("eclipse.png")
plt.show()


m = batman.TransitModel(params, t, supersample_factor = 7, exp_time = 0.001)
