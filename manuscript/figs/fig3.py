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

from __future__ import print_function
import numpy as np
import math
import matplotlib.pyplot as plt
import timeit
import batman
from pylab import *

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
#matplotlib.rcParams.update({'font.size':14})
plt.figure(figsize=(4,4))

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

flux = m.LightCurve(params)		#calculates light curve

#generates Figure FIXME: max err as a function of function call time
zs = np.linspace(0., 1., 1000)
rp = 0.1
u = [0., 0.7, 0.0, -0.3]
n = 20
ts = []
errs = []
fac = np.logspace(-3, -1, n)
for i in range(n):
	m.set_fac(fac[i])
	wrapped = wrapper(m.LightCurve, params)
	t = timeit.timeit(wrapped,number=10)/10.
	ts.append(t)
	print(t)
	err = m.calc_err() 
	errs.append(err)
plt.plot(np.array(ts), np.array(errs), color='k')
print(np.min(errs), np.max(errs))
print(np.max(ts), np.min(ts))
plt.xlim((1.0e-3, 5.0e-2))
plt.ylim((1.0e-3, 10.))
plt.yscale('log')
plt.xscale('log')
plt.xlabel("Execution time (s)")
plt.ylabel("Truncation error (ppm)")
plt.tight_layout()
plt.savefig("f3.pdf", dpi=300)
plt.show()




