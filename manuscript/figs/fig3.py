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
matplotlib.rcParams.update({'font.size':11})
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
   
t = np.linspace(-0.013, 0.013, 100)  	#times at which to calculate light curve	

params.limb_dark = "quadratic"          #limb darkening model
params.u = [0.1, 0.3]       	#limb darkening coefficients
m = batman.TransitModel(params, t)      #initializes model
f1 = m.light_curve(params)
wrapped = wrapper(m.light_curve, params)
time = timeit.timeit(wrapped,number=100)/100.
#plt.axvline(time, color='0.5', linestyle='dashed')
plt.gca().annotate("", xy=(time, .1), xycoords='data', xytext=(time, 1), textcoords='data', arrowprops=dict(arrowstyle="->")) 

print("quadratic", time)

params.limb_dark = "nonlinear"          #limb darkening model
params.u = [0.0, 0.7, 0.0, -0.3]       	#limb darkening coefficients
m = batman.TransitModel(params, t)      #initializes model
flux = m.light_curve(params)		#calculates light curve

#generates Figure 3: max err as a function of function call time
n = 100
ts = []
errs = []
fac = np.logspace(-2, 0.3, n)
for i in range(n):
	m.fac = fac[i]
	wrapped = wrapper(m.light_curve, params)
	if i<10: t = timeit.timeit(wrapped,number=200)/200.
	else: t = timeit.timeit(wrapped,number=1000)/1000.
	ts.append(t)
	err = m.calc_err() 
	print(t, err)
	errs.append(err)
plt.plot(np.array(ts), np.array(errs), color='k')
print(np.min(errs), np.max(errs))
print(np.max(ts), np.min(ts))


plt.gca().tick_params('both', length=8, width=1.2, which='major')
plt.gca().tick_params('both', length=4, width=0.8, which='minor')

plt.xlim((2.0e-5, 1.0e-3))
plt.ylim((1.0e-1, 300.))
plt.yscale('log')
plt.xscale('log')
plt.xlabel("Execution time (s)")
plt.ylabel("Truncation error (ppm)")
plt.tight_layout()
plt.savefig("f3.eps", dpi=300)
plt.show()




