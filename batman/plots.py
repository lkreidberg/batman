from __future__ import print_function
import numpy as np
import math
import matplotlib.pyplot as plt
from .transitmodel import *
import timeit

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped


def make_plots():
	"""zs = np.linspace(0., 1., 1000)
	rp = 0.1
	u = [0., 0.7, 0.0, -0.3]
	f = occultnl.occultnl(zs, rp, u[0], u[1], u[2], u[3], 1.0e-2, 4)
	fhi = occultnl.occultnl(zs, rp, u[0], u[1], u[2], u[3], 1.0e-4, 4)
	fquad = occultquad.occultquad(zs, rp, 0.1, 0.3, 4)
	#for i in range(len(f)): print "z, fnl, fquad", zs[i], f[i], fquad[i]

	for i in range(1,16):
		wrapped = wrapper(occultquad.occultquad, zs, rp, 0.1, 0.3, i)
		t = timeit.timeit(wrapped,number=1)
		print i, t

	plt.plot(zs, (f - fhi)*1.0e6)
	plt.plot(zs, (fhi - fquad)*1.0e6, color='r')
	plt.axvline(0.9)
	plt.show()"""


	#generates Figure FIXME: max err as a function of function call time
	"""zs = np.linspace(0., 1., 1000)
	rp = 0.1
	u = [0., 0.7, 0.0, -0.3]
	n = 20
	ts = []
	errs = []
	f_ref = occultnl.occultnl(zs, rp, u[0], u[1], u[2], u[3], 1.0e-4, 4)
	fac = np.logspace(-3, -1, n)
	for i in range(n):
		wrapped = wrapper(occultnl.occultnl, zs, rp, u[0], u[1], u[2], u[3], fac[i], 12)
		t = timeit.timeit(wrapped,number=10)/10.
		ts.append(t)
		print t
		f= occultnl.occultnl(zs, rp, u[0], u[1], u[2], u[3], fac[i], 12)
		err = np.max(np.abs(f - f_ref))
		errs.append(err)
	plt.plot(np.array(ts), np.array(errs)*1.0e6)
	plt.yscale('log')
	plt.xscale('log')
	plt.xlabel("Time (s)")
	plt.ylabel("Max Err (ppm)")
	plt.show()"""
	

	

