import numpy as np
import math
import matplotlib.pyplot as plt
from transitmodel import TransitParams
from transitmodel import TransitModel
from . import occultnl 
from . import occultquad
#import timeit

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped


def test():
	"""zs = np.linspace(0., 1., 1000)
	rp = 0.1
	u = [0., 0.7, 0.0, -0.3]
	f = occultnl.occultnl(zs, rp, u[0], u[1], u[2], u[3], 1.0e-2)
	fhi = occultnl.occultnl(zs, rp, u[0], u[1], u[2], u[3], 1.0e-4)
	fquad = occultquad.occultquad(zs, rp, 0.1, 0.3)
	for i in range(len(f)): print "z, fnl, fquad", zs[i], f[i], fquad[i]

	plt.plot(zs, (f - fhi)*1.0e6)
	plt.plot(zs, (fhi - fquad)*1.0e6, color='r')
	plt.axvline(0.9)
	plt.show()"""

	print "Starting tests..."
	params = TransitParams()
	params.t0 = 0.
	params.per = 1.
	params.rp = 0.115
	params.a = 15.23
	params.inc = 1.555
	params.ecc = 0.
	params.w = math.pi/2.
	params.u = np.array([0.0, 0.7, 0.0, -0.3])

	failures = 0

	t = np.linspace(params.t0-params.per/30., params.t0 + params.per/30., 1000)
	err_max = 0.1
	limb_dark = "nonlinear"

	m = TransitModel(params, t, err_max, limb_dark)
	nonlinear_lc = m.LightCurve(params)
	
	#generates Figure FIXME: max err as a function of function call time
	"""zs = np.linspace(0., 1., 1000)
	rp = 0.1
	u = [0., 0.7, 0.0, -0.3]
	n = 20
	ts = []
	errs = []
	f_ref = occultnl.occultnl(zs, rp, u[0], u[1], u[2], u[3], 1.0e-4)
	fac = np.logspace(-3, -1, n)
	for i in range(n):
		wrapped = wrapper(occultnl.occultnl, zs, rp, u[0], u[1], u[2], u[3], fac[i])
		t = timeit.timeit(wrapped,number=10)/10.
		ts.append(t)
		print t
		f= occultnl.occultnl(zs, rp, u[0], u[1], u[2], u[3], fac[i])
		err = np.max(np.abs(f - f_ref))
		errs.append(err)
	plt.plot(np.array(ts), np.array(errs)*1.0e6)
	plt.yscale('log')
	plt.xscale('log')
	plt.xlabel("Time (s)")
	plt.ylabel("Max Err (ppm)")
	plt.show()"""
	
	err = m.calc_err()
	if err > err_max: failures += 1


	limb_dark = "quadratic"
	params.u = [0.1,0.3]
	m = TransitModel(params, t, err_max, limb_dark)
	quadratic_lc = m.LightCurve(params)

	if np.max(np.abs(quadratic_lc-nonlinear_lc))> err_max: failures += 1

	print "Tests finished with " + "{0}".format(failures) + " failures"


