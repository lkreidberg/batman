from __future__ import print_function
import numpy as np
import math
import matplotlib.pyplot as plt
from .transitmodel import *

def test():
	print("Starting tests...")
	failures = 0
	

	params = TransitParams()
	params.t0 = 0.
	params.per = 1.
	params.rp = 0.115
	params.a = 15.23
	params.inc = 87.
	params.ecc = 0.
	params.w = 90.
	"""params.u = np.array([0.0, 0.7, 0.0, -0.3])
	params.limb_dark = "nonlinear"

	t = np.linspace(params.t0+params.per/60., params.t0 + params.per/30., 1000)
	err_max = 0.5
	
	m = TransitModel(params, t, err_max)
	nonlinear_lc = m.LightCurve(params)
	err = m.calc_err()
	if err > err_max: failures += 1

	params.limb_dark = "quadratic"
	params.u = [0.1,0.3]
	#m = TransitModel(params, t, err_max)
	for i in range(10000000): quadratic_lc = m.LightCurve(params)

	if np.max(np.abs(quadratic_lc-nonlinear_lc))*1.0e6 > err_max: failures += 1"""

	t = np.array([0.002])
	params.u = np.array([0.1, 0.3])
	params.limb_dark = "quadratic"
	m = TransitModel(params, t)
	lc = m.LightCurve(params)

	print("Tests finished with " + "{0}".format(failures) + " failures")

	

