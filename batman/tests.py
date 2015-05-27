import numpy as np
import math
from transitmodel import TransitParams
from transitmodel import TransitModel

def test():
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
	
	err = m.calc_err()
	if err > err_max: failures += 1


	limb_dark = "quadratic"
	params.u = [0.1,0.3]
	m = TransitModel(params, t, err_max, limb_dark)
	quadratic_lc = m.LightCurve(params)

	if np.max(np.abs(quadratic_lc-nonlinear_lc))> err_max: failures += 1

	print "Tests finished with " + "{0}".format(failures) + " failures"
