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
from .transitmodel import *

def test():
	print("Starting tests...")
	failures = 0

	params = TransitParams()
	params.t0 = 0.
	params.per = 1.0
	params.rp = 0.1
	params.a = 15.23
	params.inc = 1.555*180./math.pi
	params.ecc = 0.
	params.w = 90. 
	params.u = np.array([0.0, 0.7, 0.0, -0.3])
	params.limb_dark = "nonlinear"

	t = np.linspace(0.01, 0.05, 1000)
	err_max = 0.7
	
	m = TransitModel(params, t, err_max)
	nonlinear_lc = m.LightCurve(params)
	err = m.calc_err()
	if err > err_max: failures += 1

	params.limb_dark = "quadratic"
	params.u = [0.1,0.3]
	m = TransitModel(params, t, err_max)
	quadratic_lc = m.LightCurve(params)

	if np.max(np.abs(quadratic_lc-nonlinear_lc))*1.0e6 > err_max: failures += 1
#	print(np.max(np.abs(quadratic_lc-nonlinear_lc))*1.0e6)
#	plt.plot((quadratic_lc - nonlinear_lc)*1.0e6)
#	plt.show()

	print("Tests finished with " + "{0}".format(failures) + " failures")

	

