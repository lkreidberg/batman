import numpy as np
import math
import batman
import matplotlib.pyplot as plt

params = batman.TransitParams()
params.t0 = 0.
params.per = 1.
params.rp = 0.115
params.a = 15.23
params.inc = 1.555
params.ecc = 0.
params.w = math.pi/2.
params.u = np.array([0.0, 0.7, 0.0, -0.3])

t = np.linspace(params.t0-params.per/30., params.t0 + params.per/30., 1000)
err_max = 0.1
limb_dark = "nonlinear"

m = batman.TransitModel(params, t, err_max, limb_dark)
m.calc_err(plot = True)

ld = ["uniform", "linear", "quadratic", "nonlinear"]
ld_coeffs = [[],[0.1],[0.1,0.3],[0.,0.7,0.,-0.3]]

for i in range(4):
	plt.subplot(220+i)
	limb_dark  = ld[i]
	params.u = ld_coeffs[i]
	m = batman.TransitModel(params, t, err_max, limb_dark)
	plt.plot(t, m.LightCurve(params))
	plt.plot(t, m.LightCurve(params))
	plt.title(limb_dark+" LD")
	plt.ylim((0.985,1.001))

plt.tight_layout()
plt.show()
