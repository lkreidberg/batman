import numpy as np
import math
import batman
import matplotlib.pyplot as plt

t0 = 0.
per = 1.
rp = 0.115
a = 15.23
inc = 1.555
ecc = 0.
w = math.pi/2.
#u = np.array([0.1, 0.3])
u = np.array([0.0, 0.7, 0.0, -0.3])
t = np.linspace(t0-per/30., t0 + per/30., 1000)
err_max = 0.1
limb_dark = "nonlinear"

m = batman.TransitModel(t, t0, per, rp, a, inc, ecc, w, u, err_max, limb_dark)
m.calc_err(plot = True)

plt.subplot(221)
limb_dark = "uniform"
u = []
m = batman.TransitModel(t, t0, per, rp, a, inc, ecc, w, u, err_max, limb_dark)
plt.plot(t, m.LightCurve(t, t0, per, rp, a, inc, ecc, w, u, err_max, limb_dark))
plt.title("Uniform LD")

plt.subplot(222)
limb_dark = "linear"
u = [0.1]
m = batman.TransitModel(t, t0, per, rp, a, inc, ecc, w, u, err_max, limb_dark)
plt.plot(t, m.LightCurve(t, t0, per, rp, a, inc, ecc, w, u, err_max, limb_dark))
plt.title("Linear LD")

plt.subplot(223)
limb_dark = "quadratic"
u = [0.1, 0.3]
m = batman.TransitModel(t, t0, per, rp, a, inc, ecc, w, u, err_max, limb_dark)
plt.plot(t, m.LightCurve(t, t0, per, rp, a, inc, ecc, w, u, err_max, limb_dark))
quad = m.LightCurve(t, t0, per, rp, a, inc, ecc, w, u, err_max, limb_dark)
plt.title("Quadratic LD")

plt.subplot(224)
limb_dark = "nonlinear"
u = [0., 0.7, 0.0, -0.3]
m = batman.TransitModel(t, t0, per, rp, a, inc, ecc, w, u, err_max, limb_dark)
nl= m.LightCurve(t, t0, per, rp, a, inc, ecc, w, u, err_max, limb_dark)
plt.plot(t, m.LightCurve(t, t0, per, rp, a, inc, ecc, w, u, err_max, limb_dark))
plt.title("Nonlinear LD")


plt.tight_layout()
plt.show()
