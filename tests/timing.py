import batman
import numpy as np
import matplotlib.pyplot as plt
import timeit
import math

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped

t0 = 0.
per = 1.
rp = 0.115
a = 15.23
inc = 1.555
ecc = 0.
w = math.pi/2.
u = np.array([0.0, 0.7, 0.0, -0.3])
t = np.linspace(t0-per/30., t0 + per/30., 1000)
err_max = 1.0
limb_dark = "nonlinear"

nit = 10

m = batman.TransitModel(t, t0, per, rp, a, inc, ecc, w, u, err_max, limb_dark)

wrapped = wrapper(m.LightCurve, t, t0, per, rp, a, inc, ecc, w, u, err_max, limb_dark)
t = timeit.timeit(wrapped,number=nit)
print "Time (s) to calculate lightcurve for nonlinear LD: ", t/nit




