import batman
import numpy as np
import matplotlib.pyplot as plt
import timeit
import math

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped

params = batman.TransitParams()
params.t0 = 0.
params.per = 1.
params.rp = 0.115
params.a = 15.23
params.inc = 1.555
params.ecc = 0.
params.w = math.pi/2.
params.u = np.array([0.0, 0.7, 0.0, -0.3])

err_max = 0.1
limb_dark = "nonlinear"

t = np.linspace(params.t0-params.per/30., params.t0 + params.per/30., 1000)

m = batman.TransitModel(params, t, err_max, limb_dark)

nit = 20
wrapped = wrapper(m.LightCurve, params)
t = timeit.timeit(wrapped,number=nit)
print "Time (s) to calculate lightcurve for nonlinear LD: ", t/nit




