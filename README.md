# batman

### About
The `batman` ("Bad-Ass Transit Model cAlculatioN") package for Python provides fast, flexible calculation of exoplanet transit light curves.  The package includes a new brute-force integration scheme that can be easily generalized to any radially symmetric limb-darkening profile, and makes use of C extension modules to speed up computation.

Currently supported limb-darkening parameterizations are: 
 - uniform 
 - linear 
 - quadratic
 - nonlinear


### Installation

Clone this repository or download the source as a zip file. 

Run `$ python setup.py install` from the source root directory to install.

To test the installation, `cd` out of the root directory and run
```
$ python -c 'import batman; batman.test()'
```

### Usage
Here is an example calculation of a light curve with nonlinear limb-darkening.

```python
import batman
import numpy as np
import matplotlib.pyplot as plt

#initializes model
params = batman.Params()
params.t0 = 0.		#transit ephemeris
params.per = 1.		#orbital period	
params.rp = 0.1		#planet radius (in stellar radii)
params.a = 15.		#semi-major axis (in stellar radii)
params.inc = 1.55	#orbital inclination	
params.ecc = 0.		#eccentricity	
params.w = 1.57		#longitude of periastron

t = np.linspace(-0.05, 0.05, 1000)         #times to calculate light curve	
err_max = 0.1                              #maximum error in light curve (in ppm)
limb_dark = "nonlinear"                    #limb darkening model
u = np.array([0., 0.7, 0., -0.3])          #limb darkening coefficients

m = batman.TransitModel(params, t, err_max, limb_dark)

#calculates model light curve
flux = m.LightCurve(params)

#plots light curve
plt.plot(t, flux)
plt.show()
```
