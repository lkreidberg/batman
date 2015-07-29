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

import numpy as np
import matplotlib.pyplot as plt
from . import _nonlinear_ld
from . import _quadratic_ld
from . import _uniform_ld
from . import _logarithmic_ld
from . import _exponential_ld
from . import _custom_ld
from . import _rsky
#from . import _mandelagol_nonlinear_ld
from math import pi
import multiprocessing
from . import openmp

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped

class TransitModel:
	"""
	Class for generating model transit light curves.	

	:param params: A :attr:`TransitParams` object containing the physical parameters of the transit
	:type params: a `TransitParams` instance

	:param t: Array of times at which to calculate the model.
	:type t: ndarray 

	:param max_err: Error tolerance (in parts per million) for the model.
	:type max_err: float, optional

	:param nthreads: Number of threads to use for parallelization. 
	:type nthreads: int, optional

	:param fac: Scale factor for integration step size
	:type fac: float, optional


	:Example:
	
	>>> m = batman.TransitModel(params, max_err = 0.5, nthreads=4)
	"""

	def __init__(self, params, t, max_err=1.0, nthreads = 1, fac = None):
		#checking for invalid input
		if  (params.limb_dark == "uniform" and len(params.u) != 0) or (params.limb_dark == "linear" and len(params.u) != 1) or \
		    (params.limb_dark == "quadratic" and len(params.u) != 2) or (params.limb_dark == "logarithmic" and len(params.u) != 2) or \
		    (params.limb_dark == "exponential" and len(params.u) != 2) or (params.limb_dark == "squareroot" and len(params.u) != 2) or \
		    (params.limb_dark == "nonlinear" and len(params.u) != 4):
			raise Exception("Incorrect number of coefficients for " +params.limb_dark + " limb darkening; u should have the form:\n \
			 u = [] for uniform LD\n \
			 u = [u1] for linear LD\n \
  			 u = [u1, u2] for quadratic, logarithmic, exponential, and squareroot LD\n \
			 u = [u1, u2, u3, u4] for nonlinear LD, or\n \
		         u = [u1, ..., un] for custom LD") 
		if params.limb_dark not in ["uniform", "linear", "quadratic", "logarithmic", "exponential", "squareroot", "nonlinear", "custom"]: 
			raise Exception("\""+params.limb_dark+"\""+" limb darkening not supported; allowed options are:\n \
				uniform, linear, quadratic, logarithmic, exponential, squareroot, nonlinear, custom")
		if max_err < 0.001: raise Exception("The lowest allowed value for max_err is 0.001. For more accurate calculation, set the integration step size explicitly with the fac parameter.")

		#initializes model parameters
		self.t = t
		self.t0 = params.t0
		self.per = params.per
		self.rp = params.rp
		self.a = params.a
		self.inc = params.inc
		self.ecc = params.ecc
		self.w = params.w
		self.u = params.u
		self.max_err = max_err
		self.limb_dark = params.limb_dark
		self.ds= _rsky._rsky(t, params.t0, params.per, params.a, params.inc*pi/180., params.ecc, params.w*pi/180.)
		if fac != None: self.fac = fac
		else: self.fac = self._get_fac()
		if nthreads==None or nthreads == 1: self.nthreads=1
		else:
			if nthreads <= multiprocessing.cpu_count()and nthreads >1 and openmp.detect(): self.nthreads = nthreads
			else: 
				if nthreads > multiprocessing.cpu_count(): raise Exception("Maximum number of threads is "+'{0:d}'.format(multiprocessing.cpu_count()))
				elif nthreads <= 1: raise Exception("Number of threads must be between 2 and {0:d}".format(multiprocessing.cpu_count()))
				else: raise Exception("OpenMP not enabled: do not set the nthreads parameter")

	def calc_err(self, plot = False):
		"""

		Calculate maximum error for transit light curve calculation.
			
		:param plot: If ``True``, plots the error in the light curve model as a function of separation of centers.
		:type plot: bool

		:return: Truncation error (parts per million)
		:rtype: float

		"""
		if self.limb_dark in ["logarithmic", "exponential", "nonlinear", "squareroot", "custom"]:
			ds = np.linspace(0., 1.1, 500)
			fac_lo = 5.0e-4
			if self.limb_dark == "nonlinear":
				f0 = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], fac_lo, self.nthreads)
				f = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.fac, self.nthreads)
			elif self.limb_dark == "squareroot":
				f0 = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[1], self.u[0], 0., 0., fac_lo, self.nthreads)
				f = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[1], self.u[0], 0., 0., self.fac, self.nthreads)
			elif self.limb_dark == "exponential":
				f0 = _exponential_ld._exponential_ld(ds, self.rp, self.u[0], self.u[1], fac_lo, self.nthreads)
				f = _exponential_ld._exponential_ld(ds, self.rp, self.u[0], self.u[1], self.fac, self.nthreads)
			elif self.limb_dark == "logarithmic":
				f0 = _logarithmic_ld._logarithmic_ld(ds, self.rp, self.u[0], self.u[1], fac_lo, self.nthreads)
				f = _logarithmic_ld._logarithmic_ld(ds, self.rp, self.u[0], self.u[1], self.fac, self.nthreads)
			else:
				f0 = _custom_ld._custom_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.u[4], self.u[5], fac_lo, self.nthreads)
				f =  _custom_ld._custom_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.u[4], self.u[5], self.fac, self.nthreads)
	
			err = np.max(np.abs(f-f0))*1.0e6
			if plot == True:
				plt.plot(ds, 1.0e6*(f-f0), color='k')
				plt.xlabel("d (separation of centers)")
				plt.ylabel("Error (ppm)") 
				plt.show()

			return err
		else: raise Exception("Function calc_err not valid for " + self.limb_dark + " limb darkening")

	def _get_fac(self):
		if self.limb_dark in ["logarithmic", "exponential", "squareroot", "nonlinear", "custom"]:
			nthreads = 1
			fac_lo, fac_hi = 5.0e-4, 1.
			ds = np.linspace(0., 1.+self.rp, 1000)
			if self.limb_dark == "nonlinear": f0 = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], fac_lo, nthreads)
			elif self.limb_dark == "squareroot": f0 = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[1], self.u[0], 0., 0., fac_lo, nthreads)
			elif self.limb_dark == "exponential": f0 = _exponential_ld._exponential_ld(ds, self.rp, self.u[0], self.u[1], fac_lo, nthreads)
			elif self.limb_dark == "logarithmic": f0 = _logarithmic_ld._logarithmic_ld(ds, self.rp, self.u[0], self.u[1], fac_lo, nthreads)
			else: f0 = _custom_ld._custom_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.u[4], self.u[5], fac_lo, nthreads)

			n = 0
			err = 0.
			while(err > self.max_err or err < 0.99*self.max_err):
				fac = (fac_lo + fac_hi)/2.
				if self.limb_dark == "nonlinear": f = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], fac, nthreads)
				elif self.limb_dark == "squareroot": f = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[1], self.u[0], 0., 0., fac, nthreads)
				elif self.limb_dark == "exponential": f = _exponential_ld._exponential_ld(ds, self.rp, self.u[0], self.u[1], fac, nthreads)
				elif self.limb_dark == "logarithmic": f = _logarithmic_ld._logarithmic_ld(ds, self.rp, self.u[0], self.u[1], fac, nthreads)
				else: f = _custom_ld._custom_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.u[4], self.u[5], fac, nthreads)

				err = np.max(np.abs(f-f0))*1.0e6
				if err> self.max_err: fac_hi = fac	
				else: fac_lo = fac
				n += 1
				if n > 1e3: raise Exception("Convergence failure in calculation of scale factor for integration step size")
			return fac
		else: return None
	
	def light_curve(self, params):
		"""
		Calculate a model light curve.

		:param params: Transit parameters
		:type params: A `TransitParams` instance

		:return: Relative flux 
		:rtype: ndarray

		:Example:

		>>> flux = m.light_curve(params)
		"""
		#recalculates rsky and fac if necessary
		if params.t0 != self.t0 or params.per != self.per or params.a != self.a or params.inc != self.inc or params.ecc != self.ecc or params.w != self.w: self.ds= _rsky._rsky(self.t, params.t0, params.per, params.a, params.inc*pi/180., params.ecc, params.w*pi/180.)
		if params.limb_dark != self.limb_dark: self.fac = self._get_fac()
		#updates transit params
		self.t0 = params.t0
		self.per = params.per
		self.rp = params.rp
		self.a = params.a
		self.inc = params.inc
		self.ecc = params.ecc
		self.w = params.w
		self.u = params.u
		self.limb_dark = params.limb_dark
		
		if params.limb_dark != self.limb_dark: raise Exception("Need to reinitialize model in order to change limb darkening option")
		if self.limb_dark == "quadratic": return (_quadratic_ld._quadratic_ld(self.ds, params.rp, params.u[0], params.u[1], self.nthreads))
		elif self.limb_dark == "linear": return _quadratic_ld._quadratic_ld(self.ds, params.rp, params.u[0], 0., self.nthreads)
		elif self.limb_dark == "nonlinear": return _nonlinear_ld._nonlinear_ld(self.ds, params.rp, params.u[0], params.u[1], params.u[2], params.u[3], self.fac, self.nthreads)
		elif self.limb_dark == "squareroot": return _nonlinear_ld._nonlinear_ld(self.ds, params.rp, params.u[1], params.u[0], 0., 0., self.fac, self.nthreads)
		elif self.limb_dark == "uniform": return _uniform_ld._uniform_ld(self.ds, params.rp, self.nthreads)
		elif self.limb_dark == "logarithmic": return (_logarithmic_ld._logarithmic_ld(self.ds, params.rp, params.u[0], params.u[1], self.fac, self.nthreads))
		elif self.limb_dark == "exponential": return (_exponential_ld._exponential_ld(self.ds, params.rp, params.u[0], params.u[1], self.fac, self.nthreads))
		elif self.limb_dark == "custom": return _custom_ld._custom_ld(self.ds, params.rp, params.u[0], params.u[1], params.u[2], params.u[3], params.u[4], params.u[5], self.fac, self.nthreads)
		#elif self.limb_dark == "mandelagol": return _mandelagol_nonlinear_ld._mandelagol_nonlinear_ld(self.ds, params.u[0], params.u[1], params.u[2], params.u[3], params.rp, len(self.ds))
		else: raise Exception("Invalid limb darkening option")
			

class TransitParams(object):
	"""
	Object to store the physical parameters of the transit.

	:param t0: Time of periastron passage (for eccentric orbits) or time of central transit (for circular orbits).
	:type t0: float

	:param per: Orbital period [in days].
	:type per: float

	:param rp: Planet radius [in stellar radii].
	:type rp: float

	:param a: Semi-major axis [in stellar radii].
	:type a: float

	:param inc: Orbital inclination [in degrees].
	:type inc: float

	:param ecc: Orbital eccentricity.
	:type ecc: float

	:param w: Argument of periapse [in degrees]
	:type w: float

	:param u: List of limb darkening coefficients.
	:type u: array_like 

	:param limb_dark: Limb darkening model (choice of "nonlinear", "quadratic", "exponential", "logarithmic", "squareroot", "linear", "uniform", or "custom")
	:type limb_dark: str

	:Example:
	
	>>> import batman
	>>> params = batman.TransitParams()
	>>> params.t0 = 0. 				#time of periastron passage (for eccentric orbits), OR
	>>>						#mid-transit time (for circular orbits)
	>>> params.per = 1.				#orbital period	
	>>> params.rp = 0.1				#planet radius (in units of stellar radii)
	>>> params.a = 15.				#semi-major axis (in units of stellar radii)
	>>> params.inc = 87.				#orbital inclination (in degrees)	
	>>> params.ecc = 0.				#eccentricity	
	>>> params.w = 90.				#longitude of periastron (in degrees) 
	>>> params.u = [0.1, 0.3] 	      	        #limb darkening coefficients
	>>> params.limb_dark = "quadratic"          	#limb darkening model
	"""
	def __init__(self):
		self.t0 = None
		self.per = None
		self.rp = None
		self.a = None
		self.inc = None
		self.ecc = None
		self.w = None
		self.u = None
		self.limb_dark = None
