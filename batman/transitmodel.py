import numpy as np
import matplotlib.pyplot as plt
from . import _nonlinear_ld
from . import _quadratic_ld
from . import _uniform_ld
from . import _custom_ld
from . import _rsky
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

	:param params:
		A :attr:`TransitParams` object containing the physical parameters of the transit

	:param t:
		Array of times at which to calculate the model.

	:param max_err: (optional)
		Error tolerance (in parts per million) for the model.

	:param nthreads: (optional)
		Number of threads to use for parallelization. 

	"""
	def __init__(self, params, t, max_err=1.0, nthreads = None):
		#checking for invalid input
		if (params.limb_dark == "uniform" and len(params.u) != 0) or (params.limb_dark == "linear" and len(params.u) != 1) or \
		    (params.limb_dark == "quadratic" and len(params.u) != 2) or (params.limb_dark == "nonlinear" and len(params.u) != 4):
			raise Exception("Incorrect number of coefficients for " +params.limb_dark + " limb darkening; u should have the form:\n \
			 u = [] for uniform LD\n \
			 u = [u1] for linear LD\n \
  			 u = [u1, u2] for quadratic LD\n \
			 u = [u1, u2, u3, u4] for nonlinear LD, or\n \
		         u = [u1, ..., un] for custom LD") 
		if params.limb_dark not in ["uniform", "linear", "quadratic", "nonlinear", "custom"]: 
			raise Exception("\""+params.limb_dark+"\""+" limb darkening not supported; allowed options are:\n \
				uniform, linear, quadratic, nonlinear, custom")

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
		self.zs= _rsky._rsky(t.tolist(), params.t0, params.per, params.a, params.inc*pi/180., params.ecc, params.w*pi/180.)
		self.fac = self._get_fac()
		if nthreads==None or nthreads == 1: self.nthreads=1
		else:
			if nthreads <= multiprocessing.cpu_count()and nthreads >1 and openmp.detect(): self.nthreads = nthreads
			else: 
				if nthreads>multiprocessing.cpu_count(): raise Exception("Maximum number of threads is "+'{0:d}'.format(multiprocessing.cpu_count()))
				elif nthreads <= 1: raise Exception("Number of threads must be between 2 and {0:d}".format(multiprocessing.cpu_count()))
				else: raise Exception("OpenMP not enabled: do not set the nthreads parameter")

	def calc_err(self, plot = False):
		"""

		Returns maximum error (in parts per million) for transit light curve calculation.
			
		:param plot: 
			If set to ``True``, plots the error in the light curve model as a function of separation of centers.

		"""
		if self.limb_dark in ["nonlinear", "custom"]:
			zs = np.linspace(0., 1.1, 500).tolist()
			fac_lo = 1.0e-4
			if self.limb_dark == "nonlinear":
				f0 = np.array(_nonlinear_ld._nonlinear_ld(zs, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], fac_lo, self.nthreads))
				f = np.array(_nonlinear_ld._nonlinear_ld(zs, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.fac, self.nthreads))
			else:
				f0 = np.array(_custom_ld._custom_ld(zs, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.u[4], self.u[5], fac_lo, self.nthreads))
				f =  np.array(_custom_ld._custom_ld(zs, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.u[4], self.u[5], self.fac, self.nthreads))
	
			err = np.max(np.abs(f-f0))*1.0e6
			if plot == True:
				plt.plot(zs, 1.0e6*(f-f0), color='k')
				plt.xlabel("d (separation of centers)")
				plt.ylabel("Error (ppm)") 
				plt.show()

			return err
		else: raise Exception("Function calc_err not valid for " + self.limb_dark + " limb darkening")

	def _get_fac(self):
		if self.limb_dark in ["nonlinear", "custom"]:
			nthreads = 1
			fac_lo, fac_hi = 1.0e-4, 1.
			zs = np.linspace(0., 1.1, 500).tolist()
			if self.limb_dark == "nonlinear": f0 = np.array(_nonlinear_ld._nonlinear_ld(zs, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], fac_lo, nthreads))
			else: f0 = np.array(_custom_ld._custom_ld(zs, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.u[4], self.u[5], fac_lo, nthreads))
			n = 0
			err = 0.
			while(err > self.max_err or err < 0.99*self.max_err):
				fac = (fac_lo + fac_hi)/2.
				if self.limb_dark == "nonlinear": f = np.array(_nonlinear_ld._nonlinear_ld(zs, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], fac, nthreads))
				else: f = np.array(_custom_ld._custom_ld(zs, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.u[4], self.u[5], fac, nthreads))

				err = np.max(np.abs(np.array(f)-np.array(f0)))*1.0e6
				if err> self.max_err: fac_hi = fac	
				else: fac_lo = fac
				n += 1
				if n>1e3: raise Exception("Convergence failure in calculation of scale factor for _nonlinear_ld")
			return fac
		else: return 0.

	def LightCurve(self, params):
		"""
		Calculates a model light curve.

		:param params:
			Transit parameter object.

		"""
		#recalculates rsky and fac if necessary
		if params.t0 != self.t0 or params.per != self.per or params.a != self.a or params.inc != self.inc or params.ecc != self.ecc or params.w != self.w: self.zs= _rsky._rsky(t, params.t0, params.per, params.a, params.inc*pi/180., params.ecc, params.w*pi/180.)
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
		if self.limb_dark == "quadratic": return np.array(_quadratic_ld._quadratic_ld(self.zs, params.rp, params.u[0], params.u[1], self.nthreads))
		elif self.limb_dark == "nonlinear": return np.array(_nonlinear_ld._nonlinear_ld(self.zs, params.rp, params.u[0], params.u[1], params.u[2], params.u[3], self.fac, self.nthreads))
		elif self.limb_dark == "linear": return np.array(_quadratic_ld._quadratic_ld(self.zs, params.rp, params.u[0], 0., self.nthreads))
		elif self.limb_dark == "uniform": return np.array(_uniform_ld._uniform_ld(self.zs, params.rp, self.nthreads))
		elif self.limb_dark == "custom": return np.array(_custom_ld._custom_ld(self.zs, params.rp, params.u[0], params.u[1], params.u[2], params.u[3], params.u[4], params.u[5], self.fac, self.nthreads))
		else: raise Exception("Invalid limb darkening option")
			

class TransitParams(object):
	"""
	An object that stores the physical parameters of the transit.

	:param t0:
		Time of periastron passage (for eccentric orbits) or time of central transit (for circular orbits).

	:param per:
		Orbital period [in days].

	:param rp:
		Planet radius [in stellar radii].

	:param a:
		Semi-major axis [in stellar radii].

	:param inc:
		Orbital inclination [in degrees].

	:param ecc:
		Orbital eccentricity.

	:param w:
	 	Argument of periapse [in degrees] (FIXME ask Dan).	

	:param u:
		List of limb darkening coefficients.

	:param limb_dark:
		Limb darkening model (choice of "nonlinear", "quadratic", "linear", "uniform", or "custom")


	"""
	def __init__(self):
		self.t0 = 0.
		self.per = 0.
		self.rp = 0.
		self.a = 0.
		self.inc = 0.
		self.ecc = 0.
		self.w = 0. 
		self.u = []
		self.limb_dark = ""
