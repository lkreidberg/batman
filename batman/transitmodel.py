import numpy as np
import matplotlib.pyplot as plt
from . import occultnl
from . import occultquad
from . import occultuniform
from . import rsky

class TransitModel:
	"""
	doc
	"""
	def __init__(self, t, t0, per, rp, a, inc, ecc, w, u, max_err, limb_dark):
		if (limb_dark == "uniform" and len(u) != 0) or (limb_dark == "linear" and len(u) != 1) or \
		    (limb_dark == "quadratic" and len(u) != 2) or (limb_dark == "nonlinear" and len(u) != 4):
			raise Exception("Incorrect number of coefficients for " +limb_dark + " limb darkening; u should have the form:\n \
			 u = [] for uniform LD\n \
			 u = [u1] for linear LD\n \
  			 u = [u1, u2] for quadratic LD\n \
			 u = [u1, u2, u3, u4] for nonlinear LD") 
		self.t = t
		self.t0 = t0
		self.per = per
		self.rp = rp
		self.a = a
		self.inc = inc
		self.ecc = ecc
		self.w = w 
		self.u = u
		self.max_err = max_err
		self.limb_dark = limb_dark
		self.zs= rsky.rsky(t, t0, per, a, inc, ecc, w)
		self.fac = self._get_fac()

#	def set_fac(self,fac):				#set scale factor manually
#		self.fac = fac

	def calc_err(self, plot = False):
		if self.limb_dark == "nonlinear":
			zs = np.linspace(0., 1.1, 500)
			fac_lo = 1.0e-4
			f0 = occultnl.occultnl(zs, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], fac_lo)
			f = occultnl.occultnl(zs, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.fac)
			err = np.max(np.abs(f-f0))*1.0e6
			print "Max err in light curve is " + "{0:0.2f}".format(err), "ppm"
			if plot == True:
				plt.plot(zs, 1.0e6*(f-f0), color='k')
				plt.xlabel("z (separation of centers)")
				plt.ylabel("Error (ppm)") 
				plt.show()
		else: raise Exception("Function calc_err not valid for " + self.limb_dark + " limb darkening")
	
	def _get_fac(self):
		if self.limb_dark == "nonlinear":
			fac_lo, fac_hi = 1.0e-4, 1.
			zs = np.linspace(0., 1.1, 500)
			f0 = occultnl.occultnl(zs, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], fac_lo)
			n = 0
			err = 0.
			while(err > self.max_err or err < 0.99*self.max_err):
				fac = (fac_lo + fac_hi)/2.
				f = occultnl.occultnl(zs, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], fac)
				err = np.max(np.abs(f-f0))*1.0e6
				if err> self.max_err: fac_hi = fac	
				else: fac_lo = fac
				n += 1
				if n>1e4: raise Exception("Convergence failure in calculation of scale factor for occultnl")
			return fac
		else: return 0.

	def LightCurve(self, t, t0, per, rp, a, inc, ecc, w, u, max_err, limb_dark):
		if limb_dark == "quadratic": return occultquad.occultquad(self.zs, rp, u[0], u[1])
		elif limb_dark == "nonlinear": return occultnl.occultnl(self.zs, rp, u[0], u[1], u[2], u[3], self.fac)
		elif limb_dark == "linear": return occultquad.occultquad(self.zs, rp, u[0], 0.)
		elif limb_dark == "uniform": return occultuniform.occultuniform(self.zs, rp)
		else: raise Exception("Limb darkening \"" + limb_dark+"\" not yet implemented.") 
			

