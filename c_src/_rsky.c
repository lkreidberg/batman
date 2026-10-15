/* The batman package: fast computation of exoplanet transit light curves
 * Copyright (C) 2015 Laura Kreidberg
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define NPY_NO_DEPRECATED_API NPY_2_0_API_VERSION
#include <Python.h>
#include "numpy/arrayobject.h"

#if defined (_OPENACC) && defined(__PGI)
#  include <accelmath.h>
#else
#  include <math.h>
#endif

#if defined (_OPENMP) && !defined(_OPENACC)
#  include <omp.h>
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

static PyObject *_rsky(PyObject *self, PyObject *args);

static PyObject *_getf(PyObject *self, PyObject *args);

#define KBLK 8		//samples solved together, so that their independent iterations overlap

/* Mean anomaly in [-pi, pi] from the orbital phase; the integer number of orbits is removed
   before multiplying by 2 pi, which is exact and more accurate than reducing n*(t - tp). */
static inline double Mphase(double t, double tp, double per)
{
	double x = (t - tp)/per;
	return 2.*M_PI*(x - rint(x));
}

/* Solves Kepler's equation E - e sin(E) = M for 0 <= e < 1 for n <= KBLK mean anomalies
   M[j] in [-pi, pi] (NaN allowed) and returns sin(E) and cos(E).

   E is odd in M, so the equation is solved for |M|. The start is E = |M| + 0.85 e (Danby 1988)
   and the iteration is Halley's method (cubic convergence), which needs one sin and one cos per
   step. E stays in [0, pi], where f(E) = E - e sin(E) - |M| is monotonic with f(0) <= 0 <= f(pi).
   Two steps are done for all samples of the block (they are independent, so that the processor
   overlaps them). A sample whose second step is smaller than 1e-6 is finished: the cubic
   convergence leaves an error far below 1e-16, and sin and cos at the new E follow from the last
   sin and cos by a Taylor expansion, so no further trigonometric call is needed. The other
   samples (about two thirds at e = 0.3) continue one by one, and bisection on [0, pi] is the
   fall-back if 30 steps did not converge (it always does).
*/
static inline void kepler_block(const double *M, int n, double e, double *sinE, double *cosE)
{
	double E[KBLK], Mr[KBLK], s[KBLK], c[KBLK], dE[KBLK];
	for (int j = 0; j < n; j++)
	{
		Mr[j] = fabs(M[j]);
		E[j] = MIN(Mr[j] + 0.85*e, M_PI);
		s[j] = 0.; c[j] = 1.; dE[j] = 0.;
	}
	for (int it = 0; it < 2; it++)
	{
		for (int j = 0; j < n; j++)
		{
			s[j] = sin(E[j]); c[j] = cos(E[j]);
			double f = E[j] - e*s[j] - Mr[j];
			double fp = 1. - e*c[j];
			double d = -f*fp/(fp*fp - 0.5*f*e*s[j]);			//Halley step
			double En = E[j] + d;
			En = (En < 0.) ? 0. : ((En > M_PI) ? M_PI : En);
			dE[j] = En - E[j];
			E[j] = En;
		}
	}
	for (int j = 0; j < n; j++)
	{
		double sg = (M[j] < 0.) ? -1. : 1.;
		double Ej = E[j], sj = s[j], cj = c[j], dj = dE[j];
		if (!(Mr[j] <= M_PI))							//NaN
		{
			sinE[j] = NAN; cosE[j] = NAN;
			continue;
		}
		int it = 2;
		while (!(fabs(dj) < 1.e-6) && it < 30)
		{
			sj = sin(Ej); cj = cos(Ej);
			double f = Ej - e*sj - Mr[j];
			double fp = 1. - e*cj;
			double d = -f*fp/(fp*fp - 0.5*f*e*sj);
			double En = Ej + d;
			En = (En < 0.) ? 0. : ((En > M_PI) ? M_PI : En);
			dj = En - Ej;
			Ej = En;
			it++;
		}
		if (fabs(dj) < 1.e-6)
		{
			sinE[j] = sg*(sj + cj*dj - 0.5*sj*dj*dj);
			cosE[j] = cj - sj*dj - 0.5*cj*dj*dj;
		}
		else									//bisection
		{
			double lo = 0., hi = M_PI;
			for (int k = 0; k < 60; k++)
			{
				Ej = 0.5*(lo + hi);
				if (Ej - e*sin(Ej) - Mr[j] > 0.) hi = Ej; else lo = Ej;
			}
			sinE[j] = sg*sin(Ej);
			cosE[j] = cos(Ej);
		}
	}
}

static PyObject *_rsky_or_f(PyObject *self, PyObject *args, int f_only)
{
	/*
		This module computes the distance between the centers of the
		star and the planet in the plane of the sky.  This parameter is
		denoted r_sky = sqrt(x^2 + y^2) in the Seager Exoplanets book
		(see the section by Murray, and Winn eq. 5).  In the Mandel & Agol
		(2002) paper, this quantity is denoted d.

		If f_only is 1, this function returns the true anomaly instead of the distance.
	*/
	double ecc, inc, a, omega, per, tc, BIGD = 100.;
	int transittype, nthreads;;
	npy_intp dims[1];
	PyArrayObject *ts, *ds;

  	if(!PyArg_ParseTuple(args,"Oddddddii", &ts, &tc, &per, &a, &inc, &ecc, &omega, &transittype, &nthreads)) return NULL;

	dims[0] = PyArray_DIMS(ts)[0];
	ds = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(ts));

	double *t_array = PyArray_DATA(ts);
	double *output_array = PyArray_DATA(ds);



	#if defined (_OPENMP) && !defined(_OPENACC)
	omp_set_num_threads(nthreads);	//specifies number of threads (if OpenMP is supported)
	#endif

	/* Quantities that do not depend on the sample are computed once. */
	const double f0 = M_PI/2. - omega;						//true anomaly corresponding to time of primary transit center
	const int circular = (ecc < 1.0e-5);
	double tp;									//time of periastron
	if(circular) tp = tc - per*f0/2./M_PI;
	else
	{
		double E0 = 2.*atan(sqrt((1. - ecc)/(1. + ecc))*tan(f0/2.));		//corresponding eccentric anomaly
		double M0 = E0 - ecc*sin(E0);
		tp = tc - per*M0/2./M_PI;
	}
	const double sini = sin(inc);
	const double sini2 = sini*sini;
	const double sino = sin(omega), coso = cos(omega);
	const double sqrt_1me2 = sqrt(1. - ecc*ecc);

	const npy_intp N = dims[0];

	if(circular)
	{
		#if defined (_OPENACC)
		#pragma acc parallel loop copyin(t_array[:N]) copyout(output_array[:N])
		#elif defined (_OPENMP)
		#pragma omp parallel for
		#endif
		for(int i = 0; i < N; i++)
		{
			double t = t_array[i];
			double f = ((t - tp)/per - (int)((t - tp)/per))*2.*M_PI;		//calculates f for a circular orbit
			if (f_only) { output_array[i] = f; continue; }
			double d;
			if (transittype == 1 && sin(f + omega)*sin(inc) <= 0.) d = BIGD; //z < 0, so d is set to large value in order to not model primary transit during secondary eclipse
			else if (transittype == 2 && sin(f + omega)*sin(inc) >= 0.) d = BIGD; //z > 0, so d is set to large value in order not to model secondary eclipse during primary transit
			else d = a*(1.0 - ecc*ecc)/(1.0 + ecc*cos(f))*sqrt(1.0 - sin(omega + f)*sin(omega + f)*sin(inc)*sin(inc));	//calculates separation of centers
			output_array[i] = d;
		}
	}
	else
	{
		/* Eccentric orbit: sin and cos of the true anomaly follow from those of the eccentric
		   anomaly, so neither tan, atan, cos(f) nor sin(f + omega) is evaluated:
		   cos f = (cos E - e)/(1 - e cos E), sin f = sqrt(1 - e^2) sin E/(1 - e cos E),
		   r = a (1 - e cos E). */
		#if defined (_OPENACC)
		#pragma acc parallel loop copyin(t_array[:N]) copyout(output_array[:N])
		#elif defined (_OPENMP)
		#pragma omp parallel for
		#endif
		for(npy_intp i0 = 0; i0 < N; i0 += KBLK)
		{
			int nb = (int)MIN((npy_intp)KBLK, N - i0);
			double M[KBLK], sE[KBLK], cE[KBLK];
			for(int j = 0; j < nb; j++) M[j] = Mphase(t_array[i0 + j], tp, per);
			kepler_block(M, nb, ecc, sE, cE);
			for(int j = 0; j < nb; j++)
			{
				double w = 1.0 - ecc*cE[j];
				if (f_only) { output_array[i0 + j] = atan2(sqrt_1me2*sE[j], cE[j] - ecc); continue; }
				double sinf = sqrt_1me2*sE[j]/w, cosf = (cE[j] - ecc)/w;
				double s = sinf*coso + cosf*sino;				//sin(f + omega)
				double ss = s*sini;
				double d;
				if (transittype == 1 && ss <= 0.) d = BIGD;
				else if (transittype == 2 && ss >= 0.) d = BIGD;
				else d = a*w*sqrt(MAX(1.0 - s*s*sini2, 0.));			//separation of centers
				output_array[i0 + j] = d;
			}
		}
	}
	return PyArray_Return((PyArrayObject *)ds);
}

static PyObject *_rsky(PyObject *self, PyObject *args)
{
	return _rsky_or_f(self, args, 0);
}



static PyObject *_getf(PyObject *self, PyObject *args)
{
	return _rsky_or_f(self, args, 1);
}


static char _rsky_doc[] = """ This module computes the distance between the centers of the \
star and the planet in the plane of the sky.  This parameter is \
denoted r_sky = sqrt(x^2 + y^2) in the Seager Exoplanets book \
(see the section by Murray, and Winn eq. 5).  In the Mandel & Agol (2002) paper, \
this quantity is denoted d.\
LK 4/27/12 """;


static char _getf_doc[] = """ This module computes the true anomaly. This parameter is \
denoted f in the Seager Exoplanets book \
(see the section by Murray, and Winn eq. 44).\
BM 1/18/16 """;


static PyMethodDef _rsky_methods[] = {
  {"_rsky", _rsky,METH_VARARGS,_rsky_doc},{"_getf", _getf,METH_VARARGS,_getf_doc},{NULL}};


#if PY_MAJOR_VERSION >= 3
	static struct PyModuleDef _rsky_module = {
		PyModuleDef_HEAD_INIT,
		"_rsky",
		_rsky_doc,
		-1,
		_rsky_methods
	};

	PyMODINIT_FUNC
	PyInit__rsky(void)
	{
		PyObject* module = PyModule_Create(&_rsky_module);
		if(!module)
		{
			return NULL;
		}
		import_array();
		return module;
	}
#else
	void init_rsky(void)
	{
	  Py_InitModule("_rsky", _rsky_methods);
	  import_array();
	}
#endif
