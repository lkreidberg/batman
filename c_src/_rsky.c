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

inline double getE(double M, double e)	//calculates the eccentric anomaly (see Seager Exoplanets book:  Murray & Correia eqn. 5 -- see section 3)
{
	/* Same iteration as before, with the repeated work removed: sin(E) is evaluated once per
	   iteration instead of twice, fmod(|x|, 2pi) = |fmod(x, 2pi)| is evaluated once, and fmod is
	   skipped where it is the identity (|x| < 2pi; 1 - e*cos(E) lies in (0, 2) for e < 1). */
	double E = M, eps = 1.0e-7;
	double fe, fs, sE, r, x;
	int ii, maxiter = 1000;
	const double twopi = 2.*M_PI;

	// modification from LK 05/07/2017:
	// add fmod to ensure convergence for diabolical inputs (following Eastman et al. 2013; Section 3.1)
	ii = 0;
	while (ii < maxiter)
	{
		sE = sin(E);
		x = E - e*sE - M;
		fe = (fabs(x) < twopi) ? x : fmod(x, twopi);
		if (!(fabs(fe) > eps)) break;
		fs = 1 - e*cos(E);
		fs = (fs >= 0. && fs < twopi) ? fs : fmod(fs, twopi);
		E = E - fe/fs;
		ii++;
	}
	if (ii >= maxiter)
	{
		PyErr_SetString(PyExc_RuntimeError, "Eccentric anomaly computation reached the maximum number of interations");
		return E;
	}
	return E;
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

	const double n = 2.*M_PI/per;	// mean motion
  const double eps = 1.0e-7;


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
	const double sqrt_ecc = sqrt((1.+ecc)/(1.-ecc));
	const double one_m_e2 = 1.0 - ecc*ecc;
	const int exactly_circular = (ecc == 0.0);					//then a(1-e^2)/(1+e cos f) = a exactly

	#if defined (_OPENACC)
	#pragma acc parallel loop copyin(t_array[:dims[0]]) copyout(output_array[:dims[0]])
	#elif defined (_OPENMP)
	#pragma omp parallel for
	#endif
	for(int i = 0; i < dims[0]; i++)
	{
		double t = t_array[i];
		double f;

		if(circular)
		{
			f = ((t - tp)/per - (int)((t - tp)/per))*2.*M_PI;			//calculates f for a circular orbit
		}
		else
		{
			double M = n*(t - tp);
			double E = getE(M, ecc);
			f = 2.*atan(sqrt_ecc*tan(E/2.));
		}
		if (f_only) {
			output_array[i] = f;
		}
		else {
			double d;
			double s = sin(f + omega);
			double ss = s*sini;
			if (transittype == 1 && ss <= 0.) d = BIGD; //z < 0, so d is set to large value in order to not model primary transit during secondary eclipse
			else if (transittype == 2 && ss >= 0.) d = BIGD; //z > 0, so d is set to large value in order not to model secondary eclipse during primary transit
			else
			{
				double r = exactly_circular ? a : a*one_m_e2/(1.0 + ecc*cos(f));
				d = r*sqrt(1.0 - s*s*sini*sini);	//calculates separation of centers
			}
			output_array[i] = d;
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
