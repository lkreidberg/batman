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

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include "numpy/arrayobject.h"
#include<math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

static PyObject *_rsky(PyObject *self, PyObject *args);

static PyObject *_getf(PyObject *self, PyObject *args);

inline double getE(double M, double e)	//calculates the eccentric anomaly (see Seager Exoplanets book:  Murray & Correia eqn. 5 -- see section 3)
{
	double E = M, eps = 1.0e-7;

	while(fabs(E - e*sin(E) - M) > eps) E = E - (E - e*sin(E) - M)/(1.0 - e*cos(E));
	return E;
}

static PyObject *_rsky(PyObject *self, PyObject *args)
{
	/*
		This module computes the distance between the centers of the 
		star and the planet in the plane of the sky.  This parameter is 
		denoted r_sky = sqrt(x^2 + y^2) in the Seager Exoplanets book 
		(see the section by Murray, and Winn eq. 5).  In the Mandel & Agol 
		(2002) paper, this quantity is denoted d.
	*/
	double ecc, inc, a, omega, per, tc, BIGD = 100.;
	int transittype;
	npy_intp dims[1];
	PyArrayObject *ts, *ds;

  	if(!PyArg_ParseTuple(args,"Oddddddi", &ts, &tc, &per, &a, &inc, &ecc, &omega, &transittype)) return NULL; 

	dims[0] = PyArray_DIMS(ts)[0]; 
	ds = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(ts));
	
	double *t_array = PyArray_DATA(ts);
	double *d_array = PyArray_DATA(ds);

	const double n = 2.*M_PI/per;	// mean motion
	const double eps = 1.0e-7;

	#pragma acc parallel loop
	for(int i = 0; i < dims[0]; i++)
	{
		double d;
		double t = t_array[i];
		
		//calculates time of periastron passage from time of inferior conjunction 
		double f = M_PI/2. - omega;								//true anomaly corresponding to time of primary transit center
		double E = 2.*atan(sqrt((1. - ecc)/(1. + ecc))*tan(f/2.));				//corresponding eccentric anomaly
		double M = E - ecc*sin(E);
		double tp = tc - per*M/2./M_PI;							//time of periastron

		if(ecc < 1.0e-5)
		{
			f = ((t - tp)/per - (int)((t - tp)/per))*2.*M_PI;			//calculates f for a circular orbit
		}
		else
		{
			M = n*(t - tp);
			E = getE(M, ecc);
			f = 2.*atan(sqrt((1.+ecc)/(1.-ecc))*tan(E/2.));
		}
		if (transittype == 1 && sin(f + omega)*sin(inc) <= 0.) d = BIGD;						//z < 0, so d is set to large value in order to not model primary transit during secondary eclipse
		else if (transittype == 2 && sin(f + omega)*sin(inc) >= 0.) d = BIGD;						//z > 0, so d is set to large value in order not to model secondary eclipse during primary transit
		else d = a*(1.0 - ecc*ecc)/(1.0 + ecc*cos(f))*sqrt(1.0 - sin(omega + f)*sin(omega + f)*sin(inc)*sin(inc));	//calculates separation of centers 
		d_array[i] = d;
	}
	return PyArray_Return((PyArrayObject *)ds);
} 



static PyObject *_getf(PyObject *self, PyObject *args)
{
	/*
		This module computes the distance between the centers of the
		star and the planet in the plane of the sky.  This parameter is
		denoted r_sky = sqrt(x^2 + y^2) in the Seager Exoplanets book
		(see the section by Murray, and Winn eq. 5).  In the Mandel & Agol
		(2002) paper, this quantity is denoted d.
	*/
	double ecc, E, inc, a, f, omega, per, M, n, tp, tc, eps, t;
	int transittype;
	npy_intp i, dims[1];
	PyArrayObject *ts, *fs;

  	if(!PyArg_ParseTuple(args,"Oddddddi", &ts, &tc, &per, &a, &inc, &ecc, &omega, &transittype)) return NULL;

	dims[0] = PyArray_DIMS(ts)[0];
	fs = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(ts));

	double *t_array = PyArray_DATA(ts);
	double *f_array = PyArray_DATA(fs);

	n = 2.*M_PI/per;	// mean motion
	eps = 1.0e-7;

	for(i = 0; i < dims[0]; i++)
	{
		t = t_array[i];

		//calculates time of periastron passage from time of inferior conjunction
		f = M_PI/2. - omega;								//true anomaly corresponding to time of primary transit center
		E = 2.*atan(sqrt((1. - ecc)/(1. + ecc))*tan(f/2.));				//corresponding eccentric anomaly
		M = E - ecc*sin(E);
		tp = tc - per*M/2./M_PI;							//time of periastron

		if(ecc < 1.0e-5)
		{
			f = ((t - tp)/per - (int)((t - tp)/per))*2.*M_PI;			//calculates f for a circular orbit
		}
		else
		{
			M = n*(t - tp);
			E = getE(M, ecc);
			f = 2.*atan(sqrt((1.+ecc)/(1.-ecc))*tan(E/2.));
		}
		f_array[i] = f;
	}
	return PyArray_Return((PyArrayObject *)fs);
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

