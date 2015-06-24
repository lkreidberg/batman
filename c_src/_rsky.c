#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include "numpy/arrayobject.h"
#include<math.h>

static PyObject *_rsky(PyObject *self, PyObject *args);

double getE(double M, double e)	//calculates the eccentric anomaly (see Seager Exoplanets book:  Murray & Correia eqn. 5 -- see section 3)
{
	double E = M, eps = 1.0e-7;
	
	while(fabs(E - e*sin(E) - M) > eps)
	{
		E = E - (E - e*sin(E) - M)/(1.0 - e*cos(E));
	}
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
	double ecc, E, inc, a, r, d, f, omega, per, M, n, t0, eps, t;
	//Py_ssize_t i, dims;
	npy_intp i, dims[1];
	PyArrayObject *ts, *zs;

  	if(!PyArg_ParseTuple(args,"Odddddd", &ts, &t0, &per, &a, &inc, &ecc, &omega)) return NULL; 

	dims[0] = PyArray_DIMS(ts)[0]; 
	zs = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(ts));
	
	double *t_array = PyArray_DATA(ts);
	double *z_array = PyArray_DATA(zs);

	n = 2.*M_PI/per;	// mean motion
	eps = 1.0e-7;
	
	for(i = 0; i < dims[0]; i++)
	{
		t = t_array[i];

		if(ecc > 1.e-5)											//calculates f for eccentric orbits
		{
			M = n*(t - t0);
			E = getE(M, ecc);
			r = a*(1.0 - ecc*cos(E));
			f = acos(a*(1.0 - ecc*ecc)/(r*ecc) - 1.0/ecc);
			if(fabs((a*(1.0 - ecc*ecc)/(r*ecc) -1.0/ecc) - 1.0) < eps) f = 0.0;
		}
		else f = ((t-t0)/per - (int)((t-t0)/per))*2.*M_PI;						//calculates f for a circular orbit
		d = a*(1.0-ecc*ecc)/(1.0+ecc*cos(f))*sqrt(1.0 - sin(omega+f)*sin(omega+f)*sin(inc)*sin(inc));	//calculates separation of centers 
		z_array[i] = d;
	}
	return PyArray_Return((PyArrayObject *)zs);
} 

static char _rsky_doc[] = """ This module computes the distance between the centers of the \
star and the planet in the plane of the sky.  This parameter is \
denoted r_sky = sqrt(x^2 + y^2) in the Seager Exoplanets book \
(see the section by Murray, and Winn eq. 5).  In the Mandel & Agol (2002) paper, \
this quantity is denoted d.\
LK 4/27/12 """;

static PyMethodDef _rsky_methods[] = {
  {"_rsky", _rsky,METH_VARARGS,_rsky_doc},{NULL}};


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

