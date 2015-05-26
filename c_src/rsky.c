/*
	This code computes the separation of centers between the occulting and occulted objects.  
	It takes as arguments an array of times t,
	eccentricity e, semi-major axis a (in units of stellar radii), inclination angle i,
	stellar radius r_s, longitude of periastron w, orbital period P, time of periastron t0, and
	an error tolerance eps for computing the eccentric anomaly.

	LK 9/14/12 
*/
#include <Python.h>
#include<numpy/arrayobject.h>
#include<math.h>

#define TWOPI 6.28318531
#define PI 3.14159265
#define IND(a,i) *((double *)(a->data+i*a->strides[0]))

static PyObject *rsky(PyObject *self, PyObject *args);

double getE(double M, double e)	//determines the eccentric anomaly (Seager Exoplanets book:  Murray & Correia eqn. 5 -- see section 3)
{
	double E = M, eps = 1.0e-7;
	
	while(fabs(E - e*sin(E) - M) > eps)
	{
		E = E - (E - e*sin(E) - M)/(1.0 - e*cos(E));
	}
	return E;
}

static PyObject *rsky(PyObject *self, PyObject *args)
{
	double ecc, E, inc, a, r, d, f, omega, per, M, n, t0, eps;
	int i;
	npy_intp dims[1];
	PyArrayObject *zs, *t;
  	if(!PyArg_ParseTuple(args,"Odddddd", &t, &t0, &per, &a, &inc, &ecc, &omega)) return NULL; 
	dims[0] = t->dimensions[0];

	zs = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_DOUBLE);

	n = TWOPI/per;	// mean motion
	eps = 1.0e-7;
	
	for(i = 0; i < dims[0]; i++)
	{
		if(ecc != 0.)
		{
			M = n*(IND(t,i) - t0);
			E = getE(M, ecc);
			r = a*(1.0 - ecc*cos(E));
			f = acos(a*(1.0 - ecc*ecc)/(r*ecc) - 1.0/ecc);
			if(fabs((a*(1.0 - ecc*ecc)/(r*ecc) -1.0/ecc) - 1.0) < eps) f = 0.0;
		}
		else f = ((IND(t, i)-t0)/per - (int)((IND(t,i)-t0)/per))*TWOPI;
		d = a*(1.0-ecc*ecc)/(1.0+ecc*cos(f))*sqrt(1.0 - sin(omega+f)*sin(omega+f)*sin(inc)*sin(inc));
		IND(zs, i) = d;
	}
	return PyArray_Return(zs);
} 

static char rsky_doc[] = "\
This code computes the distance between the centers of the\n\
star and the planet in the plane of the sky.  This parameter is\n\ 
denoted r_sky = sqrt(x^2 + y^2) in the Seager Exoplanets book\n\
(see the section by Murray, and Winn eq. 5).  In the Mandel & Agol (2002) paper,\n\
this quantity is denoted d.\n\
LK 4/27/12 ";

static PyMethodDef rsky_methods[] = {
  {"rsky", rsky,METH_VARARGS,rsky_doc},{NULL}};

void initrsky(void)
{
  Py_InitModule("rsky", rsky_methods);
  import_array();
}

