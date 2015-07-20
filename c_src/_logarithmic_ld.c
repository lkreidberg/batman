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
#include <math.h>

#if defined (_OPENMP)
#  include <omp.h>
#endif

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

static PyObject *_logarithmic_ld(PyObject *self, PyObject *args);

double intensity(double r, double u1, double u2, double norm)
{
	if(r > 0.99995) r = 0.99995;
	double mu = sqrt(1. - r*r);
	return (1. - u1*(1.-mu) - u2*mu*log(mu))/norm; 
}

double area(double d, double r, double R)
{
	/*
	Returns area of overlapping circles with radii r and R; separated by a distance d
	*/
	double arg1 = (d*d + r*r - R*R)/(2.*d*r); 	
	double arg2 = (d*d + R*R - r*r)/(2.*d*R); 
	double arg3 = MAX((-d + r + R)*(d + r - R)*(d - r + R)*(d + r + R), 0.);

	if(r <= R - d) return M_PI*r*r;						//planet completely overlaps stellar circle
	else if(r >= R + d) return M_PI*R*R;					//stellar circle completely overlaps planet
	else return r*r*acos(arg1) + R*R*acos(arg2) - 0.5*sqrt(arg3);		//partial overlap
}

static PyObject *_logarithmic_ld(PyObject *self, PyObject *args)
{
	double rprs, d, fac, A_i, r, I; 
	int nthreads;
	npy_intp i, dims[1];
	double dr, A_f, r_in, r_out, delta, u1, u2;
	
	PyArrayObject *zs, *flux;
  	if(!PyArg_ParseTuple(args,"Oddddi", &zs, &rprs, &u1, &u2, &fac, &nthreads)) return NULL;

	dims[0] = PyArray_DIMS(zs)[0]; 
	flux = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(zs));	//creates numpy array to store return flux values

	double *f_array = PyArray_DATA(flux);
	double *z_array = PyArray_DATA(zs);

	/*
		NOTE:  the safest way to access numpy arrays is to use the PyArray_GETITEM and PyArray_SETITEM functions.
		Here we use a trick for faster access and more convenient access, where we set a pointer to the 
		beginning of the array with the PyArray_DATA (e.g., f_array) and access elements with e.g., f_array[i].
		Success of this operation depends on the numpy array storing data in blocks equal in size to a C double.
		If you run into trouble along these lines, I recommend changing the array access to something like:
			d = PyFloat_AsDouble(PyArray_GETITEM(zs, PyArray_GetPtr(zs, &i))); 
		where zs is a numpy array object.


		Laura Kreidberg 07/2015
	*/
	
	#if defined (_OPENMP)
	omp_set_num_threads(nthreads);	//specifies number of threads (if OpenMP is supported)
	#endif

	double norm = (-3.*u1 + 2.*u2 + 9.)*M_PI/9.0; 	//normalization for intensity profile (faster to calculate it once, rather than every time intensity is called)		

	#if defined (_OPENMP)
	#pragma omp parallel for private(d, r_in, r_out, delta, r, dr, A_i, A_f, I)
	#endif
	for(i = 0; i < dims[0]; i++)
	{
		d = z_array[i];
		r_in = MAX(d - rprs, 0.);					//lower bound for integration
		r_out = MIN(d + rprs, 1.0);					//upper bound for integration

		if(r_in >= 1.) f_array[i] = 1.0;				//flux = 1. if the planet is not transiting
		else
		{
			delta = 0.;						//variable to store the integrated intensity, \int I dA
			r = r_in;						//starting radius for integration
			dr = fac*acos(r); 					//initial step size 

			r += dr;						//first step

			A_i = 0.;						//initial area
	
			while(r < r_out)
			{
				A_f = area(d, r, rprs);				//calculates area of overlapping circles
				I = intensity(r - dr/2., u1, u2, norm); 	//intensity at the midpoint
				delta += (A_f - A_i)*I;				//increase in transit depth for this integration step
				dr = fac*acos(r);  				//updating step size
				r = r + dr;					//stepping to next element
				A_i = A_f;					//storing area
			}
			dr = r_out - r + dr;  					//calculating change in radius for last step  
			r = r_out;						//final radius for integration
			A_f = area(d, r, rprs);					//area for last integration step
			I = intensity(r - dr/2., u1, u2, norm); 		//intensity at the midpoint 
			delta += (A_f - A_i)*I;					//increase in transit depth for this integration step

			f_array[i] = 1.0 - delta;	//flux equals 1 - \int I dA 
		}
	}
	return PyArray_Return((PyArrayObject *)flux);

} 

static char _logarithmic_ld_doc[] = "This extension module returns a limb darkened light curve for a logarithmic stellar intensity profile.";

static PyMethodDef _logarithmic_ld_methods[] = {
  {"_logarithmic_ld", _logarithmic_ld, METH_VARARGS, _logarithmic_ld_doc},{NULL}};

#if PY_MAJOR_VERSION >= 3
	static struct PyModuleDef _logarithmic_ld_module = {
		PyModuleDef_HEAD_INIT,
		"_logarithmic_ld",
		_logarithmic_ld_doc,
		-1, 
		_logarithmic_ld_methods
	};

	PyMODINIT_FUNC
	PyInit__logarithmic_ld(void)
	{
		PyObject* module = PyModule_Create(&_logarithmic_ld_module);
		if(!module)
		{
			return NULL;
		}
		import_array(); 
		return module;
	}
#else

	void init_logarithmic_ld(void)
	{
	  	Py_InitModule("_logarithmic_ld", _logarithmic_ld_methods);
		import_array(); 
	}
#endif

