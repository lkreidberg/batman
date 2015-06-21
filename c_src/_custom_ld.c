#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include<numpy/arrayobject.h>
#include<math.h>
#include<stdio.h>

#if defined (_OPENMP)
#  include <omp.h>
#endif

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

static PyObject *_custom_ld(PyObject *self, PyObject *args);

double intensity(double r, double u1, double u2, double u3, double u4, double u5, double u6);

double area(double d, double r, double R)
{
	/*
	Returns area of overlapping circles with radii r and R; separated by a distance d
	*/
	double arg1 = (d*d+r*r-R*R)/(2.*d*r), arg2 = (d*d+R*R-r*r)/(2.*d*R), arg3 = MAX((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R), 0.);

	if(r<=R-d) return M_PI*r*r;						//planet completely overlaps stellar circle
	else if(r>=R+d) return M_PI*R*R;					//stellar circle completely overlaps planet
	else return r*r*acos(arg1) + R*R*acos(arg2) - 0.5*sqrt(arg3);		//partial overlap
}

static PyObject *_custom_ld(PyObject *self, PyObject *args)
{
	double rprs, d, fac, A_i, r, I, dr, A_f, r_in, r_out, delta, u1, u2, u3, u4, u5, u6;
	int i, nthreads;
	
	npy_intp dims[1], idx;
	PyArrayObject *zs, *flux;
  	if(!PyArg_ParseTuple(args,"Oddddddddi", &zs, &rprs, &u1, &u2, &u3, &u4, &u5, &u6, &fac, &nthreads)) return NULL; //parses input arguments
	
	dims[0] = PyArray_DIMS(zs)[0]; 
	flux = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(zs));	//creates numpy array to store return flux values

	#if defined (_OPENMP)
	omp_set_num_threads(nthreads);	//specifies number of threads (if OpenMP is supported)
	#endif

	#if defined (_OPENMP)
	#pragma omp parallel for private(d, r_in, r_out, delta, r, dr, A_i, A_f, I)
	#endif
	for(i = 0; i < dims[0]; i++)
	{
		idx = (npy_intp)i;
		d = *(double*)PyArray_GetPtr(zs, &idx);				//separation of centers
		r_in = MAX(d - rprs, 0.);					//lower bound for integration
		r_out = MIN(d + rprs, 1.0);					//upper bound for integration

		if(r_in >= 1.) *(double*)PyArray_GetPtr(flux, &idx) = 1.;	//flux = 1. if the planet is not transiting
		else
		{
			delta = 0.;						//variable to store the integrated intensity, \int I dA

			r = r_in;						//starting radius for integration
			dr = fac*acos(r); 					//initial step size 
			r += dr;						//first step
			A_i = 0.;						//initial area
	
			while(r<r_out)
			{
				A_f = area(d, r, rprs);				//calculates area of overlapping circles
				I = intensity(r-dr/2.,u1,u2, u3, u4, u5, u6); 	//intensity at the midpoint
				delta += (A_f - A_i)*I;				//increase in transit depth for this integration step
				dr = fac*acos(r);  				//updating step size
				r = r+dr;					//stepping to next element
				A_i = A_f;					//storing area
			}
			dr = r_out -r + dr;  					//calculating change in radius for last step
			r = r_out;						//final radius for integration
			A_f = area(d, r, rprs);					//area for last integration step
			I = intensity(r-dr/2.,u1,u2, u3, u4, u5, u6); 		//intensity at the midpoint 
			delta += (A_f - A_i)*I;					//increase in transit depth for this integration step

		 	*(double*)PyArray_GetPtr(flux, &idx) = 1. - delta;	//flux equals 1 - \int I dA
		}
	}
	return PyArray_Return(flux);
} 

static char _custom_ld_doc[] = "This extension module returns a limb darkened light curve for a custom stellar intensity profile.";

static PyMethodDef _custom_ld_methods[] = {
  {"_custom_ld", _custom_ld, METH_VARARGS, _custom_ld_doc},{NULL}};

void init_custom_ld(void)
{
  Py_InitModule("_custom_ld", _custom_ld_methods);
  import_array();
}

