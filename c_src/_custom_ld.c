#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include<numpy/arrayobject.h>
#include<math.h>
#include<stdio.h>

#if defined (_OPENMP)
#  include <omp.h>
#endif

#define TWOPI 6.28318531		//FIXME more precise!
#define PI 3.14159265
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

static PyObject *_custom_ld(PyObject *self, PyObject *args);

double intensity(double r, double u1, double u2, double u3, double u4);

double area(double d, double r, double R)
{
	double arg1 = (d*d+r*r-R*R)/(2.*d*r), arg2 = (d*d+R*R-r*r)/(2.*d*R), arg3 = MAX((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R), 0.);

	if(r<=R-d) return PI*r*r;
	else if(r>=R+d) return PI*R*R;
	else return r*r*acos(arg1) + R*R*acos(arg2) - 0.5*sqrt(arg3);
}

static PyObject *_custom_ld(PyObject *self, PyObject *args)
{
	double rprs, d, fac, A_i, r, I, u1, u2, u3, u4; 
	int i, nthreads;
	double dr, A_f, r_in, r_out, delta;
	
	npy_intp dims[1], idx;
	PyArrayObject *zs, *flux, *u;
  	if(!PyArg_ParseTuple(args,"OdOdi", &zs, &rprs, &u, &fac, &nthreads)) return NULL;
	
	dims[0] = PyArray_DIMS(zs)[0]; 

	flux = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(zs));

	/*for(i = 0; i < (int)PyArray_DIMS(us); i++)
	{
		idx0 = (npy_intp)i;
		u[i] = *(double*)PyArray_GetPtr(us, &idx0);		//stores limb darkening coefficients in an array
		printf("u = %f\n", u[i]);
	}*/

	i = 0;
	idx = (npy_intp)i;
	u1 = *(double*)PyArray_GetPtr(u, &idx);
	i = 1;
	idx = (npy_intp)i;
	u2 = *(double*)PyArray_GetPtr(u, &idx);
	i = 2;
	idx = (npy_intp)i;
	u3 = *(double*)PyArray_GetPtr(u, &idx);
	i = 3;
	idx = (npy_intp)i;
	u4 = *(double*)PyArray_GetPtr(u, &idx);



	#if defined (_OPENMP)
	omp_set_num_threads(nthreads);
	#endif

	#if defined (_OPENMP)
	#pragma omp parallel for private(d, r_in, r_out, delta, r, dr, A_i, A_f, I)
	#endif
	for(i = 0; i < dims[0]; i++)
	{
		idx = (npy_intp)i;
		d = *(double*)PyArray_GetPtr(zs, &idx);
		r_in = MAX(d - rprs, 0.);		
		r_out = MIN(d + rprs, 1.0);		

		if(r_in >= 1.) *(double*)PyArray_GetPtr(flux, &idx) = 1.;
		else
		{
			delta = 0.;

			r = r_in;
			dr = fac*acos(r);  
			r += dr;
			A_i = 0.;
	
			while(r<r_out)
			{
				A_f = area(d, r, rprs);
				I = intensity(r-dr/2.,u1, u2, u3, u4); 
				delta += (A_f - A_i)*I;
				dr = fac*acos(r);  
				r = r+dr;
				A_i = A_f;
			}
			dr = r_out -r + dr;  
			r = r_out;
			A_f = area(d, r, rprs);
			I = intensity(r-dr/2.,u1, u2, u3, u4); 
			delta += (A_f - A_i)*I;

		 	*(double*)PyArray_GetPtr(flux, &idx) = 1. - delta;	
		}
	}
	
	return PyArray_Return(flux);
} 

static char _custom_ld_doc[] = "docstring";

static PyMethodDef _custom_ld_methods[] = {
  {"_custom_ld", _custom_ld, METH_VARARGS, _custom_ld_doc},{NULL}};

void init_custom_ld(void)
{
  Py_InitModule("_custom_ld", _custom_ld_methods);
  import_array();
}

