#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>

#if defined (_OPENMP)
#  include <omp.h>
#endif

static PyObject *_uniform_ld(PyObject *self, PyObject *args)
{
	int i, nthreads;
	double z, p, kap0, kap1, pi = acos(-1.);

	PyArrayObject *zs, *flux;
	npy_intp dims[1], idx;
	
  	if(!PyArg_ParseTuple(args,"Odi", &zs, &p, &nthreads)) return NULL;		//parses function input

	dims[0] = PyArray_DIMS(zs)[0]; 
	flux = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(zs));		//creates numpy array to store output flux values
	
	if(fabs(p-0.5)<1.e-3) p = 0.5;

	#if defined (_OPENMP)
	omp_set_num_threads(nthreads);
	#endif

	#if defined (_OPENMP)
	#pragma omp parallel for private(z, kap1, kap0)
	#endif
	for(i=0; i<dims[0]; i++)
	{
		idx = (npy_intp)i;
		z = *(double*)PyArray_GetPtr(zs, &idx);
		
		if(z >= 1.+p) *(double*)PyArray_GetPtr(flux, &idx) = 1.;		//no overlap
		if(p >= 1. && z <= p - 1.) *(double*)PyArray_GetPtr(flux, &idx) = 0.;	//total eclipse of the star
		else if(z <= 1.-p) *(double*)PyArray_GetPtr(flux, &idx) = 1.-p*p;	//planet is fully in transit
		else									//planet is crossing the limb
		{
			kap1=acos(fmin((1.-p*p+z*z)/2./z,1.));
			kap0=acos(fmin((p*p+z*z-1.)/2./p/z,1.));
			*(double*)PyArray_GetPtr(flux, &idx) = 1. - (p*p*kap0+kap1 -0.5*sqrt(fmax(4.*z*z-pow(1.+z*z-p*p, 2.), 0.)))/pi;
		}
	}

	return PyArray_Return(flux);
}


static char _uniform_ld_doc[] = "This extension module returns a limb darkened light curve for a uniform stellar intensity profile.";

static PyMethodDef _uniform_ld_methods[] = {
  {"_uniform_ld", _uniform_ld, METH_VARARGS, _uniform_ld_doc},{NULL}};

void init_uniform_ld(void)
{
  Py_InitModule("_uniform_ld", _uniform_ld_methods);
  import_array();
}

