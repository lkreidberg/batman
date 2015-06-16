#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>

#if defined (_OPENMP)
#  include <omp.h>
#endif

static PyObject *occultuniform(PyObject *self, PyObject *args)
{
	int i, nthreads;
	double z, p, kap0, kap1, pi = acos(-1.);

	PyArrayObject *zs, *flux;
	npy_intp dims[1], idx;
	
  	if(!PyArg_ParseTuple(args,"Odi", &zs, &p, &nthreads)) return NULL;

	dims[0] = PyArray_DIMS(zs)[0]; 
	flux = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(zs));
	
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
		
		if(z >= 1.+p) *(double*)PyArray_GetPtr(flux, &idx) = 1.;
		if(p >= 1. && z <= p - 1.) *(double*)PyArray_GetPtr(flux, &idx) = 0.;
		else if(z <= 1.-p) *(double*)PyArray_GetPtr(flux, &idx) = 1.-p*p;
		else	
		{
			kap1=acos(fmin((1.-p*p+z*z)/2./z,1.));
			kap0=acos(fmin((p*p+z*z-1.)/2./p/z,1.));
			*(double*)PyArray_GetPtr(flux, &idx) = 1. - (p*p*kap0+kap1 -0.5*sqrt(fmax(4.*z*z-pow(1.+z*z-p*p, 2.), 0.)))/pi;
		}
	}

	return PyArray_Return(flux);
}


static char occultuniform_doc[] = "LK 05/2015";

static PyMethodDef occultuniform_methods[] = {
  {"occultuniform", occultuniform, METH_VARARGS, occultuniform_doc},{NULL}};

void initoccultuniform(void)
{
  Py_InitModule("occultuniform", occultuniform_methods);
  import_array();
}

