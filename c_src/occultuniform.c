/*
C  This routine computes the lightcurve for occultation
C  of a quadratically limb-darkened source without microlensing.

	This code is a translation of the Mandel/Agol Fortran code into C.
	LK 9/13/12
*/
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
//#include <omp.h>

#define IND(a,i) *((double *)(a->data+i*a->strides[0]))

static PyObject *occultuniform(PyObject *self, PyObject *args)
{
	int i;
	double z, p, kap0, kap1, pi = acos(-1.);

	PyArrayObject *zs, *muo1;
	npy_intp dims[1];
	
  	if(!PyArg_ParseTuple(args,"Od", &zs, &p)) return NULL;

	dims[0] = zs->dimensions[0];
	muo1 = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
	
	if(fabs(p-0.5)<1.e-3) p = 0.5;

	for(i=0; i<dims[0]; i++)
	{
		z = IND(zs, i);
		
		if(z >= 1.+p) IND(muo1,i) = 1.;
		if(p >= 1. && z <= p - 1.) IND(muo1,i) = 0.;
		else if(z <= 1.-p) IND(muo1,i) = 1.-p*p;
		else	
		{
			kap1=acos(fmin((1.-p*p+z*z)/2./z,1.));
			kap0=acos(fmin((p*p+z*z-1.)/2./p/z,1.));
			IND(muo1,i) = 1. - (p*p*kap0+kap1 -0.5*sqrt(fmax(4.*z*z-pow(1.+z*z-p*p, 2.), 0.)))/pi;
		}
	}

	return PyArray_Return(muo1);
}


static char occultuniform_doc[] = "LK 05/2015";

static PyMethodDef occultuniform_methods[] = {
  {"occultuniform", occultuniform, METH_VARARGS, occultuniform_doc},{NULL}};

void initoccultuniform(void)
{
  Py_InitModule("occultuniform", occultuniform_methods);
  import_array();
}

