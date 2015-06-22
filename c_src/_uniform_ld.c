#include <Python.h>
#include <math.h>

#if defined (_OPENMP)
#  include <omp.h>
#endif

static PyObject *_uniform_ld(PyObject *self, PyObject *args)
{
	int nthreads;
	double z, p, kap0, kap1;

	PyObject *zs, *flux;
	Py_ssize_t i, dims;
	
  	if(!PyArg_ParseTuple(args,"Odi", &zs, &p, &nthreads)) return NULL;		//parses function input

	dims = PyList_Size(zs);
	flux = PyList_New(dims);	//creates numpy array to store return flux values
	
	if(fabs(p-0.5)<1.e-3) p = 0.5;

	#if defined (_OPENMP)
	omp_set_num_threads(nthreads);
	#endif

	#if defined (_OPENMP)
	#pragma omp parallel for private(z, kap1, kap0)
	#endif
	for(i=0; i<dims; i++)
	{
		z = PyFloat_AsDouble(PyList_GetItem(zs,i)); 
		
		if(z >= 1.+p) PyList_SetItem(flux, i, Py_BuildValue("d", 1.));		//no overlap
		if(p >= 1. && z <= p - 1.) PyList_SetItem(flux, i, Py_BuildValue("d", 0.));	//total eclipse of the star
		else if(z <= 1.-p) PyList_SetItem(flux, i, Py_BuildValue("d", 1.-p*p));	//planet is fully in transit
		else									//planet is crossing the limb
		{
			kap1=acos(fmin((1.-p*p+z*z)/2./z,1.));
			kap0=acos(fmin((p*p+z*z-1.)/2./p/z,1.));
			PyList_SetItem(flux, i, Py_BuildValue("d", 1. - (p*p*kap0+kap1 -0.5*sqrt(fmax(4.*z*z-pow(1.+z*z-p*p, 2.), 0.)))/M_PI));
		}
	}

	return Py_BuildValue("O", flux); 
}


static char _uniform_ld_doc[] = "This extension module returns a limb darkened light curve for a uniform stellar intensity profile.";

static PyMethodDef _uniform_ld_methods[] = {
  {"_uniform_ld", _uniform_ld, METH_VARARGS, _uniform_ld_doc},{NULL}};

#if PY_MAJOR_VERSION >= 3
	static struct PyModuleDef _uniform_ld_module = {
		PyModuleDef_HEAD_INIT,
		"_uniform_ld",
		_uniform_ld_doc,
		-1, 
		_uniform_ld_methods
	};

	PyMODINIT_FUNC
	PyInit__uniform_ld(void)
	{
		return PyModule_Create(&_uniform_ld_module);
	}
#else

	void init_uniform_ld(void)
	{
	  Py_InitModule("_uniform_ld", _uniform_ld_methods);
	}
#endif

