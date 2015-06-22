#include <Python.h>
#include<math.h>

#if defined (_OPENMP)
#  include <omp.h>
#endif

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

static PyObject *_nonlinear_ld(PyObject *self, PyObject *args);

double intensity(double r, double u1, double u2, double u3, double u4, double norm)
{
	if(r > 0.99995) r = 0.99995;
	double sqrtmu = pow(1.-r*r,0.25);
	return (1. - u1*(1.-sqrtmu) - u2*(1. - pow(sqrtmu,2.)) - u3*(1.-pow(sqrtmu, 3.)) - u4*(1.-pow(sqrtmu,4.)))/norm; 	
}

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

static PyObject *_nonlinear_ld(PyObject *self, PyObject *args)
{
	double rprs, d, fac, A_i, r, I; 
	int nthreads;
	Py_ssize_t i, dims;
	double dr, A_f, r_in, r_out, delta, u1, u2, u3, u4;
	
	PyObject *zs, *flux;
  	if(!PyArg_ParseTuple(args,"Oddddddi", &zs, &rprs, &u1, &u2, &u3, &u4, &fac, &nthreads)) return NULL;

	dims = PyList_Size(zs);
	flux = PyList_New(dims);	//creates numpy array to store return flux values

	#if defined (_OPENMP)
	omp_set_num_threads(nthreads);	//specifies number of threads (if OpenMP is supported)
	#endif

	double norm = (-u1/10.-u2/6.-3.*u3/14.-u4/4.+0.5)*2.*M_PI; 	//normalization for intensity profile (faster to calculate it once, rather than every time intensity is called)		

	#if defined (_OPENMP)
	#pragma omp parallel for private(d, r_in, r_out, delta, r, dr, A_i, A_f, I)
	#endif
	for(i = 0; i < dims; i++)
	{
		d = PyFloat_AsDouble(PyList_GetItem(zs,i)); 
		r_in = MAX(d - rprs, 0.);					//lower bound for integration
		r_out = MIN(d + rprs, 1.0);					//upper bound for integration

		if(r_in >= 1.) PyList_SetItem(flux, i, Py_BuildValue("d", 1.0));	//flux = 1. if the planet is not transiting
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
				I = intensity(r-dr/2.,u1,u2, u3, u4, norm); 	//intensity at the midpoint
				delta += (A_f - A_i)*I;				//increase in transit depth for this integration step
				dr = fac*acos(r);  				//updating step size
				r = r+dr;					//stepping to next element
				A_i = A_f;					//storing area
			}
			dr = r_out -r + dr;  					//calculating change in radius for last step
			r = r_out;						//final radius for integration
			A_f = area(d, r, rprs);					//area for last integration step
			I = intensity(r-dr/2.,u1,u2, u3, u4, norm); 		//intensity at the midpoint 
			delta += (A_f - A_i)*I;					//increase in transit depth for this integration step

			PyList_SetItem(flux, i, Py_BuildValue("d", 1.0-delta));	//flux equals 1 - \int I dA 
		}
	}
	return Py_BuildValue("O", flux); 

} 

static char _nonlinear_ld_doc[] = "This extension module returns a limb darkened light curve for a nonlinear stellar intensity profile.";

static PyMethodDef _nonlinear_ld_methods[] = {
  {"_nonlinear_ld", _nonlinear_ld, METH_VARARGS, _nonlinear_ld_doc},{NULL}};

#if PY_MAJOR_VERSION >= 3
	static struct PyModuleDef _nonlinear_ld_module = {
		PyModuleDef_HEAD_INIT,
		"_nonlinear_ld",
		_nonlinear_ld_doc,
		-1, 
		_nonlinear_ld_methods
	};

	PyMODINIT_FUNC
	PyInit__nonlinear_ld(void)
	{
		return PyModule_Create(&_nonlinear_ld_module);
	}
#else

	void init_nonlinear_ld(void)
	{
	  Py_InitModule("_nonlinear_ld", _nonlinear_ld_methods);
	}
#endif

