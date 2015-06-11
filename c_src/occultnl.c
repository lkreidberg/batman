#include <Python.h>
#include<numpy/arrayobject.h>
#include<math.h>
#include<stdio.h>

#define TWOPI 6.28318531		//FIXME more precise!
#define PI 3.14159265
#define IND(a,i) *((double *)(a->data+i*a->strides[0]))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

static PyObject *occultnl(PyObject *self, PyObject *args);

double intensity_nl(double r, double u1, double u2, double u3, double u4)
{
	if(r > 0.99995) r = 0.99995;
	double sqrtmu = pow(1.-r*r,0.25);
	double norm = (-u1/10.-u2/6.-3.*u3/14.-u4/4.+0.5)*TWOPI; 		//calculate norm by integrating I(r)r dr dtheta
	return (1. - u1*(1.-sqrtmu) - u2*(1. - pow(sqrtmu,2.)) - u3*(1.-pow(sqrtmu, 3.)) - u4*(1.-pow(sqrtmu,4.)))/norm; 	
}

double area(double d, double r, double R)
{
	double arg1 = (d*d+r*r-R*R)/(2.*d*r), arg2 = (d*d+R*R-r*r)/(2.*d*R), arg3 = MAX((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R), 0.);

	if(r<=R-d) return PI*r*r;
	else if(r>=R+d) return PI*R*R;
	else return r*r*acos(arg1) + R*R*acos(arg2) - 0.5*sqrt(arg3);
}

static PyObject *occultnl(PyObject *self, PyObject *args)
{
	double rprs, d, fac, A_i, r, I; 
	int i, n;
	double dr, A_f, r_in, r_out, delta, u1, u2, u3, u4;
	
	npy_intp dims[1];
	PyArrayObject *zs, *flux;
  	if(!PyArg_ParseTuple(args,"Odddddd", &zs, &rprs, &u1, &u2, &u3, &u4, &fac)) return NULL;
	
	dims[0] = zs->dimensions[0];

	flux = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_DOUBLE);

	for(i = 0; i < dims[0]; i++)
	{
		d = IND(zs, i);
		r_in = MAX(d - rprs, 0.);		
		r_out = MIN(d + rprs, 1.0);		

		if(r_in >= 1.) IND(flux, i) = 1.;
		else
		{
			delta = 0.;

			r = r_in;
			dr = fac*acos(r);  
			r += dr;
			A_i = 0.;
	
			n = 0;
			while(r<r_out)
			{
				A_f = area(d, r, rprs);
				I = intensity_nl(r-dr/2.,u1,u2, u3, u4); 
				delta += (A_f - A_i)*I;
				dr = fac*acos(r);  
				r = r+dr;
				A_i = A_f;
			}
			dr = r_out -r + dr;  
			r = r_out;
			A_f = area(d, r, rprs);
			I = intensity_nl(r-dr/2.,u1,u2, u3, u4); 
			delta += (A_f - A_i)*I;

			IND(flux, i) = 1. - delta;	
		}
	}

	return PyArray_Return(flux);
} 

static char occultnl_doc[] = "docstring";

static PyMethodDef occultnl_methods[] = {
  {"occultnl", occultnl, METH_VARARGS, occultnl_doc},{NULL}};

void initoccultnl(void)
{
  Py_InitModule("occultnl", occultnl_methods);
  import_array();
}

