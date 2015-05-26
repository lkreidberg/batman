#include <Python.h>
#include<numpy/arrayobject.h>
#include<math.h>
#include<stdio.h>

#define TWOPI 6.28318531
#define PI 3.14159265
#define IND(a,i) *((double *)(a->data+i*a->strides[0]))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

static PyObject *occultnl(PyObject *self, PyObject *args);

double intensity_quad(double r, double u1, double u2)
{
	double mu = sqrt(1.-r*r), norm = TWOPI*(-2.*u1-u2+6.)/12.;
	return (1. - u1*(1.-mu) - u2*(1. - mu)*(1 - mu))/norm; 	
}

double intensity_nl(double r, double u1, double u2, double u3, double u4)
{
	double sqrtmu = pow(1.-r*r,0.25);
	double norm = (-u1/10.-u2/6.-3.*u3/14.-u4/4.+0.5)*TWOPI; 		//calculate norm by integrating I(r)r dr dtheta
	return (1. - u1*(1.-sqrtmu) - u2*(1. - pow(sqrtmu,2.)) - u3*(1.-pow(sqrtmu, 3.)) - u4*(1.-pow(sqrtmu,4.)))/norm; 	
}

double area(double d, double r, double R)
{
	double arg1 = (d*d+r*r-R*R)/(2.*d*r), arg2 = (d*d+R*R-r*r)/(2.*d*R), arg3 = MAX((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R), 0.);
	if(arg2>1.) arg2 = 1.;
	if(arg2<-1.) arg2 = -1.0;
	if(arg1>1.) arg1 = 1.;
	if(arg1<-1.) arg1 = -1.0;
	
	return r*r*acos(arg1) + R*R*acos(arg2) - 0.5*sqrt(arg3);
}

static PyObject *occultnl(PyObject *self, PyObject *args)
{
	double rprs, d, fac, temp, r, I_ave; 
	int i, n;
	double dr, dA, r_in, r_out, delta, u1, u2, u3, u4;
	
	npy_intp dims[1];
	PyArrayObject *zs, *flux;
  	if(!PyArg_ParseTuple(args,"Odddddd", &zs, &rprs, &u1, &u2, &u3, &u4, &fac)) return NULL;
	
	dims[0] = zs->dimensions[0];

	flux = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
	

	for(i = 0; i < dims[0]; i++)
	{
		d = IND(zs, i);
		r_in = d - rprs;		
		r_out = MIN(d + rprs, 1.0);		
		if(r_in >= 1.) IND(flux, i) = 1.;
		else if (d == 0.)
		{
			r_in = 0;
			r_out = rprs;
			r = r_out;
			delta = 0.;
			while(r>0.)
			{
				dr = fac*acos(r);
				if(r-dr <0.) break;
				I_ave = (intensity_nl(r,u1,u2, u3, u4) + 4.*intensity_nl(r-dr/2., u1, u2, u3, u4) + intensity_nl(r-dr, u1, u2, u3, u4))/6.;
			//	I_ave = intensity_nl(r-dr/2.,u1,u2, u3, u4); 
				delta += I_ave*(r-dr/2.)*dr*TWOPI;
				r = r-dr;
			}
			dr = r;
			I_ave = (intensity_nl(r,u1,u2, u3, u4) + 4.*intensity_nl(r-dr/2., u1, u2, u3, u4) + intensity_nl(r-dr, u1, u2, u3, u4))/6.;
			delta += I_ave*(r-dr/2.)*dr*TWOPI;
			IND(flux, i) = 1. - delta;	
		}
		else
		{
			delta = 0.;
			temp = 0.;			

			r = r_in;
			n = 0;
			while(r<r_out)
			{
				dr = fac*acos(r);  
				r = r+dr;
				if(r>1.)
				{
					r = r -dr;
					break;
				}
				dA = area(d, r, rprs);
				I_ave = (intensity_nl(r,u1,u2, u3, u4) + 4.*intensity_nl(r-dr/2., u1, u2, u3, u4) + intensity_nl(r-dr, u1, u2, u3, u4))/6.;
				//I_ave = intensity_nl(r-dr/2.,u1,u2, u3, u4); 
				delta += (dA - temp)*I_ave;
				temp = dA;
			}
			dr = r_out - r;
			r = r+dr;
			dA = area(d, r, rprs);
			I_ave = (intensity_nl(r,u1,u2, u3, u4) + 4.*intensity_nl(r-dr/2., u1, u2, u3, u4) + intensity_nl(r-dr, u1, u2, u3, u4))/6.;
			delta += (dA - temp)*I_ave;

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

