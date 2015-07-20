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

double *occultuniform(PyArrayObject *z0, double p, int nz);

static PyObject *_mandelagol_nonlinear_ld(PyObject *self, PyObject *args)
{
/*	Comments here
*/
	int i, j, nz, nr, i1, i2, nmax = 513;
	double rl, *mulimbf, *zt0, *mulimb, *mulimb0, *mulimbp, dt=0., *t, *th, *r, sig, *mulimb1, *mulimbhalf, *mulimb3half, *mulimb2, sig1, sig2, omega, dmumax, fac, *mu, f1, f2;
	PyArrayObject *z0, *muo1;
	npy_intp dims[1];

	double pi = acos(-1.);
	double c1, c2, c3, c4, p;
  	if(!PyArg_ParseTuple(args,"Odddddi", &z0, &c1, &c2, &c3, &c4, &p, &nz)) return NULL;
	rl = p;

	dims[0] = z0->dimensions[0];
	muo1 = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
	
	mulimbf = (double *)malloc(5*nz*sizeof(double));
	zt0 = (double *)malloc(nz*sizeof(double));
	mulimb0 = (double *)malloc(nz*sizeof(double));
	mulimb = (double *)malloc(nz*sizeof(double));
	mulimbp = (double *)malloc(nz*sizeof(double));
	mulimbhalf = (double *)malloc(nz*sizeof(double));
	mulimb1 = (double *)malloc(nz*sizeof(double));
	mulimb3half = (double *)malloc(nz*sizeof(double));
	mulimb2 = (double *)malloc(nz*sizeof(double));
	mu = (double *)malloc(nz*sizeof(double));
	
	t = (double *)malloc(nmax*sizeof(double));
	th = (double *)malloc(nmax*sizeof(double));
	r = (double *)malloc(nmax*sizeof(double));

	mulimb0 = occultuniform(z0, p, nz);
	i1 = nz;
	i2 = 0; //???
	fac = 0.;
	
	for(i=0; i<nz; i++)
	{
		zt0[i] = IND(z0, i);
		mulimbf[i] = 1.;
		mulimbf[nz + i] = 0.8;
		mulimbf[2*nz + i] = 2./3.;
		mulimbf[3*nz + i] = 4./7.;
		mulimbf[4*nz + i] = 0.5;
		mulimb[i] = mulimb0[i];
		if(mulimb0[i] != 1.)
		{
			i1 = fmin(i1,i);
			i2 = fmax(i2,i);
		}
		fac = fmax(fac, fabs(mulimb0[i]-1.));
	}
	
	omega=4.*((1.-c1-c2-c3-c4)/4.+c1/5.+c2/6.+c3/7.+c4/8.);
	nr=2;
	dmumax=1.;

	while(dmumax > fac*2.e-4)
	//while(dmumax > fac*3.7e-5)
	{
		//printf("%0.10f\n", dmumax);
		for(i = i1; i <= i2; i++) mulimbp[i] = mulimb[i];	
		nr *= 2;
		dt = 0.5*pi/(double)nr;
		for(j = 0; j <= nr; j++)
		{
			t[j] = dt*(double)(j);
			th[j] = t[j]+0.5*dt;
			r[j] =sin(t[j]);
		}	
		sig = sqrt(cos(th[nr-1]));
		for(i = i1; i <= i2; i++) 
		{
			mulimbhalf[i] = sig*sig*sig*mulimb0[i]/(1. - r[nr-1]);
			mulimb1[i] = sig*mulimbhalf[i];
			mulimb3half[i] = sig*mulimb1[i];
			mulimb2[i] = sig*mulimb3half[i];
		}	

		//printf("%0.9f\t%0.9f\n", sig, mulimbhalf[i1]);

		for(j=1; j<nr; j++)	//I checked these limits and I think they're right
		{
			for(i=0; i<nz; i++)
			{
				IND(z0, i) = zt0[i]/r[j];	
			}
			mu = occultuniform(z0, rl/r[j], nz);
			//printf("%d\t%0.9f\n", j, mu[0]);
			sig1 = sqrt(cos(th[j-1]));
			sig2 = sqrt(cos(th[j]));
			dmumax = 0.;
			for(i=i1; i<=i2; i++)
			{
				f1 = r[j]*r[j]*mu[i]/(r[j] - r[j-1]);
				f2 = r[j]*r[j]*mu[i]/(r[j+1] - r[j]);
				mulimbhalf[i] += f1*pow(sig1,3) - f2*pow(sig2,3);
				mulimb1[i] += f1*pow(sig1,4) - f2*pow(sig2,4);
				mulimb3half[i] += f1*pow(sig1,5) - f2*pow(sig2,5);
				mulimb2[i] += f1*pow(sig1,6) - f2*pow(sig2,6);
				mulimb[i] = ((1.-c1-c2-c3-c4)*mulimb0[i]+c1*mulimbhalf[i]*dt+c2*mulimb1[i]*dt+c3*mulimb3half[i]*dt+c4*mulimb2[i]*dt)/omega;
				if(mulimb[i]+mulimbp[i] != 0.) 
				{
					dmumax = fmax(dmumax, fabs(mulimb[i]-mulimbp[i])/(mulimb[i]+mulimbp[i]));
					//printf("%0.10f\t%0.10f\n",mulimb[i], mulimbp[i]);
				}
			}
		}
	}
	for(i=i1; i<=i2; i++)
	{
		mulimbf[i] = mulimb0[i];
		mulimbf[nz+i] = mulimbhalf[i]*dt;
		mulimbf[2*nz+i] = mulimb1[i]*dt;
		mulimbf[3*nz+i] = mulimb3half[i]*dt;
		mulimb0[i] = mulimb[i];
	}	

	for(i=0; i<nz; i++)
	{
		IND(z0, i) = zt0[i];
		IND(muo1, i) = mulimb0[i];
	}

	free(mulimbf);
	free(zt0);
	free(mulimb0);
	free(mulimb);
	free(mulimbp);
	free(mulimbhalf);
	free(mulimb1);
	free(mulimb3half);
	free(mulimb2);
	free(mu);
	free(t);
	free(th);
	free(r);

	//printf("%0.6f\t%0.6f\n", IND(muo1,0), IND(muo1,1));
	return PyArray_Return(muo1);
}


double *occultuniform(PyArrayObject *z0, double p, int nz)
{
	int i;
	double *muo1, z, kap0, kap1, pi = acos(-1.);

	muo1 = (double *)malloc(nz*sizeof(double));

	if(fabs(p-0.5)<1.e-3) p = 0.5;
	for(i=0; i < nz; i++)
	{
		z = IND(z0, i);
		
		if(z >= 1.+p) muo1[i] = 1.;
		if(p >= 1. && z <= p - 1.) muo1[i] = 0.;
		else if(z <= 1.-p) muo1[i] = 1.-p*p;
		else	
		{
			kap1=acos(fmin((1.-p*p+z*z)/2./z,1.));
			kap0=acos(fmin((p*p+z*z-1.)/2./p/z,1.));
			muo1[i] = 1. - (p*p*kap0+kap1 -0.5*sqrt(fmax(4.*z*z-pow(1.+z*z-p*p, 2.), 0.)))/pi;
		}
	}
	return muo1;
}

static char _mandelagol_nonlinear_ld_doc[] = "LK 05/2015";

static PyMethodDef _mandelagol_nonlinear_ld_methods[] = {
  {"_mandelagol_nonlinear_ld", _mandelagol_nonlinear_ld, METH_VARARGS, _mandelagol_nonlinear_ld_doc},{NULL}};

void init_mandelagol_nonlinear_ld(void)
{
  Py_InitModule("_mandelagol_nonlinear_ld", _mandelagol_nonlinear_ld_methods);
  import_array();
}

