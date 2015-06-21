#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>

#if defined (_OPENMP)
#  include <omp.h>
#endif

double rc(double x, double y);
double rj(double x, double y, double z, double p);
double ellec(double k);
double ellk(double k);
double rf(double x, double y, double z);
double max(double a, double b);
double min(double a, double b);
static PyObject *_quadratic_ld(PyObject *self, PyObject *args);

static PyObject *_quadratic_ld(PyObject *self, PyObject *args)
{
/*	Input: *************************************
 
	 zs   	impact parameter in units of rs
	 u1   	linear    limb-darkening coefficient (gamma_1 in Mandel & Agol 2002)
	 u2   	quadratic limb-darkening coefficient (gamma_2)
	 p    	occulting star size in units of rs
	
	 Output: ***********************************
	
	 flux 	fraction of flux at each zs for a limb-darkened source
	
	 Limb darkening has the form:
	 I(r)=[1-u1*(1-sqrt(1-(r/rs)^2))-u2*(1-sqrt(1-(r/rs)^2))^2]/(1-u1/3-u2/6)/pi
*/
	int i, nz, nthreads;
	double u1, u2, p, *mu, *mu0,  *lambdad, *etad, \
		*lambdae, lam, pi, x1, x2, x3, z, omega, kap0 = 0.0, kap1 = 0.0, \
		q, Kk, Ek, Pk, n;
	PyArrayObject *zs, *flux;
	npy_intp dims[1], idx;

  	if(!PyArg_ParseTuple(args,"Odddi", &zs, &p, &u1, &u2, &nthreads)) return NULL;

	dims[0] = PyArray_DIMS(zs)[0];
	nz = dims[0];

	lambdad = (double *)malloc(nz*sizeof(double));
	lambdae = (double *)malloc(nz*sizeof(double));
	etad = (double *)malloc(nz*sizeof(double));
	mu = (double *)malloc(nz*sizeof(double));
	mu0 = (double *)malloc(nz*sizeof(double));

	if(fabs(p - 0.5) < 1.0e-3) p = 0.5;
	flux = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(zs));

	omega=1.0-u1/3.0-u2/6.0;
	pi=acos(-1.0);

	#if defined (_OPENMP)
	omp_set_num_threads(nthreads);
	#endif

	#if defined (_OPENMP)
	#pragma omp parallel for private(z, x1,x2,x3,n,q,Kk,Ek,Pk,kap0,kap1)
	#endif
	for(i = 0; i < nz; i++)
	{	
		idx = (npy_intp)i;
		z = *(double*)PyArray_GetPtr(zs, &idx);
		x1 = pow((p - z), 2.0);
		x2 = pow((p + z), 2.0);
		x3 = p*p - z*z;

		//source is unocculted:
		if(z >= 1.0 + p)					
		{
			//printf("zone 1\n");
			lambdad[i] = 0.0;
			etad[i] = 0.0;
			lambdae[i] = 0.0;
			*(double*)PyArray_GetPtr(flux, &idx) = 1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*lambdad[i]+u2*etad[i])/omega;
			mu0[i] = 1.0-lambdae[i];
			continue;
		}
		//source is completely occulted:
		if(p >= 1.0 && z <= p - 1.0)			
		{
			//printf("zone 2\n");
			lambdad[i] = 0.0;
			etad[i] = 0.5;		//error in Fortran code corrected here, following Jason Eastman's python code
			lambdae[i] = 1.0;
			 *(double*)PyArray_GetPtr(flux, &idx)= 1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*(lambdad[i] + 2.0/3.0)+u2*etad[i])/omega;
			mu0[i] = 1.0-lambdae[i];
			continue;
		}
		//source is partly occulted and occulting object crosses the limb:
		if(z >= fabs(1.0 - p) && z <= 1.0 + p)	
		{				
			//printf("zone 3\n");
			kap1 = acos(min((1.0-p*p+z*z)/2.0/z,1.0));
			kap0 = acos(min((p*p+z*z-1.0)/2.0/p/z,1.0));
			lambdae[i] = p*p*kap0+kap1;
			lambdae[i] = (lambdae[i] - 0.50*sqrt(max(4.0*z*z-pow((1.0+z*z-p*p), 2.0),0.0)))/pi;
		}
		//occulting object transits the source but doesn't completely cover it:
		if(z <= 1.0 - p)				
		{					
			//printf("zone 4\n");
			lambdae[i] = p*p;
		}
		//edge of the occulting star lies at the origin
		if(fabs(z-p) < 1.0e-3*(z+p))		
		{
			z = p;
		//	printf("zone 5\n");
			if(p == 0.5)	
			{
				//printf("zone 6\n");
				lambdad[i] = 1.0/3.0-4.0/pi/9.0;
				etad[i] = 3.0/32.0;
				*(double*)PyArray_GetPtr(flux, &idx) = 1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*lambdad[i]+u2*etad[i])/omega;
				mu0[i] = 1.0-lambdae[i];
				continue;
			}
			else if(z >= 0.5)
			{
				//printf("zone 5.1\n");
				lam = 0.5*pi;
				q = 0.5/p;
				Kk = ellk(q);
				Ek = ellec(q);
				lambdad[i] = 1.0/3.0+16.0*p/9.0/pi*(2.0*p*p-1.0)*Ek- \
				 	 	(32.0*pow(p, 4.0)-20.0*p*p+3.0)/9.0/pi/p*Kk;
				etad[i] = 1.0/2.0/pi*(kap1+p*p*(p*p+2.0*z*z)*kap0- \
				              	(1.0+5.0*p*p+z*z)/4.0*sqrt((1.0-x1)*(x2-1.0)));
				continue;
			}
			else if(z<0.5)	
			{
			//	printf("zone 5.2\n");
				lam = 0.50*pi;
				q = 2.0*p;
				Kk = ellk(q);
				Ek = ellec(q);
				lambdad[i] = 1.0/3.0+2.0/9.0/pi*(4.0*(2.0*p*p-1.0)*Ek + (1.0-4.0*p*p)*Kk);
				etad[i] = p*p/2.0*(p*p+2.0*z*z);
				*(double*)PyArray_GetPtr(flux, &idx) = 1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*lambdad[i]+u2*etad[i])/omega;
				mu0[i] = 1.0-lambdae[i];
				continue;
			}
			*(double*)PyArray_GetPtr(flux, &idx) = 1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*lambdad[i]+u2*etad[i])/omega;
			mu0[i] = 1.0-lambdae[i];
			continue;
		}
		//occulting star partly occults the source and crosses the limb:
		if((z > 0.5 + fabs(p -0.5) && z < 1.0 + p) || (p > 0.5 && z > fabs(1.0-p)*1.0001 \
			&& z < p))
		{
			//printf("zone 3.1\n");
			lam = 0.50*pi;
			q = sqrt((1.0-pow((p-z), 2.0))/4.0/z/p);
			Kk = ellk(q);
			Ek = ellec(q);
			n = 1.0/x1-1.0;
			Pk = Kk-n/3.0*rj(0.0,1.0-q*q,1.0,1.0+n);
			lambdad[i] = 1.0/9.0/pi/sqrt(p*z)*(((1.0-x2)*(2.0*x2+ \
			        x1-3.0)-3.0*x3*(x2-2.0))*Kk+4.0*p*z*(z*z+ \
			        7.0*p*p-4.0)*Ek-3.0*x3/x1*Pk);
			if(z < p) lambdad[i] += 2.0/3.0;
			etad[i] = 1.0/2.0/pi*(kap1+p*p*(p*p+2.0*z*z)*kap0- \
				(1.0+5.0*p*p+z*z)/4.0*sqrt((1.0-x1)*(x2-1.0)));
			*(double*)PyArray_GetPtr(flux, &idx) = 1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*lambdad[i]+u2*etad[i])/omega;
			mu0[i] = 1.0-lambdae[i];
			continue;
		}
		//occulting star transits the source:
		if(p <= 1.0  && z <= (1.0 - p)*1.0001)	
		{
			//printf("zone 4.1\n");
			lam = 0.50*pi;
			q = sqrt((x2-x1)/(1.0-x1));
			Kk = ellk(q);
			Ek = ellec(q);
			n = x2/x1-1.0;
			Pk = Kk-n/3.0*rj(0.0,1.0-q*q,1.0,1.0+n);
			lambdad[i] = 2.0/9.0/pi/sqrt(1.0-x1)*((1.0-5.0*z*z+p*p+ \
			         x3*x3)*Kk+(1.0-x1)*(z*z+7.0*p*p-4.0)*Ek-3.0*x3/x1*Pk);
			if(z < p) lambdad[i] += 2.0/3.0;
			if(fabs(p+z-1.0) <= 1.0e-4)
			{
				lambdad[i] = 2.0/3.0/pi*acos(1.0-2.0*p)-4.0/9.0/pi* \
				            sqrt(p*(1.0-p))*(3.0+2.0*p-8.0*p*p);
			}
			etad[i] = p*p/2.0*(p*p+2.0*z*z);
		}
		*(double*)PyArray_GetPtr(flux, &idx) = 1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*lambdad[i]+u2*etad[i])/omega;
		mu0[i] = 1.0-lambdae[i];
	}
	free(lambdae);
	free(lambdad);
	free(etad);
	free(mu);
	free(mu0);

	return PyArray_Return(flux);
}
double min(double a, double b)
{
	if(a < b) return a;
	else return b;
}

double max(double a, double b)
{
	if(a > b) return a;
	else return b;
}
double rc(double x, double y)
{
	double rc, ERRTOL,TINY,SQRTNY,BIG,TNBG,COMP1,COMP2,THIRD,C1,C2, C3,C4;
	ERRTOL=0.04; TINY=1.69e-38; SQRTNY=1.3e-19; BIG=3.0e37;
	TNBG=TINY*BIG; COMP1=2.236/SQRTNY; COMP2=TNBG*TNBG/25.0;
	THIRD=1.0/3.0; C1=0.3; C2=1.0/7.0; C3=0.375; C4=9.0/22.0;

	double alamb,ave,s,w,xt,yt;
	if(x < 0.0 || y == 0.0 || (x+fabs(y)) < TINY || (x+fabs(y)) > BIG || (y < -COMP1 && x > 0 && x < COMP2)){
		printf("Invalid argument(s) in rc\n");
		return 0;
	}
	if(y > 0.0)
	{
		xt=x;
		yt=y;
		w=1.0;
	}
	else
	{
		xt=x-y;
		yt=-y;
		w=sqrt(x)/sqrt(xt);
	}
	s = ERRTOL*10.0;
	while(fabs(s) > ERRTOL)
	{
		alamb = 2.0*sqrt(xt)*sqrt(yt)+yt;
		xt = 0.25*(xt+alamb);
		yt  =0.25*(yt+alamb);
		ave = THIRD*(xt+yt+yt);
		s = (yt-ave)/ave;
	}
	rc = w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave);
	return rc;
}

double rj(double x, double y, double z, double p)
{
	double rj, ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6,C7,C8, tempmax;
	
	ERRTOL=0.05; TINY=2.5e-13; BIG=9.0e11; C1=3.0/14.0;
	C2=1.0/3.0; C3=3.0/22.0; C4=3.0/26.0; C5=.750*C3;
     	C6=1.50*C4; C7=.50*C2; C8=C3+C3;
	
	double  a = 0.0,alamb,alpha,ave,b = 0.0,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee, \
     		fac,pt,rcx = 0.0,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt;
      
	if(x < 0.0 || y < 0.0 || z < 0.0 || (x+y) < TINY || (x+z) < TINY || (y+z) < TINY || fabs(p) < TINY \
		|| x > BIG || y > BIG || z > BIG || fabs(p) > BIG)
	{
		return 0;
	}
	sum=0.0;
	fac=1.0;
	if(p > 0.0)
	{
		xt=x;
		yt=y;
		zt=z;
		pt=p;
	}
	else
	{
		xt = min(x, y);
		xt = min(xt,z);
		zt = max(x, y);
		zt = max(zt, z);
		yt = x+y+z-xt-zt;
		a = 1.0/(yt-p);
		b = a*(zt-yt)*(yt-xt);
		pt = yt+b;
		rho = xt*zt/yt;
		tau = p*pt/yt;
		rcx = rc(rho,tau);
	}
	tempmax = ERRTOL*10.0;
	while(tempmax > ERRTOL)
	{
		sqrtx = sqrt(xt);
		sqrty = sqrt(yt);
		sqrtz = sqrt(zt);
		alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		alpha = pow((pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz), 2.0);
		beta = pt*(pt+alamb)*(pt + alamb);
		sum = sum+fac*rc(alpha,beta);
		fac = 0.25*fac;
		xt = 0.25*(xt+alamb);
		yt = 0.25*(yt+alamb);
		zt = 0.250*(zt+alamb);
		pt = 0.25*(pt+alamb);
		ave = 0.2*(xt+yt+zt+pt+pt);
		delx = (ave-xt)/ave;
		dely = (ave-yt)/ave;
		delz = (ave-zt)/ave;
		delp = (ave-pt)/ave;
		tempmax = max(fabs(delx), fabs(dely));
		tempmax = max(tempmax, fabs(delz));
		tempmax = max(tempmax, fabs(delp));
	}
	ea = delx*(dely+delz)+dely*delz;
	eb = delx*dely*delz;
	ec = delp*delp;
	ed = ea-3.0*ec;
	ee = eb+2.0*delp*(ea-ec);
	rj = 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4)) +\
		delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave));
	if(p < 0.0) rj=a*(b*rj+3.0*(rcx-rf(xt,yt,zt)));
	return rj;  
}
	
double ellec(double k)
{
	double m1,a1,a2,a3,a4,b1,b2,b3,b4,ee1,ee2,ellec;
	// Computes polynomial approximation for the complete elliptic
	// integral of the second kind (Hasting's approximation):
	m1=1.0-k*k;
	a1=0.44325141463;
	a2=0.06260601220;
	a3=0.04757383546;
	a4=0.01736506451;
	b1=0.24998368310;
	b2=0.09200180037;
	b3=0.04069697526;
	b4=0.00526449639;
	ee1=1.0+m1*(a1+m1*(a2+m1*(a3+m1*a4)));
	ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*log(1.0/m1);
	ellec=ee1+ee2;
	return ellec;
}

double ellk(double k)
{
	double a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,ellk, ek1,ek2,m1;
	// Computes polynomial approximation for the complete elliptic
	// integral of the first kind (Hasting's approximation):
	m1=1.0-k*k;
	a0=1.38629436112;
	a1=0.09666344259;
	a2=0.03590092383;
	a3=0.03742563713;
	a4=0.01451196212;
	b0=0.5;
	b1=0.12498593597;
	b2=0.06880248576;
	b3=0.03328355346;
	b4=0.00441787012;
	ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)));
	ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*log(m1);
	ellk=ek1-ek2;
	return ellk;
}

double rf(double x, double y, double z)
{
	double rf, ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4, tempmax;
	
	ERRTOL=0.08; TINY=1.5e-38; BIG=3.0e37; THIRD=1.0/3.0;
	C1=1.0/24.0; C2=0.1; C3=3.0/44.0; C4=1.0/14.0;
	
	double alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;

	if(min(x,y) < 0.0 || z < 0.0 || min(x+y, x+z) < TINY || y+z < TINY || max(x,y) > BIG || z > BIG)
	{
		printf("Invalid argument(s) in rf\n");
		return 0;
	}
	xt=x;
	yt=y;
	zt=z;
	tempmax = ERRTOL*10.0;
	while(tempmax > ERRTOL)
	{
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=THIRD*(xt+yt+zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
		tempmax = max(fabs(delx), fabs(dely));
		tempmax = max(tempmax, fabs(delz));
	}
	e2=delx*dely-delz*delz;
	e3=delx*dely*delz;
	rf=(1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave);
	return rf;
}

static char _quadratic_ld_doc[] = "This extension module returns a limb darkened light curve for a quadratic stellar intensity profile.";

static PyMethodDef _quadratic_ld_methods[] = {
  {"_quadratic_ld", _quadratic_ld,METH_VARARGS,_quadratic_ld_doc},{NULL}};

void init_quadratic_ld(void)
{
  Py_InitModule("_quadratic_ld", _quadratic_ld_methods);
  import_array();
}

