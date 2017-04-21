#if defined (_OPENACC) && defined(__PGI)
#  include <accelmath.h>
#else
#  include <math.h>
#endif

#if defined (_OPENMP) && !defined(_OPENACC)
#  include <omp.h>
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

inline double intensity(double x, double* args); //must be defined in the C file including this header

inline double area(double d, double x, double R)
{
	/*
	Returns area of overlapping circles with radii x and R; separated by a distance d
	*/
	double arg1 = (d*d + x*x - R*R)/(2.*d*x);
	double arg2 = (d*d + R*R - x*x)/(2.*d*R);
	double arg3 = MAX((-d + x + R)*(d + x - R)*(d - x + R)*(d + x + R), 0.);

	if(x <= R - d) return M_PI*x*x;							//planet completely overlaps stellar circle
	else if(x >= R + d) return M_PI*R*R;						//stellar circle completely overlaps planet
	else return x*x*acos(arg1) + R*R*acos(arg2) - 0.5*sqrt(arg3);			//partial overlap
}

void calc_limb_darkening(double* f_array, double* d_array, int N, double rprs, double fac, int nthreads, double* intensity_args)
{
	/*
		NOTE:  the safest way to access numpy arrays is to use the PyArray_GETITEM and PyArray_SETITEM functions.
		Here we use a trick for faster access and more convenient access, where we set a pointer to the
		beginning of the array with the PyArray_DATA (e.g., f_array) and access elements with e.g., f_array[i].
		Success of this operation depends on the numpy array storing data in blocks equal in size to a C double.
		If you run into trouble along these lines, I recommend changing the array access to something like:
			d = PyFloat_AsDouble(PyArray_GETITEM(ds, PyArray_GetPtr(ds, &i)));
		where ds is a numpy array object.


		Laura Kreidberg 07/2015
	*/

	#if defined (_OPENMP) && !defined(_OPENACC)
	omp_set_num_threads(nthreads);	//specifies number of threads (if OpenMP is supported)
	#endif

	#if defined (_OPENACC)
	#pragma acc parallel loop copyout(f_array[:N]) present(intensity_args)
	#elif defined (_OPENMP)
	#pragma omp parallel for
	#endif
	for(int i = 0; i < N; i++)
	{
		double d = d_array[i];
		double x_in = MAX(d - rprs, 0.);					//lower bound for integration
		double x_out = MIN(d + rprs, 1.0);					//upper bound for integration

		if(x_in >= 1.) f_array[i] = 1.0;				//flux = 1. if the planet is not transiting
		else
		{
			double delta = 0.;						//variable to store the integrated intensity, \int I dA
			double x = x_in;						//starting radius for integration
			double dx = fac*acos(x); 					//initial step size

			x += dx;						//first step

			double A_i = 0.;						//initial area

			while(x < x_out)
			{
				double A_f = area(d, x, rprs);				//calculates area of overlapping circles
				double I = intensity(x - dx/2., intensity_args); 	//intensity at the midpoint
				delta += (A_f - A_i)*I;				//increase in transit depth for this integration step
				dx = fac*acos(x);  				//updating step size
				x = x + dx;					//stepping to next element
				A_i = A_f;					//storing area
			}
			dx = x_out - x + dx;  					//calculating change in radius for last step  FIXME
			x = x_out;						//final radius for integration
			double A_f = area(d, x, rprs);					//area for last integration step
			double I = intensity(x - dx/2., intensity_args); 		//intensity at the midpoint
			delta += (A_f - A_i)*I;					//increase in transit depth for this integration step

			f_array[i] = 1.0 - delta;	//flux equals 1 - \int I dA
		}
	}
}
