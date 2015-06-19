#include <math.h>

#define TWOPI 6.28318531		//FIXME more precise!
/*
	Here is a dummy example with I(r) \propto 1 - u1*r + u2*r^2
	The normalization constant is calculated by constraining the integrated intensity to equal 1:
		\int_r \int_theta {I(r)*r*dr*dtheta}/norm = 1
*/

double intensity(double r, double u1, double u2, double u3, double u4)
{
	if(r > 0.99995) r = 0.99995;
	double sqrtmu = pow(1.-r*r,0.25);
	double norm = (-u1/10.-u2/6.-3.*u3/14.-u4/4.+0.5)*TWOPI; 		//calculate norm by integrating I(r)r dr dtheta
	return (1. - u1*(1.-sqrtmu) - u2*(1. - pow(sqrtmu,2.)) - u3*(1.-pow(sqrtmu, 3.)) - u4*(1.-pow(sqrtmu,4.)))/norm; 	
}
