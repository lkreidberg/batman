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

#include <math.h>

/*
	- The intensity function returns the stellar intensity at a radius r; where 0 <= r <= 1
	- The function arguments are the normalized radius (r), and limb darkening coefficients c1, ..., un
	- see below for an example intensity profile for I(r) \propto 1 - c1*(1-sqrt(1-r^2)) - c2*ln((sqrt(1-r^2)+c)/(1+c))
	- The normalization constant is calculated by constraining the integrated intensity to equal 1:
		\int_r \int_theta {I(r)*r*dr*dtheta}/norm = 1
*/

double intensity(double r, double c1, double c2, double c3, double c4, double c5, double c6) 
{
	if(r > 0.99995) r = 0.99995;
	double mu = sqrt(1. - r*r);
	double norm = 2.*M_PI*(-c1/6. - c2*c3/2. + c2/4. + 0.5 + c2*c3*c3*log(1. + 1./c3)/2.); 
	return (1. - c1*(1. - mu) - c2*log((mu + c3)/(1. + c3)))/norm;
}
