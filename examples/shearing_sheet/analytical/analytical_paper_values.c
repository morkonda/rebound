/**
 * @file 	problem.c
 * @brief 	Example problem: shearing sheet.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This problem uses shearing sheet boundary
 * conditions. Particle properties resemble those found in 
 * Saturn's rings. 
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>
#include <unistd.h>

#define __USE_C99_MATH


double OMEGA;
int number_of_particles = 0;
FILE *fp;
FILE *f_p;
double a;
double M_saturn;
double a_0;
double J_m;
double r_h;
double toomre_wavelength;
double total_mass;
double total_mass_1;
double total_mass_2;
double total_mass_3;
double total_mass_4;
double total_mass_5;
double x_left_1;
double x_left_2;
double x_right_1;
double x_right_2;
double slope_left;
double slope_right;
double cut_off_boundary_left;
double cut_off_boundary_right;
double surface_density_lower;
double surface_density_upper;
double shift;
double moonlet_x;
double moonlet_radius;
double moonlet_mass;
extern double minimum_collision_velocity;

double result_left_1();
double result_left_2();
double result_left_3();
double result_right_1();
double result_right_2();
double result_right_3();
double torque_function();
double migration_rate_function();


void main(int argc, char* argv[])
{
	//remove("position_moonlet.txt");
	//FILE *fp_moonlet;
	//fp_moonlet = fopen("position_moonlet.txt", "a");
	//fclose(fp_moonlet);
        //for(int i=0; i<argc;i++)
	//{
	//	printf("%d\t%s\n",i,argv[i]);
	//}
	
	// Setup constants

	OMEGA 					= 0.00013143527;				// 1/s  
	double G 				= 6.67428e-11;					// N / (1e-5 kg)^2 m^2
	M_saturn				= 5.e26;					// kg
	//a_0 					= pow(((G*M_saturn)/(OMEGA*OMEGA)),(1./3.));	// m
	//a_0 					= 134912000. ;					// m


	double particle_density			= 400.;						// kg/m^3
	double moonlet_density			= 300.;						// kg/m^3
	double particle_radius_min		= 1;						// m
	double particle_radius_max		= 4;						// m
	double particle_radius_slope 		= -3;	
	double increment 			= 120.;
	double boxsize 				= 100.*increment;				// m
	double boxsize_x			= boxsize*2.;					// m
	double boxsize_y			= boxsize*2.;					// m

	cut_off_boundary_left			= -fabs(boxsize - (0.2*boxsize));
	cut_off_boundary_right			= fabs(boxsize - (0.2*boxsize));
	//double cut_off_boundary_left		= -1700.;
	//double cut_off_boundary_right		= 1700.;
	
	moonlet_radius				= 150./2.;		//scaled moonlet radius 				// m
	double temp_moonlet_radius		= 150./2.;
	//moonlet_radius			= 1000.;		//Bleriot's radius					// m
	a_0					= pow(((G*(M_saturn))/(OMEGA*OMEGA)),(1./3.));			// m
	moonlet_mass				= (3.*M_saturn)*pow((moonlet_radius/a_0),3.);
	r_h 					= a_0*pow((moonlet_mass/(3*M_saturn)),(1./3.));

	// Initial conditions

	bool result 				= false;
	double shift_implemented		= 0.;
	double mass 				= 0.;
	double temp_x 				= 0.;
	double temp_sigma 			= 0.;
	double number_density 			= 0.;
	double surface_area 			= 0.;

	//double x_left_1		 	= -32.4*increment;	//real x_left_1
	x_left_1				= -3.*r_h; 

	//double slope_left 			= -0.048;		//real slope						// kg m^-2 m^-1
	//double slope_right	 		= 0.0144;		//real slope						// kg m^-2 m^-1
	slope_left				= -4.8;		//scaled from the real slope
	slope_right				= 1.44;		//scaled from the real slope

	//double surface_density_upper		= 338.;			//real density upper
	//double surface_density_lower		= 300.;			//real density lower
	surface_density_upper			= 338./20.;	//scaled density upper
	surface_density_lower	 		= 0.;		//scaled density lower
	
	double alpha				= 2.46;
	double K_0				= 0.69676999;
	double K_1				= 0.45731925;
	double n_0				= OMEGA;
									// 1/s
	x_left_2	 			= (1./slope_left)*(surface_density_lower-surface_density_upper+(slope_left*x_left_1));
	x_right_1				= -x_left_2;
	x_right_2 				= (1./slope_right)*(surface_density_upper-surface_density_lower+(slope_right*x_right_1));


	moonlet_x 				= 0.;										// m


        double a                                = a_0 + moonlet_x;      //actual semi-major axis                                // m
        double J_m                              = a_0*a_0*n_0;


        total_mass_1 = (surface_density_upper)*((boxsize_x/2.)+x_left_1)*boxsize_y;
        total_mass_2 = ((x_left_2-x_left_1)*surface_density_lower)+(0.5*(x_left_2-x_left_1)*(surface_density_upper-surface_density_lower));
        total_mass_3 = (surface_density_lower)*(x_right_1-x_left_2)*boxsize_y;
        total_mass_4 = ((x_right_2-x_right_1)*surface_density_lower)+(0.5*(x_right_2-x_right_1)*(surface_density_upper-surface_density_lower));
        total_mass_5 = (surface_density_upper)*((boxsize_x/2.)-x_right_2)*boxsize_y;
        total_mass = total_mass_1+total_mass_2+total_mass_3+total_mass_4+total_mass_5;

printf("Total mass particles = %f\n", total_mass);

printf("M_ring/M_moonlet = %f\n", total_mass/moonlet_mass);



double result_left_1()
{
	double lower_limit = -boxsize_x/2.;
	double upper_limit = -2.5*r_h;
        double constant_neg_mdr = 3.*n_0*r_h*r_h*((64.*pow(G*moonlet_mass,2.)*a)/(243.*pow(OMEGA,3.)))*(pow(((2.*K_0)+(K_1)),2.));

	double integral_left_1_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double integral_lower_value = ((1./(3.*pow(integral_lower_limit,3.)))-(alpha/(2.*a*pow(integral_lower_limit,2.))));	
		double integral_upper_value = ((1./(3.*pow(integral_upper_limit,3.)))-(alpha/(2.*a*pow(integral_upper_limit,2.))));

		if(integral_lower_limit == lower_limit)
		 {
			integral_lower_value = 0.;
		 }

		double temp_result = (1./(r_h*r_h))*constant_neg_mdr*integral_sigma*(integral_upper_value-integral_lower_value);
		return temp_result;
	 }
	
		
	double integral_left_1_for_function_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		double integral_lower_value = (+((slope_left)/(2.*pow(integral_lower_limit,2.)))+((-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x))/(3.*pow(integral_lower_limit,3.)))-((alpha*slope_left)/(a*integral_lower_limit))-((alpha*(-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x)))/(2.*a*pow(integral_lower_limit,2.))));
		double integral_upper_value = (+((slope_left)/(2.*pow(integral_upper_limit,2.)))+((-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x))/(3.*pow(integral_upper_limit,3.)))-((alpha*slope_left)/(a*integral_upper_limit))-((alpha*(-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x)))/(2.*a*pow(integral_upper_limit,2.))));

		if(integral_lower_limit == lower_limit)
		 {
			integral_lower_value = 0.;
		 }		

		double temp_result = (1./(r_h*r_h))*constant_neg_mdr*(integral_upper_value-integral_lower_value);
		return temp_result;
	 }


//x_left_1 is to the right of (-2.5*r_h)
	if(upper_limit <= x_left_1)
	 {
		return integral_left_1_for_constant_sigma(lower_limit, upper_limit, surface_density_upper);
	 }
	
//x_left_1 is to the left of (-2.5*r_h) and x_left_2 is to the right of (-2.5*r_h)
	else if(upper_limit > x_left_1 && upper_limit < x_left_2) 
	 {
		double temp_result_1 = integral_left_1_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);		
		double temp_result_2 = integral_left_1_for_function_sigma(x_left_1, upper_limit);	
		return temp_result_1 + temp_result_2;
	 }
	
//x_left_1 and x_left_2 are to the left of (-2.5*r_h)
	else if(upper_limit > x_left_2)
	 {
		double temp_result_1 = integral_left_1_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);
		double temp_result_2 = integral_left_1_for_function_sigma(x_left_1, x_left_2);
		double temp_result_3 = integral_left_1_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);
		return temp_result_1 + temp_result_2 + temp_result_3;
	 }

//none of the above is true
	else
	 {
		printf("Missing in result_left_1");
		return 0.;
	 }
}

double result_left_2()
{
	double lower_limit = -2.5*r_h;
	double upper_limit = -1.8*r_h;
	double constant_neg_cdr = 3.*n_0*r_h*r_h*(r_h/(2.*a_0))*J_m;

	double integral_left_2_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double integral_upper_value = -(1./3.)*pow(integral_upper_limit,3.);
		double integral_lower_value = -(1./3.)*pow(integral_lower_limit,3.);

		double temp_result = (1./(r_h*r_h*r_h))*constant_neg_cdr*integral_sigma*(integral_upper_value - integral_lower_value);
		return temp_result;
	 }

	double integral_left_2_for_function_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		double integral_upper_value = -(((slope_left*pow(integral_upper_limit,4.))/4.)+(((-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x))*pow(integral_upper_limit,3.))/3.));
		double integral_lower_value = -(((slope_left*pow(integral_lower_limit,4.))/4.)+(((-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x))*pow(integral_lower_limit,3.))/3.));

		double temp_result = (1./(r_h*r_h*r_h))*constant_neg_cdr*(integral_upper_value - integral_lower_value);
		return temp_result;
	 }

//x_left_1 is to the right of (-1.8*r_h)
	if(upper_limit <= x_left_1)
	 {
		return integral_left_2_for_constant_sigma(lower_limit, upper_limit, surface_density_upper);
	 }

//x_left_1 is in between (-2.5*r_h) & (-1.8*r_h) and x_left_2 is to the right of (1.8*r_h)
	else if(lower_limit < x_left_1 && upper_limit > x_left_1 && upper_limit < x_left_2)
	 {
		double temp_result_1 = integral_left_2_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);
		double temp_result_2 = integral_left_2_for_function_sigma(x_left_1, upper_limit);
		return temp_result_1 + temp_result_2;
	 }

//x_left_1 and x_left_2 are in between (-2.5*r_h) and (-1.8*r_h)
	else if(lower_limit < x_left_1 && upper_limit > x_left_2)
	 {
		double temp_result_1 = integral_left_2_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);
		double temp_result_2 = integral_left_2_for_function_sigma(x_left_1, x_left_2);
		double temp_result_3 = integral_left_2_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);
		return temp_result_1 + temp_result_2 + temp_result_3;
	 } 
	
//x_left_1 is to the left of (-2.5*r_h) and x_left_2 is in between (-2.5*r_h) & (-1.8*r_h)
	else if(lower_limit > x_left_1 && lower_limit < x_left_2 && upper_limit > x_left_2)
	 {
		double temp_result_1 = integral_left_2_for_function_sigma(lower_limit, x_left_2);
		double temp_result_2 = integral_left_2_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);	
		return temp_result_1 + temp_result_2;
	 }

//both x_left_1 and x_left_2 are to the left of (-2.5*r_h)
	else if(lower_limit > x_left_2)
	 {
		return integral_left_2_for_constant_sigma(lower_limit, upper_limit, surface_density_lower);
 	 }

//x_left_1 is to the left of (-2.5*r_h) and x_left_2 is to the right of (-1.8*r_h)
	else if(lower_limit > x_left_1 && upper_limit < x_left_2)
	 {
		return integral_left_2_for_function_sigma(lower_limit, upper_limit);
	 }

//none of the above is true
	else
	 {
		printf("Missing in result_left_2");
		return 0.;
	 }	
}


double result_left_3()
{
	double lower_limit = -1.8*r_h;
	double upper_limit = 0.;
	double constant_neg_hr = 3.*n_0*r_h*r_h*(r_h/a_0)*J_m;

	double integral_left_3_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double integral_upper_value = -(1./3.)*pow(integral_upper_limit,3.);
		double integral_lower_value = -(1./3.)*pow(integral_lower_limit,3.);

		double temp_result = (1./(r_h*r_h*r_h))*constant_neg_hr*integral_sigma*(integral_upper_value - integral_lower_value);	
		return temp_result; 
	 }
	
	double integral_left_3_for_function_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		double integral_upper_value = -(((slope_left*pow(integral_upper_limit,4.))/4.)+(((-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x))*pow(integral_upper_limit,3.))/3.));
		double integral_lower_value = -(((slope_left*pow(integral_lower_limit,4.))/4.)+(((-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x))*pow(integral_lower_limit,3.))/3.));

		double temp_result = (1./(r_h*r_h*r_h))*constant_neg_hr*(integral_upper_value - integral_lower_value);	
		return temp_result; 
	 } 

//both x_left_1 and x_left_2 are to the right of (-1.8*r_h)	
	if(x_left_1 > lower_limit)
	 {
		double temp_result_1 = integral_left_3_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);		
		double temp_result_2 = integral_left_3_for_function_sigma(x_left_1, x_left_2);
		double temp_result_3 = integral_left_3_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);
		return temp_result_1 + temp_result_2 + temp_result_3;
	 }

//x_left_1 is to the left of (-1.8*r_h) and x_left_2 is to the right of (-1.8*r_h)
	else if(x_left_1 < lower_limit && x_left_2 > lower_limit)
	 {
		double temp_result_1 = integral_left_3_for_function_sigma(lower_limit, x_left_2);		
		double temp_result_2 = integral_left_3_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);
		return temp_result_1 + temp_result_2;
	 }

//both x_left_1 and x_left_2 are to the left of (-1.8*r_h)	
	else if(x_left_2 < lower_limit)
	 {
		double temp_result = integral_left_3_for_constant_sigma(lower_limit, upper_limit, surface_density_lower);
		return temp_result;
	 }

//none of the above is true
	else
	 {
		printf("Missing in result_left_3");
		return 0.;
	 }
}


//Clockwise: -ve torque
double result_right_1()
{
	double lower_limit = 0.;
	double upper_limit = 1.8*r_h;
	double constant_pos_hr = 3.*n_0*r_h*r_h*(r_h/a_0)*J_m;

	double integral_right_1_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double integral_upper_value = (1./3.)*pow(integral_upper_limit,3.);
		double integral_lower_value = (1./3.)*pow(integral_lower_limit,3.);

		double temp_result = (1./(r_h*r_h*r_h))*constant_pos_hr*integral_sigma*(integral_upper_value - integral_lower_value);
		return temp_result; 
	 }

	double integral_right_1_for_function_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		double integral_upper_value = (((slope_right*pow(integral_upper_limit,4.))/4.)+(((-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x))*pow(integral_upper_limit,3.))/3.));
		double integral_lower_value = (((slope_right*pow(integral_lower_limit,4.))/4.)+(((-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x))*pow(integral_lower_limit,3.))/3.));

		double temp_result = (1./(r_h*r_h*r_h))*constant_pos_hr*(integral_upper_value - integral_lower_value);
		return temp_result; 
	}

//x_right_1 is to the right of (1.8*r_h)
	if(x_right_1 > upper_limit)
	 {

		double temp_result = integral_right_1_for_constant_sigma(lower_limit, upper_limit, surface_density_lower);
		return temp_result;
	 }

//x_right_1 is to the left of (1.8*r_h) and x_right_2 is to the right of (1.8*r_h)
	else if(x_right_1 < upper_limit && x_right_2 > upper_limit)
	 {
		double temp_result_1 = integral_right_1_for_constant_sigma(lower_limit, x_right_1, surface_density_lower);
		double temp_result_2 = integral_right_1_for_function_sigma(x_right_1, upper_limit);
		return temp_result_1 + temp_result_2;
	 }

//both x_right_1 and x_rigth_2 are to the left of (1.8*r_h)
	else if(x_right_2 < upper_limit)
	 {
		double temp_result_1 = integral_right_1_for_constant_sigma(lower_limit, x_right_1, surface_density_lower);
		double temp_result_2 = integral_right_1_for_function_sigma(x_right_1, x_right_2);
		double temp_result_3 = integral_right_1_for_constant_sigma(x_right_2, upper_limit, surface_density_upper);
		return temp_result_1 + temp_result_2 + temp_result_3;
	 }

//none of the above is true
	else
	 {
		printf("Missing in result_right_1");
		return 0.;
	 }	
}


double result_right_2()
{
	double lower_limit = 1.8*r_h;
	double upper_limit = 2.5*r_h;
	double constant_pos_cdr = 3.*n_0*r_h*r_h*(r_h/(2.*a_0))*J_m;

	double integral_right_2_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double integral_upper_value = (1./3.)*pow(integral_upper_limit,3.);
		double integral_lower_value = (1./3.)*pow(integral_lower_limit,3.);

		double temp_result = (1./(r_h*r_h*r_h))*constant_pos_cdr*integral_sigma*(integral_upper_value - integral_lower_value);	
		return temp_result; 
	 }
	
	double integral_right_2_for_function_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		double integral_upper_value = (((slope_right*pow(integral_upper_limit,4.))/4.)+(((-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x))*pow(integral_upper_limit,3.))/3.));
		double integral_lower_value = (((slope_right*pow(integral_lower_limit,4.))/4.)+(((-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x))*pow(integral_lower_limit,3.))/3.));

		double temp_result = (1./(r_h*r_h*r_h))*constant_pos_cdr*(integral_upper_value - integral_lower_value);
		return temp_result;
	 }

//both x_right_1 and x_right_2 are to the right of (2.5*r_h)
	if(x_right_1 > upper_limit)
	 {
		return integral_right_2_for_constant_sigma(lower_limit, upper_limit, surface_density_lower);
	 }

//x_right_1 is in between (1.8*r_h) & (2.5*r_h) and x_right_2 is to the right of (2.5*r_h)
	else if(x_right_1 > lower_limit && x_right_1 < upper_limit && x_right_2 > upper_limit)
	 {
		double temp_result_1 = integral_right_2_for_constant_sigma(lower_limit, x_right_1, surface_density_lower);
		double temp_result_2 = integral_right_2_for_function_sigma(x_right_1, upper_limit);
		return temp_result_1 + temp_result_2;
	 }

//both x_right_1 and x_right_2 are in between (1.8*r_h) & (2.5*r_h)
	else if(x_right_1 > lower_limit && x_right_2 < upper_limit)
	 {
		double temp_result_1 = integral_right_2_for_constant_sigma(lower_limit, x_right_1, surface_density_lower);
		double temp_result_2 = integral_right_2_for_function_sigma(x_right_1, x_right_2);
		double temp_result_3 = integral_right_2_for_constant_sigma(x_right_2, upper_limit, surface_density_upper);
		return temp_result_1 + temp_result_2 + temp_result_3;
	 }

//x_right_1 is to the left of (1.8*r_h) and x_right_2 is in between (1.8*r_h) & (2.5*r_h)
	else if(x_right_1 < lower_limit && x_right_2 > lower_limit && x_right_2 < upper_limit)
	 {
		double temp_result_1 = integral_right_2_for_function_sigma(lower_limit, x_right_2);
		double temp_result_2 = integral_right_2_for_constant_sigma(x_right_2, upper_limit, surface_density_upper);
		return temp_result_1 + temp_result_2;
	 }
//both x_right_1 and x_right_2 are to the left of (1.8*r_h)
	else if(x_right_2 < lower_limit)
	 {
		return integral_right_2_for_constant_sigma(lower_limit, upper_limit, surface_density_upper);
	 }

//x_right_1 is to the left of (1.8*r_h) and x_right_2 is to the right of (2.5*r_h)
	else if(x_right_1 < lower_limit && x_right_2 > upper_limit)
	 {
		return integral_right_2_for_function_sigma(lower_limit, upper_limit);
	 }

//none of the above is true
	else
	 {
		printf("Missing in result_right_2");
		return 0.;
	 }
}


double result_right_3()
{
	double lower_limit = 2.5*r_h;
	double upper_limit = boxsize_x/2.;
	double constant_pos_mdr = 3.*n_0*r_h*r_h*((64.*pow(G*moonlet_mass,2.)*a)/(243.*pow(OMEGA,3.)))*(pow(((2.*K_0)+(K_1)),2.));

	double integral_right_3_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double integral_lower_value = (-(1./(3.*pow(integral_lower_limit,3.)))-(alpha/(2.*a*pow(integral_lower_limit,2.))));
		double integral_upper_value = (-(1./(3.*pow(integral_upper_limit,3.)))-(alpha/(2.*a*pow(integral_upper_limit,2.))));

		if(integral_upper_limit == upper_limit)
		 {
			integral_upper_value = 0.;
		 }

		double temp_result = (1./(r_h*r_h))*constant_pos_mdr*integral_sigma*(integral_upper_value-integral_lower_value);
	 	return temp_result;
	 }

	double integral_right_3_for_function_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		double integral_lower_value = (-((slope_right)/(2.*pow(integral_lower_limit,2.)))-((-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x))/(3.*pow(integral_lower_limit,3.)))-((alpha*slope_right)/(a*integral_lower_limit))-((alpha*(-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x)))/(2.*a*pow(integral_lower_limit,2.))));
		double integral_upper_value = (-((slope_right)/(2.*pow(integral_upper_limit,2.)))-((-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x))/(3.*pow(integral_upper_limit,3.)))-((alpha*slope_right)/(a*integral_upper_limit))-((alpha*(-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x)))/(2.*a*pow(integral_upper_limit,2.))));

		if(integral_upper_limit == upper_limit)
		 {
			integral_upper_value = 0.;
		 }

		double temp_result = (1./(r_h*r_h))*constant_pos_mdr*(integral_upper_value - integral_lower_value);
		return temp_result; 
	}

//both x_right_1 and x_right_2 are to the left of (2.5*r_h)
	if(x_right_2 <= lower_limit)
	 {
		double temp_1 = integral_right_3_for_constant_sigma(lower_limit, upper_limit, surface_density_upper);
	 	return temp_1;
	 }
	
//x_right_1 is to the left of (2.5*r_h) and x_right_2 is to the right of (2.5*r_h)
	else if(x_right_1 < lower_limit && x_right_2 > lower_limit)
	 {
		double temp_result_1 = integral_right_3_for_function_sigma(lower_limit, x_right_2);
		double temp_result_2 = integral_right_3_for_constant_sigma(x_right_2, upper_limit, surface_density_upper);
		return temp_result_1 + temp_result_2;
	 }

//both x_right_1 and x_right_2 are to the right of (2.5*r_h)
	else if(x_right_1 > lower_limit)
	 {
		double temp_result_1 = integral_right_3_for_constant_sigma(lower_limit, x_right_1, surface_density_lower);
		double temp_result_2 = integral_right_3_for_function_sigma(x_right_1, x_right_2);
		double temp_result_3 = integral_right_3_for_constant_sigma(x_right_2, upper_limit, surface_density_upper);
		return temp_result_1 + temp_result_2 + temp_result_3;
	 }

//none of the above is true
	else
	 {
		printf("Missing in result_right_3");
		return 0.;
	 }
}


double torque_function()
{
	double temp_result_left_1 = fabs(result_left_1());
	double temp_result_left_2 = -fabs(result_left_2());
	double temp_result_left_3 = -fabs(result_left_3());
	double result_torque_left = temp_result_left_1 + temp_result_left_2 + temp_result_left_3;

	double temp_result_right_1 = fabs(result_right_1());
	double temp_result_right_2 = fabs(result_right_2());
	double temp_result_right_3 = -fabs(result_right_3());
	double result_torque_right = temp_result_right_1 + temp_result_right_2 + temp_result_right_3;

	double result_torque = result_torque_left + result_torque_right;	

	printf("Torque_1 = %16.16e Nm \n", temp_result_left_1);
	printf("Torque_2 = %16.16e Nm \n", temp_result_left_2);
	printf("Torque_3 = %16.16e Nm \n", temp_result_left_3);
	printf("Torque_4 = %16.16e Nm \n", temp_result_right_1);
	printf("Torque_5 = %16.16e Nm \n", temp_result_right_2);
	printf("Torque_6 = %16.16e Nm \n", temp_result_right_3);
	printf("Torque (-ve x) = %16.16e Nm\n", result_torque_left);
	printf("Torque (+ve x) = %16.16e Nm\n", result_torque_right);
	return result_torque;		
}


double migration_rate_function();
{
	double result_torque = torque_function();
	double result_migration_rate = (2.*result_torque)/(moonlet_mass*n_0*a_0);		// m/s
	result_migration_rate = (result_migration_rate/1000.)*(365.*24.*60.*60.);		// km/yr
	printf("Calculated Torque = %e Nm \n", result_torque);
	printf("Analytical Migration rate = %16.16f km/yr\n", result_migration_rate);
}

}

