
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>
#include <unistd.h>

#define __USE_C99_MATH

FILE *fp;
FILE *f_p;

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
	// Setup constants
	double alpha					= 2.46;
	double K_0					= 0.69676999;
	double K_1					= 0.45731925;
	double OMEGA 					= 0.00013143527;					// 1/s  
	double n_0					= OMEGA;
	double G 					= 6.67428e-11;						// N / (1e-5 kg)^2 m^2
	double M_saturn					= 5.e26;						// kg
	double increment 				= 120.;
	double boxsize 					= 100.*increment;					// m
	double boxsize_x				= boxsize*2.;						// m
	double boxsize_y				= boxsize*2.;						// m
	double moonlet_radius				= 150./2.;				 		// m
	double temp_moonlet_radius			= 150./2.;						// m
	double a_0					= pow(((G*(M_saturn))/(OMEGA*OMEGA)),(1./3.));		// m
	double moonlet_mass				= (3.*M_saturn)*pow((moonlet_radius/a_0),3.);
	double r_h 					= a_0*pow((moonlet_mass/(3.*M_saturn)),(1./3.));

	// Initial conditions
	double x_left_1					= -3.*r_h; 
	double slope_left				= -4.8;
	double slope_right				= 1.44;
	double surface_density_upper			= 338./20.;
	double surface_density_lower	 		= 0.;
	double x_left_2	 				= (1./slope_left)*(surface_density_lower-surface_density_upper+(slope_left*x_left_1));
	double x_right_1				= -x_left_2;
	double x_right_2 				= (1./slope_right)*(surface_density_upper-surface_density_lower+(slope_right*x_right_1));
	double moonlet_x 				= 0.;							// m
        double a					= a_0 + moonlet_x;      //actual semi-major axis        // m
        double J_m					= a_0*a_0*n_0;


        double total_mass_1 = (surface_density_upper)*((boxsize_x/2.)+x_left_1)*boxsize_y;
        double total_mass_2 = ((x_left_2-x_left_1)*surface_density_lower)+(0.5*(x_left_2-x_left_1)*(surface_density_upper-surface_density_lower));
        double total_mass_3 = (surface_density_lower)*(x_right_1-x_left_2)*boxsize_y;
        double total_mass_4 = ((x_right_2-x_right_1)*surface_density_lower)+(0.5*(x_right_2-x_right_1)*(surface_density_upper-surface_density_lower));
        double total_mass_5 = (surface_density_upper)*((boxsize_x/2.)-x_right_2)*boxsize_y;
        double total_mass = total_mass_1+total_mass_2+total_mass_3+total_mass_4+total_mass_5;

	printf("\nM_ring/M_moonlet = %f\n", total_mass/moonlet_mass);



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

//	printf("Torque (-ve x) = %16.16e Nm\n", result_torque_left);
//	printf("Torque (+ve x) = %16.16e Nm\n", result_torque_right);
//      printf("Calculated Torque = %e Nm \n", result_torque);
	return result_torque;		
}


double migration_rate_function();
{
	double result_torque = torque_function();
	double result_migration_rate = (2.*result_torque)/(moonlet_mass*n_0*a_0);		// m/s
	result_migration_rate = (result_migration_rate/1000.)*(365.*24.*60.*60.);		// km/yr
	printf("Analytical Migration rate = %16.16f km/yr\n", result_migration_rate);
}

}

