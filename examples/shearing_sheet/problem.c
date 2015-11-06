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
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "input.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"
#include "display.h"

#define __USE_C99_MATH


extern double OMEGA;
int number_of_particles = 0;
FILE *fp;
FILE *f_p;
double a;
double M_saturn;
double alpha;
double K_0;
double K_1;
double n_0;
double a_0;
double J_m;
double r_h;
double toomre_wavelength;
double x_left_1;
double x_left_2;
double x_right_1;
double x_right_2;
double slope_left;
double slope_right;
double increment;
double cut_off_boundary_left;
double cut_off_boundary_right;
double surface_density_lower;
double surface_density_upper;
double shift;
double moonlet_x;
double moonlet_radius;
double temp_moonlet_radius;
double moonlet_mass;
extern double minimum_collision_velocity;

extern double (*coefficient_of_restitution_for_velocity)(double); 
double coefficient_of_restitution_bridges(double v);
bool check(double x, double y);

double sigma_function(double x);
double toomre_wavelength_function(double x);
extern double opening_angle2;


bool check(double x, double y)
{
	if(x <= cut_off_boundary_left || x >= cut_off_boundary_right)
		return false;

	if(x < x_left_1)
		if(y <= surface_density_upper)
			return true;

	if(x >= x_left_1 && x < x_left_2)
		if(y <= ((slope_left*x)+(surface_density_upper)-(slope_left*x_left_1)))
			return true;

	if(x >= x_left_2 && x < x_right_1)
		if(y <= surface_density_lower)
			return true;

	if(x >= x_right_1 && x < x_right_2)
		if(y <= ((slope_right*x)+(surface_density_lower)-(slope_right*x_right_1)))
			return true;

	if(x >= x_right_2)
		if(y <= surface_density_upper)
			return true;

	return false;
}

double sigma_function(double x)
{
	if(x < x_left_1)
		return surface_density_upper;

	if(x >= x_left_1 && x < x_left_2)
		return ((slope_left*x)-(slope_left*x_left_1)+surface_density_upper);

	if(x >= x_left_2 && x < x_right_1)
		return surface_density_lower;

	if(x >= x_right_1 && x < x_right_2)
		return ((slope_right*x)-(slope_right*x_right_1)+surface_density_lower);

	if(x >= x_right_2)
		return surface_density_upper;
}

void problem_init(int argc, char* argv[])
{
	remove("position_moonlet.txt");
        for(int i=0; i<argc;i++)
	{
		printf("%d\t%s\n",i,argv[i]);
	}

	// Setup constants
	#ifdef GRAVITY_TREE
	opening_angle2		= .4;
	#endif // GRAVITY_TREE
	OMEGA 					= 0.00013143527;				// 1/s
	//Simulation runs for 30 orbits
	tmax                           		= 30.*2.*M_PI/OMEGA;  
	G 					= 6.67428e-11;					// N / (1e-5 kg)^2 m^2
	M_saturn				= 5.e26;					// kg
	alpha					= 2.46;
	K_0					= 0.69676999;
	K_1					= 0.45731925;
	n_0					= OMEGA;	
	softening 				= 0.1;						// m
	dt 					= 1e-3*2.*M_PI/OMEGA;				// s
	#ifdef OPENGL
	display_rotate_z			= 20;						// Rotate the box by 20 around the z axis, then 
	display_rotate_x			= 60;						// rotate the box by 60 degrees around the x axis	
	#ifdef LIBPNG
	system("mkdir png");
	#endif // LIBPNG
	#endif // OPENGL
	root_nx = 2; root_ny = 2; root_nz = 1;
	//nghostx = 2; nghosty = 2; nghostz = 0; 							// Use two ghost rings
	nghostx = 0; nghosty = 0; nghostz = 0;

	double particle_density			= 400.;										// kg/m^3
	double moonlet_density			= 300.;										// kg/m^3
	double particle_radius_min		= 1;										// m
	double particle_radius_max		= 4;										// m
	double particle_radius_slope 		= -3;	

	increment 				= 15.;
	boxsize 				= 100.*increment;								// m
	cut_off_boundary_left			= -fabs(boxsize - (0.2*boxsize));						// m
	cut_off_boundary_right			= fabs(boxsize - (0.2*boxsize));						// m

	//printf("Cut_off_boundary: Left = %f\tRight = %f\n", cut_off_boundary_left, cut_off_boundary_right);	
	//printf("Box: Left = %f\tRight = %f\n", -boxsize, boxsize);

	init_box();
	
	moonlet_radius				= 150./2.;	 								// m
	temp_moonlet_radius			= 150./2.;									// m
	a_0					= pow(((G*(M_saturn))/(OMEGA*OMEGA)),(1./3.));					// m
	moonlet_mass                         	= moonlet_density*(4./3.)*M_PI*pow(temp_moonlet_radius, 3.);                    // kg
	r_h 					= a_0*pow((moonlet_mass/(3.*M_saturn)),(1./3.));				// m

	// Use Bridges et al coefficient of restitution.
	coefficient_of_restitution_for_velocity = coefficient_of_restitution_bridges;
	minimum_collision_velocity = particle_radius_min*OMEGA*0.001;  				// small fraction of the shear

	bool result 				= false;
	double shift_implemented		= 0.;
	double mass 				= 0.;
	double temp_x 				= 0.;
	double temp_sigma 			= 0.;
	double number_density 			= 0.;
	double surface_area 			= 0.;

	// Initial Conditions

	x_left_1				= -3.*r_h; 									// m
	x_left_2				= -2.*r_h;									// m
	surface_density_upper			= 338./10.;									// kg m^-2
	surface_density_lower			= 0.;										// kg m^-2
	x_right_1				= 1.75*r_h;									// m
	x_right_2				= 3.*r_h;									// m	

	slope_left				= (surface_density_lower-surface_density_upper)/(x_left_2-x_left_1);
	slope_right				= (surface_density_upper-surface_density_lower)/(x_right_2-x_right_1);

	moonlet_x				= 0.;										// m

/*
	slope_left				= -4.8;		//scaled from the real slope					// kg m^-2 m^-1
	slope_right				= 1.44;		//scaled from the real slope					// kg m^-2 m^-1
	surface_density_upper			= 338./20.;	//scaled density upper						// kg m^-2
	surface_density_lower	 		= 0.;		//scaled density lower						// kg m^-2
	x_left_2	 			= (1./slope_left)*(surface_density_lower-surface_density_upper+(slope_left*x_left_1));
	x_right_1				= -x_left_2;
	x_right_2 				= (1./slope_right)*(surface_density_upper-surface_density_lower+(slope_right*x_right_1));
*/

	printf("Left Coordinates = (%f, %f) m\t surface_density = (%f, %f) kg m^-2\t Right coordinates = (%f, %f) m\n", x_left_1, x_left_2, surface_density_lower, surface_density_upper, x_right_1, x_right_2);

	printf("slopes = (%f, %f)\n", slope_left, slope_right);

	fp                                      = fopen("position_x.txt", "w");
        fprintf(fp, "%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f\n", x_left_1, x_left_2, x_right_1, x_right_2, surface_density_lower, surface_density_upper, slope_left, slope_right, r_h);


        double total_mass_1 = (surface_density_upper)*((boxsize_x/2.)+x_left_1)*boxsize_y;
        double total_mass_2 = ((x_left_2-x_left_1)*surface_density_lower)+(0.5*(x_left_2-x_left_1)*(surface_density_upper-surface_density_lower));
        double total_mass_3 = (surface_density_lower)*(x_right_1-x_left_2)*boxsize_y;
        double total_mass_4 = ((x_right_2-x_right_1)*surface_density_lower)+(0.5*(x_right_2-x_right_1)*(surface_density_upper-surface_density_lower));
        double total_mass_5 = (surface_density_upper)*((boxsize_x/2.)-x_right_2)*boxsize_y;
        double total_mass = total_mass_1+total_mass_2+total_mass_3+total_mass_4+total_mass_5;

	printf("M_ring/M_moonlet = %f\n", total_mass/moonlet_mass);

	while(mass<total_mass)
	{
        	temp_x = tools_uniform(-boxsize_x/2.,boxsize_x/2.);
        	temp_sigma = tools_uniform(0, (surface_density_upper+50.));
        	result = check(temp_x, temp_sigma);
        	if(result)
         	{
                	struct particle pt;
                	pt.x                            = temp_x;
                	pt.y                            = tools_uniform(-boxsize_y/2.,boxsize_y/2.);
                	pt.z                            = tools_normal(1.);
                	pt.vx                           = 0;
                	pt.vy                           = -1.5*pt.x*OMEGA;
                	pt.vz                           = 0;
                	pt.ax                           = 0;
                	pt.ay                           = 0;
                	pt.az                           = 0;
                	double radius                   = tools_powerlaw(particle_radius_min,particle_radius_max,particle_radius_slope);
                	//#ifndef COLLISIONS_NONE
                	pt.r                            = radius;
                	//#endif
                	double particle_mass            = particle_density*4./3.*M_PI*radius*radius*radius;
                	pt.m                            = particle_mass;
                	particles_add(pt);
                	mass                           += particle_mass;
                	number_of_particles++;
                	number_density                  = total_mass/(particle_mass*surface_area);
                	fprintf(fp, "%f \t %f \t %f \n", pt.x, temp_sigma, number_density);
         	}
	}

        double shift_left = -800;
        double shift_right = 800;
        double force_left;
        double force_right;
        int particles_left = 0;
        int particles_right = 0;

        for(int b=0; b<20; b++)
         {
                double force = 0;
                force_left = 0;
                force_right = 0;
                particles_left = 0;
                particles_right = 0;
                shift = (shift_left+shift_right)/2.;

                for(int i=0;i<N;i++)
                 {	
			for(int k=-nghosty; k <= nghosty; k++)
		 	 {	
                        	struct particle p = particles[i];
                        	double dx = p.x - shift;
                        	double dy = p.y + ((double)k)*boxsize_y;
                        	force = (G*p.m*moonlet_mass)/((dx*dx)+(dy*dy));

                        	if(dx < moonlet_x)
                         	 {
                                	force_left += force;
                                	particles_left += 1;
                         	 }
                        	else
                         	 {
                                	force_right += force;
                                	particles_right += 1;
                         	 }
		 	 }
                 }

                if(fabs(force_left) > fabs(force_right))
                        shift_left = shift;
                else
                        shift_right = shift;
        }

	for(int i = 0; i<N; i++)
	 {
		particles[i].x -= shift;
		particles[i].vy = -1.5*particles[i].x*OMEGA;

	 } 

	boundaries_check();
	tree_update();

        //printf("Difference in Particles = %d - %d = %d\n", particles_left, particles_right, particles_left - particles_right);
        //printf("Shift = %f\n", shift);

        struct particle pt;
        pt.x                                    = moonlet_x;
        pt.y                                    = 0.;
        pt.z                                    = 0.;
        pt.vx                                   = 0.;
        pt.vy                                   = -1.5*pt.x*OMEGA;
        //pt.vy                                 = 0;
        pt.vz                                   = 0;
        pt.ax                                   = 0;
        pt.ay                                   = 0;
        pt.az                                   = 0;
        pt.r                                    = moonlet_radius;
        pt.m                                    = moonlet_mass;
        particles_add(pt);
        number_of_particles++;
        fclose(fp);


	a				        = a_0 + moonlet_x;	//actual semi-major axis				// m
	J_m 					= a_0*a_0*n_0;
	//r_h                                   = a_0*(pow((moonlet_mass/(3.*M_saturn)),(1./3.)));				// m
	surface_area 				= (boxsize_x/2.)*(boxsize_y/2);

	x_left_1 -= shift;
	x_left_2 -= shift;
	x_right_1 -= shift;
	x_right_2 -= shift;
	cut_off_boundary_left -= shift;
	cut_off_boundary_right -= shift;

        //fprintf(fp, "%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f\n", x_left_1, x_left_2, x_right_1, x_right_2, surface_density_lower, surface_density_upper, slope_left, slope_right, r_h);


	double sigma_left(double x)
	{
		return (slope_left*x)-(slope_left*x_left_1)+surface_density_upper;
	}

	double sigma_right(double x)
	{
		return (slope_right*x)-(slope_right*x_right_1)+surface_density_lower;
	}

	f_p = fopen("initial_profile.txt","w");
	fprintf(f_p,"%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f %f \t %f \t %f \n", r_h, cut_off_boundary_left, surface_density_upper, x_left_1, slope_left, x_left_2, surface_density_lower, moonlet_x, x_right_1, slope_right, x_right_2, cut_off_boundary_right, moonlet_mass, a, a_0);
	fclose(f_p);
	
	double toomre_wavelength_function(double x)
	{
        	return (2.*M_PI*M_PI*(sigma_function(x))/OMEGA/OMEGA*G);
	}
}                

                 
double coefficient_of_restitution_bridges(double v)
{
	// assumes v in units of [m/s]
	double eps = 0.32*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
}


void problem_inloop()
{
}


void problem_output()
{
	if (output_check(1e-3*2.*M_PI/OMEGA))
	{
		output_timing();
		//output_append_velocity_dispersion("veldisp.txt");
	}

	//output the particles profile every orbit
	if (output_check((2.*M_PI/OMEGA)))
	{
		char fp_position_data[1024];
		sprintf(fp_position_data,"Data/position_%2.1f.txt",t/(2.*M_PI/OMEGA));
		output_ascii(fp_position_data);
	}

	double temp_x;
	double temp_sigma;
	bool temp_result;	

	for(int i=0;i<number_of_particles;i++)
	 {
		if(particles[i].x > x_left_2 && particles[i].x < 0 && particles[i].r != moonlet_radius)
		 {
			temp_result = false;
			while(temp_result == false)
			 {
                        	temp_x 			= tools_uniform(-boxsize_x/2.,x_left_2);
                        	temp_sigma 		= tools_uniform(0, (surface_density_upper+50.));
                        	temp_result 		= check(temp_x, temp_sigma);

			 }
			particles[i].x 			= temp_x;
			particles[i].vy                 = -1.5*particles[i].x*OMEGA;
		 }

		else if(particles[i].x > 0 && particles[i].x < x_right_1 && particles[i].r != moonlet_radius)
		 {
			temp_result = false;
			while(temp_result == false)
			 {
				temp_x			= tools_uniform(x_right_1, boxsize_x/2.);
				temp_sigma		= tools_uniform(0, (surface_density_upper+50.));
				temp_result		= check(temp_x, temp_sigma);
			 }
			particles[i].x			= temp_x;
			particles[i].vy			= -1.5*particles[i].x*OMEGA;
		 }
	 }

	//output the moonlet profile 10 times per orbit 
	if (output_check((2.*M_PI/OMEGA)/10.))
	//if(0)
	{
		FILE* fp_moonlet_data = fopen("position_moonlet.txt","a+");
		for(int i=0;i<number_of_particles;i++)
		 {
			if (particles[i].r == moonlet_radius)
		         { 
		     		struct particle ml = particles[i];
				fprintf(fp_moonlet_data,"%f\t%e\t%e\t%e\t%e\t%e\t%e\n", (t/(2.*M_PI/OMEGA)), ml.x, ml.y, ml.z, ml.vx, ml.vy, ml.vz);
				fclose(fp_moonlet_data);
			 }
		 }
	}
}


double initial_density_profile(int side)
{
	double temp_x;
	double temp_sigma;
	bool temp_result = false;
	if(side == -1)
	 {
        	while(temp_result == false)
                {
                	temp_x                  = tools_uniform(-boxsize_x/2.,x_left_2);
			//printf("\ntemp_x_left = %f\t", temp_x);
                        temp_sigma              = tools_uniform(0, (surface_density_upper+50.));
                        temp_result             = check(temp_x, temp_sigma);
				
                }
		//printf("\ntemp_x_left = %f",temp_x);
		return temp_x;
	 }

	else if(side == 1)
	 {
        	while(temp_result == false)
                  {      
                	temp_x                  = tools_uniform(x_right_1, boxsize_x/2.);
                        temp_sigma              = tools_uniform(0, (surface_density_upper+50.));
                      	temp_result             = check(temp_x, temp_sigma);
                  }
		return temp_x;
	 }
}


void problem_finish()
{
}
