import math
import numpy
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

with open('position_moonlet.txt','r') as particle_data_file:
	data = particle_data_file.readlines()

with open('position_x.txt','r') as file_to_find_number_of_particles:
	number_data = file_to_find_number_of_particles.readlines()

with open('initial_profile.txt','r') as initial_profile_file:
	profile_data_lines = initial_profile_file.readline()

profile_data = profile_data_lines.split();
#print(profile_data)

display_initial_profile		= 0
display_output_migration_rate	= 1
display_output_position		= 0
display_output_velocity		= 0
display_output_delta_a		= 0
OMEGA				= 0.00013143527
G                               = 6.67428e-11
alpha                           = 2.46
K_0                             = 0.69676999
K_1                             = 0.45731925
n_0                             = OMEGA

n				= len(number_data)
x_data 				= []
y_data 				= []
z_data 				= []
v_x_data			= []
v_y_data			= []
v_z_data			= []
time_orbits_data 		= []
time_years_data			= []
delta_a_m_data			= []
delta_a_km_data			= []


r_h				= float(profile_data[0])
cut_off_boundary_left		= float(profile_data[1])
surface_density_upper		= float(profile_data[2])
x_left_1			= float(profile_data[3])
slope_left			= float(profile_data[4])
x_left_2			= float(profile_data[5])
surface_density_lower		= float(profile_data[6])
moonlet_x			= float(profile_data[7])
x_right_1			= float(profile_data[8])
slope_right			= float(profile_data[9])
x_right_2			= float(profile_data[10])
cut_off_boundary_right		= float(profile_data[11])
moonlet_mass			= float(profile_data[12])
a				= float(profile_data[13])
a_0				= float(profile_data[14])

left_limit			= cut_off_boundary_left + (0.1*cut_off_boundary_left) 
right_limit			= cut_off_boundary_right + (0.1*cut_off_boundary_right)

J_m                             = a_0*a_0*n_0
constant 			= ((64.*pow(G*moonlet_mass,2.)*a)/(243.*pow(OMEGA,3.)))*(pow(((2.*K_0)+(K_1)),2.))


for line in data:
	time_position_velocity 	= line.split()
	time_orbits		= float(time_position_velocity[0])				# orbits
	x 			= float(time_position_velocity[1])				# m
	y 			= float(time_position_velocity[2])				# m
	z 			= float(time_position_velocity[3])				# m
	v_x 			= float(time_position_velocity[4])				# m/s
	v_y			= float(time_position_velocity[5])				# m/s
	v_z			= float(time_position_velocity[6])				# m/s
	time_years		= ((2.*math.pi*time_orbits)/OMEGA)/(60.*60.*24.*365.)		# years
	delta_a_m		= ((2.*v_y)/OMEGA) + (4.*x)					# m
	delta_a_km		= delta_a_m/1000.						# km
	time_orbits_data.append(time_orbits)
	time_years_data.append(time_years)
	x_data.append(x)
	y_data.append(y)
	z_data.append(y)
	v_x_data.append(v_x)
	v_y_data.append(v_y)
	v_z_data.append(v_z)
	delta_a_m_data.append(delta_a_m)
	delta_a_km_data.append(delta_a_km)

def sigma_function(b_hat):
	
	x = b_hat*r_h	

	if x < cut_off_boundary_left:
		return 0.
	if x >= cut_off_boundary_left and x <= x_left_1:
		return surface_density_upper
	if x > x_left_1 and x <= x_left_2:
		return (slope_left*x)-(slope_left*x_left_1)+surface_density_upper
	if x > x_left_2 and x <= x_right_1:
		return surface_density_lower
	if x > x_right_1 and x <= x_right_2:
		return (slope_right*x)-(slope_right*x_right_1)+surface_density_lower
	if x > x_right_2 and x <= cut_off_boundary_right:
		return surface_density_upper
	if x > cut_off_boundary_right:
		return 0. 

def delta_j_function(b_hat):

	b = b_hat*r_h
	
	if b_hat <= -2.5:
		return constant*(1/pow(abs(b),5.))*(1+((alpha*abs(b))/a))
	
	elif b_hat > -2.5 and b_hat <= -1.8:
		return ((abs(b_hat))*r_h*J_m)/(2.*a_0)
	
	elif b_hat > -1.8 and b_hat <= 0.:
		return ((abs(b_hat))*r_h*J_m)/a_0;

	elif b_hat > 0 and b_hat <= 1.8:
		return (b_hat*r_h*J_m)/a_0;

	elif b_hat > 1.8 and b_hat <= 2.5:
		return (b_hat*r_h*J_m)/(2.*a_0);

	elif b_hat > 2.5:
		return constant*(1/pow(b,5.))*(1+((alpha*b)/a))
	
	else:
		print("Error")

if display_initial_profile == 1:
	x_rh_data 			= []
	sigma_data 			= []
	delta_j_data 			= []
	torque_data 			= []
	integrated_torque_data		= []
	for i in numpy.arange(left_limit, right_limit, 10.):
		x_rh_data.append(i/r_h)
		sigma 			= sigma_function(i/r_h)
		delta_j 		= ((delta_j_function(i/r_h))/(J_m*r_h))
		#print delta_j
		torque 			= sigma*delta_j
		sigma_data.append(sigma)
		delta_j_data.append(delta_j*pow(10.,3.))
		torque_data.append(torque)

	zero_position = len(x_rh_data)/2	

	for i in x_rh_data:
		temp_torque = 0.
		if i < 0.:
			for j in x_rh_data[:zero_position+1]:
				if j > i:
					temp_torque += delta_j_function(j)*sigma_function(j)
		if i > 0.:
			for j in x_rh_data[zero_position:]:
				if j < i:
					temp_torque += delta_j_function(j)*sigma_function(j)

		integrated_torque_data.append(temp_torque)

	f, axarr = plt.subplots(3, sharex=True)
	axarr[0].plot(x_rh_data, sigma_data, 'r')
	axarr[0].set_ylabel("Surface Density [kg/m^2]")
	axarr[0].set_title('Surface Density profile')

	axarr[1].plot(x_rh_data, delta_j_data)
	axarr[1].set_yscale('log')
	axarr[1].set_ylabel('Delta_J/(J_m*r_h) [km^-1]')
	axarr[1].set_title('Normalised Angular Momentum Transfer')

	axarr[2].plot(x_rh_data, integrated_torque_data, 'c')
	axarr[2].set_xlabel("Normalised Impact Parameter b/r_h")
	axarr[2].set_ylabel("Integrated Torque value [Nm]")
	axarr[2].set_title("Integrated Torque profile")

	"""axarr[2].plot(x_rh_data, torque_data, 'c')
	axarr[2].set_xlabel("Normalised Impact Parameter b/r_h")
	axarr[2].set_ylabel('Torque [Nm]')
	axarr[2].set_title('Torque profile')
	"""
	plt.show()	

if display_output_migration_rate == 1:
	coeff = numpy.polyfit(time_years_data[:90], delta_a_km_data[:90], 1)
	polynomial = numpy.poly1d(coeff)
	migration_rate = coeff[0]
	print "Numerical Migration rate =", migration_rate, "km/yr"
	plt.plot(time_years_data, delta_a_km_data, 'ro', label = "Delta a")
	plt.plot(time_years_data, polynomial(time_years_data), 'r')
	plt.xlabel('Time (years)')
	plt.ylabel('Delta a (km)')
	plt.legend()
	plt.show()

if display_output_position == 1:
	plt.plot(time_orbits_data, x_data, 'r', label = "x position")
	plt.plot(time_orbits_data, y_data, 'g', label = "y position")
	plt.title("Position of the moonlet versus time with "+str(n)+" particles in the simulation")
	plt.xlabel("Time (orbits)")
	plt.ylabel("Position (m)")
	plt.legend(loc='lower left')
	plt.show()

if display_output_velocity == 1:
	plt.plot(time_orbits_data, v_x_data, 'g', label="x velocity")
	plt.title("Velocity of the moonlet versus time with "+str(n)+" particles in the simulation")
	plt.xlabel("Time (orbits)")
	plt.ylabel("Velocity (m/s)")
	plt.legend(loc = 'lower left')
	plt.show()

if display_output_delta_a == 1:
	plt.plot(time_years_data, delta_a_km_data, 'r', label="delta_a plot")
	plt.title("Change in semimajor axis of the moonlet versus time with "+str(n)+" particles in the simulation")
	plt.xlabel("Time (years)")
	plt.ylabel("Delta a (km)")
	plt.legend(loc = 'lower left')
	plt.show()
