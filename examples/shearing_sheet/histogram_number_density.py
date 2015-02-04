import math
import numpy as np
import scipy as sp
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit

with open('position_x.txt','r') as particle_data_file:
        first_line = particle_data_file.readline()
        data = particle_data_file.readlines()

number_of_orbits			= 30
step					= 5
display_histogram_1_orbit		= 0
display_histogram_n_orbits		= 1
display_surface_density			= 1
display_viscosity_fit			= 0
x_data 					= []
sigma_data 				= []
number_density 				= []
number_density_for_fit			= []
n 					= 0
co 					= ['cyan','pink','lightblue','silver','violet','brown','black','green','darkblue','red']
temp_first_line                 	= first_line.split()
x_left_1                        	= float(temp_first_line[0])
x_left_2                        	= float(temp_first_line[1])
x_right_1                       	= float(temp_first_line[2])
x_right_2                       	= float(temp_first_line[3])
surface_density_lower           	= float(temp_first_line[4])
surface_density_upper           	= float(temp_first_line[5])
slope_left                      	= float(temp_first_line[6])
slope_right                     	= float(temp_first_line[7])
x_limit_lower                   	= x_left_1 + (0.2*x_left_1)
x_limit_upper                   	= x_right_2 + (0.2*x_right_2)

for line in data:
	temp = line.split()
        x_data.append(float(temp[0]))
        sigma_data.append(float(temp[1]))
        number_density.append(float(temp[2]))
        n += 1

def get_amp(a):
	#returns the average of the first a number of elements in number_density_for_fit
	temp_sum = 0
	for i in number_density_for_fit[:a]:
		temp_sum +=i
	return temp_sum/len(number_density_for_fit[:a])

def fit_func(x, dx, shift):
	#returns the value of the function that is been fit for a given x, dx and shift 
	#shift = x_left_1 + ((x_left_1 - x_left_2)/2.)
	#shift = [-370]
        return get_amp(5)/(1.+np.exp((x-shift)/dx))

if display_surface_density == 1:
	title_str = "Plot of Surface density versus x for " + str(n) + " particles"
	plt.plot(x_data, sigma_data, 'g^')
	plt.xlabel("x (m)")
	plt.ylabel("surface density (kg m^-2)")
	plt.title(title_str)
	#plt.xlim(x_limit_lower, x_limit_upper)
	plt.show()

if display_histogram_1_orbit == 1:
	x_data				= []

	title_str = "plot of Number density versus x for " + str(n) + " particles"
	plt.title(title_str)
	plt.xlabel("x (m)")
	plt.ylabel("Number density (m^-1)")
	with open('Data/position_0.0.txt', 'r') as particle_data_1_orbit:
		data_1_orbit = particle_data_1_orbit.readlines()
	
	for line in data_1_orbit:
		temp = line.split()
		x_data.append(float(temp[0]))

	label_str 			= "After adding particles"
	plt.hist(x_data, 50, histtype = 'step', color = 'red', label = label_str)
	plt.show()

if display_histogram_n_orbits == 1:
	title_str = "Plot of Number density versus x for " + str(n) + " particles with a fit"
	plt.title(title_str)
	plt.xlabel("x (m)")
	plt.ylabel("Number density (m^-1)")
	dx 				= []
	dshift				= []
	viscosity_data 			= []
	sum_dx 				= 0.
	dshift_sum			= 0.
	counter 			= 0
	viscosity_sum			= 0.
	OMEGA                           = 0.00013143527		# 1/s
	
	for i in range(0,number_of_orbits, step):
	#for i in range(10,100):
		j 			= i*1.0
		#j 			= i/10.0
		x_data 			= []
		x 			= []
		number_density_for_fit	= []
		with open('Data/position_%2.1f.txt' %j, 'r') as particle_data_file:
			data = particle_data_file.readlines()

 		for line in data:
			temp = line.split()
			x_data.append(float(temp[0]))
	
		if i==0:
			label_str 	= "After adding particles"
		elif i==1:
			label_str	= "After 1 orbit"
		else:
			label_str	= "After %d orbits" %i

		H = plt.hist(x_data, 50, histtype = 'step',color = co[i/5], label = label_str)
		for i in range(len(H[1])-1):
			if H[1][i] < 200.:
				x.append(H[1][i])
				number_density_for_fit.append(H[0][i])

		if display_viscosity_fit == 1:
			try:
				coeff, covar = curve_fit(fit_func, x, number_density_for_fit, p0 = (60., -370.))
			except:
				print("Could not fit data")
			dx.append(coeff[0])
			dshift.append(coeff[1])
			plt.plot(x, fit_func(x, coeff[0], coeff[1]), 'r--')
	
	if display_viscosity_fit == 1:
		for i in dx:
			sum_dx += i
			counter += 1
			t = counter*2.*math.pi/OMEGA
			viscosity_data.append((10000.*i*i)/t)
		
		for i in viscosity_data:
			viscosity_sum += i

		for i in dshift:
			dshift_sum += i
	
		shift_average = dshift_sum/len(dshift)
		viscosity = viscosity_sum/len(viscosity_data)	
		avg_dx = float(sum_dx/counter)
		#print "dx = ", avg_dx, " m"
		print "Shift =", shift_average*100., "cm"
		print "Viscosity =",viscosity, "cm^2/s"
	plt.legend(loc='upper center')
	plt.show()
