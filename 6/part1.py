import random
import numpy as np
import os,sys
import matplotlib.pyplot as plt
import json

def definite_integral_show(start, end, N):	
	x = np.arange(start, end, 0.01)
	y = 3*x*x
	x_rand = start + np.random.random(N)*(end-start)
	maxF = max(y)
	y_rand = 0 +  np.random.random(N)*maxF

	ind_below = np.where(y_rand < (3*x_rand*x_rand))
	ind_above = np.where(y_rand >= (3*x_rand*x_rand))
	v = maxF*(end-start)*len(ind_below[0])/N
	
	return v



def main():
	# PART A
	I_values = []
	num_points = np.array(range(1,1000))
	I_values1 = [definite_integral_show(0,1,i) for i in num_points]
	I_values1 = np.array(I_values1)

	plt.plot(num_points,I_values1)
	plt.xlabel("Number of Points")
	plt.ylabel("Integral Value")
	plt.title("a. I vs Number of Points")
	plt.show()
	print("a. Plotted")
	# plt.savefig("1_a.png")
	
	# PART B

	N = 20
	num_of_trials = 100
	I_values2 = [definite_integral_show(0,1,N) for i in range(num_of_trials)]
	I_values2 = np.array(I_values2)

	plt.plot(np.array(range(num_of_trials)),I_values2)
	plt.title('b. I vs trials(N=20)')
	plt.xlabel("Number of Trials")
	plt.ylabel("Integral Value for N = 20")
	plt.show()
	# plt.savefig("1_b.png")
	# plt.figure()

	print("b. Standard Deviation = %d  Num_trials = %d is %f"%(N,num_of_trials,np.std(I_values2)))
	# PART C
	N = 1000
	num_of_trials = N
	I_values3 = [definite_integral_show(0,1,N) for i in range(num_of_trials)]
	I_values3 = np.array(I_values3)

	plt.plot(np.array(range(num_of_trials)),I_values3)
	plt.title('c. I vs trials (N=1000)')
	plt.xlabel("Number of Trials")
	plt.ylabel("Integral Value for N = 1000")
	plt.show()
	# plt.savefig("1_c.png")
	print("c. Standard Deviation = %d Num_trials = %d is %f"%(N,num_of_trials,np.std(I_values3)))
	# plt.figure()
	# PART D
	std_vs_N = []
	N = np.array(range(1,1000))
	num_of_trials = 100
	print("d. In progress:")
	for n in N:
		print('\r' + str( round((n*100)/len(N),2)) + "%", end = '')
		I_values = [definite_integral_show(0,1,n) for i in range(num_of_trials)]
		I_values = np.array(I_values)
		std_vs_N.append(np.std(I_values))
	std_vs_N = np.array(std_vs_N)
	plt.plot(N,std_vs_N)
	plt.plot(N,1/np.sqrt(N))
	plt.title('d. Comparison of Monte Carlo with 1/sqrt(N)')
	plt.legend(('Standard Deviation of 100 trials','1/sqrt(N)'))
	plt.xlabel("Number of Points")
	# plt.savefig("1_d.png")
	plt.show()
	# plt.figure()
	print()
main()