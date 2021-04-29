import matplotlib.pyplot as plt
import numpy as np

def definite_integral_show(start, end, N):

	x = np.arange(start, end, 0.01)
	y = np.arange(start, end, 0.01)
	z = np.multiply(x**2,y)
	maxF = max(z)

	x_rand = start + np.random.random(N)*(end-start)
	y_rand = start + np.random.random(N)*(end-start)
	z_rand = 0 +  np.random.random(N)*maxF

	ind_below = np.where(z_rand < np.multiply(x_rand**2,y_rand))
	ind_above = np.where(z_rand >= np.multiply(x_rand**2,y_rand))

	v = maxF*(end-start)*len(ind_below[0])/N
	return v



def main() :
	# PART A
	num_points = np.array(range(1,1000))
	I_values1 = [definite_integral_show(0,1,i) for i in num_points]
	I_values1 = np.array(I_values1)

	plt.plot(num_points,I_values1)
	plt.xlabel("Number of Points (N)")
	plt.ylabel("Integral Value")
	plt.title("a. I vs N")
	plt.show()
	# plt.savefig("2_a.png")
	print("a. Plotted")

	# PART B
	N = 20
	num_trails = 100
	I_values2 = [definite_integral_show(0,1,N) for i in range(num_trails)]
	I_values2 = np.array(I_values2)

	plt.plot(np.array(range(num_trails)),I_values2)
	plt.xlabel("Number of Trials")
	plt.ylabel("Integral Value for N = 20")
	plt.title('b. I vs trials (N=20)')
	plt.show()
	# plt.savefig("2_b.png")

	print("b. STANDARD DEVIATION FOR N = %d AND NUM_TRAILS = %d is %f"%(N,num_trails,np.std(I_values2)))

	# PART C
	N = 1000
	num_trails = 100
	I_values3 = [definite_integral_show(0,1,N) for i in range(num_trails)]
	I_values3 = np.array(I_values3)

	plt.plot(np.array(range(num_trails)),I_values3)
	plt.xlabel("Number of Trials")
	plt.ylabel("Integral Value for N = 1000")
	plt.title('c. I vs trials (N=1000)')
	plt.show()
	# plt.savefig("2_c.png")

	print("c. STANDARD DEVIATION FOR N = %d AND NUM_TRAILS = %d is %f"%(N,num_trails,np.std(I_values3)))

	# PART D
	N = np.array(range(1,1000))
	num_trails = 500
	std_vs_N = []
	print("d. In progress: ")
	for n in N:
		print('\r' + str( round((n*100)/len(N),2)) + "%", end = '')
		I_values = [definite_integral_show(0,1,n) for i in range(num_trails)]
		I_values = np.array(I_values)
		std_vs_N.append(np.std(I_values))
	print()

	std_vs_N = np.array(std_vs_N)
	plt.plot(N,std_vs_N)
	plt.plot(N,1/np.sqrt(N))
	plt.xlabel("Number of Points (N)")
	plt.legend(('Standard Deviation (500 trials)','1/sqrt(N)'))
	plt.title('d. Comparison of Monte Carlo with 1/sqrt(N)')
	plt.show()
	# plt.savefig("2_d.png")


main()