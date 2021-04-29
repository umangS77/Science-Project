import random 
import matplotlib.pyplot as plt

xr = []
yr = []
for INTERVAL in range(1,200,1):
	circle_points= 0
	pi=0
	for i in range(INTERVAL**2): 

		rand_x= random.uniform(-1, 1) 
		rand_y= random.uniform(-1, 1) 

		origin_dist= rand_x**2 + rand_y**2

		if origin_dist<= 1: 
			circle_points+= 1

	pi = 4* circle_points/ (INTERVAL**2) 

	xr.append(INTERVAL)
	yr.append(pi)


plt.title('Estimated value of Pi using N inputs vs N')
plt.xlabel('N')
plt.ylabel('Pi')
plt.plot(xr,yr)

plt.show()
