import random 
import matplotlib.pyplot as plt
from math import *

################################################################################3

# function to calculate nCr

def ncr(n, r):
	p = 1
	k = 1
	if (n - r < r):
		r = n - r
	if (r != 0):
		while (r>0):
			p = p * n
			k = k * r
			n = int(n)
			r = int(r)
			p = int(p)
			k = int(k)
			m = gcd(p, k)
			p = p // m
			k = k // m
			n = n - 1
			r = r - 1


	else:
		p = 1
	return p

################################################################################3


# for meeting after N steps:

N=50 # number of steps

#computational graph

plot1 = plt.figure(1)
xr = [0]
yr1 = [1]
yr2 = [1]
yr3 = [0]
yr4 = [0]
zero = [0]
for INTERVAL in range(1,N,1):
	c=0
	d=0
	e=0
	zero.append(0)
	sqsum=0
	for j in range(INTERVAL**2):
		a=0
		b=0
		for i in range(INTERVAL):
			p1= random.randint(0,1)
			p2= random.randint(0,1)
			if(p1 == 0):
				p1=-1
			if(p2 == 0):
				p2=-1
			a = a + p1
			b = b + p2

		if(a==b):
			c=c+1
		if(a==b and a==0):
			d=d+1

		e=e+a;
		sqsum += a*a

	xr.append(INTERVAL)
	yr1.append(c/(INTERVAL**2))
	yr2.append(d/(INTERVAL**2))
	yr3.append(e/(INTERVAL**2))
	yr4.append((sqsum)/INTERVAL**2)

plt.title('Probability of meeting after N steps vs N')
plt.xlabel('N')
plt.ylabel('Probability')


#analytical graph
anay1 = [1]
anay2 = [1]
for i in range(1,N,1):
	c=ncr(2*i,i)/(4**i)
	if(i%2==1):
		d=0
	else:
		d=ncr(i,i/2)/(2**i)
		d=d*d
	anay1.append(c)
	anay2.append(d)

plt.plot(xr, anay1, label = "analytical", color = "orange")

#computational graph
plt.plot(xr, yr1, label = "computational", color = "#00b3ff")

plt.legend()

################################################################################

# for meeting at origin after N steps

plot2 = plt.figure(2)
plt.title('Probability of meeting at origin after N steps vs N')
plt.xlabel('N')
plt.ylabel('Probability')

#analytical graph
plt.plot(xr, anay2, label = "analytical",color = "orange")

# computational graph
plt.plot(xr, yr2, label = "computational",color = "#00b3ff")

plt.legend()

################################################################################

# mean displacement

plot2 = plt.figure(3)
plt.title('Mean Displacement of drunk after N steps vs N')
plt.xlabel('N')
plt.ylabel('Mean Displacement')
plt.plot(xr, zero, label = "analytical",color = "orange")
plt.plot(xr, yr3, label = "computational",color = "#00b3ff")
plt.legend()


################################################################################

# mean square displacement

plot2 = plt.figure(4)
plt.title('Mean Square Displacement of drunk after N steps vs N')
plt.xlabel('N')
plt.ylabel('Mean Square Displacement')
plt.plot(xr, xr, label = "analytical", color = "orange")
plt.plot(xr, yr4, label = "computational",color = "#00b3ff")
plt.legend()

################################################################################

plt.show()
