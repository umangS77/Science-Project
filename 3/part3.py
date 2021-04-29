import random
import numpy as np
import os,sys
import matplotlib.pyplot as plt
import json

T,kb,alpha,tau,m = (np.random.rand() for i in range(5))
tau = tau*0.1
msd = []
print("T- ",T,end="\n")
print("kb - ",kb,end="\n")
print("alpha -",alpha,end="\n")
print("tau -",tau,end="\n")
print("m -",m,end="\n")

t = np.linspace(0, 100, 10000)

for i in t:

    if i*tau < 0.3:
        msd.append((kb*T*i**2)/m) 
    else:
        msd.append((2*kb*T*i)/alpha)

msd = np.array(msd)
plt.ylabel("MSD")
plt.xlabel("time")
plt.plot( t,msd )
plt.grid(color='lightgray',linestyle='-',linewidth = 0.5)
plt.show()
# plt.savefig("msd.png")