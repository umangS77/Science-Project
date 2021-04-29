import random
import numpy as np
import os,sys
import matplotlib.pyplot as plt
import json
M ,K,H = (np.random.rand() for i in range(3))
H = abs(H)

print("M is ",M,end="\n")
print("K is ",K ,end="\n")
print("H is ",H,end="\n")
t = np.linspace(0, 2*np.pi, 10000)

a = np.sqrt((2*H)/K)
b = np.sqrt(2*H*M)
x,p,p1 = ([] for i in range(3))

plt.plot( a*np.cos(t) , b*np.sin(t) )
plt.grid(color='lightgray',linestyle='-',linewidth = 0.5)
plt.ylabel("P")
plt.xlabel("X")
plt.show()

# plt.savefig("figure.png")