import math
import numpy as np
import matplotlib.pyplot as plt
import plotting

h=0.5
omega = 2 * np.pi
omega = 1
trange = np.arange(0, 8, h)

def gfunc(t):
	return np.cos(omega*t)
	
def truedg1(t):
	return - omega * np.sin(omega * t)
	
def dg1(y, olddydx):
	d1 = olddydx + (dg2(y) * h)
	return d1

def dg2(y):
	return (-omega**2)* y
	
def hfunc(t):
	return np.exp(-omega*t)
	
def truedh1(y):
	return dh1(y)
	
def dh1(y, olddydx=0):
	d1 = - omega * y
	return d1
	
def dh2(y):
	d2 = (omega **2)*y
	return d2 
	
functions = [[gfunc, truedg1, dg1, dg2], [hfunc, truedh1, dh1, dh2]]
funclabels = ["y = cos (kx)", "y = exp(-kx)"]

for j in range (0, len(functions)):
	plt.subplot(1, len(functions), j+1)
	[f, truedf1, df1, df2] = functions[j]
	functionname = funclabels[j]
	plotting.run(f, truedf1, df1, df2, trange, h, functionname)

	
plt.suptitle('Fits with h = ' + str(h))
	
plt.legend()
plt.tight_layout()
plt.savefig('graphs/singlependulum.pdf')
plt.close()
