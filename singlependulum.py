import math
import numpy as np
import matplotlib.pyplot as plt
import numericalmethods as nm

h=0.1
omega = 2 * np.pi
trange = np.arange(0, 20, h)

def f(t):
	return np.cos(omega*t)

def d2ydx2(y):
	return -(omega**2)* y
	
truey = f(trange)
fity=[1.]
fitdydx=[0.]

for i in range(1, trange.size):
	oldy = fity[i-1]
	olddydx = fitdydx[i-1]
	
	fit, dydx = nm.leapfrog(d2ydx2, oldy, olddydx, h)
	
	fity.append(fit)
	fitdydx.append(dydx)
	
print len(truey), len(trange), len(fity), len(fitdydx)

plt.plot(trange,truey,  label="True Y")
plt.plot(trange,fity, label="Fity")
plt.legend()
plt.show()
