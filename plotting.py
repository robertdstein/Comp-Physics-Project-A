import math
import numpy as np
import matplotlib.pyplot as plt
import numericalmethods as nm	

def run(f, truedf1, df1, df2, trange, h, functionname):
	
	truey = f(trange)
	
	labels = ["Explicit Euler", "Leapfrog", "Implicit Euler", "Runge-Kutta 4"]
	methods = [nm.expliciteuler, nm.leapfrog, nm.impliciteuler, nm.rk4]
	
	for k in range(0, len(methods)):
		
		label = labels[k]
		method = methods[k]
		
		fity=[truey[0],]
		fitdydx=[truedf1(truey[0]),]
		
		if method == nm.leapfrog:		
			fity.append(truey[1])
			fitdydx.append(truedf1(truey[1]))
		
		for i in range(len(fity), trange.size):
			oldy = fity[i-1]
			olddydx = fitdydx[i-1]
			
			oldoldy = fity[i-2]
			oldolddydx = fitdydx[i-2]
			
			fit, dydx = method(df1=df1, df2=df2, yi=oldy, olddydx=olddydx, h=h, oldoldy=oldoldy, oldolddydx=oldolddydx)
			
			fity.append(fit)
			fitdydx.append(dydx)
		
		plt.plot(trange,fity, label=label)
	
	plt.plot(trange,truey, '--', color='black', label="True Y")
	plt.ylim(-2,2)
	plt.title(functionname)
