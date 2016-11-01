import math
import numpy as np
import matplotlib.pyplot as plt
import numericalmethods as nm	

def run(y0, dy0, df1, df2, trange, h, functionname,truey=None):
	
	labels = ["Explicit Euler", "Leapfrog", "Implicit Euler", "Runge-Kutta 4"]
	methods = [nm.expliciteuler, nm.leapfrog, nm.impliciteuler, nm.rk4]
	
	data = []
	
	for k in range(0, len(methods)):
		
		label = labels[k]
		method = methods[k]
		
		fity=[y0]
		fitdydx=[dy0]
		
		if method == nm.leapfrog:
			fit, dydx = nm.rk4(df1=df1, df2=df2, yi=y0, olddydx=dy0, h=h)		
			fity.append(fit)
			fitdydx.append(dydx)
		
		print label, fity, fitdydx
		
		for i in range(len(fity), trange.size):
			oldy = fity[i-1]
			olddydx = fitdydx[i-1]
			
			oldoldy = fity[i-2]
			oldolddydx = fitdydx[i-2]
			
			fit, dydx = method(df1=df1, df2=df2, yi=oldy, olddydx=olddydx, h=h, oldoldy=oldoldy, oldolddydx=oldolddydx)
			
			fity.append(fit)
			fitdydx.append(dydx)
		
		plt.plot(trange,fity, label=label)
		
		data.append([label, trange, fity])
		
	if truey != None:
		if truey[0] == y0:
			plt.plot(trange,truey, '--', color='black', label="True Y")
		else:
			raise Exception("y0 does not match first entry in true y")
		
	plt.title(functionname)
	
	return data
	
