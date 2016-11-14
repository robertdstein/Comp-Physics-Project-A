import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from classes import *

#Set the initial stepsize
h=0.02

#Set the time range in seconds
tmax = 20

#Set the initial values for theta and d(theta)/dt
y0 = [0.1]
dy0 = [0.0]

v0 = np.matrix([y0, dy0])

#Set the alternate initial values for theta and d(theta)/dt
vbig = np.matrix([[0.75 * math.pi], [0.0]])

#Sets the pendulum length and mass
l= [1.0]
m= [1.0]

#calculates the scaling variables a and b
a = math.sqrt(g/(l[0]))
b= a*m[0]*l[0]

trange = np.arange(0, tmax*a, h)

#Toggle on to find threshold values for h
#Toggle off for much faster code running
find_thresholds=False

#Iterates over damping values and initial values

for D in [0.0, 0.2]:
	for v in [v0, vbig]:
		print "Running with Damping Constant D =", D
	
		def f(vector):
			d0 = vector.item(1)
			d1 = -np.sin(vector.item(0)) - (D *vector.item(1))
			d = np.matrix([[d0],[d1]])
			return d	
		
		#Generate and print a dynamic title 
		if D > 0.0:	
			title = "Single_Pendulum_Damped"
		else:
			title = "Single_Pendulum_Undamped"
		
		if v[0] == vbig[0]:
			title += "_Large"
		
		print title
		
		#Simulates solutions using all FDMs
		
		sim = simulation(v, f=f, trange=trange, h=h, a=a)
		sim.addallmethods()
		
		#Produce graphs for theta the change in energy, and the components of total energy	
		sim.plottheta(title=title)
		sim.plottotalenergy(title=title)
		sim.plotcomponentenergy(title=title)
		
		#Optionally find Stability Thresholds
		if find_thresholds:
			sim.findstability()
		
		#Optional line to produce animated mp4 files illustrating the simulation
		#Uncomment line to test
		#~ sim.makeanimation(title=title)
		
		
	
