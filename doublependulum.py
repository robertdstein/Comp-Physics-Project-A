import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from classes import *

#Set the initial stepsize
h=0.02

#Set the time range in seconds
tmax = 20

#Set the initial values for theta, d(theta)/dt, phi and d(phi)/dt
y0 = [0.1]
dy0 = [0.0]
z0 = [0.0]
dz0 = [0.0]

v0 = np.matrix([y0, dy0, z0, dz0])

#Set the value for the lengths (must be equal in this case)
l= np.matrix([[1.0], [1.0]])

#Create a range of ratios mass2/mass1
msmall= np.matrix([[1.0], [0.01]])
meq = np.matrix([[1.0], [1.0]])
mbig = np.matrix([[1.0], [100.0]])

ms = [msmall, meq, mbig]

#Iterate over masses and undamped/damped simulations
for m in ms:
	for D in [0.0, (m[0]*math.sqrt(g*l.item(0)))]:
		
		#Calculate the scaling variables a, R, G
		a = math.sqrt(g/(l[0]))
		R = m[1]/m[0]
		G = D/(m[0]*math.sqrt(g*l.item(0)))
		
		print "Running with Damping Constant G =", G, "and ratio R =", R
		
		trange = np.arange(0, tmax*a, h)
		
		def f(vector):
			"""The Differential Equations for evolving the 4D vector.
			
			For a double pendulum, vector format is a matrix of:
			
			[theta]
			[d(theta)/dt]
			[phi]
			[d(phi)/dt]
			
			Returns a 4D vector representing the change in each variable over the step.
			"""
			d0 = vector.item(1)
			d1 = -((R+1) * vector.item(0)) - (G*vector.item(1)) + (R * vector.item(2))
			d2 = vector.item(3)
			d3 = ((R+1) * vector.item(0)) + (G*(1-(1./R))*vector.item(1)) - ((R+1) * vector.item(2)) - (vector.item(3)*G/R)
			d = np.matrix([[d0],[d1], [d2], [d3]])
			return d
			
		#Generate and print a dynamic title 	
		title = "Double_Pendulum"
		if R<1.0:
			title+="_Small"
		elif R > 1.0:
			title += "_Big"
		else:
			title += "_Eq"
			
		if G > 0.0:
			title += "_Damped"
		else:
			title += "_Undamped"
			
		print title
		
		#Create a simulation with the Runge-Kutta 4 method
		sim = simulation(v0, f=f, trange=trange, h=h, lengths=l, masses=m, a=a)
		sim.addfdm("Runge-Kutta 4")
		
		#Check whether the solution is stable
		for method in sim.methods:
			print method.name, "is", method.stable
			
		#Produce graphs for theta, both angles, the change in energy, the components of total energy, and velocity	
		sim.plottheta(title=title)
		sim.plotangles(title=title)
		sim.plottotalenergy(title=title)
		sim.plotcomponentenergy(title=title)
		sim.plotvelocity(title=title)
		
		#Optional line to produce animated mp4 files illustrating the simulation
		#Uncomment line to test
		#~ sim.makeanimation(title=title)
