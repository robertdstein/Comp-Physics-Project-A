import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from plotting import *

h=0.02
omega = 2 * np.pi
trange = np.arange(0, 10, h)

y0 = [1.0]
dy0 = [0.0]
z0 = [0.0]
dz0 = [0.0]

v0 = np.matrix([y0, dy0, z0, dz0])

print "v0", v0


l= np.matrix([[1.0], [1.0]])
m= np.matrix([[1.0], [3.0]])

a = 1.0
b= 1.0

D=0.2

R = m[0]/m[1]
G = D/(m[0]*math.sqrt(g*l.item(0)))

def f(vector):
	d0 = vector.item(1)
	d1 = -(R+1) * vector.item(0) - G*vector.item(1) + R * vector.item(2)
	d2 = vector.item(3)
	d3 = (R+1) * vector.item(0) - G*(1-(1./R))*vector.item(1) - (R+1) * vector.item(2) - vector.item(3)*G/R
	d = np.matrix([[d0],[d1], [d2], [d3]])
	return d
	

#~ fig = plt.figure()		
title = "Double_Pendulum"
sim = simulation(v0, f=f, trange=trange, h=h, lengths=l, masses=m)
#~ sim.addfdm("Implicit Euler")
sim.addallmethods()
sim.plottheta(title=title)
sim.plottotalenergy(title=title)
sim.makeanimation(title=title)
	
