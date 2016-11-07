import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from plotting import *

h=0.002
omega = 2 * np.pi
trange = np.arange(0, 100, h)

y0 = [1.0]
dy0 = [0.0]

v0 = np.matrix([y0, dy0])

print "v0", v0


l= [1.0]
m= [1.0]

a = 1.0
b= 1.0

D=0.2

def f(vector):
	d0 = vector.item(1)
	c1 = g/(l[0]*a**2)
	c2 = (b*D)/(l[0]*m[0]*a)
	d1 = -c1* np.sin(vector.item(0)) - (c2 * vector.item(1))
	d = np.matrix([[d0],[d1]])
	return d
	

#~ fig = plt.figure()		
title = "Single_Pendulum"
sim = simulation(v0, f=f, trange=trange, h=h)
#~ sim.addfdm("Implicit Euler")
sim.addallmethods()
sim.plottheta(title=title)
sim.plottotalenergy(title=title)
#~ sim.makeanimation(title=title)
	
