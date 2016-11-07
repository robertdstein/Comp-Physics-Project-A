import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from plotting import *

h=0.02
omega = 2 * np.pi
trange = np.arange(0, 1000, h)

y0 = 1.0
dy0 = 0.0


l= 1.0
m= 1.0

a = 1.0
b= 1.0

D=0.2

c1 = g/(l*a**2)
c2 = (b*D)/(l*m*a)

print "y0", y0, "dy0", dy0
print "c1", c1, "c2", c2
	
def df1(y, olddydx):
	d1 = olddydx + (df2(y, olddydx) * h)
	return d1

def df2(y, olddydx):
	return -c1* np.sin(y) - (c2 * olddydx)

#~ fig = plt.figure()		
title = "Single_Pendulum"
sim = simulation(theta0=y0, dtheta0=dy0, df1=df1, df2=df2, trange=trange, h=h)
sim.addallmethods()
sim.plottheta(title=title)
sim.plotenergy(title=title)
sim.makeanimation(title=title)

#~ plt.title('Single Pendulum with h = ' + str(h))
#~ plt.ylabel(r"$\Theta$[deg]")
#~ plt.xlabel('time (s)')
	
#~ plt.legend()
#~ plt.tight_layout()

#~ fig.set_size_inches(15, 10)
#~ plt.savefig('graphs/singlependulum.pdf')
#~ plt.close()


#~ for entry in data:
	
	#~ label = entry[0]
	#~ times = entry[1]
	#~ thetas = entry[2]
	#~ dthetas = entry[3]
	
	#~ xvals=[]
	#~ yvals=[]
	#~ gpevals = []
	#~ kevals = []
	
	
	#~ for i in range(0, len(thetas)):
		#~ val = thetas[i]
		#~ dval = dthetas[i]
		
		
		#~ xpos = l*np.sin(val)
		#~ ypos = -l*np.cos(val)
		#~ xvals.append(xpos)
		#~ yvals.append(ypos)
		#~ if ypos > 0:
			#~ print ypos
	
	#~ plotwidth = l + 0.5
	
	#~ fig = plt.figure()
	#~ ax = fig.add_subplot(111, autoscale_on=False, xlim=(-plotwidth, plotwidth), ylim=(-plotwidth, plotwidth))
	#~ line, = ax.plot([], [], 'o-', color='red', lw=2)
	#~ time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
	
	#~ def init():
		#~ line.set_data([], [])
		#~ time_text.set_text('')
		#~ return line, time_text
		
	#~ fps = 5
		
	#~ frameskip=1.0/(h*fps)
	
	#~ nframes = int(float(len(trange))/float(frameskip))
	#~ print len(times), frameskip, nframes
	
	#~ def animate(i):
		#~ intval = int(i*frameskip)
		#~ x = [0.0, xvals[intval]]
		#~ y = [0.0, yvals[intval]]
		#~ time = "Time = " + str(times[intval]) + " s"
		#~ line.set_data(x, y)
		#~ time_text.set_text(time)
		#~ return line, time_text
	
	#~ plt.title('Single Pendulum with ' + label + " and h = " + str(h))	
	#~ anim = animation.FuncAnimation(fig, animate, nframes, blit=True, init_func=init)
	#~ anim.save('animations/single_pendulum' + label+ '.mp4')
	#~ plt.close()
	
