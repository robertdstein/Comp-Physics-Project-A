import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import plotting

h=0.02
omega = 2 * np.pi
trange = np.arange(0, 40, h)

y0 = 1.5
dy0 = 0.0

g = 9.81
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

fig = plt.figure()		
functionname = "Single Pendulum"
data = plotting.run(y0, dy0, df1, df2, trange, h, functionname)

plt.title('Single Pendulum with h = ' + str(h))
plt.ylabel(r"$\Theta$[deg]")
plt.xlabel('time (s)')
	
plt.legend()
#~ plt.tight_layout()

fig.set_size_inches(15, 10)
plt.savefig('graphs/singlependulum.pdf')
plt.close()


for entry in data:
	
	label = entry[0]
	times = entry[1]
	thetas = entry[2]
	xvals=[]
	yvals=[]
	for val in thetas:
		xpos = l*np.sin(val)
		ypos = -l*np.cos(val)
		xvals.append(xpos)
		yvals.append(ypos)
		if ypos > 0:
			print ypos
		
	fig = plt.figure()
	ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
	line, = ax.plot([], [], 'o-', lw=2)
	
	def init():
	    line.set_data([], [])
	    return line
	    
	    
	fps = 5
	    
	frameskip=1.0/(h*fps)
	
	nframes = int(float(len(trange))/float(frameskip))
	print len(times), frameskip, nframes
	
	def animate(i):
		x = [0.0, xvals[int(i*frameskip)]]
		y = [0.0, yvals[int(i*frameskip)]]
		line.set_data(x, y)
		return line
		
	anim = animation.FuncAnimation(fig, animate, nframes, blit=True, init_func=init)
	anim.save('animations/single_pendulum' + label+ '.mp4')
	plt.close()
	
