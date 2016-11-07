import math
import numpy as np
import matplotlib.pyplot as plt
import numericalmethods as nm
from matplotlib import animation

g = 9.81

class simulation:
	"""A container for a simulation
	"""
	def __init__(self, theta0, dtheta0, df1, df2, h, trange, fps=5, mass=1.0, length=1.0):
		self.methods =[]
		self.theta0 = theta0
		self.dtheta0 = dtheta0
		self.df1 = df1
		self.df2 = df2
		self.h = h
		self.trange = trange
		self.fps = fps
		self.mass = mass
		self.length = length

	def addfdm(self, methodname):
		newmethod = fdm(functionname=methodname, theta0=self.theta0, dtheta0=self.dtheta0, df1 = self.df1, df2 = self.df2, h = self.h, fps = self.fps, mass=self.mass, length=self.length, trange = self.trange)
		self.methods.append(newmethod)
		
	def addallmethods(self):
		for name in ["Explicit Euler", "Leapfrog", "Implicit Euler", "Runge-Kutta 4"]:
			self.addfdm(name)
		
	def plottheta(self, title, exactsolution=None):
		fig = plt.figure()
		
		for method in self.methods:
			x, y, label = method.toplot()
			plt.plot(x, y, label=label)
			
		if exactsolution != None:
			y = exactsolution(self.trange)
			plt.plot(self.trange, y, label="Exact Solution")
		
		plt.title(str(title) + ' with h = ' + str(self.h))
		plt.ylabel(r"$\Theta$[deg]")
		plt.xlabel('time (s)')
			
		plt.legend()
		
		fig.set_size_inches(15, 10)
		plt.savefig('graphs/' + str(title) + '_Theta.pdf')
		plt.close()
		
	def plotenergy(self, title):
		fig = plt.figure()
		
		
		for i in range(0, len(self.methods)):
			method = self.methods[i]
			nrows = 2
			ncolumns = float(len(self.methods))/2.
			plt.subplot(2, int(ncolumns) + int(ncolumns-int(ncolumns)), i+1)
			x, y, label = method.deltaEplot()
			plt.plot(x, y, label=label)
		
			plt.title(method.name)
			plt.ylabel(r"$\Delta$E")
			plt.xlabel('time (s)')

		plt.suptitle(str(title) + ' with h = ' + str(self.h))	

		fig.set_size_inches(20, 20)
		plt.savefig('graphs/' + str(title) + '_Energy.pdf')
		plt.close()
		
	def makeanimation(self,title):
		for method in self.methods:
			method.animation(title)
			
class fdm:
	"""A class for one finite differential method
	"""
	def __init__(self, functionname, theta0, dtheta0, df1, df2, h, fps, mass, length, trange):
		print "Initialising FDM", functionname
		self.name = functionname
		self.trange=trange
		self.yieldfunction(functionname)
		self.fits=[]
		fit0 = frame(trange[0], theta0, dtheta0, length, mass)
		self.fits.append(fit0)
		self.df1 = df1
		self.df2 = df2
		self.h = h
		self.fps = fps
		self.frameskip = 1.0/(h*fps)
		self.mass=mass
		self.length=length
		self.stable="Absolutely Stable"
		self.simulate()
		
	def yieldfunction(self, name):
		if name == "Explicit Euler":
			self.function = nm.expliciteuler
		elif name == "Leapfrog":
			self.function = nm.leapfrog
		elif name == "Implicit Euler":
			self.function = nm.impliciteuler
		elif name == "Runge-Kutta 4":
			self.function = nm.rk4
		else:
			raise Exception("Function "+ str(functionname) + " is not found!")
		
	def simulate(self):
		if self.function == nm.leapfrog:
			oldframe = self.fits[0]
			theta, dtheta = nm.impliciteuler(df1=self.df1, df2=self.df2, oldframe=oldframe, h=self.h)
			fit1 = frame(self.trange[1], theta, dtheta, self.length, self.mass)	
			self.fits.append(fit1)
		
		oldE = self.fits[0].totalenergy
		
		for i in range(len(self.fits), self.trange.size):
			t = self.trange[i]
			oldframe = self.fits[i-1]
			#~ oldE = oldframe.totalenergy
			if i < 2:
				oldoldframe = None
			else:
				oldoldframe = self.fits[i-2]
			
			theta, dtheta = self.function(df1=self.df1, df2=self.df2, oldframe=oldframe, h=self.h, oldoldframe = oldoldframe)
			
			fit = frame(t, theta, dtheta, self.length, self.mass)
			
			self.fits.append(fit)
			
			newE = fit.totalenergy
			
			if newE > oldE:
				if newE > (oldE + (10**-10)):
					self.stable = "Unstable"
				else:
					self.stable = "Stable to threshold"
			else:
				pass
				
			oldE = newE
				
	def toplot(self):
		y = []
		for oneframe in self.fits:
			yval = oneframe.theta 
			y.append(yval)
			
		label = self.name + " (" + str(self.stable) + ")"
		return self.trange, y, label
		
	def deltaEplot(self):
		deltaEs = []
		frame0 = self.fits[0]
		E0 = frame0.totalenergy
		
		for oneframe in self.fits:
			E1 = oneframe.totalenergy
			deltaE = E1-E0
			deltaEs.append(deltaE)
			E0 = E1
					
		label = self.name + " (" + str(self.stable) + ")"	
		return self.trange, deltaEs, label
		
	def animation(self, title):
		plotwidth = self.length + 0.5
		
		fig = plt.figure()
		ax = fig.add_subplot(111, autoscale_on=False, xlim=(-plotwidth, plotwidth), ylim=(-plotwidth, plotwidth))
		line, = ax.plot([], [], 'o-', color='red', lw=2)
		time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
		energy_text = ax.text(0.05, 0.8, '', transform=ax.transAxes)
		
		def animinit():
			line.set_data([], [])
			time_text.set_text('')
			energy_text.set_text('')
			return line, time_text, energy_text
			
		global baseE
		baseE = self.fits[0].totalenergy
		
		#~ def updatebase(E):
			#~ global baseE
			#~ baseE = E
		
		def animate(i):
			intval = int(i*self.frameskip)
			currentframe = self.fits[intval]
			x = [0.0, currentframe.xpos]
			y = [0.0, currentframe.ypos]
			time = "Time = " + str(currentframe.time) + " s"
			E = currentframe.totalenergy
			global baseE
			if E > baseE:
				energy_text.set_color('red')
			else:
				energy_text.set_color('black')
			energy = "Energy = " + str(E) + " J"
			line.set_data(x, y)
			time_text.set_text(time)
			energy_text.set_text(energy)
			
			baseE = E
			return line, time_text, energy_text
		
		nframes = int(float(len(self.trange))/float(self.frameskip))
		
		plt.title(str(title) + ' with ' + self.name + " and h = " + str(self.h))	
		anim = animation.FuncAnimation(fig, animate, nframes, blit=True, init_func=animinit)
		anim.save('animations/' + str(title) + self.name+ '.mp4')
		plt.close()
				
class frame:
	"""A single snapshot in a simulation
	"""
	def __init__(self, t, theta, dtheta, length, mass):
		self.time = t
		self.mass = mass
		self.theta = theta
		self.dtheta =  dtheta
		self.xpos = length*np.sin(theta)
		self.ypos = -length*np.cos(theta)
		self.velocity = length * dtheta
		self.findenergy(length)
		
	def findenergy(self, length):
		"Finds energy of system"
		self.kpe = 0.5 * (self.velocity**2) * self.mass
		self.gpe = g * (length+self.ypos) * self.mass
		self.totalenergy = self.kpe + self.gpe
		
#~ def run(y0, dy0, df1, df2, trange, h, functionname,truey=None):
	
	#~ labels = ["Explicit Euler", "Leapfrog", "Implicit Euler", "Runge-Kutta 4"]
	#~ methods = [nm.expliciteuler, nm.leapfrog, nm.impliciteuler, nm.rk4]
	
	#~ data = []
	
	#~ for k in range(0, len(methods)):
		
		#~ label = labels[k]
		#~ method = methods[k]
		
		#~ fity=[y0]
		#~ fitdydx=[dy0]
		
		#~ if method == nm.leapfrog:
			#~ fit, dydx = nm.rk4(df1=df1, df2=df2, yi=y0, olddydx=dy0, h=h)		
			#~ fity.append(fit)
			#~ fitdydx.append(dydx)
		
		#~ print label, fity, fitdydx
		
		#~ for i in range(len(fity), trange.size):
			#~ oldy = fity[i-1]
			#~ olddydx = fitdydx[i-1]
			
			#~ oldoldy = fity[i-2]
			#~ oldolddydx = fitdydx[i-2]
			
			#~ fit, dydx = method(df1=df1, df2=df2, yi=oldy, olddydx=olddydx, h=h, oldoldy=oldoldy, oldolddydx=oldolddydx)
			
			#~ fity.append(fit)
			#~ fitdydx.append(dydx)
		
		#~ plt.plot(trange,fity, label=label)
		
		#~ data.append([label, trange, fity, fitdydx])
		
	#~ if truey != None:
		#~ if truey[0] == y0:
			#~ plt.plot(trange,truey, '--', color='black', label="True Y")
		#~ else:
			#~ raise Exception("y0 does not match first entry in true y")
		
	#~ plt.title(functionname)
	
	#~ return data
	
