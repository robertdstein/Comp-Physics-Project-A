import math
import numpy as np
import matplotlib.pyplot as plt
import numericalmethods as nm
from matplotlib import animation
import cmath

g = 9.81
names = ["Explicit Euler", "Leapfrog", "Implicit Euler", "Runge-Kutta 4"]

class simulation:
	"""A container for a simulation
	"""
	def __init__(self, v0, f, h, trange, fps=10, masses=np.matrix([[1.0]]), lengths=np.matrix([[1.0]]), a=1.0):
		"""Args:
			v0:	Vector containing starting values for variables
			f: Function acting on vectors v to evolve system
			h: Stepsize for system evolution
			fps: Frames per second for generations of animations
			trange: Range of scaled t values to be simulated
			lengths: Matrix containing each pendulum length
			masses: Matrix containing each pendulum mass
			a: Time-scaling variable a
		"""
		print "Initial Conditions Vector:", v0
		self.methods =[]
		self.v0 = v0
		self.f = f
		self.h = h
		self.trange = trange
		self.fps = fps
		self.masses = masses
		self.lengths = lengths
		self.a = a

	def addfdm(self, methodname):
		"""Add a new FDM to the simulation
		Args:
			functionname: 	Name of fdm to be used.
							Can be one of:
								"Explicit Euler"
								"Leapfrog"
								"Implicit Euler"
								"Runge-Kutta 4"
		"""
		print "Simulating with FDM", methodname
		newmethod = fdm(functionname=methodname, v0=self.v0, f=self.f, h = self.h, fps = self.fps, lengths=self.lengths, masses=self.masses, trange = self.trange, a = self.a)
		self.methods.append(newmethod)
		
	def addallmethods(self):
		"""Add all available FDMs to simulation
		""" 
		for name in names:
			self.addfdm(name)
			
	def findstability(self):
		"""Find the stability threshold for each method
		"""
		vals = []
		#Iterate over methods
		for name in names:
			print "Finding Threshold for", name
			hval = 0.1
			hlimit = 0.001
			newmethod = fdm(functionname=name, v0=self.v0, f=self.f, h = hval, fps = self.fps, lengths=self.lengths, masses=self.masses, trange = self.trange, a = self.a)
			
			if newmethod.stable == "Stable":
				#Iterate over increasing h
				while (newmethod.stable == "Stable") and hval > hlimit:
					hval *= 1.02
					newmethod = fdm(functionname=name, v0=self.v0, f=self.f, h = hval, fps = self.fps, lengths=self.lengths, masses=self.masses, trange = self.trange, a = self.a)
			else:
				#Iterate over decreasing h	
				while (newmethod.stable != "Stable") and hval > hlimit:
					hval *= 0.98
					newmethod = fdm(functionname=name, v0=self.v0, f=self.f, h = hval, fps = self.fps, lengths=self.lengths, masses=self.masses, trange = self.trange, a = self.a)
			
			#Checks if minumum step has been reached
			if hval < hlimit:
				print "Minimum step", hlimit, "reached, true threshold unknown!"
				hval = None
			else:
				print "Threshold is h = ", hval
			vals.append(hval)
			
		#print identified thresholds
		print " \n Thresholds:"
		print names
		print vals
		
	def plottheta(self, title):
		"""Plots the theta values for all methods
		"""
		fig = plt.figure()
		
		#Loop over simulated FDMs
		for method in self.methods:
			x, y, label = method.toplot()
			plt.plot(x, y, label=label)
		
		plt.title(str(title) + ' with h = ' + str(self.h))
		plt.ylabel(r"$\Theta$[deg]")
		plt.xlabel('time (s)')
			
		plt.legend()
		
		fig.set_size_inches(15, 10)
		plt.savefig('graphs/' + str(title) + '_Theta.pdf')
		plt.close()
		
	def plotangles(self, title):
		"""Plots the angle values for all methods
		"""
		fig = plt.figure()
		
		#Loop over simulated FDMs
		for method in self.methods:
			x, ys, labels = method.anglesplot()
			#Loop over angles
			for i in range(0, len(ys)):
				y = ys[i]
				label = labels[i]
				plt.plot(x, y, label=label)
		
		plt.title(str(title) + ' with h = ' + str(self.h))
		plt.ylabel(r"$\Theta$[deg]")
		plt.xlabel('time (s)')
			
		plt.legend()
		
		fig.set_size_inches(15, 10)
		plt.savefig('graphs/' + str(title) + '_Angle.pdf')
		plt.close()
		
	def plottotalenergy(self, title):
		"""Plot the change in energy between each step
		"""
		fig = plt.figure()
		#Loop over simulated FDMs
		for i in range(0, len(self.methods)):
			method = self.methods[i]
			nrows = 2
			ncolumns = float(len(self.methods) + 1)/2.
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
		
	def plotcomponentenergy(self, title):
		"""Plot the energy components for each step
		"""
		fig = plt.figure()
		
		for i in range(0, len(self.methods)):
			method = self.methods[i]
			
			x, ys, labels = method.splitEplot()
			nrows = len(self.methods)
			
			#Plots GPE, KE and system energy in all cases
			#If there are 2 pendulums, plots Pendulum 1 and 2 Energies
			if len(ys) > 5:
				ax1 = plt.subplot(nrows,2,  (2*i)+1)
				ax2 = plt.subplot(nrows,2,  (i*2)+2)
			else:
				ax1 = None
				ax2 = plt.subplot(nrows,1,  i+1)	
			
				 
			for j in range(0, len(ys)):
				if j > 2:						
					y = ys[j]
					label = labels[j]
					ax2.plot(x, y, label=label)
				elif len(ys) > 5:
					y = ys[j]
					label = labels[j]
					ax1.plot(x, y, label=label)
				elif j == 2:
					y = ys[j]
					label = labels[j]
					ax2.plot(x, y, label=label)
				
				axs = [ax2]
				
				if len(ys) > 5:
					axs.append(ax1)
					
				for ax in axs:
					ax.set_title(method.name)
					ax.set_ylabel("E")
					ax.set_xlabel('time (s)')
					ax.legend()

		plt.suptitle(str(title) + ' with h = ' + str(self.h))	
		
		fig.set_size_inches(15, 20)
		
		if len(self.methods) == int(1):
			fig.set_size_inches(20, 8)
		plt.savefig('graphs/' + str(title) + '_Energy_Components.pdf')
		plt.close()
		
	def plotvelocity(self, title):
		"""Plots the velocity for each method
		"""
		fig = plt.figure()
		
		#Iterate over simulated FDMs
		for i in range(0, len(self.methods)):
			method = self.methods[i]
			nrows = 2
			ncolumns = float(len(self.methods) + 1)/2.
			plt.subplot(2, int(ncolumns) + int(ncolumns-int(ncolumns)), i+1)
			x, ys, labels = method.velplot()
			for i in range(0, len(ys)):
				y = ys[i]
				label = labels[i]
				plt.plot(x, y, label=label)
		
			plt.title(method.name)
			plt.ylabel("Vel m/s")
			plt.xlabel('time (s)')
			plt.legend()

		plt.suptitle(str(title) + ' with h = ' + str(self.h))	
		
		fig.set_size_inches(20, 20)
		plt.savefig('graphs/' + str(title) + '_Velocity.pdf')
		plt.close()
		
	def makeanimation(self,title):
		"""Make Animations for each method
		"""
		for method in self.methods:
			method.animation(title)
			
class fdm:
	"""A class for one finite differential method.
	"""
	def __init__(self, functionname, v0, f, h, fps, lengths, masses, trange, a):
		"""Args:
			functionname: 	Name of fdm to be used.
							Can be one of:
								"Explicit Euler"
								"Leapfrog"
								"Implicit Euler"
								"Runge-Kutta 4"
			v0:	Vector containing starting values for variables
			h: Stepsize for system evolution
			f: Function acting on vectors v to evolve system
			fps: Frames per second for generations of animations
			lengths: Matrix containing each pendulum length
			masses: Matrix containing each pendulum mass
			trange: Range of scaled t values to be simulated
			a: Time-scaling variable a
		"""
		self.name = functionname
		self.trange=trange
		self.yieldfunction(functionname)
		self.f = f
		self.h = h
		self.fps = fps
		self.frameskip = a/(h*fps)
		self.masses=masses
		self.lengths=lengths
		self.stable="Stable"
		self.a = a
		
		#Creates the first frame, then simulates solution
		self.fits=[]
		fit0 = frame(trange[0], v0, lengths, masses, a)
		self.fits.append(fit0)

		self.simulate()
		
	def yieldfunction(self, name):
		""" Matches a given name to a corresponding numerical method function.
		Functions are defined in numericalmethods.py
		"""
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
		"""Simulates a solution for the function
		"""
		
		#Adds a second frame if required for Leapfrog
		if self.function == nm.leapfrog:
			oldframe = self.fits[0]
			v = nm.rk4( oldframe=oldframe,f=self.f, h=self.h)
			fit1 = frame(self.trange[1], v, self.lengths, self.masses, a = self.a)	
			self.fits.append(fit1)
		
		#Finds the initial system energy
		oldE = self.fits[0].systemenergy
		
		#Loops over all scaled-t values to simulate each frame
		for i in range(len(self.fits), self.trange.size):
			t = self.trange[i]
			oldframe = self.fits[i-1]
			if i < 2:
				oldoldframe = None
			else:
				oldoldframe = self.fits[i-2]
			
			#Finds the variable values for the step, then creates a new frame
			v = self.function(oldframe=oldframe,f=self.f, h=self.h, oldoldframe = oldoldframe)
			
			fit = frame(t, v, self.lengths, self.masses, a = self.a)
			
			self.fits.append(fit)
			
			#Checks that the system has not exceeded initial energy
			newE = fit.systemenergy
			
			if newE > (oldE * 1.001):
				self.stable = "Unstable"
				
	def toplot(self):
		"""Returns an array containing each theta value
		Also a time value array and a label
		"""
		y = []
		times = []
		for oneframe in self.fits:
			yval = oneframe.vector.item(0)
			times.append(oneframe.time)	
			y.append(yval)
			
		label = self.name + " (" + str(self.stable) + ")"
		return times, y, label
		
	def anglesplot(self):
		"""Returns arrays containing each angle value
		Also a time value array and labels
		"""
		npend = len(self.fits[0].pendulums)
		if npend < 2:
			raise Exception("Insufficient Pendulums!")
		
		y = [[], []]
		labels = ["Theta", "Phi"]
		times = []
		
		for oneframe in self.fits:
			for i in range(0, len(oneframe.pendulums)):
				pen = oneframe.pendulums[i]
				angle = pen.theta
				y[i].append(angle)
			times.append(oneframe.time)	

		return times, y, labels
		
	def deltaEplot(self):
		"""Returns arrays containing change in energy
		Also a time value array and a label
		"""
		deltaEs = []
		times = []
		frame0 = self.fits[0]
		E0 = frame0.systemenergy
		
		for oneframe in self.fits:
			E1 = oneframe.systemenergy
			deltaE = E1-E0
			deltaEs.append(deltaE)
			E0 = E1
			times.append(oneframe.time)		
					
		label = self.name + " (" + str(self.stable) + ")"	
		return times, deltaEs, label
		
	def splitEplot(self):
		"""Returns arrays containing each energy component value
		Also a time value array and labels
		"""
		Es = [[], [], [], []]
		times = []
		npend = len(self.fits[0].pendulums)
		
		#Produces an array/label pair for each pendulum energy and their sum
		#Produces an array/label pair for system GPE/KE and their sum
		labels = []
		for i in range(0, npend):
			Es.append([])
			label = "Pendulum " + str(i+1)
			labels.append(label)
		labels.append("Sum of Pendulums")
		labels.append("GPE")
		labels.append("KE")
		labels.append("Sum of Components")
			
		for oneframe in self.fits:
			gpe = 0.0
			ke = 0.0
			E = 0.0
			times.append(oneframe.time)
			for i in range(0, npend):
				pen = oneframe.pendulums[i]
				Ei = pen.totalenergy
				Es[i].append(Ei)
				E += Ei
				gpe += pen.gpe
				ke += pen.ke
			Es[npend].append(E)
			Es[npend+1].append(gpe)
			Es[npend+2].append(ke)
			Es[npend+3].append(gpe+ke)
					
		return times, Es, labels
		
	def velplot(self):
		"""Returns arrays containing each velocity value
		Also a time value array and labels
		"""
		vels = []
		times=[]
		npend = len(self.fits[0].pendulums)
		
		#Produces an array/label pair for each velocity and their sum
		labels = []
		for i in range(0, npend):
			vels.append([])
			label = "Pendulum " + str(i+1)
			labels.append(label)
		vels.append([])
		labels.append("Sum")
			
		for oneframe in self.fits:
			sumvel = 0
			times.append(oneframe.time)
			for i in range(0, npend):
				pen = oneframe.pendulums[i]
				vel = abs(pen.velocity)
				vels[i].append(vel)
				sumvel += vel
			vels[(npend+1)-1].append(sumvel)
				
		return times, vels, labels
		
	def animation(self, title):
		"""Produces an animation of the simulation
		"""
		#Sets a figure larger than maximum pendulum radius
		plotwidth = np.sum(self.lengths) + 0.5
		
		fig = plt.figure()
		ax = fig.add_subplot(111, autoscale_on=False, xlim=(-plotwidth, plotwidth), ylim=(-plotwidth, plotwidth))
		
		#Define Initial values for pendulum, time, energy and velocity
		line, = ax.plot([], [], 'o-', color='red', lw=2)
		time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
		energy_text = ax.text(0.05, 0.8, '', transform=ax.transAxes)
		vel_text = ax.text(0.05, 0.7, '', transform=ax.transAxes)
		
		def animinit():
			line.set_data([], [])
			time_text.set_text('')
			energy_text.set_text('')
			vel_text.set_text('')
			return line, time_text, energy_text, vel_text
		
		#Sets system energy to indicate if E exceeds initial energy
		global baseE
		baseE = self.fits[0].systemenergy
		
		#Defines the animation function to be iterated over
		def animate(i):
			intval = int(i*self.frameskip)
			currentframe = self.fits[intval]
			x = [0.0]
			y = [0.0]
			vel = ""
			for pen in currentframe.pendulums:
				x.append(pen.xpos)
				y.append(pen.ypos)
				vel+= "Vel =" + str.format('{0:.3f}', abs(pen.velocity)) + " "
			time = "Time = " + str.format('{0:.1f}', currentframe.time) + " s"
			
			#Changes energy font colour if it exceeds initial energy
			E = currentframe.systemenergy
			global baseE
			if E > baseE:
				energy_text.set_color('red')
			else:
				energy_text.set_color('black')
			energy = "Energy = " + str.format('{0:.3f}', E) + " J"
			line.set_data(x, y)
			time_text.set_text(time)
			vel_text.set_text(vel)
			energy_text.set_text(energy)
			
			return line, time_text, energy_text, vel_text
			
		#Produce mp4 file and  close animation	
		
		nframes = int(float(len(self.trange))/float(self.frameskip))
		
		plt.title(str(title) + ' with ' + self.name + " and h = " + str(self.h))	
		anim = animation.FuncAnimation(fig,  animate, nframes, blit=True, init_func=animinit)
		anim.save('animations/' + str(title) + self.name+ '.mp4',fps=self.fps)
		plt.close()
		
				
class frame:
	"""A single snapshot of a simulation.
	Each frame contains one or more pendulums.
	"""
	def __init__(self, t, v, lengths, masses, a):
		"""	Args:
			t: Time of frame
			v: Vector containing all angles and derivatives
			lengths: Matrix containing each pendulum length
			masses: Matrix containing each pendulum mass
			a: Time-scaling variable a
		"""
		self.pendulums = []
		self.systemenergy = 0.0
		self.vector = v
		self.a = a
		
		#Calculate the True Time
		self.time=t/a
		
		xpos = 0.0
		ypos = 0.0

		equilibrium = 0.0
		vel = 0.0
		
		#Iterate over pendulums
		for i in range(0, len(lengths)):
			angle, dangle, length, mass = v.item(2*i), v.item((2*i)+1), lengths.item(i), masses.item(i)
			
			#Calculate the new x, y and equilibrium positions, as well as the unscaled velocity
			xpos += length*np.sin(angle)
			ypos += -length*np.cos(angle)
			
			equilibrium += -length
			
			vel += cmath.rect(length*dangle, angle)
			
			#Creates a new pendulum, adds it to the frame, and updates system energy
			
			newpendulum = pendulum(t, angle, dangle, xpos, ypos, length, mass, equilibrium, vel, a)
			
			self.pendulums.append(newpendulum)
			
			self.systemenergy += newpendulum.totalenergy
			
		
class pendulum:
	"""One pendulum within a frame.
	"""
	def __init__(self, t, angle, dangle, xpos, ypos, length, mass, equilibrium, vel, a):
		"""	Args:
			t: Time of frame
			angle: Angle between given pendulum and vertical
			dangle: Derivative of this angle w.r.t. time
			xpos: X position of pendulum
			ypos: Y position of pendulum
			length: Length of pendulum string
			mass: Mass of Pendulum
			equilibrium: Equilibrium height of pendulum at rest
			vel: Unscaled velocity of pendulum
			a: Time-scaling variable a

		"""
		self.mass = mass
		self.theta = angle
		self.dtheta =  dangle
		self.xpos = xpos
		self.ypos = ypos
		self.a = a
		
		#Calculate the unscaled time and velocity, as well as the energy
		self.time = t/a
		self.velocity = vel*a
				
		self.findenergy(equilibrium)
		
	def findenergy(self, equilibrium):
		"""Finds energy of system
		total energy = GPE + KE
		"""
		self.ke = 0.5 * self.mass * (abs(self.velocity)**2)
		self.gpe = g * (self.ypos-equilibrium) * self.mass
		self.totalenergy = self.ke + self.gpe
	
