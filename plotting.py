import math
import numpy as np
import matplotlib.pyplot as plt
import numericalmethods as nm
from matplotlib import animation

g = 9.81

class simulation:
	"""A container for a simulation
	"""
	def __init__(self, v0, f, h, trange, fps=5, masses=np.matrix([[1.0]]), lengths=np.matrix([[1.0]])):
		self.methods =[]
		self.v0 = v0
		self.f = f
		self.h = h
		self.trange = trange
		self.fps = fps
		self.masses = masses
		self.lengths = lengths

	def addfdm(self, methodname):
		newmethod = fdm(functionname=methodname, v0=self.v0, f=self.f, h = self.h, fps = self.fps, masses=self.masses, lengths=self.lengths, trange = self.trange)
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
		
	def plottotalenergy(self, title):
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
		
	def plotcomponentenergy(self, title):
		fig = plt.figure()
		
		for i in range(0, len(self.methods)):
			method = self.methods[i]
			nrows = 2
			ncolumns = float(len(self.methods))/2.
			plt.subplot(2, int(ncolumns) + int(ncolumns-int(ncolumns)), i+1)
			x, ys, labels = method.splitEplot()
			for i in range(0, len(ys)):
				y = ys[i]
				label = labels[i]
				plt.plot(x, y, label=label)
		
			plt.title(method.name)
			plt.ylabel("E")
			plt.xlabel('time (s)')
			plt.legend()

		plt.suptitle(str(title) + ' with h = ' + str(self.h))	
		
		fig.set_size_inches(20, 20)
		plt.savefig('graphs/' + str(title) + '_Energy_Components.pdf')
		plt.close()
		
	def makeanimation(self,title):
		for method in self.methods:
			method.animation(title)
			
class fdm:
	"""A class for one finite differential method
	"""
	def __init__(self, functionname, v0, f, h, fps, masses, lengths, trange):
		print "Initialising FDM", functionname
		self.name = functionname
		self.trange=trange
		self.yieldfunction(functionname)
		self.fits=[]
		fit0 = frame(trange[0], v0, lengths, masses)
		self.fits.append(fit0)
		self.f = f
		self.h = h
		self.fps = fps
		self.frameskip = 1.0/(h*fps)
		self.masses=masses
		self.lengths=lengths
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
			v = nm.rk4(f=self.f, oldframe=oldframe, h=self.h)
			fit1 = frame(self.trange[1], v, self.lengths, self.masses)	
			self.fits.append(fit1)
		
		oldE = self.fits[0].systemenergy
		
		for i in range(len(self.fits), self.trange.size):
			t = self.trange[i]
			oldframe = self.fits[i-1]
			if i < 2:
				oldoldframe = None
			else:
				oldoldframe = self.fits[i-2]
			
			v = self.function(f = self.f, oldframe=oldframe, h=self.h, oldoldframe = oldoldframe)
			
			fit = frame(t, v, self.lengths, self.masses)
			
			
			for j in range(0, len(fit.pendulums)):
				pen = fit.pendulums[j]
				oldpen = oldframe.pendulums[j]
				dist = math.sqrt((pen.xpos- oldpen.xpos)**2 + (pen.ypos - oldpen.ypos)**2)
				velestimate = dist/self.h
				
				#~ print "Vel estimate", velestimate, "v", pen.velocity, "oldv", oldpen.velocity
			
			self.fits.append(fit)
			
			newE = fit.systemenergy
			
			if newE > oldE:
				if newE > (oldE + (10**-10)):
					self.stable = "Unstable"
				else:
					self.stable = "Stable to threshold"
			else:
				pass
				
			#~ oldE = newE
				
	def toplot(self):
		y = []
		for oneframe in self.fits:
			yval = oneframe.vector.item(0)
			y.append(yval)
			
		label = self.name + " (" + str(self.stable) + ")"
		return self.trange, y, label
		
	def deltaEplot(self):
		deltaEs = []
		frame0 = self.fits[0]
		E0 = frame0.systemenergy
		
		for oneframe in self.fits:
			E1 = oneframe.systemenergy
			deltaE = E1-E0
			deltaEs.append(deltaE)
			E0 = E1
					
		label = self.name + " (" + str(self.stable) + ")"	
		return self.trange, deltaEs, label
		
	def splitEplot(self):
		Es = [[], [], [], []]
		npend = len(self.fits[0].pendulums)
		
		labels = []
		for i in range(0, npend):
			Es.append([])
			label = "Pendulum " + str(i+1)
			labels.append(label)
		labels.append("GPE")
		labels.append("KE")
		labels.append("Sum Pendulums")
		labels.append("Sum components")
			
		for oneframe in self.fits:
			gpe = 0.0
			ke = 0.0
			E = 0.0
			for i in range(0, npend):
				pen = oneframe.pendulums[i]
				Ei = pen.totalenergy
				Es[i].append(Ei)
				E += Ei
				gpe += pen.gpe
				ke += pen.ke
			Es[npend].append(gpe)
			Es[npend+1].append(ke)
			Es[npend+2].append(E)
			Es[npend+3].append(ke + gpe)
					
		return self.trange, Es, labels
		
	def animation(self, title):
		plotwidth = np.sum(self.lengths) + 0.5
		
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
		baseE = self.fits[0].systemenergy
		
		def animate(i):
			intval = int(i*self.frameskip)
			currentframe = self.fits[intval]
			x = [0.0]
			y = [0.0]
			for pen in currentframe.pendulums:
				x.append(pen.xpos)
				y.append(pen.ypos)
			time = "Time = " + str(currentframe.time) + " s"
			E = currentframe.systemenergy
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
	def __init__(self, t, v, lengths, masses):
		self.pendulums = []
		self.systemenergy = 0.0
		self.vector = v
		self.time=t
		allxpos = [0.0]
		allypos = [0.0]

		intlength = 0.0
		sumlength = 0.0
		for val in lengths:
			sumlength += val.item(0)
			intlength += sumlength
		avlength = 0.0
		
		oldthetas = []
		olddthetas = []
		oldlengths=[]
		
		for i in range(0, len(lengths)):
			theta, dtheta, length, mass = v.item(2*i), v.item((2*i)+1), lengths.item(i), masses.item(i)
			xpos = allxpos[i] + length*np.sin(theta)
			allxpos.append(xpos)
			ypos = allypos[i] -length*np.cos(theta)
			allypos.append(ypos)
			
			oldthetas.append(theta)
			olddthetas.append(dtheta)
			oldlengths.append(length)
			avlength += length
			vx=0.0
			vy = 0.0
			lastdt=0.0
			vel=0.0
			
			for j in range(0, len(olddthetas)):
				oldt = oldthetas[j]
				olddt = olddthetas[j]
				oldl = oldlengths[j]
				oldx = allxpos[j]
				oldy = allypos[j]
				dist = math.sqrt((xpos-oldx)**2 + (ypos-oldy)**2)
				vx += oldl * np.sin(oldt)*(olddt)
				vy += oldl * np.cos(oldt)*(olddt)
				#vel += np.absolute(dist * olddt)
				lastdt = olddt
				
			#~ if vel > 1.00:
				#~ print "Theta", theta, "Vel", vel, "Dtheta", dtheta
			
			newpendulum = pendulum(t, theta, dtheta, xpos, ypos, vx, vy, length, mass, avlength, vel)
			self.systemenergy += newpendulum.totalenergy
			
			self.pendulums.append(newpendulum)
			
		
class pendulum:
	"""One pendulum within a frame
	"""
	def __init__(self, t, theta, dtheta, xpos, ypos, vx, vy, length, mass, avlength, vel):
		self.time = t
		self.mass = mass
		self.theta = theta
		self.dtheta =  dtheta
		self.xpos = xpos
		self.ypos = ypos
		self.vx= vx
		self.vy = vy
		self.velocity = math.sqrt((vx**2) + (vy**2))
		self.findenergy(avlength)
		
	def findenergy(self, avlength):
		"Finds energy of system"
		self.ke = 0.5 * (self.velocity**2) * self.mass
		self.gpe = g * (self.ypos + avlength) * self.mass
		self.totalenergy = self.ke + self.gpe
	
