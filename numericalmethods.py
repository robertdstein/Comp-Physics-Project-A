
def expliciteuler(oldframe, f, h, oldoldframe=None):
	"""Uses the Explicit Euler method to evolve the variable vector.
	
	Args:
		oldframe: previous state to provide starting vector.
		f: the fuction to yield each variable derivative.
		h: stepsize to evolve variables
	
	Returns:
		New variable vector
	"""
	
	v = oldframe.vector
	deltav =  f(v)
	newv = v + (deltav*h)
	return  newv
	
def leapfrog(oldframe, f, h, oldoldframe):
	"""Uses the Leapfrog method to evolve the variable vector.
	
	Args:
		oldframe: most recent state to provide n-1 vector.
		oldoldframe: penultimate state to provide n-2 vector
		f: the fuction to yield each variable derivative.
		h: stepsize to evolve variables
	
	Returns:
		New variable vector
	"""
	if oldoldframe == None:
		raise Exception("No n-2 frame found!")
	else:
		v = oldframe.vector
		oldv = oldoldframe.vector
		deltav =  f(v)
		newv = oldv + (2*deltav*h)
		return newv
	
def rk4(oldframe, f, h, oldoldframe=None):
	"""Uses the Runge-Kutta 4 method to evolve the variable vector.
	
	Args:
		oldframe: previous state to provide starting vector.
		f: the fuction to yield each variable derivative.
		h: stepsize to evolve variables
	
	Returns:
		New variable vector
	"""
	v = oldframe.vector
	k1 = f(v)
	k2 = f(v + k1*h/2.)
	k3 = f(v + k2*h/2.)
	k4 = f(v + k3*h)
	newv = v + (h*(k1 + 2*(k2+k3) + k4)/6.)
	return  newv
	
def impliciteuler(oldframe, f, h, oldoldframe=None):
	"""Uses the Implicit Euler method to evolve the variable vector.
	
	Args:
		oldframe: previous state to provide starting vector.
		f: the fuction to yield each variable derivative.
		h: stepsize to evolve variables
	
	Returns:
		New variable vector
	"""
	v = oldframe.vector
	guessv = (h*f(v)) + v
	deltav = f(guessv)
	newv = v + (h*deltav)
	return newv
