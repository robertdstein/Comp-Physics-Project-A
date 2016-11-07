
def expliciteuler(f, oldframe, h, oldoldframe=None):
	v = oldframe.vector
	deltav =  f(v)
	newv = v + (deltav*h)
	return  newv
	
def leapfrog(f, oldframe, h, oldoldframe):
	if oldoldframe == None:
		raise Exception("No n-2 frame found!")
	else:
		v = oldframe.vector
		oldv = oldoldframe.vector
		deltav =  f(v)
		newv = oldv + (2*deltav*h)
		return newv
	
def rk4(f, oldframe, h, oldoldframe=None):
	v = oldframe.vector
	k1 = f(v)
	k2 = f(v + k1*h/2.)
	k3 = f(v + k2*h/2.)
	k4 = f(v + k3*h)
	newv = v + (h*(k1 + 2*(k2+k3) + k4)/6.)
	return  newv
	
def impliciteuler(f, oldframe, h, oldoldframe=None):
	v = oldframe.vector
	guessv = (h*f(v)) + v
	deltav = f(guessv)
	newv = v + (h*deltav)
	return newv
