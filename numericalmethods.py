

def expliciteuler(df1,yi, olddydx, h):
	yf = yi + (olddydx*h)
	dydx =  olddydx + df1(yi)
	return  yf, dydx
	
def leapfrog(df2, yi, olddydx,  h):
	yf = yi + (olddydx*h) + (0.5 * df2(yi)*(h**2))
	a0 = df2(yi)
	a1 = df2(yf)
	dydx = olddydx + 0.5*(a0+a1)*h
	return  yf, dydx
	
def rk4(df1, yi, olddydx,  h):
	k1 = df1(yi)
	k2 = df1(yi + k1*h/2)
	k3 = df1(yi + k2*h/2)
	k4 = df1(yi + k3*h)
	dydx = (k1 + 2*(k2+k3) + k4)/6
	yf = yi + (h*dydx)
	print yi, yf, k1, k2, k3, k4, dydx
	return  yf, dydx
	
def impliciteuler(df1,yi, olddydx, h):
	guessdydx =  olddydx + df1(yi)
	guessy = yi + (guessdydx*h)
	dydx = olddydx + df1(guessy)
	yf = yi + (dydx*h)
	return  yf, dydx
	

