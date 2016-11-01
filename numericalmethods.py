
def expliciteuler(df1, yi, olddydx, h, df2=0, oldoldy=0, oldolddydx=0):
	yf = yi + (olddydx*h)
	dydx =  df1(yi, olddydx)
	return  yf, dydx
	
def leapfrog(df2, yi, olddydx,  h, df1=0, oldoldy=0, oldolddydx=0):
	#~ yf = yi + (olddydx*h) + (0.5 * df2(yi)*(h**2))
	#~ a0 = df2(yi)
	#~ a1 = df2(yf)
	#~ dydx = olddydx + 0.5*(a0+a1)*h
	
	yf = oldoldy + (2*olddydx*h)
	dydx = oldolddydx + (2*df2(yi, olddydx=0)*h)
	
	return  yf, dydx
	
def rk4(df1, yi, olddydx,  h, df2=0, oldoldy=0, oldolddydx=0):
	k1 = df1(yi, olddydx)
	k2 = df1(yi + k1*h/2., olddydx)
	k3 = df1(yi + k2*h/2., olddydx)
	k4 = df1(yi + k3*h, olddydx)
	dydx = (k1 + 2*(k2+k3) + k4)/6.
	yf = yi + (h*dydx)
	return  yf, dydx
	
def impliciteuler(df1, yi, olddydx, h, df2=0, oldoldy=0, oldolddydx=0):
	guessdydx =  df1(yi, olddydx)
	guessy = yi + (guessdydx*h)
	dydx = df1(guessy, olddydx)
	yf = yi + (dydx*h)
	return  yf, dydx
	

