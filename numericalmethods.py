

def euler(d2ydx2,yi, olddydx, h):
	yf = yi + (olddydx*h)
	dydx =  olddydx + (h*d2ydx2(yi))
	return  yf, dydx
	
def leapfrog(d2ydx2, yi, olddydx,  h):
	yf = yi + (olddydx*h) + (0.5 * d2ydx2(yi)*(h**2))
	a0 = d2ydx2(yi)
	a1 = d2ydx2(yf)
	dydx = olddydx + 0.5*(a0+a1)*h
	return  yf, dydx
