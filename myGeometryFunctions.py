import numpy as np

	
def seg_intersect(x1,y1,x2,y2,x3,y3,x4,y4):
	
	x13 = x1 - x3
	y13 = y1 - y3
	x21 = x2 - x1
	y21 = y2 - y1
	x43 = x4 - x3
	y43 = y4 - y3
	
	ts = (x13*y21-x21*y13)/(x43*y21-x21*y43)
	us = (x13*y43-x43*y13)/(x43*y21-x21*y43)
	
	return ((ts>=0) and (ts<=1) and (us>=0) and (us<=1))
	
	
def segCrossPoly(polyX,polyY,p1,p2):
	crossed = []
	NPoints = len(polyX)-1
	
	for i in range(NPoints):
		if seg_intersect(polyX[i],polyY[i],polyX[i+1],polyY[i+1],p1[0],p1[1],p2[0],p2[1]):
			crossed.append(i)
		
	return crossed

def polyLengths(polyX,polyY):
	NPoints = len(polyX)-1
	Lengths = np.zeros(NPoints)
	
	for i in range(NPoints):
		dX = polyX[i+1]-polyX[i]
		dY = polyY[i+1]-polyY[i]
		
		Lengths[i] = np.sqrt(np.sum(dX**2+dY**2))
		
	return Lengths
	
def polyNorms(polyX,polyY):
	NPoints = len(polyX)-1
	normX = np.zeros(NPoints)
	normY = np.zeros(NPoints)
	
	for i in range(NPoints):
		dX = polyX[i+1]-polyX[i]
		dY = polyY[i+1]-polyY[i]
		L = np.sqrt(np.sum(dX**2+dY**2))
		
		normX[i] = -dY/L
		normY[i] = dX/L
		
	return normX,normY

def mirrorNorm(normX,normY,vX_in,vY_in):
	
	vecProj = normX*vX_in+normY*vY_in
	vX_out = vX_in - 2 * normX*vecProj
	vY_out = vY_in - 2 * normY*vecProj
	
	return vX_out, vY_out
	