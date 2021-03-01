import numpy as np
import relativisticTransforms as rT
import matplotlib.path as mpltPath
from myGeometryFunctions import *	
from driftDiffusionSimulatorBase import driftDiffusionSimulatorBase

class driftDiffusionSimulatorDome(driftDiffusionSimulatorBase):
	'''
	driftDiffusionSimulatorDome simulates hydrodynamic flow of dirac electrons
	through a constriction and into a semi-circular "dome". The arc of the dome
	is the drain and a source injects through the constriction at the radial 
	center of the arc.
	'''
	def __init__(self):
	
		self.constrictionWidth = 1.
		self.sourceDrainRatio = 1.
		self.injectorWidth = 1.
		self.injectorHeight = 0.5
		self.fieldResolution = 0.1
		self.domeD = 3.
	
		self.setConstrictionWidth(1.)
		
		super().__init__()
		
		self.domeD=self.boxL
		self.setDiameter(self.boxL)
		self.setFieldResolution(self.fieldResolution)
        
	def updateBody(self):
		rDome = self.domeD/2.
		xc = self.constrictionWidth/2.
		xi = self.injectorWidth/2.
		yi = self.injectorHeight
		yw = self.domeD/2.
		
		thetas = np.linspace(np.pi,0.,15)
		
		self.borderX = np.cos(thetas)*self.domeD/2.
		self.borderY = np.sin(thetas)*self.domeD/2.
		self.borderX = np.append(self.borderX,
									np.array([xc,xi,-xi,-xc,self.borderX[0]]))
		self.borderY = np.append(self.borderY,
									np.array([0,-yi,-yi,0,self.borderY[0]]))
			
		s,d,m,r,f = 2,1,0,-1,-2 #'source','drain','mirror','rough'
		
		self.edgeStyle = [d,d,d,d,d,d,d,d,d,d,d,d,d,d,m,m,s,m,m]	
		self.boxRange = [[-self.domeD/2.,self.domeD/2.],
										[-self.injectorHeight,self.domeD/2.]]
        
		self.setUpEdges()
		self.setFieldResolution(self.fieldResolution)
		
	def diffusiveWalls(self):
		s,d,m,r = 2,1,0,-1 #'source','drain','mirror','rough'
		
		self.edgeStyle = [d,d,d,d,d,d,d,d,d,d,d,d,d,d,r,r,s,r,r]
		self.setUpEdges()
		
	def mirrorWalls(self):
		s,d,m,r = 2,1,0,-1 #'source','drain','mirror','rough'
		
		self.edgeStyle = [d,d,d,d,d,d,d,d,d,d,d,d,d,d,m,m,s,m,m]
		self.setUpEdges()

	def setInjectorShape(self,width_in, height_in):
		self.injectorWidth = width_in
		self.injectorHeight = height_in
		self.updateBody()
	
	def setConstrictionWidth(self, width_in):
		self.constrictionWidth = width_in
		self.updateBody()
		
		
	def setDiameter(self,D):
		self.boxL = D
		self.domeD = D
		self.updateBody()
		self.addTip()
		
		
	def setFieldResolution(self,dX):
		self.fieldResolution = dX
		Nr = int(np.round(self.domeD/dX))
		
		self.Nx = Nr
		self.Ny = int(np.ceil(Nr/2*(1+self.injectorHeight/self.domeD*2.)))
		
		self.rho = np.zeros((self.Nx,self.Ny))
		self.Px = np.zeros((self.Nx,self.Ny))
		self.Py = np.zeros((self.Nx,self.Ny))
		self.Erho = np.zeros((self.Nx,self.Ny))
		
		_,self.histX,self.histY = np.histogram2d(np.zeros(10),np.zeros(10),
				bins = [self.Nx,self.Ny],range = self.boxRange)
		
