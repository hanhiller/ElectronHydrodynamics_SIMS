import numpy as np
import relativisticTransforms as rT
import matplotlib.path as mpltPath
from myGeometryFunctions import *	
from driftDiffusionSimulatorBase import driftDiffusionSimulatorBase

class driftDiffusionSimulatorFullCircle(driftDiffusionSimulatorBase):
	'''
	driftDiffusionSimulatorFullCircle simulates hydrodynamic flow of dirac electrons
	through a constriction and into a circular device. The enitre arc of the circle
	is the drain and a source injects through a constriction at the arc's radial center
	'''
	def __init__(self):
	
		self.constrictionWidth = 0.5
		self.sourceDrainRatio = 1.
		self.injectorWidth = .6
		self.injectorHeight = 1.5
		self.fieldResolution = 0.1
		self.diameter = 3.
	 
		self.setConstrictionWidth(1.)
		
		super().__init__()
		
		self.diameter=self.boxL
		self.setDiameter(self.boxL)
		self.setFieldResolution(self.fieldResolution)
        
	def updateBody(self):
		R = self.diameter/2.
		xc = self.constrictionWidth/2.
		xi = self.injectorWidth/2.
		yi = self.injectorHeight
		yw = R
		
		thetas = np.linspace(1.39*np.pi,-.39*np.pi,25)
		
		self.borderX = np.cos(thetas)*R
		self.borderY = np.sin(thetas)*R
		self.borderX = np.append(self.borderX,
									np.array([xc,xi,-xi,-xc,self.borderX[0]]))
		self.borderY = np.append(self.borderY,
									np.array([0,-yi,-yi,0,self.borderY[0]]))
			
		s,d,m,r,f = 2,1,0,-1,-2 #'source','drain','mirror','rough'
		
		self.edgeStyle = [d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,m,m,s,m,m]	
		self.boxRange = [[-R, R], [-R, R]]
        
		self.setUpEdges()
		self.setFieldResolution(self.fieldResolution)
		
	def diffusiveWalls(self):
		s,d,m,r = 2,1,0,-1 #'source','drain','mirror','rough'
		
		self.edgeStyle = [d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,r,r,s,r,r]
		self.setUpEdges()
		
	def mirrorWalls(self):
		s,d,m,r = 2,1,0,-1 #'source','drain','mirror','rough'
		
		self.edgeStyle = [d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,m,m,s,m,m]
		self.setUpEdges()

	def setInjectorShape(self,width_in, height_in):
		self.injectorWidth = width_in
		self.injectorHeight = height_in
		self.updateBody()
	
	def setConstrictionWidth(self, width_in):
		self.constrictionWidth = width_in
		self.updateBody()
		
		
	def setDiameter(self, D):
		self.boxL = D
		self.diameter = D
		self.updateBody()
		self.addTip()
		
		
	def setFieldResolution(self,dX):
		self.fieldResolution = dX
		Nr = int(np.round(self.diameter/dX))
		
		self.Nx = Nr
		self.Ny = int(np.ceil(Nr/2*(1+self.injectorHeight/self.diameter*2.)))
		
		self.rho = np.zeros((self.Nx,self.Ny))
		self.Px = np.zeros((self.Nx,self.Ny))
		self.Py = np.zeros((self.Nx,self.Ny))
		self.Erho = np.zeros((self.Nx,self.Ny))
		
		_,self.histX,self.histY = np.histogram2d(np.zeros(10),np.zeros(10),
				bins = [self.Nx,self.Ny],range = self.boxRange)
		
