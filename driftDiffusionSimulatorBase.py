import numpy as np
import relativisticTransforms as rT
from scipy.spatial.distance import pdist,squareform
import matplotlib.path as mpltPath
from myGeometryFunctions import *	

import time
	
class driftDiffusionSimulatorBase:
	'''
	driftDiffusionSimulatorBase is a simulation object that serves to simulate
	electron flow in 2D. Specifically, it is designed to simulate dirac-like 
	quasiparticles (conical dispersion) and allows them to scatter off of each-
	other while conserving momentum and energy. An effective temperature is 
	included by setting a minimum energy for available scattering that is small-
	er (slightly) than the average energy of particles (roughly the Fermi ener-
	gy). Region boundaries are defined with edge types either being a "source",
	"drain","mirror", or "rough surface". Net currents are sourced by preferen-
	tially sourcing particles from drain compared with source electrode. 
	'''

	def __init__(self):
		self.rho=np.array([])
		self.Temp=.9
		self.boxL= 3.
		self.partDensity = 100
		self.rOverlap = 0.05
		self.sourceDrainRatio = 1.
		self.p_scatter = 0.01
		self.probeTipDiameter = 1
		self.NumberDrainBins = 10
		self.probeCenterX = 0
		self.probeCenterY = 0
		
		self.setBoundaryLength(self.boxL)
		self.setOverlapRadius(self.rOverlap)
		self.updateBody()
		self.addTip()
		self.calcArea()
		self.setNpart(self.partDensity, self.Area)
		self.setEmin(self.partDensity, self.Temp)
		self.updateNparticles()
		self.findOverlaps()
		self.updateScatterProb(0.01)
		self.setStepSize(.01)
		self.setFieldResolution(.1)
		
		self.timeCount = 0
		self.initializationCount = 100000000
		
	def updateBody(self):
		#builds the border and defines the edge properties
		xw = self.boxL/2.
		yw = self.boxL/2.
		
        # left border
		self.borderX = np.array([-xw,-xw])
		self.borderY = np.array([-yw,yw])
        
		# create binned drain edge (top)
		BinLength = self.boxL/self.NumberDrainBins
		for bin in np.arange(1, self.NumberDrainBins):
 		    self.borderX = np.append(self.borderX, self.borderX[bin]+BinLength )
 		    self.borderY = np.append(self.borderY, yw)
        
        # right and bottom borders
		self.borderX = np.append(self.borderX, [[xw,xw,-xw]] )
		self.borderY = np.append(self.borderY, [[yw,-yw,-yw]] )
        
		s = 2#'source'
		d = 1#'drain'
		m = 0#'mirror'
		r =-1#'rough'
		f =-2#'false'
		
		self.edgeStyle = [m]
		for bin in np.arange(self.NumberDrainBins):
 		    self.edgeStyle += [d]
                          
		self.edgeStyle += [m, s]
		
		self.setUpEdges()
		

	def addTip(self):
        
		# add circular probe tip
        
		s = 2#'source'
		d = 1#'drain'
		m = 0#'mirror'
		r =-1#'rough'
		f =-2#'false'
        
		thetas = np.linspace(np.pi, 3*np.pi, 100)
		self.borderX = np.append( self.borderX, self.probeCenterX + np.cos(thetas)*self.probeTipDiameter/2)
		self.borderX = np.append( self.borderX, self.borderX[0])
		self.borderY = np.append( self.borderY, self.probeCenterY + np.sin(thetas)*self.probeTipDiameter/2)
		self.borderY = np.append( self.borderY, self.borderY[0])
		
		# add edgeStyles for edges that form the probe
		self.edgeStyle = np.append(self.edgeStyle, f)
		for angle in thetas:
		    if angle == 3*np.pi:  # don't add edge for last point on circle
		        break
		    self.edgeStyle = np.append(self.edgeStyle, m)

		self.edgeStyle = np.append(self.edgeStyle, f)  # add last fictisous edge
		self.setUpEdges()
        
	def loadBody(self,borderX_in,borderY_in,edgeStyle_in):
		#allows arbitrary boundry 
		self.borderX = borderX_in
		self.borderY = borderY_in
		self.edgeStyle = edgeStyle_in
		
		self.setUpEdges()
		
	def setUpEdges(self):
		#updates parameters based on the borders and edge properties
		self.edgeStyle = np.array(self.edgeStyle)
		Ncorners = len(self.edgeStyle)
		
		#calculate norms and edge lengths
		self.normX,self.normY = polyNorms(self.borderX,self.borderY)
		self.lengths = polyLengths(self.borderX,self.borderY)
		
		s = 2#'source'
		d = 1#'drain'
		m = 0#'mirror'
		r =-1#'rough'

		#find lengths of source and drain edges 
		self.sources = [i for i,x in enumerate(self.edgeStyle) if x == s]
		self.sourceLength = np.sum(self.lengths[self.sources])
		self.drains = [i for i,x in enumerate(self.edgeStyle) if x == d]
		self.drainLength = np.sum(self.lengths[self.drains])
		
		#calculate probability of injecting from a source edge
		self.sourceProbs = (self.lengths[self.sources]*self.sourceDrainRatio/
					(self.drainLength+self.sourceLength*self.sourceDrainRatio))
					
		#calculate probability of injecting from a drain edge
		self.drainProbs = (self.lengths[self.drains]/
					(self.drainLength+self.sourceLength*self.sourceDrainRatio))
		
		self.Probs = np.zeros(Ncorners)
		self.Probs[self.sources] = self.sourceProbs
		self.Probs[self.drains] = self.drainProbs

		self.ProbIdx = (np.where(self.Probs>0))[0]
		self.cumProb = np.cumsum(self.Probs[self.ProbIdx])
		
		#initialize Nabsorbed and Ninjected for counting statistics
		self.Nabsorbed = np.zeros(Ncorners)
		self.Ninjected = np.zeros(Ncorners)
		
		#sets border path for edge crossing functions
		self.borderPath = mpltPath.Path(
									np.vstack((self.borderX,self.borderY)).T)
	
	def calcArea(self):
        # calculates area of the device
		edgeCount = 0
		for borderEdge in self.edgeStyle:
			if borderEdge == -2:
				break
			edgeCount += 1
         
		Area = 0
		for i in range(edgeCount):
			Area += abs(self.borderX[i]*self.borderY[i+1] - self.borderX[i+1]*self.borderY[i])
		Area = (Area/2 - np.pi*(self.probeTipDiameter/2)**2)
		self.Area = Area
		
	def setNpart(self, partDensity, Area):
        # ptc. density and area set the number of particles
		self.Npart = int(partDensity*Area)
    
	def setEmin(self, partDensity, Temp):
        # min scattering energy as a function of temp and density
		self.Emin = 1 - (3.41475*10**-3)*( (Temp**1.5)/partDensity)**1.25
        
		assert self.Emin > 0.8, 'Non-physcial minimum scattering energy; raise particle density for this temperature'
    
	def setSourceDrainRatio(self,Ratio_in):
		#sets relative probability of injection from the source and drain edges
		self.sourceDrainRatio = Ratio_in
		self.updateBody()
		self.addTip()
		
	def injectPosition(self):
		#calculates the position to inject particle based on set probabilities
		injectIdx = self.ProbIdx[np.where(self.cumProb>np.random.rand())[0][0]]
		injectFraction = np.random.rand()
		x1 = self.borderX[injectIdx+1]
		x0 = self.borderX[injectIdx]
		y1 = self.borderY[injectIdx+1]
		y0 = self.borderY[injectIdx]
		
		injectX = x0+injectFraction*(x1-x0)
		injectY = y0+injectFraction*(y1-y0)
		
		return injectX,injectY,injectIdx

	def injectFromContact(self,contactIdx):
		#calculates position to inject along a specific contact
		injectFraction = np.random.rand()
		x1 = self.borderX[contactIdx+1]
		x0 = self.borderX[contactIdx]
		y1 = self.borderY[contactIdx+1]
		y0 = self.borderY[contactIdx]
		
		injectX = x0+injectFraction*(x1-x0)
		injectY = y0+injectFraction*(y1-y0)
		
		return injectX,injectY

				
	def findNewInject(self):
		#searches for a new place to inject so there is no overlap with existing
		#particles
		injectX,injectY,injectIdx = self.injectPosition()
		while self.getRegionOccupied(injectX,injectY):
			injectX,injectY = self.injectFromContact(injectIdx)
		return injectX,injectY,injectIdx
		
	def setOverlapRadius(self, R_in):
		#sets radius of particle interactions
		self.rOverlap = R_in
		
	def updateNparticles(self):
		#updates the number of simultaneous particles in the simulation
		Npart_in = self.Npart
		
		# randomly assign position
		self.Xpos = self.boxL*(np.random.rand(Npart_in)-0.5)
		self.Ypos = self.boxL*(np.random.rand(Npart_in)-0.5)
		
		# find which particles are out of bounds
		outOfBounds = np.invert(self.borderPath.contains_points(
											np.vstack((self.Xpos,self.Ypos)).T))
		
		#keep randomly searching until all particles are in bounds
		while any(outOfBounds):
			for idx in np.where(outOfBounds)[0].tolist():
				self.Xpos[idx] = self.boxL*(np.random.rand()-0.5)
				self.Ypos[idx] = self.boxL*(np.random.rand()-0.5)
			outOfBounds = np.invert(self.borderPath.contains_points(
											np.vstack((self.Xpos,self.Ypos)).T))
		self.findOverlaps()

		#assign random velocity direction
		thetas = np.random.rand(Npart_in)*2.*np.pi
		self.vX = np.cos(thetas)
		self.vY = np.sin(thetas)
		#assign fixed momentum amplitude
		self.pR = np.ones(np.shape(thetas))
		
		#intializes particular statistics 
		self.L_NoScatter = np.zeros(Npart_in)
		self.traces = [[] for i in range(Npart_in)]
		
		#set up indices for quick retreival 
		self.i_lookup = np.zeros(int(self.Npart*(self.Npart-1)/2),dtype = int)
		self.j_lookup = np.zeros(int(self.Npart*(self.Npart-1)/2),dtype = int)
		_idx = 0
		for i in range(self.Npart):
			for j in range(i+1,self.Npart):
				self.i_lookup[_idx] = i
				self.j_lookup[_idx] = j
				_idx += 1

		
		
	def findOverlaps(self):
		#finds particle pairs that are within a the rOverlap radius. Successive-
		#ly checks the separation in X and then Y before calculating the cartes-
		#ian distance to save computational resource
		
		self.overlaps = squareform(self.checkOverlaps())
	
	def checkOverlaps(self):
		r = np.vstack((self.Xpos,self.Ypos)).T
		return pdist(r)< self.rOverlap
	

	def getRegionOccupied(self,X_in,Y_in):
		#checks if a given X_in and Y_in are in an already occupied region
		occupied = False
		for i in range(self.Npart):
			absDx = np.abs(self.Xpos[i]-X_in)
			if absDx < self.rOverlap:
				absDy = np.abs(self.Ypos[i]-Y_in)
				if absDy < self.rOverlap:
					if absDx**2+absDy**2 < self.rOverlap**2:
						occupied = True
		return occupied
		
	
	def randomPointInBounds(self):
		#finds a single coordinate-pair that is inbounds
		x,y = self.boxL*(np.random.rand()-0.5),self.boxL*(np.random.rand()-0.5)
		inBounds = self.borderPath.contains_points(np.vstack((x,y)).T)
		while not inBounds:
			x = self.boxL*(np.random.rand()-0.5)
			y = self.boxL*(np.random.rand()-0.5)
			inBounds = self.borderPath.contains_points(np.vstack((x,y)).T)
		return x,y
		
		
	def updateScatterProb(self,p_in):
		#sets the probability of scattering in the bulk
		self.p_scatter = p_in
		
	def setStepSize(self,dX):
		#sets the distance traversed in a single timestep
		self.DeltaX = dX
	
	def setBoundaryLength(self,L):
		#sets the overall boundary of the system
		self.boxL = L
		self.boxRange = [[-L/2.,L/2.],[-L/2.,L/2.]]
	
	def setFieldResolution(self,dX):
		#sets several spatial histograms with box size "dX"
		Nr = int(np.round((self.boxRange[0][1]-self.boxRange[0][0])/dX))
		
		self.Nx = Nr
		self.Ny = int(np.round((self.boxRange[1][1]-self.boxRange[1][0])/dX))

		self.rho = np.zeros((self.Nx,self.Ny))  # particle density
		self.Px = np.zeros((self.Nx,self.Ny))   # X velocity density
		self.Py = np.zeros((self.Nx,self.Ny))   # Y velocity density
		self.Erho = np.zeros((self.Nx,self.Ny)) # Energy Density
		
		#sets the X,Y mesh
		_,self.histX,self.histY = np.histogram2d(
							self.Xpos,self.Ypos,bins=[self.Nx,self.Ny],
							range = self.boxRange,weights = self.vX*self.pR)
		
	def randCos(self):
		#returns an angle between -pi/2 and pi/2 following a cosine distribution
		theta = np.pi*(np.random.rand()-.5)
		while np.random.rand()>np.cos(theta):
			theta = np.pi*(np.random.rand()-.5)
		return theta
		
	def randCosNorm(self,normX,normY):
		#returns a vector selected from a cosine distribution relative to an
		#input normalvector.
		theta = self.randCos()
		vXNew = -normX*np.cos(theta)-normY*np.sin(theta)
		vYNew = -normY*np.cos(theta)+normX*np.sin(theta)
		return vXNew, vYNew

	
	def timestep(self):
		# executes a single timestep of the simulation
		
		#propagate all particles in the direction of their velocity vector
		self.Xpos += self.vX*self.DeltaX
		self.Ypos += self.vY*self.DeltaX
		#update mean free path statistics
		self.L_NoScatter+=self.DeltaX
		dt1 = 0
		dt2 = 0
		dt3 = 0
		dt4 = 0
		
		#find which particles ended up out of bounds
		outOfBounds = np.invert(
			self.borderPath.contains_points(np.vstack((self.Xpos,self.Ypos)).T))
		_t2 = time.time()

		for i in range(self.Npart):
			
			if outOfBounds[i]:
				# backtracks if out of bounds
				X0 = self.Xpos[i] - self.vX[i]*self.DeltaX
				Y0 = self.Ypos[i] - self.vY[i]*self.DeltaX

				while not self.borderPath.contains_point((X0,Y0)):
					theta = np.random.rand()*2*np.pi
					X0 = self.Xpos[i] - np.sin(theta)*self.DeltaX
					Y0 = self.Ypos[i] - np.cos(theta)*self.DeltaX
				
				#finds index of the line(s) crossed
				lineCrossed = segCrossPoly(self.borderX,self.borderY,
								np.array((X0,Y0)),
								np.array((self.Xpos[i],self.Ypos[i])))
								
				self.Xpos[i],self.Ypos[i] = X0, Y0
				
				lenCrossed = len(lineCrossed)
				if lenCrossed == 0:  # no segments crossed; exactly betweeen two
					self.vX[i] *= -1
					self.vY[i] *= -1
					#	print('error!')
				
				elif lenCrossed == 1:  # 1 segment crossed; the expected case
					#select the first in the list of lines crossed (convenience)
					idxCross = lineCrossed[0]
					if self.edgeStyle[idxCross] == -1:  #diffusive edge
						self.scatterFromDiffusiveEdge(i,idxCross)
						
					if self.edgeStyle[idxCross] == 0:	#mirror edge
						self.reflectFromMirrorEdge(i,idxCross)
						
					if self.edgeStyle[idxCross] > 0:		#source or drain
						#re-inject particle and take statistics
						self.consumeAndReinject(i,idxCross)
				
				else:
                    
                    
					r1 = np.where(self.edgeStyle[lineCrossed]>0)[0]
					if len(r1)>0: #if crossed a source/drain
						#consume by random selection among the source and drain
						#edges in the lineCrossed list
						self.consumeAndReinject(
											i,r1[np.random.randint(len(r1))])
					else:
						#scatter from an edge whose norm has a component in the
						#direction of the velocity
						r2 = np.random.randint(lenCrossed)
						idxCross = lineCrossed[r2]  # randomly select which edge is actually crossed
						while (self.vX[i]*self.normX[idxCross] +
										self.vY[i]*self.normY[idxCross] < 0):
							r2 = np.random.randint(lenCrossed)
							idxCross = lineCrossed[r2]
						
						if self.edgeStyle[idxCross] == -1:  #diffusive edge
							self.scatterFromDiffusiveEdge(i,idxCross)
						
						if self.edgeStyle[idxCross] == 0:	#mirror edge
							self.reflectFromMirrorEdge(i,idxCross)
					
						
			elif self.p_scatter > np.random.rand(1):
				#randomize direction in bulk
				theta = np.random.rand(1)*2.*np.pi
	
				self.vX[i] = np.cos(theta)
				self.vY[i] = np.sin(theta)
		
		
		_overlaps = self.checkOverlaps()
		
		_idx_i = self.i_lookup[_overlaps]
		_idx_j = self.j_lookup[_overlaps]
		
		for i,j in zip(_idx_i,_idx_j):
			if not self.overlaps[i,j]:
				#perform relativistic scattering
				p1 = np.array([self.vX[i]*self.pR[i],
										self.vY[i]*self.pR[i]])
				p2 = np.array([self.vX[j]*self.pR[j],
										self.vY[j]*self.pR[j]])
				#calculate outgoing momenta
				p3,p4 = rT.scatterMasslessParticles(
												p1,p2,self.Emin)
				if not (np.array_equal(p1,p3) and 
										np.array_equal(p2,p4)):
				
					#measure momentum amplitude
					self.pR[i] = np.sqrt(np.sum(p3**2))
					self.pR[j] = np.sqrt(np.sum(p4**2))
					#set velocities
					self.vX[i],self.vY[i] = p3/self.pR[i]
					self.vX[j],self.vY[j] = p4/self.pR[j]
					
					self.L_NoScatter[i]=0
					self.L_NoScatter[j]=0
				
			self.overlaps[i,j] = True
			
				
		#update global time index
		self.timeCount+=1

		if self.timeCount == self.initializationCount:
			#reset statistics after initalizationCount number of timesteps
			self.Px*=0
			self.Py*=0
			self.rho*=0
			self.Erho*=0
			self.Ninjected*=0
			self.Nabsorbed*=0

		if np.mod(self.timeCount,5)==0:
			#update histograms every 5 timesteps
			h,_,_ = np.histogram2d(self.Xpos,self.Ypos,bins=[self.Nx,self.Ny],
										range = self.boxRange,weights = self.vX)
			self.Px+=h
			
			h,_,_ = np.histogram2d(self.Xpos,self.Ypos,bins=[self.Nx,self.Ny],
										range = self.boxRange,weights = self.vY)
			self.Py+=h
			
			h,_,_ = np.histogram2d(self.Xpos,self.Ypos,bins=[self.Nx,self.Ny],
										range = self.boxRange,weights = self.pR)
			self.Erho+=h
			
			h,_,_ = np.histogram2d(self.Xpos,self.Ypos,bins=[self.Nx,self.Ny],
										range = self.boxRange)
			self.rho+=h
		
				
		
	def scatterFromDiffusiveEdge(self,partIdx,edgeIdx):
		self.vX[partIdx],self.vY[partIdx] = self.randCosNorm(
										self.normX[edgeIdx],self.normY[edgeIdx])
						
	def reflectFromMirrorEdge(self,partIdx,edgeIdx):
		self.vX[partIdx],self.vY[partIdx] = mirrorNorm(self.normX[edgeIdx],
						self.normY[edgeIdx],self.vX[partIdx],self.vY[partIdx])
						
	def consumeAndReinject(self,partIdx,edgeIdx):
		xNew,yNew,idxNew = self.injectPosition()
		self.Nabsorbed[edgeIdx]+=1		
		self.Ninjected[idxNew]+=1
		
		vXNew,vYNew = self.randCosNorm(
						self.normX[idxNew],self.normY[idxNew])
		
		self.vX[partIdx],self.vY[partIdx] = vXNew,vYNew
		#small offset to make sure particle is in bounds
		self.Xpos[partIdx] = xNew+0*vXNew*self.DeltaX*.0001
		self.Ypos[partIdx] = yNew+0*vYNew*self.DeltaX*.0001
	
	def saveState(self,fname):#, fnameREDUCED):
		#saves all object data into an npz file
		varNames = list(self.__dict__.keys())
		varNames.remove('i_lookup')
		varNames.remove('j_lookup')
		varNames.remove('overlaps')
		self.varNames = varNames
		
		saveFunction = "np.savez(fname, "
		for name in varNames:
			saveFunction = saveFunction+name+' = self.'+name+', '
			
		saveFunction = saveFunction[:-2] + ')'
		exec(saveFunction)
        
		# redcuedVarNames = [self.Erho, self.Px, self.Py, self.Nabsorbed, self.Ninjected,
						# self.histX, self.histY, self.borderX, self.borderY, self.edgeStyle]
		# reducedSaveFunction = "np.savez(fnameREDUCED, "
		# for name in reducedVarNames:
			# reducedSaveFunction = reducedSaveFunction+name+' = self.'+name+', '
		# reducedSaveFunction = reducedSaveFunction[:-2] + ')'
		# exec(reducedSaveFunction)
		
		
	def runAndSave(self,Nsteps,dNsave,fname):#,fnameREDUCED):
		
		t0 = time.time()
		for i in range(Nsteps):
			
			self.timestep()
			
			if i%dNsave == 0:
				t1 = time.time()
				print('%d:  %g'%(i,(t1-t0)))
				t0 = t1
			
				self.saveState(fname)#,fnameREDUCED)
