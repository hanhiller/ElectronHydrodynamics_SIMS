import sys
import configparser as cp
import numpy as np
import scipy.io
import os
import time
import multiprocessing
from driftDiffusionSimulatorDome import driftDiffusionSimulatorDome
from myGeometryFunctions import *

'''
runDomeScript.py reads in a script text file and runs a number of simulations of
a "driftDiffusionSimulatorDome" simulation object. It enables starting from a 
fresh simulation, or taking the results of a previous simulation as initial cond-
itions. 
'''

def main():
	
	initFile = sys.argv[1]
	
	config = cp.ConfigParser(interpolation=cp.ExtendedInterpolation())
	config.read(initFile)

	constrictionWidth = config['DEFAULT'].getfloat('Constriction Width', 0.3)
	pScatter = config['DEFAULT'].getfloat('Scattering Probability', 0)
	Temp = config['DEFAULT'].getfloat('Temperature', 100)
	sourceDrainRatio = config['DEFAULT'].getfloat('Source Drain Ratio', 1.2)
	partDensity = config['DEFAULT'].getint('Particle Density', 100)
	Nsteps = config['DEFAULT'].getint('time steps', 10000)
	dNsave = config['DEFAULT'].getint('save interval', 1000)
	domeD = config['DEFAULT'].getfloat('dome diameter', 3)
	collisionDist = config ['DEFAULT'].getfloat('collision overlap', 0.05)
	stepSize = config['DEFAULT'].getfloat('step size', 0.01)
	injectH,injectW = [float(s) for s in config['DEFAULT'].get(
					'injector height and width','0.5,1').split(',')]
    
	probeCenterX = config['DEFAULT'].getfloat('probeCenterX', 0)
	probeCenterY = config['DEFAULT'].getfloat('probeCenterY', 0)

	Ncpu = config['DEFAULT'].getint('Number of CPUs', 1)
	outPath = config['DEFAULT'].get('base output path', './SIM_data')
	
	diffusiveEdges = config['DEFAULT'].getboolean('diffusive edges?',False)
	
	initCondFile = config['DEFAULT'].get('inital conditions','')
	
	splitName = initFile.split('.')
	if len(splitName) == 2:
		fnameBase = initFile .split('.')[0]
	elif len(splitName) == 3:
		fnameBase = '.'.join([initFile .split('.')[0], initFile .split('.')[1]])
	elif len(splitName) == 4:
		fnameBase = '.'.join([initFile .split('.')[0], initFile .split('.')[1], initFile .split('.')[2]])
	
	# fnameBaseREDUCED = fnameBase+'_REDUCED'
    
	if not os.path.isdir(outPath+fnameBase):
		os.mkdir(outPath+fnameBase)
		# os.mkdir(outPath+fnameBaseREDUCED)

	print("Simulation temperature:", Temp)

	dSims = []
	for iterationNum in range(Ncpu):
		
		#set up the simulation from script file
		dSim = driftDiffusionSimulatorDome()
		dSim.probeCenterX = probeCenterX
		dSim.probeCenterY = probeCenterY
		dSim.Temp = Temp
		dSim.partDensity = partDensity
		dSim.setEmin(dSim.partDensity, dSim.Temp)
		dSim.setConstrictionWidth(constrictionWidth)
		dSim.setInjectorShape(injectW,injectH)
		dSim.setDiameter(domeD)
		dSim.updateBody()
		dSim.addTip()
		dSim.calcArea()
		dSim.setNpart(dSim.partDensity, dSim.Area)
		dSim.updateNparticles()
		dSim.updateScatterProb(pScatter)
		dSim.setSourceDrainRatio(sourceDrainRatio)
		dSim.setOverlapRadius(collisionDist)
		dSim.setStepSize(stepSize)
		dSim.updateNparticles()
		if diffusiveEdges:
			dSim.diffusiveWalls()		

		for i in range(int(dSim.Npart/2)):
			thetas = np.random.rand()*2.*np.pi
			dSim.vX[i] = np.cos(thetas)
			dSim.vY[i] = np.sin(thetas)
			dSim.vX[i+int(dSim.Npart/2)] = -np.cos(thetas)
			dSim.vY[i+int(dSim.Npart/2)] = -np.sin(thetas)
		
		if initCondFile:
			if os.path.isfile(outPath+initCondFile+'/'+initCondFile+("_%03d"%iterationNum)+".npz"):
				mat = np.load(outPath+initCondFile+'/'+initCondFile+("_%03d"%iterationNum)+".npz")
				dSim.Xpos = mat['Xpos']
				dSim.Ypos = mat['Ypos']
				dSim.vX = mat['vX']
				dSim.vY = mat['vY']
				dSim.pR = mat['pR']
			
		#tic = time.time()
		dSims.append(dSim)

	jobs=[]
	for iterationNum,dSim in enumerate(dSims):
		fnameOut = outPath+fnameBase+"/"+fnameBase+("_%03d"%iterationNum)+".npz"
		#fnameOutREDUCED = outPath+fnameBaseREDUCED+"/"+fnameBaseREDUCED+("_%03d"%iterationNum)+".npz"
	
		p = multiprocessing.Process(target=dSim.runAndSave, args=(Nsteps,dNsave,fnameOut)#,fnameOutREDUCED))
		jobs.append(p)
		p.start()
	p.join()

	


if __name__== "__main__":
  main()