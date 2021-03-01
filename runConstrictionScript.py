import sys
import configparser as cp
import numpy as np
import scipy.io
import os
import time
import multiprocessing
from driftDiffusionSimulatorConstriction import driftDiffusionSimulatorConstriction
from myGeometryFunctions import *


def main():
	
	initFile = sys.argv[1]
	
	config = cp.ConfigParser(interpolation=cp.ExtendedInterpolation())
	config.read(initFile)

	pScatter = config['DEFAULT'].getfloat('Scattering Probability', 0)
	Emin = config['DEFAULT'].getfloat('Minimum Energy', 1)
	sourceDrainRatio = config['DEFAULT'].getfloat('Source Drain Ratio', 1.2)
	Npart = config['DEFAULT'].getint('Number of Particles', 100)
	outPath = config['DEFAULT'].get('base output path', './')
	Nsteps = config['DEFAULT'].getint('time steps', 10000)
	dNsave = config['DEFAULT'].getint('save interval', 1000)
	Ncpu = config['DEFAULT'].getint('Number of CPUs', 1)
	boxL = config['DEFAULT'].getfloat('bounding box', 3)
	fieldRes = config['DEFAULT'].getfloat('histogram resolution',0.1)
	
	constrWidth,constrThick = [float(s) for s in config['DEFAULT'].get(
					'Constriction Width and Thickness','0.3,.05').split(',')]

	diffusiveEdges = config['DEFAULT'].getboolean('diffusive edges?',False)
	
	fnameBase = initFile.split('.')[0]
	
	if not os.path.isdir(outPath+fnameBase):
		os.mkdir(outPath+fnameBase)
	
	dSims = []
	for iterationNum in range(Ncpu):
		
	
		dSim = driftDiffusionSimulatorConstriction()
		dSim.setConstrictionWidth(constrWidth)
		dSim.setWallThickness(constrThick)
		dSim.setBoundaryLength(boxL)
		dSim.updateBody()
		dSim.setFieldResolution(fieldRes)
		dSim.updateScatterProb(pScatter)
		dSim.Emin=Emin
		dSim.setSourceDrainRatio(sourceDrainRatio)
		dSim.updateNparticles_NoOverlap(int(Npart))
		if diffusiveEdges:
			dSim.diffusiveWalls()		

		for i in range(int(Npart/2)):
			thetas = np.random.rand()*2.*np.pi
			dSim.vX[i] = np.cos(thetas)
			dSim.vY[i] = np.sin(thetas)
			dSim.vX[i+int(Npart/2)] = -np.cos(thetas)
			dSim.vY[i+int(Npart/2)] = -np.sin(thetas)

		#tic = time.time()
		dSims.append(dSim)

	jobs=[]
	for iterationNum,dSim in enumerate(dSims):
		fnameOut = outPath+fnameBase+"/"+fnameBase+("_%03d"%iterationNum)+".npz"	
	
		p = multiprocessing.Process(target=dSim.runAndSave, args=(Nsteps,dNsave,fnameOut))
		jobs.append(p)
		p.start()
	p.join()

		


if __name__== "__main__":
  main()