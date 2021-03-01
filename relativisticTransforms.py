import numpy as np



def boost(pBar,E,vBar):
		'''
		This function takes a particle with vector momentum pBar and energy E and performs a lorentz boost with vector velocity
		vBar. It returns a new momentum pout and energy Eout. Velocities are in units of c
		'''
		
		vSqr = np.sum(vBar**2)
		gamma = (1.-vSqr)**-0.5
		pDotV = np.sum(pBar*vBar);
		
		if vSqr>0:   
			#perform boost
			Eout = gamma*(E-pDotV)
			pout = pBar+(gamma - 1.)*pDotV*vBar/vSqr - gamma*E*vBar
		else:
			#do this if there is no boost
			Eout = E
			pout = pBar
		return pout,Eout
		
		
		
def findCollinearBoost(pBar1,pBar2):
	'''
	This function computes the velocity vector that one should perform a lorentz boost along to get the two incoming momentum
	vectors to be collinear. This current version assumes massless particles: velocity is c along the momentum direction.
	'''
	#calculate momentum directions
	theta1 = np.arctan2(pBar1[1],pBar1[0])
	theta2 = np.arctan2(pBar2[1],pBar2[0])
	
	# find direction along which both vectors have the same projected amplitude. 
	thetaV = np.arctan2((np.cos(theta1)-np.cos(theta2)),(np.sin(theta2)-np.sin(theta1)))
	
	# calculate amplitude of the projected amplitude
	vAmp = np.cos(theta1)*np.cos(thetaV)+np.sin(theta1)*np.sin(thetaV)
	
	# output velocity vector to boost along
	vOut = np.array([vAmp*np.cos(thetaV),vAmp*np.sin(thetaV)])
	
	return vOut

def findCenterOfMassBoost(pBar1,E1,pBar2,E2):

	'''
	This function takes two particles momenta and energies and boosts them into a center of mass frame. It provides the output
	momenta/energies as well as the two sequental boost velocities to get to the COM frame. the velocity v1 defines the boost
	required to get the vectors to be collinear and the velocity v2 is the boost to get both particles to have the same energy.
	This currently assumes massless particles, since the COM is more generally defined as the case of equal and opposite momenta. 
	'''
	#boost to bring momenta to be collinear
	vBoost1 = findCollinearBoost(pBar1,pBar2)
	pBar1mid,E1mid = boost(pBar1,E1,vBoost1)
	pBar2mid,E2mid = boost(pBar2,E2,vBoost1)
	
	#calculate boost to get equal energies (and thus momenta for massless particles)
	p1 = np.sqrt(np.sum(pBar1mid**2))
	p2 = np.sqrt(np.sum(pBar2mid**2))
	vBoost2 = pBar1mid/p1*(E1mid-E2mid)/(p1+p2)
	#perform the second boost
	pBar1out,E1out = boost(pBar1mid,E1mid,vBoost2)
	pBar2out,E2out = boost(pBar2mid,E2mid,vBoost2)
	
	return pBar1out,E1out,pBar2out,E2out,vBoost1,vBoost2



def scatterMasslessParticles(p1,p2,Emin):
	# boost to get to COM frame
	E1 = np.sqrt(np.sum(p1**2))
	E2 = np.sqrt(np.sum(p2**2))
	
	p10,E10,p20,E20,v1,v2 = findCenterOfMassBoost(p1,E1,p2,E2)

	# randomize output angles
	theta0 = np.random.rand()*2.*np.pi
	c0 = np.cos(theta0)
	s0 = np.sin(theta0)
	rotMat = np.array([[c0,s0],[-s0,c0]])
	p10 = np.matmul(p10,rotMat)
	p20 = np.matmul(p20,rotMat)

	# reverse 2nd boost
	p11,E11 = boost(p10,E10,-v2)
	p21,E21 = boost(p20,E20,-v2)
	# reverse 1st boost
	p12,E12 = boost(p11,E11,-v1)
	p22,E22 = boost(p21,E21,-v1)
	
	if E12 <Emin or E22 <Emin:
		return p1,p2
	else:
		if np.random.rand()<0.5:
			return p12,p22
		else:
			return p22,p12