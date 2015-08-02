""" 
blocks.py contains the Block Class

Each Block has:
	(.name) string identifying the block's name (or ID)
	(.m) material, defined the problem statement
			 materials are just a dictionary
	(.state) physical state
	(.F) list of fluxes, which return a dictionary
	(.S) list of sources, which return a dictionary
		Sources and Fluxes do not need to be ordered 
		since they are never explicitly globally unwrapped
	(.t) time
	(.T) Time function, returns dict of coefficient on time terms

The equation for the block is
R(state) = Sum(Fluxes(state)) + Sum(Sources(state)) = 0
In the unsteady case
d/d(state) = R(state)

All blocks are connected through fluxes, defined in flux.py

These are blocks rather than volumes as they have no geometric info
All geometric info is in the fluxes, which are effectively
boundary conditions on the blocks
------------------------------------
function tests are run by doctest
python blocks.py

"""
from collections import OrderedDict

class Block(object):
	""" 
	Block Class

	__init__: 	Object Constructor

	input(s):   (s) string corresponding to block name
							(m) material (Pass in None if not needed)
							(t) time
							(initialStates) Optional key-value pairs for
								 			 initial conditions
	output(s):	None
	"""

	def __init__(self,s,material,t = 0,**initialStates):
		self.name = s
		self.m = material
		self.state = OrderedDict(initialStates)
		self.F = []
		self.S = []	
		self.t = t
		self.T = lambda B : dict([(s,1) for s in B.state])

	"""
	These overload the [] operator, such that the 
	states can be set much more easily

	"""
	def __getitem__(self,key):
		return self.state[key]

	def __setitem__(self,key,val):
		self.state[key] = val
	"""
	addFlux:
	addSource: Block Setup Functions

	input(s):  (F,S) Flux objects, Source objects
	output(s): None

	these functions don't add anything new, but make
	building blocks easier using addFlux/addSource instead
	of directly modifying the lists
	"""
	def addFlux(self,F):
		self.F.append(F)
		F.B = self

	def removeFlux(self,name):
		self.F[:] = [F for F in self.F if F.name != name]

	def addSource(self,S):
		self.S.append(S)

	def removeSource(self,name):
		self.S[:] = [S for S in self.S if S.name != name]

	def updateSource(self,name,key,value):
		for S in self.S:
			if(S.name == name):
				S.updateParameters(key,value)

	"""
	R: Residual function

	input(s):  None
	output(s): dict of states and residuals corresponding to state

	sums over sources and fluxes to calculate the residual
	"""
	def R(self):
		R = OrderedDict([(s,0) for s in self.state])
		for d in [F.flux() for F in self.F]+[S.S(self) for S in self.S]:
			for s in d:
				R[s] += d[s]
		return R

	"""
	__repr__

	input(s):  None
	output(s): str representation of the object. 

	allows for using print direction on an object
	"""	
	def __repr__(self):
		return "block named " +self.name+" with "+str(len(self.state)) \
			+ " state variables of material "+self.m['name']+"\n"  \
			+ ",".join([s + '=' + str(self[s]) for s in self.state])

if __name__ == "__main__":
    import doctest
    doctest.testmod()

