"""
problem.py contains the Problem Class

Each Problem has:
	(.b) the blocks to solve for in the problem
	(.bc) the blocks used as boundary blocks
	(.mapping) the mapping between the local and global systems


------------------------------------
function tests are run by doctest
python problem.py
------------------------------------
"""

"""
fsolve is used for solving steady state
ode is used for solving transient
"""
from scipy.optimize import fsolve
from scipy.integrate import odeint
from collections import OrderedDict
import numpy as np
import numdifftools as nd
import sys

class Problem(object):
	""" 
	Problem Class

	__init__:		Problem Constructor

	input(s):   (blocks) relevant blocks
							(boundaries) blocks on the boundaries
								these are needed for unsteady problems
	output(s):	None
	"""
	def __init__(self,blocks,boundaries = [],**parameters):
		self.b = blocks
		namelist = set([b.name for b in blocks])
		if(len(namelist) != len(blocks)):
			sys.exit("multiple blocks have the same name")
		self.bc = boundaries
		self.mapping = [(i, k) for i, b in enumerate(blocks) for k in b.state.keys()]


	"""
	__repr__ 		overloading for print command

	input(s):   None
	output(s):	other information
	"""
	def __repr__(self):
		return "problem with "+str(len(b))+" blocks"
	"""
	__getitem__, __setitem__:	get/set blocks by integer key in list

	input(s):   access blocks by index
	output(s):	as shown
	"""
	def __getitem__(self, key):
		return self.b[key]
	def __setitem__(self, key, value):
		self.b[key] = value


	"""
	update:			Updates the blocks by unwrapping the new solution

	input(s):   (solution) global array of floats corresponding to mapping
	output(s):	None
	"""     
	def update(self,solution):
		for ix, (i,k) in enumerate(self.mapping):
			self.b[i][k] = solution[ix]
		for bc in self.bc:
			for s in bc.state:
				bc[s] = bc.S[0].S(bc)[s]

	def updateUnst(self,t):
		for b in self.b + self.bc:	
			b.t = t

	def getSolutionVec(self):
		solution = [None]*len(self.mapping)
		for ix, (i,k) in enumerate(self.mapping):
			solution[ix] = self.b[i][k]
		return solution

	"""
	r:					Global residual function r(solution) 

	input(s):    (solution) global array of floats corresponding to mapping
	output(s):	R(solution) global array of floats corresponding to residual

	updates solution first, then computes
	should be passed into another function
	"""
	def r(self,solution):
		self.update(solution)
		return [self.b[i].R()[v] for i,v in self.mapping]

	def rVec(self,solution):
		self.update(solution.tolist())
		return np.array([self.b[i].R()[v] for i,v in self.mapping])

	def rUnst(self,solution,t):
		self.updateUnst(t)
		self.update(solution)
		return [self.b[i].R()[v]/self.b[i].T(self.b[i])[v] for i,v in self.mapping]

	"""
	solve:			wrapper for chosen (non)linear solver

	input(s):   None
	output(s):	None

	unwraps blocks, passes into solver, finishes by updating blocks one last time
	""",
	def solve(self,t=0):
		# self.updateUnst(t)
		solution = fsolve(self.r, self.getSolutionVec())

	def jacobian(self):
		Jfun = nd.Jacobian(self.rVec)
		return Jfun(self.getSolutionVec())
	"""
	solveUnst:	solve the transient problem

	input(s):   (ti) initial time
							(tf) final time
							(n)  number of timesteps
	output(s):	None

	unwraps blocks, passes into solver, finishes by updating blocks one last time
	"""

	def solveUnst(self,t):
		solution = [None]*len(self.mapping)
		# This has the unsteady part
		# Solver, just live and let live	
		for ix, (i,k) in enumerate(self.mapping):
			solution[ix] = self.b[i][k]
		soln = odeint(self.rUnst, solution, t,hmax=(t[-1]-t[0])/len(t),rtol = 1e-4, atol = 1e-4)

		# final update
		self.updateUnst(t[-1])
		self.update(soln[-1,:])

		# Lets output all the steps
		fullSolution = dict([(b.name + '_'+s,[]) for b in self.b for s in b.state])
		for j in range(0,len(t)):
			for ix, (i,k) in enumerate(self.mapping):
				fullSolution[self.b[i].name+'_'+k].append(soln[j,ix])
		fullSolution['t'] = t
		return fullSolution
		
if __name__ == "__main__":
    import doctest
    doctest.testmod()

