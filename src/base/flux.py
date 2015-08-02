""" 
flux.py contains a list of flux functions 

This is meant to serve as a wrapper
for all flux functions

Each flux computes the heat transfer rate between two blocks
Each flux contains:
	The block it belongs to (.B)
	The block its connected to (.N)
	The type (.f), this is a pointer to the function
	The geometry between the left and right state (.G) (or parameters)

Rather than define flux functions in here, define them in the problem file
and pass them into here. 
------------------------------------
function tests are run by doctest
python flux.py
------------------------------------

"""
class Flux(object):
	"""
	Flux object. Each flux is effectively a boundary condition

	__init__:		Flux Constructor

	input(s):   (N) Neighboring block
							(f) Flux function name
							(P) Parameters
							(name) identifying name

	output(s):	None
	"""
	def __init__(self,N,f,P=None,name=''):
		self.B = None # this will be set when its added to the block
		self.N = N
		self.F = f 
		self.P = P
		self.name = name

	"""
	The following is a wrapper for flux function choices defined with
	Fluxes should be define as such
	input(s):   (B) Block
							(N) Neighboring block
							(P) Parameters (optional) 
	output(s):	dict with entries for the flux for each state contributed to

	Fluxes can have blocks with different physical states
	provided the flux function passes contributions to some of the states

	For example, L.state could have 'T','rho' and R.state could have just 'T'
	as long as the flux between them only affects 'T'


	The bulk of the hard-coding and material calls are in here
	"""

	def flux(self):
		return self.F(self.B,self.N,self.P)

	"""
	__repr__ 		overloading for print command

	input(s):   None
	output(s):	name
	"""
	def __repr__(self):
		return name+" "+self.F.__name__
		
if __name__ == "__main__":
	import doctest
	doctest.testmod()

