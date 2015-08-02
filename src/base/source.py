"""
source.py contains the Source class

Each Source has:
	(.S) source function, evaluated 
  (.p) parameters


------------------------------------
function tests are run by doctest
python blocks.py
------------------------------------
>>> b = blocks.Block('test','water')
>>> b.state['T'] = 20.0
>>> Sa = Source('const',**{'T':0.5})
>>> Sa.S(b)['T']
0.5
>>> b.state['T'] = 10.0
>>> Sa.S(b)['T']
0.5

"""
class Source(object):
	""" 
	Source Class

	input(s):   (s) string corresponding to function name
							(parameters) optional dictionary with arguments for the source functions
	output(s):	None
	"""
	def __init__(self,s,name='',**parameters):
		self.S = eval('self.'+s)
		self.name = name
		self.p = parameters

	"""
	input(s):   (b) Block
	output(s):	dict corresponding to b.state variables
	"""
	def const(self,b):
		return dict([(state,self.p.get(state,0.0)) for state in b.state])

	def time(self,b):
		return dict([(state,self.p[state](b.t)) for state in b.state])

	def linear(self,b):
		return dict([(state,self.p.get(state,0.0)*b.state[state]) for state in b.state])

	def updateParameters(self,key,value):
		self.p[key] = value

	"""
	__repr__ 		overloading for print command

	input(s):   None
	output(s):	name
	"""
	def __repr__(self):
		return name
if __name__ == "__main__":
  import doctest
  import blocks
  doctest.testmod()

