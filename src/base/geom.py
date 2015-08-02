""" 
geom.py contains the Geometry and GeemetrySet classes

Each geometry represents a planar, rectangular facade,
given by a set of coordinates in counter-clockwise order,
such that normals point outwards and the first point is on 
the bottom right.

The class itself acts as a wrapper around the coordinates,
and coordinates can be accessed simply by iterating over the
object itself, and standard get/set/len type operations.

Quantities that are used often are stored, and quantities that
are not required often, are calculations for now.

Requires the shapely library to function properly

Each Geometry Object has:
  (.n) normal vector, as an np.array
  (.orient) orientation angle, in radians
  (.tilt) tilt angle, in radians
  (.coords) coordinates
  (._index) a number representing its location in the GeometrySet
            let the GeometrySet deal with this, dont access directly
  (.nX) the number of blocks in the X direction
  (.nY) the number of blocks in the Y direction
  (.dir) (NSEW/Roof)
  (.data) a dictionary of data
  (.R) rotation matrix, transposed for ease of use,
       transpose(R(tilt)*R(orient))
  (.poly) Polygon Representation of the surface in x-z plane

Each GeometrySet Object is a collection of Geometry Objects,
in an easily accessible way, allowing for queries
it is more or less, an overloaded list

geometry objects do not know what 'index' they have
or what other geometry objects shade them, this is all
stored in the GeometrySet

Each GeometrySet Object has:
  (.g) a list of geometry objects
  (.s) a list of subsets for each geometry object
  (.m) a list of matches, similar geometries, 
       for each geometry object
  (.name) a name
  (.dataNames) a list of data stored in the geometry

There is also a SimpleGeometry object,
which only has orientation and tilt, and really acts
as a way to use functions that only need these two things, 
without needing coordinates. Rather than have a base class
for geometry that inherits, this, its just easier to have this trivial
class created, and pass it in the place of geometry objects when
coordinates and such is not needed
"""

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
try:
  from shapely.geometry import Polygon
  SHAPELY = True
except ImportError:
  SHAPELY = False

class SimpleGeometry(object):
  """ 
  SimpleGeometry Class

  __init__:   Object Constructor

  input(s):   orientation (radians)
              tilt (radians)
  output(s):  None
  """
  def __init__(self,orient,tilt):
    self.orient = orient
    self.tilt = tilt
  def __repr__(self):
    return 'SimpleGeometry Object with orient '+str(self.orient) \
      +' and tilt '+str(self.tilt)
class Geometry(object):
  """ 
  Geometry Class

  __init__:   Object Constructor

  input(s):   coordinates
  output(s):  None
  """

  def __init__(self,coords):
    if(len(coords) != 4):
      sys.exit("surface with more than 4 coordinates.\n \
        only rectangular facades supported for now\n")

    self.coords = coords

    self.n = np.cross(coords[2]-coords[0],coords[1]-coords[0])
    self.n = self.n/np.linalg.norm(self.n)
    # correction for ordering based on z   
    if(np.dot(self.n,[0,0,1]) < 0.):
      self.n = -self.n
    self.orient = np.pi/2.+np.arctan2(self.n[1],self.n[0])
    self.tilt = np.pi/2.- np.arctan2(np.sqrt(self.n[0]*self.n[0]
      +self.n[1]*self.n[1]),self.n[2])

    self.nX = 1
    self.nY = 1
    self._index = 0

    # its a roof if the tilt is greater than 51 degrees
    if(abs(self.orient) <= np.pi/4. and self.tilt < 2.*np.pi/7.):
      self.dir = 'south'
    elif(abs(self.orient-np.pi/2.) <= np.pi/4. and self.tilt < 2.*np.pi/7.):
      self.dir = 'east'
    elif(abs(self.orient+np.pi) <= np.pi/4. and self.tilt < 2.*np.pi/7.):
      self.dir = 'north'
    elif(abs(self.orient+np.pi/2.) <= np.pi/4. and self.tilt < 2.*np.pi/7.):
      self.dir = 'west'
    elif(self.tilt > 2.*np.pi/7.):
      self.dir = 'roof'
    else:
      print self.orient
      sys.exit('invalid orientation')

    Rx = np.matrix([[1,0,0],[0,np.cos(self.tilt), -np.sin(self.tilt)],
      [0,np.sin(self.tilt), np.cos(self.tilt)]],float)
    Rz = np.matrix([[np.cos(self.orient), -np.sin(self.orient), 0], 
      [np.sin(self.orient), np.cos(self.orient), 0], [0,0,1]],float)
    self.R = np.transpose(Rx*Rz)

    if SHAPELY:
      self.poly = Polygon([(c[:,0],c[:,2]) for c in 
        [np.asmatrix(self[i])*self.R
        for i in range(len(self))]])

    self.data = {}

  """
              overloading of container functions

  input(s):   None
  output(s):  name
  """

  def __len__(self):
    return 4
  def __iter__(self):
    return iter(self.coords)
  def __getitem__(self,key):
    return self.coords[key]
  def __setitem__(self,key,val):
    if(len(val) != 3):
      sys.exit("coordinates must arrays of length 3")
    self.coords[key] = val
  def __repr__(self):
    return "Facade with height "+str(self.height())+" width " \
      + str(self.width()) + " in direction "+self.dir

  """
  height:     computes height and width of block
  width:
  
  input(s):   None
  output(s):  height, width
  """
  def height(self):
    return np.linalg.norm(self.coords[1]-self.coords[0])

  def width(self):
    return np.linalg.norm(self.coords[2]-self.coords[1])

class GeometrySet(object):
  """ 
  GeometrySet Class

  __init__:   Object Constructor

  input(s):   
  output(s):  None
  """
  def __init__(self,name):
    self.g = []
    self.s = None
    self.m = None
    self.name = name
    self.dataNames = None
    if not SHAPELY:
      print "\n!!! Warning: cannot compute shading \
      without shapely package!!! \n \
      https://pypi.python.org/pypi/Shapely\n"

  """
  len,iter,       these allow the GeometrySet to
  get/set item    be treated as a container, and
  append          iterated over

  input(s):       key, integer
                  val, Geometry Object
  output(s):      as shown
  """

  def __len__(self):
    return len(self.g)
  def __iter__(self):
    return iter(self.g)
  def __getitem__(self,key):
    return self.g[key]
  def __setitem__(self,key,val):
    self.g[key] = val

  """
  index:          This is a useful function, 
                  returning the index of a Geometry Object
                  in the container

  input(s):       Geometry Object
  output(s):      integer indicating location
  """  
  def index(self,geom):
    return geom._index

  def append(self,geom):
    geom._index = len(self)
    self.g.append(geom)

  """
  __repr__    overloading for print command

  input(s):   None
  output(s):  name
  """
  def __repr__(self):
    return "GeometrySet " + name +" with "+str(len(self))+" surfaces"

  """
  getMatch:   returns the match or subset list for particular
  getSubset   surface

  input(s):   Geometry Object
  output(s):  List of Geometry Objects
  """
  def getMatches(self,A):
    return self.m[self.index(A)]

  def getSubset(self,A):
    return self.s[self.index(A)]

  """
  computeBlockCounts:   sets nX, nY for all geometries

  input(s):   height and width
  output(s):  none
  """
  def computeBlockCounts(self,blockHeight,blockWidth):
    for g in self.g:
      g.nX = int(g.width()/blockWidth)
      g.nY = int(g.height()/blockHeight)

  """
  collectDataNames:   sets the list of strings corresponding
                      to all the data stored in the geometry
                      really should only be used for output

  input(s):           none
  output(s):           none
  """
  def collectDataNames(self):
    if self.dataNames is None:
        # collect DataNames
      dataNames = set()
      for g in self.g:
        for k in g.data.keys():
          dataNames.add(k)
      self.dataNames = list(dataNames)
  """
  initializeData:   initialize data arrays as arrays of length nY

  input(s):         names, list of data to initialize
  output(s):        none
  """      
  def initializeData(self,names):
    for g in self.g:
      for name in names:
        g.data[name] = np.zeros(g.nY)
  """
  computeSubsets:   computes the subset of surfaces that
                    could shade another surface, based
                    on solar rays. See doc/shaded.tex

  input(s):   none
  output(s):  none
  """
  
  def computeSubsets(self):
    if self.s is not None:
      return
    self.s = []
    for i in range(len(self.g)):
      self.s.append([])  
    if not SHAPELY:
      return  
    for A in self.g:

      min_z = min([c[2] for c in A])
      firstpass = [B for B in self.g if (np.dot(A.n,B.n) > 1e-5
      and self.g.index(A) != self.g.index(B) 
      and max([c[2] for c in B]) >= min_z)]

      secondpass = []
      for B in firstpass:
        for i in range(4):
          p = B[i]-A[i]   
          if np.dot(p,A.n) > 1e-10:
            secondpass.append(B)
            break
      # check if C blocks B from shading A
      for B in secondpass:
        polygonB = B.poly
        for C in secondpass:
          if(self.g.index(B) == self.g.index(C)): continue

          pC = [] # projected coordinates
          for i in range(4):
            n = B[i]-A[i]
            # project point from A through C onto B in direction of B-A
            s = np.dot(B.n,(B[0]-C[i]))/np.dot(B.n,n)
            if s < 1e-10:
              continue
            pC.append(C[i]+n*s) # projected points

          polygonB = self.getIntersection(polygonB,pC,B.R)

        if(polygonB.area/B.poly.area > 1e-5):
          self.s[self.index(A)].append(B)
  """
  similar:    determines whether two objects are geometrically similar
              same normal, height, width, and shading
  input(s):   two geometry objects
  output(s):  True/False
  """
  def similar(self,A,B):
    tol = 1e-4
    if(abs(A.height() - B.height()) > tol):
      return False
    if(abs(A.width() - B.width()) > tol):
      return False
    if(np.dot(A.n,B.n) < 1.-tol):
      return False
    return True

  """
  computeMatches:   computes the subset of surfaces that are
                    geometrically similar, see above

  input(s):         none
  output(s):        none
  """
  
  def computeMatches(self):
    if self.m is not None:
      return
    self.m = []
    for i in range(len(self.g)):
      self.m.append([])
    if not SHAPELY:
      return      
    for A in self.g:
      # B is in the match of A if:
      # it hasnt been added, its similar, its not exactly the same one
      # and for now, its subset is either empty or size 1
      similar = [B for B in self.g 
        if (B not in self.m[self.index(A)] and 
          self.similar(A,B) and self.index(A) != self.index(B)
          and len(self.s[self.index(A)]) == len(self.s[self.index(B)])
          and len(self.s[self.index(A)]) <= 1)]

      # now if they both have sizes of 1 or less
      for B in similar:
        if len(self.s[self.index(A)]) == 1:
          AA = self.s[self.index(A)][0]
          BB = self.s[self.index(B)][0]
          if self.similar(AA,BB) and \
            sum([np.dot(AA[i]-A[i],BB[i]-B[i]) for i in range(4)]):
            self.m[self.index(A)].append(B)
            self.m[self.index(B)].append(A)             
        else:
          self.m[self.index(A)].append(B)
          self.m[self.index(B)].append(A)

  """
  getIntersection:  returns the remaining part of polygonA
                    not intersected by the coordinates in B
                    A - intersection(A,B)

  input(s):         polygonA, Polygon from geometry object A
                    B, list of coordinates
                    R, rotation matrix
                    plot (optional argument to plot these)
  output(s):        Polygon object

  This input format is used since we know what A is
  on the x-z plane, and then make a bunch of comparisons with other 
  geometry objects, B, and this reduces the number of times we create
  polygons, while simplifying the user API.
  """
  def getIntersection(self,polygonA,B,R,plot=False):
    if not SHAPELY:
      sys.exit("unable to compute intersections without shapely package")
    if(len(B) < 3):
      return polygonA
    RB = (B*R)
    polygonB = Polygon([(RB[i,0],RB[i,2]) for i in range(4)])
    intersection = polygonA-polygonB.intersection(polygonA)
    
    if(plot):
      plotIntersections(polygonA,polygonB,intersection)

    return intersection

  """
  getUnshadedFraction:  gets the sunlit fraction of facade A

  inputs(s):            A, Geometry Object
                        n, vector
  outputs(s):           float, between 0 and 1
  """
  def getUnshadedFraction(self,A,n):
    if not SHAPELY:
      return 1.0
    # is the vector pointing into the surface
    dotAn = np.dot(A.n,n);
    if(dotAn) > 1e-10:
      return 0.

    polygonA = A.poly
    for B in self.getSubset(A):
      pPB = []
      for c in B:
        s = np.dot(A.n,(A[0]-c))/dotAn
        if s < 1e-10:
          continue          
        pPB.append(c+n*s)
      polygonA = self.getIntersection(polygonA,pPB,A.R)
    return polygonA.area/A.poly.area

"""
plotIntersections:    support for getIntersections
                      plots the intersection in the x-z plane

input(s):       pA,pB,pI (Polygon Objects)
output(s):      none
"""      
def plotIntersections(pA,pB,pI,iA,iB):
  plt.plot([x for (x,y) in pI.exterior.coords],
    [y for (x,y) in pI.exterior.coords],
    label='intersection',linewidth=5.0)
  plt.plot([x for (x,y) in pA.exterior.coords],
    [y for (x,y) in pA.exterior.coords],
    label='A',linewidth=2.0,color='red')
  plt.plot([x for (x,y) in pB.exterior.coords],
    [y for (x,y) in pB.exterior.coords],
    label='B',linewidth=2.0,color='black')
  plt.legend(loc=0)
  plt.axis('equal')
  plt.show()

"""
plot:           support for GeometrySet

input(s):       GeometrySet Object
output(s):      none
""" 
def plot(GeometrySet):
  ax = plt.gcf().gca(projection='3d')
  ax.view_init(60,210)
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  ax.set_zlabel('z')
  for g in GeometrySet:
    xyz = []
    xyz_ = []
    for i in range(3):
      xyz.append([c[i] for c in g])
      xyz_.append(np.mean(xyz[i]))
      xyz[i].append(xyz[i][0])

    ax.plot(xyz[0],xyz[1],zs=xyz[2],color='black')
    ax.quiver(xyz_[0],xyz_[1],xyz_[2],g.n[0],g.n[1],g.n[2],
      length=5.0,linewidth=2.0,color='red')
    ax.text(xyz_[0],xyz_[1],xyz_[2],str(GeometrySet.index(g)))
  plt.show()