from numpy import *
  
import copy

LOG_ERROR=1
LOG_INFO=2
LOG_WARNING=3
LOG_ALLINFO=4

class loggable:
  """
  Simple class allowing for logging messages through it's debug() method.
  """
  log_level=LOG_WARNING
  _last_message = None
  _last_count = 0
  
  @classmethod
  def debug(cls,text,level):
    if loggable.log_level >= level:
      if text == cls._last_message:
        cls._last_count += 1
      elif cls._last_count > 1:
        print 'repeated (%d times)' % cls._last_count
      else:
        print text
        cls._last_message = text
        cls._last_count = 1

# Conversion factors:
  # Conversion from kbar to ev/A^3. From units:
  # You have: kbar
  # You want: ev/angstrom**3
  #  * 0.00062415097
  #  / 1602.1765
KBAR_TO_EVA3 = 0.00062415097

EVA3_TO_GPA = 1.6021765e+02

  # electron/angstrom^2 to coulomb/m^2
  # You have: e/angstrom**2
  # You want: coulomb/m**2
  #   * 16.021765
  #   / 0.062415097
EA2_TO_CM2 = 16.021765

  # You have: (e * bohrmagneton / (eV * angstrom^2)) * 4 * pi * 1e-7 V * s/(A m)*c
  # You want:
  #         Definition: 0.34938003
  #
COUPLING_TO_GAUSS = 0.34938003

  # You have: bohrmagneton / angstrom**3
  # You want: A /m
  #         * 9274009
BMA3_TO_AM = 9274009

  # You have: bohrmagneton /eV
  # You want: tesla^-1
  #         * 5.7883817e-05
BMEV_TO_T = 5.7883817e-05

  # You have: e*e / angstrom  /eV
  # You want: epsilon0
  #      * 180.95126
EEAEV_TO_EPSILON0 = 180.95126

  # You have: bohrmagneton / angstrom^4 * bohrmagneton /angstrom^4 * 1/(eV/angstrom^5 )
  # You want: 1/mu0
  #      * 0.00067458168
BM2A3EV_TO_1MU0 = 0.00067458168


def metro_dist( v1, v2 ):
  """
  Helper function used in fixing "jumped ions" and problematic polarizations.
  """
  return abs(v1 - v2).max()

def mat2str(m,mode=' % 7.4f'):
  if len( m.shape ) == 1:
    m = array([m])
  if len( m.shape ) > 2:
    m = m.reshape( m.shape[0], -1 )
      
  return '\n'.join([' '.join([mode%col for col in row]) for row in array(m) ] )+'\n'


def tensor2voigt( tens, transform_directions = [[1,2]] ):
  #return mat( [ [ a[p,0,0], a[p,1,1], a[p,2,2], (a[p,1,2]+a[p,2,1])/2., (a[p,0,2]+a[p,2,0])/2., (a[p,1,0]+a[p,0,1])/2. ] for p in range(a.shape[0]) ] )
  a = array(tens)
  for t_d in transform_directions:
    if t_d is list:
      # Assume two indices (real tensor) 
      x11 = [slice(0,dim) for dim in a.shape]
      x11[t_d[0]] = 0
      x11[t_d[1]] = 0
      x22 = [slice(0,dim) for dim in a.shape]
      x22[t_d[0]] = 1
      x22[t_d[1]] = 1
      x33 = [slice(0,dim) for dim in a.shape]
      x33[t_d[0]] = 2
      x33[t_d[1]] = 2
      x12 = [slice(0,dim) for dim in a.shape]
      x12[t_d[0]] = 0
      x12[t_d[1]] = 1
      x13 = [slice(0,dim) for dim in a.shape]
      x13[t_d[0]] = 0
      x13[t_d[1]] = 2
      x23 = [slice(0,dim) for dim in a.shape]
      x23[t_d[0]] = 1
      x23[t_d[1]] = 2
      x21 = [slice(0,dim) for dim in a.shape]
      x21[t_d[0]] = 1
      x21[t_d[1]] = 0
      x31 = [slice(0,dim) for dim in a.shape]
      x31[t_d[0]] = 2
      x31[t_d[1]] = 0
      x32 = [slice(0,dim) for dim in a.shape]
      x32[t_d[0]] = 2
      x32[t_d[1]] = 1
      position = t_d[1]
    else:
      # Assume flat setting (unrolled tensor)
      x11 = [slice(0,dim) for dim in a.shape]
      x11[t_d] = 0
      x22 = [slice(0,dim) for dim in a.shape]
      x22[t_d] = 4
      x33 = [slice(0,dim) for dim in a.shape]
      x33[t_d] = 8
      x12 = [slice(0,dim) for dim in a.shape]
      x12[t_d] = 1
      x13 = [slice(0,dim) for dim in a.shape]
      x13[t_d] = 2
      x23 = [slice(0,dim) for dim in a.shape]
      x23[t_d] = 5
      x21 = [slice(0,dim) for dim in a.shape]
      x21[t_d] = 3
      x31 = [slice(0,dim) for dim in a.shape]
      x31[t_d] = 6
      x32 = [slice(0,dim) for dim in a.shape]
      x32[t_d] = 7
      position = t_d+1
    
    a = array([a[x11], a[x22], a[x33], (a[x32]+a[x23])/2, (a[x13]+a[x31])/2, (a[x21]+a[x12])/2])
    a = rollaxis(a, 0, position)

  if isinstance(tens, matrix):
    return mat(a)
  else:
    return a    

def iterate_all_indices(shape,shape_min=None):
  """
  Generator function for looping through all possible combinations of the values in given ranges.
  """
  if shape_min is None:
    shape_min = [0]*len(shape)

  for i in range(int(shape_min[0]),int(shape[0])):
    if len(shape) > 1:
      for k in iterate_all_indices(shape[1:], shape_min[1:]):
        yield [i]+k
    else:
      yield [i]    

def iterate_shifts():
  """
  Generator function for iterating over all neighbours of 3D cube.
  """
  for dx,dy,dz in iterate_all_indices([3,3,3]):
    yield (dx-1,dy-1,dz-1)
  
def rotatetensor(tens, rot):
  """
  Rotate (general - matrix) third-rank tensor. 
  FIXME: Make it more general.
  """
  res = zeros(tens.shape)
  
  #for i,j,k,r,s,t in iterate_all_indices(tens.shape + res.shape):
  #  res[r,s,t] += rot[r,i]*rot[s,j]*rot[t,k]*tens[i,j,k]
  #return res

  for indices in iterate_all_indices(tens.shape + res.shape):
    mult = 1.0
    left_indices = indices[len(indices)/2:]
    right_indices = indices[0:len(indices)/2]
    
    for i in range(len(left_indices)):
      mult *= rot[left_indices[i], right_indices[i] ]
    
    res[tuple(left_indices)] +=  mult*tens[tuple(right_indices)]
  return res

def inv_saxis_rotation(saxis):
  """
  Conversion from cartesian to internal (reported by VASP) coordinates for magnetic moment.
  See: http://cms.mpi.univie.ac.at/vasp/vasp/node159.html
  """
  alpha = math.atan2(saxis[1], saxis[0])
  beta = math.atan2(sqrt(saxis[0]*saxis[0]+saxis[1]*saxis[1]), saxis[2])
  sa = sin(alpha)
  ca = cos(alpha)
  sb = sin(beta)
  cb = cos(beta)
  return mat([[cb*ca, cb*sa, -sb],
              [  -sa,    ca,   0],
              [sb*ca, sb*sa,  cb]])
  
def saxis_rotation(saxis):
  """
  Conversion from internal (reported by VASP) to cartesian coordinates for magnetic moment.
  See: http://cms.mpi.univie.ac.at/vasp/vasp/node159.html
  """ 
  return linalg.inv(inv_saxis_rotation(saxis))

def safe_inv(matrix, threshold = 0, iterate_components = None):
  """
  In a lot of places we need to safely invert a (almost) Hermitian matrix. This is the way 
  we want to do that: by symmetrizing the matrix and skipping eigenvalues close to 0,   
  """
  (ee,vv) = linalg.eig( (matrix+matrix.T)/2. )
  vv=mat(vv)

  s_ee = sorted(ee)
  if threshold < 0:
    threshold = (s_ee[-threshold-1]+s_ee[-threshold])/2.

  if iterate_components is None:
    iterate_components = range(matrix.shape[0])
    
  new_ee = zeros(ee.shape)
  skipcount = 0
  for i in iterate_components:
    if abs(ee[i]) < threshold:
      new_ee[i] = 0.
      skipcount += 1
    else:
      new_ee[i]=1./ee[i]

  return vv*mat(diag(new_ee))*vv.T, skipcount

def invert_with_warning(matrix, threshold, warning_msg, expected_skip):
  inv_m, skip = safe_inv(matrix, threshold)
  if skip != expected_skip:
    loggable.debug(warning_msg % skip, LOG_WARNING)
    loggable.debug('Near-zero threshold was: %f'%threshold, LOG_ALLINFO)
    loggable.debug('Eigenvalues of the matrix were:', LOG_ALLINFO)
    loggable.debug(mat2str( linalg.eig((matrix+matrix.T)/2.)[0]), LOG_ALLINFO)
  
  return inv_m                 
    

