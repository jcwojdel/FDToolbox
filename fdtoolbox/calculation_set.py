from numpy import *


from fdtoolbox.utility import *

import os
import atexit
import re
try:
  import gzip
  have_gzip = True
except:
  have_gzip = False

class calculation(loggable):
  """
  A class for keeping all data read from OUTCAR for a single calculation. It reads the following
  properties of the calculation and stores them in the following attributes:
  - self.unit_cell - unit cell (3x3 A)
  - self.recip_cell - reciprocal unit cell (3x3 A**-1)
  - self.volume - volume calculated from the unit cell (1 A**3)
  - self.energy - total energy reported by VASP (1 eV)
  - self.name - name of the system as put in POSCAR
  - self.species - list of species numbers 
  - self.num_atoms - number of atoms
  - self.atoms - atomic positions internal coordinates (self.num_atomsx3 1)
  - self.forces - forces acting on ions as energy per internal coordinate derivative (self.num_atomsx3 eV)
  - self.stress - stress on the cell in cartesian coordinates (3x3 eV/A**3)
  - self.magnetization - total cell magnetization in cartesian coordinates (3 e)
  - self.polarization - polarization as read by load_polarization() method (3 e/A**2)
  For other polarization related attributes see load_polarization() method documentation.
  """

  saxis = [ 0., 0., 1. ]
  
  def get__volume(self):
    if not hasattr(self,'__volume'):
      self.__volume = linalg.det(self.unit_cell)
    return self.__volume    
  volume = property( get__volume, None )
  
  def load_polarization(self, filename):
    """
    This function tries to read Berry phase calculations from three files:
    filename+'_berry_1' filename+'_berry_2' filename+'_berry_3'
    If the files are not readable, no polarization attribute is set.
    On successful run this method fills the following attributes:
    - self.polarization - polarization as calculated from the three files, in cartesian coordinates (3 e/A**2)
    - self.ev_phase - expectation value reported by VASP for each of the three directions (3x3 -e*A)
    - self.ev_recip - the mean of the expectation value projected to internal coordinates (3 -e)
    - self.berry_phase - berry phases of all strings in each of the directions (3xN 1)
    - self.polarization_quant - polarization quant for the calculation (3x3 e/A**2)
    - self.i_polarization_quant - inverse of polarization quant
    """
    self.berry_phase=3*[None]
    self.berry_ev=3*[None]
    self.berry_ion=3*[None]
    self.species = []
    for kdirection in range(1,4):
      fname = ''.join( [ filename, '_berry_%d' % kdirection ] )
      if os.access( fname, os.R_OK ):
        ifile = file( fname )
      elif have_gzip and os.access(fname+'.gz', os.R_OK ):
        ifile = gzip.open(fname+'.gz')
      else:
        self.berry_phase=3*[array([[0.,0.]])]
        self.berry_ev=3*[array([0., 0., 0.])]
        self.berry_ion=3*[[[0., 0., 0.]]]
        self.polarization_quant = self.unit_cell.copy() / self.volume
        self.berry_multiplier = 1.
        self.polarization_quant *= self.berry_multiplier
        self.i_polarization_quant = linalg.inv(self.polarization_quant)
        self.ev_recip = dot( mean(self.berry_ev,axis=0), self.recip_cell )
        self.recalculate_polarization()
        return      
      line = " "
      self.berry_phase[kdirection-1] = []
      self.berry_ev[kdirection-1] = []
      while line:
        line = ifile.readline()
        if line.startswith("   ISPIN  ="):
          data = line.split()
          self.numspin = int(data[2])
        elif line.startswith("   LNONCOLLINEAR ="):
          # This is a bit tricky - noncollinear calculations are officially numspin==1, but
          # they do contain separate electrons, not electron pairs, and therefore should be 
          # treated as numspin==2 for calculation of polarizations.
          data = line.split()
          if data[2] == "T":
            self.numspin = 2
        elif line.startswith("Expectation value"):
          line = ifile.readline()
          data = line.split() 
          self.berry_ev[kdirection-1].append(  [ float(data[3][:-1]), float(data[4][:-1]), float(data[5]) ] )
        elif line.startswith("K-point string"):
          data = line.split()
          last_kpoint_weight = float(data[5])
        elif line.startswith("Berry-Phase term:  "):
          data = line.split()
          self.berry_phase[kdirection-1].append( [float(data[2]), last_kpoint_weight] )
        elif line.startswith("ionic term:"):
          line = ifile.readline()
          data = line.split() 
          self.berry_ion[kdirection-1] = [ [ float(data[3][:-1]), float(data[4][:-1]), float(data[5]) ] ]
      
      self.berry_phase[kdirection-1] = array(self.berry_phase[kdirection-1])
      self.berry_ev[kdirection-1] = mean(array(self.berry_ev[kdirection-1]), axis=0)
          
    self.polarization_quant = self.unit_cell.copy() / self.volume
            
    # In case of nonpolarised calculations, we have two electrons per each Berry string
    # So the quantum needs to be doubled.
    self.berry_multiplier = 3-self.numspin
    self.polarization_quant *= self.berry_multiplier
    self.i_polarization_quant = linalg.inv(self.polarization_quant)
    
    self.ev_recip = dot( mean(self.berry_ev,axis=0), self.recip_cell )
    
    self.recalculate_polarization()
    
  def berry_term(self):
    """
    Returns Berry term part of the polarization for each of the IGPAR directions (i.e. 
    k-point weighted sum of Berry terms for each k-point). The term is multiplied by proper 
    multiplier to make sure that closed-shell calculations sum up to proper number of electrons. 
    """
    bp = zeros((1,3))
    for i in range(3):
      # berry_phase[i][:,0] -> berry term, self.berry_phase[i][:,1] -> k-point weight
      bp[0,i] = sum(self.berry_phase[i][:,0]*self.berry_phase[i][:,1])
    return self.berry_multiplier*mat(bp)
    
  def ionic_term(self):
    """
    Returns ionic part of polarization. Contrary to VASP it calculates this with 0 as a 
    reference point (VASP uses fract(1/2,1/2,1/2)).
    """
    ion=zeros((1,3))
    for atom, zval in zip(self.atoms, self.zvals):
      ion += -zval*atom
    return mat(ion)
    
      
  def recalculate_polarization(self):
    """
    Recalculates the polarization from the data gathered from the OUTCAR_berry_[123] files.
    This function should be called whenever we tamper with Berry phases (to fix problems with
    polarization quantum).
    """    
    # Gather data from all three directions
    
    self.polarization = [self.berry_term()*self.unit_cell] + mean(self.ionic_term(),axis=0) + mean(self.berry_ev,axis=0)

    # Note that VASP reports data in 'electrons * A', meaning the lack of minus sign
    # and not scaling by cell volume.
    self.polarization /= -self.volume

    return self.polarization
    # At the end, the polarization is expressed in e/A**2
    


  def load_from_outcar(self, filename):
    """
    This function reads atom positions, unit cell, forces and stress tensor from the VASP
    OUTCAR file
    """
    if os.access(filename, os.R_OK):
      F=file(filename)
    elif have_gzip and os.access(filename+'.gz', os.R_OK):
      F=gzip.open(filename+'.gz')
    else:
      raise Exception("Missing file")
      self.debug('Could not read %s' % filename, LOG_ERROR)
    
    line=" "
    while line:
      line = F.readline()
      if line.startswith("   number of dos"): 
        #   number of dos      NEDOS =    301   number of ions     NIONS =     10
        data = line.split()
        self.num_atoms = int(data[-1])
      elif line.startswith(" POSITION"):
        # POSITION                                       TOTAL-FORCE (eV/Angst)
        # -----------------------------------------------------------------------------------
        #     -0.00003     -0.00002      0.07380         0.000005      0.000005     -0.000031
        #      0.00003     -0.00002      6.85381         0.000003      0.000006     -0.000031
        F.readline()
        #data = fromfile( F, dtype('float'), 6*self.num_atoms, ' ' ).reshape((self.num_atoms,6))
        data = zeros((self.num_atoms,6))
        for i in range(self.num_atoms):
          data[i,:] = [float(val) for val in F.readline().split() ]
        
        self.atoms  = mat(data[:,0:3])
        self.forces = mat(data[:,3:6])
      elif line.startswith(" number of electron"):
        # number of electron   94.0000001 magnetization    0.0000060  -0.0364612   0.0000000
        data=line.split()
        self.__magnetization = mat(zeros((1,3)))
        if len(data) == 8:
          self.__magnetization[0,0]=float(data[5])
          self.__magnetization[0,1]=float(data[6])
          self.__magnetization[0,2]=float(data[7])
        elif len(data) == 6:
          self.debug("Only collinear magnetization available", LOG_WARNING)
          self.__magnetization[0,2]=float(data[5])
        else:
          self.debug("No magnetization available", LOG_WARNING)

      elif line.startswith(" VOLUME and BASIS"):
        # VOLUME and BASIS-vectors are now :
        # -----------------------------------------------------------------------------
        #  energy-cutoff  :      500.00
        #  volume of cell :      118.77
        #      direct lattice vectors                 reciprocal lattice vectors
        #     0.000000641  3.180190650  4.520060784    -0.000000095  0.209629759  0.073746221
        #    -2.754138309 -1.590096578  4.519978046    -0.181544893 -0.104816800  0.073746248
        #     2.754139668 -1.590094072  4.519978030     0.181544989 -0.104816797  0.073746194
        F.readline()
        F.readline()
        F.readline()
        F.readline()
        data = zeros((3,6))
        for i in range(3):
          data[i,:] = [float(val) for val in F.readline().split() ]
        
        self.unit_cell = mat(data[:,0:3])
      elif line.startswith("  in kB"):
        # Note that we are reading the line expressed in kBar, because due to the units it prints out 
        # more significant digits. It is then converted to eV/A**3
        #  FORCE on cell =-STRESS in cart. coord.  units (eV/reduce length):
        #  Direction    X        Y        Z        XY       YZ       ZX
        #  --------------------------------------------------------------------------------------
        #  Alpha Z    504.77    504.77    504.77
        #  Ewald    -2811.56  -2811.69  -2888.52      0.00      0.00      0.00
        #  Hartree    969.00    968.94    916.13      0.00      0.00      0.00
        #  E(xc)     -495.50   -495.50   -495.50      0.00      0.00      0.00
        #  Local       97.66     97.85    230.25      0.00      0.00      0.00
        #  n-local    -61.49    -61.50    -58.22      0.00      0.00      0.00
        #  augment    450.97    450.97    448.01      0.00      0.00      0.00
        #  Kinetic   1346.16   1346.17   1343.09      0.00      0.00      0.00
        #  -------------------------------------------------------------------------------------
        #  Total        0.01      0.01      0.00      0.00      0.00      0.00
        #  in kB        0.07      0.07      0.05      0.00      0.00      0.00
        data=line.split()
        self.stress = mat(zeros((3,3)))
        self.stress[0,0] = float(data[2])
        self.stress[1,1] = float(data[3])
        self.stress[2,2] = float(data[4])
        self.stress[0,1] = self.stress[1,0] = float(data[5])
        self.stress[2,1] = self.stress[1,2] = float(data[6])
        self.stress[0,2] = self.stress[2,0] = float(data[7])
        # Conversion from kbar to ev/A^3.
        self.stress *= KBAR_TO_EVA3
      elif line.startswith("  free  energy"):
        #  free  energy   TOTEN  =       -71.184108 eV
        data = line.split()
        self.energy = float(data[-2])
      elif line.startswith(" POSCAR ="):
        # Name of the system
        self.name = line[line.find("=")+1:]
      elif line.startswith("   ions per type ="):
        # Ionic species information
        self.species = line.split()[4:]
      elif line.startswith("   ZVAL   ="):
        # Effective Z of each nuclei
        # Note that it might be multiline 
        data = []
        while line.startswith("   ZVAL   ="):
          data.extend( line.split()[2:] )
          line = F.readline()
          
        self.zvals = []
        for num,sp in enumerate(self.species):
          self.zvals.extend(int(sp)*[float(data[num])])
        
      elif line.startswith(" total charge "):
        F.readline()
        F.readline()
        F.readline()
        self.charges = [0.]*self.num_atoms
        for i in range(self.num_atoms):
          line = F.readline()
          data = line.split()
          self.charges[i] = float(data[-1])
      elif line.startswith(" magnetization ("):
        data = line.split()
        if data[1]=="(x)":
          mindx=0
        if data[1]=="(y)":
          mindx=1
        if data[1]=="(z)":
          mindx=2
        F.readline()
        F.readline()
        F.readline()
        if not hasattr(self,"proj_magn"):
          self.proj_magn = zeros((self.num_atoms, 3))
        for i in range(self.num_atoms):
          line = F.readline()
          data = line.split()
          self.proj_magn[i,mindx] = float(data[-1])
        
        

    self.load_polarization( filename )
    self.fileID = filename
    F.close()
    
  def get__magnetization(self):
    """
    Returns the magnetization vector in Cartesian coordinates. The self.__magnetization always
    stores the magnetization as read from OUTCAR. Note that for optimization purposes the rotated
    magnetization is stored between calls. If self.saxis is changed, the invalidate_magnetization 
    method _must_ be called on all Calculation class objects.    
    """
    # Rotation of reported magnetic moment - according to:
    # http://cms.mpi.univie.ac.at/vasp/vasp/node159.html
    if not hasattr(self, "__rotated_magnetization"):
      self.__rotated_magnetization = dot( self.__magnetization, inv_saxis_rotation(self.saxis) )
      
    return self.__rotated_magnetization 
  
  def invalidate_magnetization(self):
    """
    Clear magnetization expressed in cartesian coordinates.
    """
    del self.__rotated_magnetization
    
  magnetization = property( get__magnetization, None )
    
    
  def load_from_poscar(self, filename):
    """
    Loads the system from POSCAR. Obviously, everything except for atomic positions is missing.
    """
    F = open( filename, 'r' )
    self.fileID = filename
    self.name = F.readline()
    scale = float(F.readline())
    self.unit_cell = scale*mat( fromfile( F, dtype('float'), 9, ' ' ).reshape((3,3)) )
    
    self.species = F.readline().split()
    self.num_atoms = 0
    for i in self.species:
      self.num_atoms += int(i)
        
    mode = F.readline()
  
    self.atoms = mat(fromfile( F, dtype('float'),  self.num_atoms*3, ' ' ).reshape(( self.num_atoms,3)))
    self.forces = zeros(self.atoms.shape)
  
    if re.search('^[cCkK]',mode):
      pass
    else:
      self.atoms = self.atoms*self.unit_cell
  
    F.close()
    
    if self.name.split()[0] == "SUPERCELL":
      self.is_supercell = True
      self.supercell_repetitions = self.name.split()[1].split('x')
      self.supercell_repetitions = [int(i) for i in self.supercell_repetitions]
    
  def save_to_poscar(self, filename):
    """
    Save the system to POSCAR
    """ 
    F = open( filename, 'w' )
    F.write( self.name )
    F.write( "    1.0\n" )
    F.write( mat2str( self.unit_cell, "%16.10f" ) )
    F.write(' '.join(self.species) )
    F.write('\n')
    F.write("Cart\n")
    F.write( mat2str( self.atoms, "%16.10f" ) )
    F.close()

  def save_to_xyz(self, filename, species = None):
    """
    Save the system to XYZ
    """ 
    F = open( filename, 'a' )
    F.write( '%d\n'%self.num_atoms )
    F.write( "XYZ\n" )
    for num,row in enumerate(self.atoms):
      if species is not None:
        F.write('%s  '%species[num])
      else:
        F.write('X%d  '%num)
      F.write( mat2str( row, "%16.10f" ) )
    F.write( "\n" )
    F.close()
    
  def save_to_arc(self, filename, species = None, header = True):
    """
    Save the system in ARC file
    """
    if header:
      F = open( filename, 'w' )
      F.write( "!BIOSYM archive 2\n" )
      F.write( "PBC=ON\n" )
    else:
      F = open( filename, 'a' )

    if species is None:
      # Prepare a list of species in the calculation if none specified
      # ARC does not allow for dummy atoms, so we fill everything with carbon
      species=self.num_atoms*['C']
            
    #FIXME: If you think this is the ugliest python code you've ever seen,
    # you are quite right! It is literal translation of some old AWK script.
    # But it works for now, so... 

    unit_cell = self.unit_cell
    a=sqrt(unit_cell[0,0]*unit_cell[0,0]+
         unit_cell[0,1]*unit_cell[0,1]+
         unit_cell[0,2]*unit_cell[0,2])
    b=sqrt(unit_cell[1,0]*unit_cell[1,0]+
         unit_cell[1,1]*unit_cell[1,1]+
         unit_cell[1,2]*unit_cell[1,2])
    c=sqrt(unit_cell[2,0]*unit_cell[2,0]+
         unit_cell[2,1]*unit_cell[2,1]+
         unit_cell[2,2]*unit_cell[2,2])
    alpha=(unit_cell[1,0]*unit_cell[2,0]+
         unit_cell[1,1]*unit_cell[2,1]+
         unit_cell[1,2]*unit_cell[2,2])/(b*c)
    beta =(unit_cell[0,0]*unit_cell[2,0]+
         unit_cell[0,1]*unit_cell[2,1]+
         unit_cell[0,2]*unit_cell[2,2])/(a*c)
    gamma=(unit_cell[0,0]*unit_cell[1,0]+
         unit_cell[0,1]*unit_cell[1,1]+
         unit_cell[0,2]*unit_cell[1,2])/(a*b)
    alpha=math.atan2(sqrt(1-alpha*alpha),alpha)
    beta =math.atan2(sqrt(1-beta *beta ),beta )
    gamma=math.atan2(sqrt(1-gamma*gamma),gamma)

    transf=zeros((3,3))
    transf[0,0]=a
    transf[1,0]=0.0
    transf[2,0]=0.0
    transf[0,1]=b*cos(gamma)
    transf[1,1]=b*sin(gamma)
    transf[2,1]=0.0
    transf[0,2]=c*cos(beta)
    transf[1,2]=c*(cos(alpha)-(cos(gamma)*cos(beta)))/sin(gamma)
    transf[2,2]=sqrt(c*c-transf[0,2]*transf[0,2]-transf[1,2]*transf[1,2])

    alpha=180*alpha/pi
    beta =180* beta/pi
    gamma=180*gamma/pi

    recip_cell = self.recip_cell.T
    frac_pos = zeros(self.atoms.shape)
    positions= zeros(self.atoms.shape)
    for i in range(self.num_atoms):
      for j in range(3):
        frac_pos[i,j]=0.
        for k in range(3):
          frac_pos[i,j]+=self.atoms[i,k]*recip_cell[j,k]
      for j in range(3):
        positions[i,j] = 0.
        for k in range(3):
          positions[i,j]+=frac_pos[i,k]*transf[j,k]

    try:
      F.write( '%80.6f\n'%self.energy )
    except:
      F.write( '\n' )
    F.write( '!DATE\n' )
    F.write( 'PBC   %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n'%( a, b, c, alpha, beta, gamma ) )
    
    for i in range(self.num_atoms):
      F.write( '%2s %-13.9f  %-13.9f  %-13.9f CORE %4d %2s %2s %6.4f %4d\n'%(
           species[i], positions[i,0], positions[i,1], positions[i,2],i,
           species[i], species[i], 0, i ) )
    F.write( 'end\n' )
    F.write( 'end\n' )

    
    
    
  def position(self):
    """
    Return atomic positions flattened to a single long vector.
    """
    return self.atoms.reshape((1,-1))
  
  def force(self, mode = 'ion'):
    """
    Return forces acting on 'ion's, 'lat'tice, or 'all' of them, flattened to one long vector.
    Note that the forces on lattice are expressed in terms of change of total energy of the 
    system per (cartesian) strain - which brings them to the same units as forces acting on ions (eV/1).
    """
    if mode == 'ion':
      return self.forces.reshape((1,-1)) / self.volume
    elif mode == 'lat':
      return self.stress.flatten(0)
    else:
      ion = self.force('ion')
      lat = self.force('lat')
      return bmat('ion, lat')

  def projected_charges(self):
    if hasattr(self, "charges"):
      return array(self.zvals)-array(self.charges)
    else:
      return None
  
  def projected_magnetizations(self):
    if hasattr(self, "proj_magn"):
      return dot( self.proj_magn, inv_saxis_rotation(self.saxis) )
    else:
      return None
    
    
  def supercell(self, shape):
    """
    Generate a supercell of the given shape
    """
    l,m,n = shape
    mult = l*m*n

    supercell = calculation()

    supercell.name = "SUPERCELL %dx%dx%d %s"%(l,m,n,self.name)
    supercell.unit_cell = multiply(self.unit_cell,array([[l,l,l],[m,m,m],[n,n,n]]))
    supercell.recip_cell = linalg.inv(supercell.unit_cell)
    supercell.num_atoms = mult*self.num_atoms
    supercell.species=["%d"%(mult*int(s)) for s in self.species]

    supercell.atoms = []
 
    for atom in array(self.atoms):
      for displ in iterate_all_indices([l,m,n]):
        supercell.atoms.append( atom + dot( array(displ), array(self.unit_cell) ) )

    supercell.atoms = mat(supercell.atoms)

    return supercell
  
  @property
  def recip_cell(self):
    return self.unit_cell.I

  @property
  def cell_a(self):
    return linalg.norm(self.unit_cell[0,:])

  @property
  def cell_b(self):
    return linalg.norm(self.unit_cell[1,:])

  @property
  def cell_c(self):
    return linalg.norm(self.unit_cell[2,:])
  
  @property
  def cell_alpha(self, unit_conv=180./pi):
    return unit_conv*arccos(dot(self.unit_cell[1,:],self.unit_cell[2,:].T)[0,0]/(self.cell_b*self.cell_c))

  @property
  def cell_beta(self, unit_conv=180./pi):
    return unit_conv*arccos(dot(self.unit_cell[0,:],self.unit_cell[2,:].T)[0,0]/(self.cell_a*self.cell_c))

  @property
  def cell_gamma(self, unit_conv=180./pi):
    return unit_conv*arccos(dot(self.unit_cell[0,:],self.unit_cell[1,:].T)[0,0]/(self.cell_a*self.cell_b))

class calculation_set(loggable):
  """
  A class representing a full set of calculations read from multiple subdirectories.
  """

  """
  These two class attributes define the modulo factor by which the berry phases are regulated when 
  aligning the phases of different calculations. Normally they should be left alone. However, in case of 
  mixed set of calculations (ISPIN=1 and ISPIN=2), the way the different spins are treated might lead to 
  a shift by 0,5 of the berry phases reported for one of the spins. The standard way of aligning to the 
  first phase of groundstate calculation will likely fail (as beta-spin electrons will be contributing 
  at exactly the noncontinguous +/-0.5 region). A quick fix to this is to divide both values by half.
  A better one would be to remember the different spin components...
  """
  BPH_MAXDIST = 0.5
  BPH_CORRECTION = 1.0

  
  # Used as a threshold for seeking if the positions of polarizations need to be fixed.
  MAX_DISPLACEMENT = 0.1
  MAX_POLARIZATION_SHIFT = 0.1

  # Threshold used to recognize translational modes
  TRANSLATIONAL_MODE_THRESHOLD = 0.001

  # The stiffness of the constraint on lattice distortions used in add_lattice_constraint()
  # and ionic distortions in me_coupling. 
  CONSTRAINT_STIFFNESS = 1e8
  
  # Mapping between descriptive names and lists of calculations read.
  mode_mapping={'all': '_calculation_list',
                'ion': '_ionic_displacements_list',
                'lat': '_lattice_displacements_list'}
  
  def __init__(self):
    self._calculation_list=[]
    self._ionic_displacements_list=[]
    self._lattice_displacements_list=[]
  
  @classmethod
  def read_directory(cls, dirname, saxis, usecache=False):
    """
    Read the whole directory containing subdirectories with OUTCARS.
    This function fills in the self._calculation_list attribute.
    At the end the self.groundstate attribute is chosen from the list 
    based on the energy.
    """
    import sys
    import pickle

    calculation.saxis = saxis
   
    if usecache:
      dataname = dirname+'_pickled.dat'
      if os.access( dataname, os.R_OK ):
        F = open(dataname,'rb')
        cs = pickle.load(F)
        F.close()
        return cs

    cs = calculation_set()

    for f in os.listdir(dirname):
      path = os.path.join(dirname,f)
      if os.path.isdir( path ):   
        fpath = os.path.join( path, 'OUTCAR' )
        if os.access( fpath, os.R_OK ) or os.access( fpath+'.gz', os.R_OK ):
          c = calculation()
          cs.debug('Reading calculation from %s' % fpath, LOG_INFO)
          c.load_from_outcar( fpath )
          cs._calculation_list.append( c )

    me = 0.0
    for c in cs._calculation_list:
      if c.energy < me:
        me = c.energy
        cs.groundstate = c
        
    if usecache:
      F = open(dataname,'wb')
      pickle.dump(cs,F,-1)
      F.close()

    return cs

  def add_calculation(self, c ):
    self._calculation_list.append( c )
    
  def try_fix_displacements(self):
    """
    Fix displacements within the calculation set. Ions are moved around by 
    lattice displacements in order to fit within self.MAX_DISPLACEMENT 
    from the groundstate (so the groundstate must be properly set before 
    calling this function).
    """
              
    def fix_one_calculation( item ):
      for idx,atom in enumerate(item.atoms):
        gs_atom = self.groundstate.atoms[idx]
        for dx,dy,dz in iterate_shifts():
          new_atom = atom+[[dx,dy,dz]]*item.unit_cell
          if metro_dist(new_atom, gs_atom) <= self.MAX_DISPLACEMENT:
            self.debug('    moving atom %d by %d,%d,%d' % (idx,dx,dy,dz),LOG_ALLINFO)
            item.atoms[idx] = new_atom[0]
      

    gs_position = self.groundstate.position()
    for item in self.calculations('ion'):
      if metro_dist(item.position(), gs_position) > self.MAX_DISPLACEMENT:
        self.debug("   fixing %s" % item.fileID, LOG_ALLINFO)
        fix_one_calculation( item )
        item.recalculate_polarization()
          
    self.debug("Corrected displacements", LOG_ALLINFO)
    self.debug(mat2str(self.displacements('ion')), LOG_ALLINFO)


  def try_fix_polarizations(self):
    """
    Try to fix the polarization if it is screwed by polarization quantum.
    """
    if self.groundstate.berry_phase[0] is None:
      return
    
    def fix_berry_phases( calc ):
      changed = False
      for kdir in range(3):
        reference_phase = self.groundstate.berry_phase[kdir][0,0]
        for ph in calc.berry_phase[kdir]:
          while ph[0]-reference_phase > calculation_set.BPH_MAXDIST:
            ph[0] -= calculation_set.BPH_CORRECTION
            changed = True
          while reference_phase-ph[0] > calculation_set.BPH_MAXDIST:
            ph[0] += calculation_set.BPH_CORRECTION
            changed = True
      return changed
        
    for calc in self._calculation_list:
      #For some reason this does not work as expected and it's safer to just fix/recalculate all calculations
      #if fix_berry_phases(calc):
        #p = calc.recalculate_polarization()
        #self.debug('Corrected polarization %s = %f %f %f'%(calc.fileID,p[0,0],p[0,1],p[0,2]), LOG_ALLINFO)
      fix_berry_phases(calc)
      p = calc.recalculate_polarization()
      self.debug('Corrected polarization %s = %f %f %f'%(calc.fileID,p[0,0],p[0,1],p[0,2]), LOG_ALLINFO)
    
  def try_polarizations_chain(self):
    def align_berry_phases( ref_calc, dest_calc ):
      for kdir in range(3):
        if ref_calc.berry_phase[kdir].shape != dest_calc.berry_phase[kdir].shape:
          self.debug("Cannot align berry phases", LOG_WARNING)
        else:
          print "IN"
          print ref_calc.berry_phase[kdir]
          print dest_calc.berry_phase[kdir]
          for i in range(ref_calc.berry_phase[kdir].shape[0]):
            while ref_calc.berry_phase[kdir][i,0] - dest_calc.berry_phase[kdir][i,0] > 0.5:
              dest_calc.berry_phase[kdir][i,0] += 1.0
            while ref_calc.berry_phase[kdir][i,0] - dest_calc.berry_phase[kdir][i,0] < -0.5:
              dest_calc.berry_phase[kdir][i,0] -= 1.0
          print "OUT"
          print dest_calc.berry_phase[kdir]
          
    
    for i in range(1,len(self._calculation_list)):
      print "align", self._calculation_list[i-1].fileID, self._calculation_list[i].fileID
      align_berry_phases(self._calculation_list[i-1],self._calculation_list[i])
      self._calculation_list[i].recalculate_polarization()
    

  def set_groundstate(self, pattern):
    """
    Set the self.groundstate attribute to a given file pattern. This is usefull for constant 
    volume calculations, where calculated energy is not be minimal with respect to lattice deformations 
    at the optimised structure.
    """
    for calc in self._calculation_list:
      if re.match(pattern,calc.fileID):
        self.groundstate = calc
        return
      
    self.debug("Cannot set groundstate to %s" % pattern, LOG_WARNING)
      
  def calculations(self,mode='all'):
    """
    Returns only a subset of calculations. In this way one can easily work on 'ion'ic or 'lat'tice
    distortions and treat them separately where neccessary.
    """
    if mode == 'ion':
      return self._ionic_displacements_list
    elif mode == 'lat':
      return self._lattice_displacements_list
    else:
      return concatenate( [ self._ionic_displacements_list, self._lattice_displacements_list ] )
  
  def set_subset(self,mode,pattern):
    """
    Set the subset of calculations based on the file pattern.
    """
    if mode == 'all':
      return
    
    current_attr = self.calculations(mode)
    del current_attr[:]
    for c in self._calculation_list:
      if re.match(pattern,c.fileID):
        current_attr.append(c)
    
  def set_ionic(self, pattern):
    self.set_subset('ion',pattern)
    
  def set_lattice(self, pattern):
    self.set_subset('lat',pattern)
      
  def displacements(self,mode='ion', calcset='ion'):
    """
    Returns displacements from the groundstate for all the calculations from the list.
    This method returns fractional coordinate displacements of ions in case of mode='ion' 
    and strain in case of mode='lat' for a specified calculation subset. In case of mode='all' 
    it returns full generalised displacement matrix.
    """
    if mode == 'ion':
      return mat(concatenate([ calc.position() - self.groundstate.position() for calc in self.calculations(calcset) ]))
    elif mode == 'lat':
      return mat(concatenate([ (self.groundstate.recip_cell*calc.unit_cell - eye(3)).flatten(1) for calc in self.calculations(calcset) ]))
    else:
      return hstack([self.displacements('ion', calcset), self.displacements('lat', calcset)])
    
  def forces(self, mode='ion', calcset='ion'):
    """
    Returns forces exerted on 'mode' subset of coordinates, as read from 'calcset' subset of all calculations. 
    In new formalism, all global properties are defined as per unit volume and that stresses are already represented like this.
    """
    return mat(concatenate([ calc.force(mode) for calc in self.calculations(calcset)]))
  
  def magnetizations(self,calcset='all'):
    """
    Returns magnetizations read from 'calcset' subset of all calculations. 
    In new formalism, all global properties are defined as per unit volume
    """
    return mat(concatenate([ c.magnetization for c in self.calculations(calcset)])) / self.groundstate.volume
  
  def magnetization_differences(self,calcset='all'):
    """
    Returns magnetization differences between calcset of calculation and the groundstate.
    """
    #if calcset == 'ion':
    return mat(concatenate([ c.magnetization-self.groundstate.magnetization for c in self.calculations(calcset)]))
    #elif calcset == 'lat':
    #  return mat(concatenate([ c.magnetization*c.recip_cell*self.groundstate.unit_cell-self.groundstate.magnetization for c in self.calculations(calcset)]))
  
  def polarizations(self,calcset='all'):
    """
    Returns polarizations read from 'calcset' subset of all calculations. 
    """
    return concatenate([ c.polarization for c in self.calculations(calcset)])
  
  def polarization_differences(self,calcset='all'):
    """
    Returns polarization differences between calcset of calculation and the groundstate.
    Note that in case of lattice distortions, comparing the polarizations directly does not 
    make much sense as they differ by the factor induced by volume cell differences and cell
    rotations. Therefore this function returns difference of berry-phases multiplied by 
    groundstate unit cell for lattice distortion subset of calculations.
    This implements the derivative of berry phases over strain formula (24) from 
    D. Vanderbilt, J Phys Chem Solid 61 (2000) 147-151   
    """
    if calcset=='ion':
      return concatenate( [ calc.polarization-self.groundstate.polarization for calc in self.calculations(calcset) ] )
    elif calcset=='lat':
      return concatenate( [ -(calc.berry_term()+calc.ev_recip-self.groundstate.berry_term()-self.groundstate.ev_recip)*self.groundstate.unit_cell/self.groundstate.volume for calc in self.calculations('lat') ] )
    else:
      return vstack([self.polarization_differences('ion'),self.polarization_differences('lat')])
  
  def force_constant_matrix(self, mode='ion', calcset='ion'):
    """
    Calculate force constant matrix (A) and residual forces (F0) at groundstate from solving set of linear
    equations involving forces (F) and displacements (d):
    F = d*A + F0
    The set of equations might be overdefined - the found solution is optimal in the leastsquare sense.
    
    Second returned value is the residual of the force at groundstate (which should be 0)
    """
    
    # FIXME: this is a workaround for quite strange behaviour...
    # When given all the calculations, the lattice part of the force constant
    # matrix becomes asymmetric and generally noisy. This does not happen when
    # only subsets are given. Even off-diagonal (ion-lattice) coefficients behave
    # very well, when only calculated on separate calculation sets.
    if calcset == 'all':
      return vstack( [ self.force_constant_matrix(mode, 'ion')[0],
                       self.force_constant_matrix(mode, 'lat')[0] ]), 0
    
    displacements = self.displacements(calcset, calcset)
    forces = self.forces(mode, calcset)
      
    (fc_matrix, rs, rk, s) = linalg.lstsq( 
      hstack( [ displacements, ones((displacements.shape[0],1)) ] ), forces )
    residual = fc_matrix[-1,:]

    # remember that when we talk about force constant matrix, we think actually second derivative
    # of the energy... which is exactly opposite sign to derivative of forces 
    fc_matrix = -fc_matrix[:-1,:]
    if hasattr(self,'_lattice_constraint') and (mode=='all' or calcset=='all' or (mode=='lat' and calcset=='lat')):
      fc_matrix[-9:,-9:] += self.CONSTRAINT_STIFFNESS*self._lattice_constraint 
    return fc_matrix, residual
         
  def magnetic_polarization_matrix(self, mode='ion'):
    """
    Returns numerically computed derivative of total magnetization vector over general coordinates.
    """
    if mode == 'ion' or mode == 'lat':
      displacements = self.displacements(mode, mode)
      mag_dif = self.magnetization_differences(mode)
      (mag_pol_mat, rs, rk, s) = linalg.lstsq( 
                                 hstack( [ displacements, ones((displacements.shape[0],1)) ] ), mag_dif )
      residual = mag_pol_mat[-1,:]
      mag_pol_mat = mag_pol_mat[:-1.:]
      return mag_pol_mat/self.groundstate.volume, residual
    else:
      impl,ir = self.magnetic_polarization_matrix('ion')
      lmpl,lr = self.magnetic_polarization_matrix('lat')
      return vstack([impl, lmpl]), ir+lr

  def electric_polarization_matrix(self, mode='ion'):
    """
    Returns numerically computed derivative of polarization vector over general coordinates.
    """
    if mode == 'ion' or mode == 'lat':
      displacements = self.displacements(mode, mode)
      pol_dif = self.polarization_differences(mode)
      (el_pol_mat, rs, rk, s) = linalg.lstsq( 
                                hstack( [ displacements, ones((displacements.shape[0],1)) ] ), pol_dif )
      residual = el_pol_mat[-1,:]
      el_pol_mat = el_pol_mat[:-1,:]
    
      if (mode == 'ion'):
        return el_pol_mat, residual
      else:
        return el_pol_mat, residual
    else:
      iepl,ir = self.electric_polarization_matrix('ion')
      lepl,lr = self.electric_polarization_matrix('lat')
      return vstack([iepl, lepl]), ir+lr
  
  def clear_lattice_constraints(self):
    """
    Clear all imposed lattice distortion constraints.
    """
    if hasattr(self,'_lattice_constraint'):
      del self._lattice_constraint
      
  def add_lattice_constraint(self, dir_plane_mat, flags=['xyshear','oopshear']):
    """
    Constrain the lattice from deforming in a given cartesian direction.
    """
    direction = mat(dir_plane_mat)
    
    if direction.shape[0] == 1:
      # We have a vector, which defines uncompressible direction
      direction /= sqrt(direction*direction.T)
      strain = direction.T*direction
    elif direction.shape[0] == 2:
      # We have two vectors which define rigid plane
      direction[0] /= sqrt(direction[0]*direction[0].T)
      direction[1] /= sqrt(direction[1]*direction[1].T)
      self.add_lattice_constraint(direction[0])                # XX
      self.add_lattice_constraint(direction[1])                # YY
      if 'xyshear' in flags:
        self.add_lattice_constraint(direction[1].T*direction[0]) # YX
        self.add_lattice_constraint(direction[0].T*direction[1]) # XY
      out_of_plane = cross(direction[0],direction[1])
      self.add_lattice_constraint(direction[0].T*out_of_plane) # XZ
      self.add_lattice_constraint(direction[1].T*out_of_plane) # YZ
      if 'oopshear' in flags:
        self.add_lattice_constraint(out_of_plane.T*direction[0]) # ZX
        self.add_lattice_constraint(out_of_plane.T*direction[1]) # ZY
      return
    elif direction.shape[0] == 3:
      # We have directly strain component generalised direction
      strain = direction
      
    strain = mat(strain.flatten(0))
    constraint = strain.T*strain
    
    if not hasattr(self,'_lattice_constraint'):
      self._lattice_constraint = mat(zeros((9,9)))
      
    self._lattice_constraint += constraint
    
    