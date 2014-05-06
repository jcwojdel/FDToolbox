import re
from numpy import *
from utility import *

class atomic_set(loggable):
  """
  A class for keeping the information about set of atoms in a unit cell.
  - self.unit_cell - unit cell (3x3 A)
  - self.recip_cell - reciprocal unit cell (3x3 A**-1)
  - self.volume - volume calculated from the unit cell (1 A**3)
  - self.name - name of the system as put in POSCAR
  - self.num_per_type - list of species numbers 
  - self.num_atoms - number of atoms
  - self.atoms - atomic positions internal coordinates (self.num_atomsx3 1)
  """  
  def get__volume(self):
    if not hasattr(self,'__volume'):
      self.__volume = linalg.det(self.unit_cell)
    return self.__volume    
  volume = property( get__volume, None )        
    
  def load_from_poscar(self, filename):
    """
    Loads the system from POSCAR. Obviously, everything except for atomic positions is missing.
    """
    with open( filename, 'r' ) as F:
      F = open( filename, 'r' )
      self.fileID = filename
      self.name = F.readline()
      scale = float(F.readline())
      self.unit_cell = mat( fromfile( F, dtype('float'), 9, ' ' ).reshape((3,3)) )
      
      # Scale < 0 means that it gives the volume we want to have
      if scale < 0.0:
        scale=(-scale/self.volume)**(1./3.)
      self.unit_cell *= scale

      # If the next line does not contain just numbers, then it is treated as a list of species
      line = F.readline()
      self.species = None
      try:
        self.num_per_type = [int(n) for n in line.split()]
      except:
        species = line.split()
        line = F.readline()
        self.num_per_type = [int(n) for n in line.split()]
        self.species = []
        for n in self.num_per_type:
          self.species.extend(n*[species[0]])
          species = species[1:]
        
      self.num_atoms = 0
      self.num_atoms=sum(self.num_per_type)
        
      mode = F.readline()
  
      self.atoms = mat(fromfile( F, dtype('float'),  self.num_atoms*3, ' ' ).reshape(( self.num_atoms,3)))
  
      if re.search('^[cCkK]',mode):
        pass
      else:
        self.atoms = self.atoms*self.unit_cell
      
    if self.name.split()[0] == "SUPERCELL":
      self.is_supercell = True
      self.supercell_repetitions = self.name.split()[1].split('x')
      self.supercell_repetitions = [int(i) for i in self.supercell_repetitions]
    
  def save_to_poscar(self, filename,direct=False,species_line=False):
    """
    Save the system to POSCAR
    """ 
    with open( filename, 'w' ) as F:
      F.write( self.name )
      F.write( "    1.0\n" )
      F.write( mat2str( self.unit_cell, "%16.10f" ) )
      if species_line:
        pos = 0
        for n in self.num_per_type: 
          F.write('%s '%self.species[pos])
          pos += n
        F.write('\n')
      F.write(' '.join([str(n) for n in self.num_per_type]) )
      F.write('\n')
      if not direct:
        F.write("Cart\n")
        F.write( mat2str( self.atoms, "%16.10f" ) )
      else:
        F.write("Direct\n")
        F.write( mat2str( dot(self.atoms,self.recip_cell), "%16.10f" ) )

  def save_to_xyz(self, filename):
    """
    Save the system to XYZ
    """ 
    with open( filename, 'a' ) as F:
      F = open( filename, 'a' )
      F.write( '%d\n'%self.num_atoms )
      F.write( "XYZ\n" )
      for num,row in enumerate(self.atoms):
        try:
          F.write('%s  '%self.species[num])
        except:
          F.write('X%d '%num)
        F.write( mat2str( row, "%16.10f" ) )
      F.write( "\n" )
    
  def save_to_arc(self, filename, header = True, comment = None):
    """
    Save the system in ARC file
    """
    if header:
      F = open( filename, 'w' )
      F.write( "!BIOSYM archive 2\n" )
      if comment is not None:
        F.write( '!%s\n'%comment )
      F.write( "PBC=ON\n" )
    else:
      F = open( filename, 'a' )
            
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
           self.species[i], positions[i,0], positions[i,1], positions[i,2],i,
           self.species[i], self.species[i], 0, i ) )
    F.write( 'end\n' )
    F.write( 'end\n' )    
    
  def position(self):
    """
    Return atomic positions flattened to a single long vector.
    """
    return self.atoms.reshape((1,-1))    
    
  def supercell(self, shape):
    """
    Generate a supercell of the given shape
    """
    l,m,n = shape
    mult = l*m*n

    supercell = atomic_set()

    supercell.name = "SUPERCELL %dx%dx%d %s"%(l,m,n,self.name)
    supercell.unit_cell = multiply(self.unit_cell,array([[l,l,l],[m,m,m],[n,n,n]]))
    supercell.recip_cell = linalg.inv(supercell.unit_cell)
    supercell.num_atoms = mult*self.num_atoms
    supercell.num_per_type=["%d"%(mult*int(s)) for s in self.num_per_type]
    supercell.species=[]
    for i,n in enumerate(self.num_per_type):
      supercell.species.extend(n*self.species[i])

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
  
  @property
  def fractional_atoms(self):
    return dot(self.atoms,self.recip_cell)
    
  def align_to(self,reference):    
    shift = (self.fractional_atoms - reference.fractional_atoms)
    ishift= shift.round()
    if abs(shift-ishift).max() > 0.4:
      self.debug('Atomic displacement larger than .4 - possible reshuffling of atoms', LOG_WARNING)
      print shift
    
    self.atoms -= dot(ishift,self.unit_cell)
    
  def save_to_res(self,filename):
    # Open output file
    with open(filename,"w") as outfile:    
      # Title
      outfile.write('TITL %s'%self.name)
      
      # Cell
      outfile.write("CELL %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n" %\
                    (0.0,self.cell_a,self.cell_b,self.cell_c,\
                    self.cell_alpha,self.cell_beta,self.cell_gamma))
      outfile.write("LATT -1\n")
      
      outfile.write("SFAC ")
      pos = 0
      for n in self.num_per_type: 
        outfile.write(' %2s'%self.species[pos])
        pos += n
      outfile.write('\n')
      
      last_sp=''
      atom_type=0
      for i in range(self.num_atoms):
        outfile.write('%2s  ' % self.species[i])
        if (self.species[i]!=last_sp):
          atom_type = atom_type+1
          last_sp   = self.species[i]
        outfile.write(' %2i' % atom_type)
        outfile.write(' %10.5f %10.5f %10.5f %10.5f %9.5f\n' % \
                            (self.fractional_atoms[i,0],self.fractional_atoms[i,1],self.fractional_atoms[i,2],1.0,0.0) )
      
      outfile.write("END\n")


    